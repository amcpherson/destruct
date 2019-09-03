import collections
import itertools
import numpy as np
import pandas as pd

import networkx
import blossomv.blossomv


BalancedRearrangment = collections.namedtuple(
    'BalancedRearrangment',
    [
        'deleted_length',
        'duplicated_length',
        'prediction_ids',
    ],
)


def detect_balanced_rearrangements(
    breakpoints,
    dec_nt_per_break=2000.,
    inc_nt_per_break=500.,
    cost_resolution=1000.,
):
    if len(breakpoints.index) == 0:
        return []

    G = networkx.Graph()

    breakpoints = breakpoints.copy()

    # Algorithm cannot handle loops created by perfect foldbacks
    # for these events, move one breakpoint by one nucleotide
    is_perfect_foldback = (
        (breakpoints['chromosome_1'] == breakpoints['chromosome_2']) &
        (breakpoints['position_1'] == breakpoints['position_2']) &
        (breakpoints['strand_1'] == breakpoints['strand_2']))
    breakpoints.loc[is_perfect_foldback, 'position_2'] += 1

    break_ends = pd.DataFrame(
        {
            'chromosome': np.concatenate([breakpoints['chromosome_1'].values, breakpoints['chromosome_2'].values]),
            'position': np.concatenate([breakpoints['position_1'].values, breakpoints['position_2'].values]),
            'strand': np.concatenate([breakpoints['strand_1'].values, breakpoints['strand_2'].values]),
        }
    )
    break_ends.drop_duplicates(inplace=True)

    # Add break end nodes
    for idx in break_ends.index:
        break_end = (
            break_ends.loc[idx, 'chromosome'],
            break_ends.loc[idx, 'position'],
            break_ends.loc[idx, 'strand'],
        )
        G.add_node(break_end)

    # Add breakpoint edges
    for idx in breakpoints.index:
        break_end_1 = (
            breakpoints.loc[idx, 'chromosome_1'],
            breakpoints.loc[idx, 'position_1'],
            breakpoints.loc[idx, 'strand_1'],
        )
        break_end_2 = (
            breakpoints.loc[idx, 'chromosome_2'],
            breakpoints.loc[idx, 'position_2'],
            breakpoints.loc[idx, 'strand_2'],
        )
        prediction_id = breakpoints.loc[idx, 'prediction_id']
        G.add_edge(break_end_1, break_end_2, edge_type='breakpoint', prediction_id=prediction_id)

    # Add reference and segment edges
    for chromosome, chrom_break_ends in break_ends.groupby('chromosome'):
        chrom_break_ends = chrom_break_ends.sort_values('position')

        # Add reference edges
        for position in chrom_break_ends['position'].values:
            break_end_1 = (chromosome, position, '+')
            break_end_2 = (chromosome, position, '-')
            G.add_edge(break_end_1, break_end_2, edge_type='reference')

        # Add segment edges
        for start, end in zip(chrom_break_ends['position'].values[:-1], chrom_break_ends['position'].values[1:]):
            break_end_1 = (chromosome, start, '-')
            break_end_2 = (chromosome, end, '+')
            length = end - start
            G.add_edge(break_end_1, break_end_2, edge_type='segment', length=length)

    # Create a genome modification graph
    #  - identical node set as genome graph
    #  - create a +1 and -1 signed edge for each original edge
    #  - color each edge:
    #     - segment edges: +1 -> red, -1 -> blue
    #     - bond edges: +1 -> blue, -1 -> red
    #  - add a cost for each edge:
    #     - segment edges: cost for increasing / decreasing per nt
    #     - breakpoint edges: decreasing inf, increasing -1

    # Create graph with duplicated edges for each sign
    H = networkx.MultiGraph()
    H.add_nodes_from(G)
    for edge in G.edges():
        H.add_edge(*edge, attr_dict=G.get_edge_data(*edge), sign=1)
        H.add_edge(*edge, attr_dict=G.get_edge_data(*edge), sign=-1)

    # Color each edge, +1 for red, -1 for blue
    # add cost for each edge
    for edge in H.edges():
        for multi_edge_idx, edge_attr in H[edge[0]][edge[1]].items():
            edge_type = edge_attr['attr_dict']['edge_type']
            edge_attr['cost'] = 0
            if edge_type == 'segment':
                edge_attr['color'] = edge_attr['sign']
                if edge_attr['sign'] == 1:
                    edge_attr['cost'] = float(edge_attr['attr_dict']['length'] / inc_nt_per_break)
                elif edge_attr['sign'] == -1:
                    edge_attr['cost'] = float(edge_attr['attr_dict']['length'] / dec_nt_per_break)
            elif edge_type == 'breakpoint':
                edge_attr['color'] = -edge_attr['sign']
                if edge_attr['sign'] == 1:
                    edge_attr['cost'] = -1.
                elif edge_attr['sign'] == -1:
                    edge_attr['cost'] = np.inf
            elif edge_type == 'reference':
                edge_attr['color'] = -edge_attr['sign']
                edge_attr['cost'] = 1. / cost_resolution
            else:
                raise ValueError('unknown edge type {}'.format(edge_type))

    # Create matching graph
    #  - duplicate nodes, one set red, one set blue
    #  - add transverse edges (v_red, v_blue)
    #  - for each original edge:
    #    - add (u_red, v_red) for red edges
    #    - add (u_blue, v_blue) for blue edges
    #    - replacate edge costs

    transverse_edge_cost = 1. / cost_resolution

    M = networkx.Graph()
    for node in H.nodes():
        transverse_edge = []
        for color in (1, -1):
            colored_node = node + (color,)
            M.add_node(colored_node)
            transverse_edge.append(colored_node)
        M.add_edge(*transverse_edge, cost=transverse_edge_cost)

    for edge in H.edges():
        for multi_edge_idx, edge_attr in H[edge[0]][edge[1]].items():
            cost = edge_attr['cost']
            if np.isinf(cost):
                continue
            color = edge_attr['color']
            colored_node_1 = edge[0] + (color,)
            colored_node_2 = edge[1] + (color,)
            M.add_edge(colored_node_1, colored_node_2, attr_dict=edge_attr, cost=cost)

    M1 = networkx.convert_node_labels_to_integers(M, label_attribute='node_tuple')

    # Min cost perfect matching
    edges = networkx.get_edge_attributes(M1, 'cost')
    for edge in edges.keys():
        if edge[0] == edge[1]:
            raise Exception('self loop {}'.format(M1[edge[0]][edge[1]]))
        edges[edge] = int(edges[edge] * cost_resolution)
    min_cost_edges = blossomv.blossomv.min_weight_perfect_matching(edges)

    # Remove unselected edges
    assert set(min_cost_edges).issubset(edges.keys())
    remove_edges = set(edges.keys()).difference(min_cost_edges)
    M2 = M1.copy()
    M2.remove_edges_from(remove_edges)

    # Re-create original graph with matched edges
    M3 = networkx.relabel_nodes(M2, mapping=networkx.get_node_attributes(M2, 'node_tuple'))

    # Create subgraph of H with only selected edges
    H1 = networkx.Graph()
    for edge in M3.edges():
        edge_attr = M3[edge[0]][edge[1]]
        node_1 = edge[0][:-1]
        node_2 = edge[1][:-1]
        if node_1 == node_2:
            continue
        if H1.has_edge(node_1, node_2):
            H1.remove_edge(node_1, node_2)
        else:
            H1.add_edge(node_1, node_2, attr_dict=edge_attr)

    # Get individual events
    rearrangements = []
    for C in networkx.connected_component_subgraphs(H1, copy=False):
        duplicated_length = 0
        deleted_length = 0
        prediction_ids = []

        for edge in C.edges():
            edge_type = C[edge[0]][edge[1]]['attr_dict']['attr_dict']['attr_dict']['edge_type']

            if edge_type == 'segment':
                sign = C[edge[0]][edge[1]]['attr_dict']['attr_dict']['sign']
                length = C[edge[0]][edge[1]]['attr_dict']['attr_dict']['attr_dict']['length']
                if sign == 1:
                    duplicated_length += length
                else:
                    deleted_length += length

            elif edge_type == 'breakpoint':
                prediction_id = C[edge[0]][edge[1]]['attr_dict']['attr_dict']['attr_dict']['prediction_id']
                prediction_ids.append(prediction_id)

        if len(prediction_ids) <= 1:
            continue

        rearrangements.append(
            BalancedRearrangment(
                deleted_length,
                duplicated_length,
                prediction_ids,
            )
        )

    return rearrangements

