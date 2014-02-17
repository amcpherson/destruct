from collections import *
import pandas as pd
import numpy as np

class Node(object):
    def __init__(self, graph, id):
        self.graph = graph
        self.id = id
        edges = self.graph.node_edge.loc[self.id]
        self.in_edges = edges[edges.in_edge == True].edge_idx
        self.out_edges = edges[edges.in_edge == False].edge_idx
    def get_in_edges(self):
        return self.in_edges
    def get_out_edges(self):
        return self.out_edges

class SequenceEdge(object):
    def __init__(self, graph, edge_idx, read_count_idx, read_count_data):
        self.graph = graph
        self.edge_idx = edge_idx
        self.read_count_idx = read_count_idx
        self.read_count_data = read_count_data
    @property
    def id(self):
        return self.read_count_data.id[self.read_count_idx]
    @property
    def chromosome1(self):
        return self.read_count_data.chromosome1[self.read_count_idx]
    @property
    def chromosome2(self):
        return self.read_count_data.chromosome2[self.read_count_idx]
    @property
    def position1(self):
        return self.read_count_data.position1[self.read_count_idx]
    @property
    def position2(self):
        return self.read_count_data.position2[self.read_count_idx]
    @property
    def readcount(self):
        return self.read_count_data.readcount[self.read_count_idx].astype(float)
    @property
    def is_disabled(self):
        return self.edge_idx in self.graph.disabled_edge_idxs
    @property
    def unregularized(self):
        return True

class IntervalEdge(SequenceEdge):
    def __init__(self, graph, edge_idx, read_count_idx, read_count_data):
        SequenceEdge.__init__(self, graph, edge_idx, read_count_idx, read_count_data)
    @property
    def edge_type(self):
        return 'interval'
    @property
    def region_length(self):
        return self.read_count_data.length[self.read_count_idx].astype(float)

class ReferenceEdge(SequenceEdge):
    def __init__(self, graph, edge_idx, read_count_idx, read_count_data):
        SequenceEdge.__init__(self, graph, edge_idx, read_count_idx, read_count_data)
    @property
    def edge_type(self):
        return 'reference'
    @property
    def region_length(self):
        return 0.0
    @property
    def readcount(self):
        return 0.0

class VariantEdge(SequenceEdge):
    def __init__(self, graph, edge_idx, read_count_idx, read_count_data):
        SequenceEdge.__init__(self, graph, edge_idx, read_count_idx, read_count_data)
    @property
    def edge_type(self):
        return 'variant'
    @property
    def region_length(self):
        return 0.0
    @property
    def readcount(self):
        return 0.0
    @property
    def unregularized(self):
        return False

class TelomereEdge(object):
    def __init__(self, graph, edge_idx, telomere_idx, telomere_data):
        self.graph = graph
        self.edge_idx = edge_idx
        self.telomere_idx = telomere_idx
        self.telomere_data = telomere_data
    @property
    def id(self):
        return self.telomere_data.id[self.telomere_idx]
    @property
    def chromosome(self):
        return self.telomere_data.chromosome[self.telomere_idx]
    @property
    def chromosome1(self):
        return self.telomere_data.chromosome[self.telomere_idx]
    @property
    def chromosome2(self):
        return self.telomere_data.chromosome[self.telomere_idx]
    @property
    def position(self):
        return self.telomere_data.position[self.telomere_idx]
    @property
    def position1(self):
        return self.telomere_data.position[self.telomere_idx]
    @property
    def position2(self):
        return self.telomere_data.position[self.telomere_idx]
    @property
    def edge_type(self):
        return 'telomere'
    @property
    def readcount(self):
        return 0.0
    @property
    def region_length(self):
        return 0.0
    @property
    def is_disabled(self):
        return self.edge_idx in self.graph.disabled_edge_idxs
    @property
    def unregularized(self):
        return self.telomere_data.is_normal[self.telomere_idx]

class BreakpointGraph(object):
    def __init__(self, interval_data, reference_data, variant_data, telomere_data):
        # Load edge data tables
        self.interval_data = interval_data
        self.reference_data = reference_data
        self.variant_data = variant_data
        self.telomere_data = telomere_data
        # Create edge list and node edge list
        self.edges = list()
        self.node_edge = list()
        for idx, row in self.interval_data.iterrows():
            edge_idx = len(self.edges)
            self.edges.append(IntervalEdge(self, edge_idx, idx, self.interval_data))
            for s in ('1', '2'):
                self.node_edge.append((row['chromosome'+s], row['position'+s], row['strand'+s], edge_idx, True))
        for idx, row in self.reference_data.iterrows():
            edge_idx = len(self.edges)
            self.edges.append(ReferenceEdge(self, edge_idx, idx, self.reference_data))
            for s in ('1', '2'):
                self.node_edge.append((row['chromosome'+s], row['position'+s], row['strand'+s], edge_idx, False))
        for idx, row in self.variant_data.iterrows():
            edge_idx = len(self.edges)
            self.edges.append(VariantEdge(self, edge_idx, idx, self.variant_data))
            self.node_edge.append((row['chromosome1'], row['position1'], row['strand1'], edge_idx, False))
            self.node_edge.append((row['chromosome2'], row['position2'], row['strand2'], edge_idx, False))
        for idx, row in self.telomere_data.iterrows():
            edge_idx = len(self.edges)
            self.edges.append(TelomereEdge(self, edge_idx, idx, self.telomere_data))
            self.node_edge.append((row['chromosome'], row['position'], row['strand'], edge_idx, False))
        # Create node edge dataframe for lookups
        self.node_edge = pd.DataFrame(self.node_edge, columns=['chromosome', 'position', 'strand', 'edge_idx', 'in_edge'])
        self.node_edge.set_index(['chromosome', 'position', 'strand'], inplace=True)
        self.node_edge.sort_index(inplace=True)
        # Create node list
        self.nodes = list()
        for idx, row in self.node_edge.reset_index()[['chromosome', 'position', 'strand']].drop_duplicates().iterrows():
            node_id = (row['chromosome'], row['position'], row['strand'])
            self.nodes.append(Node(self, node_id))
        # Default empty disabled edges
        self.disabled_edge_idxs = set()
    def set_disabled_edge(self, idx):
        self.disabled_edge_idxs.add(idx)
    def clear_disabled_edges(self):
        self.disabled_edge_idxs = set()
