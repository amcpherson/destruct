import math
from collections import *
import pypeliner

import pybreakcopies



class EulerianLP(object):

    def __init__(self):
        self.lambda_variant = 0.0
        self.lambda_telomere = 0.0

    def solve(self, bg):

    # edge_info.read_count = python::extract<double>(edge_tuple[0]);
    # edge_info.region_length = python::extract<double>(edge_tuple[1]);
    # edge_info.obj_coeff = python::extract<double>(edge_tuple[2]);
    # edge_info.upper_bound = python::extract<double>(edge_tuple[3]);

        edge_infos = list()

        for edge in bg.edges:

            # Calculate objective coefficient based on regularization
            obj_coeff = 0.0
            if edge.edge_type == 'variant' and not edge.unregularized:
                obj_coeff = self.lambda_variant
            elif edge.edge_type == 'telomere' and not edge.unregularized:
                obj_coeff = self.lambda_telomere

            # Calculate upper bound for enabled/disabled edges
            upper_bound = 1000.0
            if edge.is_disabled:
                upper_bound = 0.0

            edge_infos.append((edge.read_count, edge.region_length, obj_coeff, upper_bound))

    # int node_idx = python::extract<int>(node_tuple[0]);

    # NodeInfo node_info;
    # node_info.edge_idx = python::extract<int>(node_tuple[1]);
    # node_info.edge_coeff = python::extract<double>(node_tuple[2]);

        node_infos = list()

        for node_idx, node in enumerate(bg.nodes):

            for edge_idx in node.get_in_edges():
                node_infos.append((node_idx, edge_idx, 1.0))

            for edge_idx in node.get_out_edges():
                node_infos.append((node_idx, edge_idx, -1.0))
                
        self.edge_copies = pybreakcopies.solve(edge_infos, node_infos)


SolutionInfo = namedtuple('SolutionInfo', ['edge_copies',
                                           'edge_types',
                                           'edge_ids',
                                           'lambda_variant',
                                           'lambda_telomere'])
SolutionInfo.__repr__ = lambda self: 'SolutionInfo'

def solve_lasso(bg, lambda_variant, lambda_telomere):

    lp = EulerianLP()

    lp.lambda_variant = lambda_variant
    lp.lambda_telomere = lambda_telomere
    
    lp.solve(bg)

    return SolutionInfo([copies for copies in lp.edge_copies],
                        [edge.edge_type for edge in bg.edges],
                        [edge.id for edge in bg.edges],
                        lambda_variant,
                        lambda_telomere)



