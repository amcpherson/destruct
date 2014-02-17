import math
from collections import *
import pypeliner

import pycplex


class cplex_job(object):
    def __init__(self, cplex_variables, cplex_constraints, verbose, threads):
        self.cplex_variables = cplex_variables
        self.cplex_constraints = cplex_constraints
        self.verbose = verbose
        self.threads = threads
    def __call__(self):
        self.objective_value, self.solution_values = pycplex.solve_min(self.cplex_variables, self.cplex_constraints, self.verbose, self.threads)

cplex_threads = 8

def cplex_solve_min(cplex_variables, cplex_constraints, verbose=True):
    return pycplex.solve_min(cplex_variables, cplex_constraints, verbose, cplex_threads)
    result = pypeliner.delegator.call_external(cplex_job(cplex_variables, cplex_constraints, verbose, cplex_threads))
    return result.objective_value, result.solution_values

def create_piecewise_linear(dispersion, region_count, region_length):

    numLeftBreakpoints = 7
    midBreakpoint = numLeftBreakpoints
    numBreakpoints = numLeftBreakpoints * 2 + 1

    if region_count == 0.0:
        return [(region_length, 0.0)]

    piecewise = list()

    # Build x values
    xVals = [0.0] * (numBreakpoints + 2)
    xVals[midBreakpoint + 1] = region_count / region_length
    xVals[midBreakpoint + 1] = max(0.001, xVals[midBreakpoint + 1])
    for k in range(midBreakpoint, -1, -1):
        xVals[k] = xVals[k + 1] / 1.2
    for k in range(midBreakpoint + 2, numBreakpoints + 2):
        xVals[k] = xVals[k - 1] * 1.2

    # Build y values
    yVals = [0.0] * (numBreakpoints + 2)
    for k in range(0, numBreakpoints + 2):
        if dispersion == 0.0:
            yVals[k] = region_length * xVals[k] - region_count * math.log(xVals[k])
        else:
            kappa = 1.0 / dispersion
            yVals[k] = -region_count * math.log(kappa * region_length * xVals[k] / (1.0 + kappa * region_length * xVals[k])) - (1.0 / kappa) * math.log(1.0 / (1.0 + kappa * region_length * xVals[k]))

    # Build slopes between x values
    for k in range(0, numBreakpoints + 1):
        x1 = xVals[k]
        y1 = yVals[k]
        x2 = xVals[k+1]
        y2 = yVals[k+1]
        slope = (y2 - y1)/(x2 - x1)
        intercept = y1 - slope * x1
        piecewise.append((slope, intercept))

    return piecewise


def calculate_median_read_count_density(read_count_densities):

    # Calculate total length
    intervals_length = 0.0
    for density, length in read_count_densities:
        intervals_length += length

    read_count_densities = sorted(read_count_densities)

    # Calculate median read count density across intervals
    median_read_count_density = 0.0
    length_accum = 0
    for density, length in read_count_densities:
        if length_accum + length > intervals_length / 2.0:
            median_read_count_density = density
            break
        length_accum += length

    return median_read_count_density


class EulerianLP(object):

    def __init__(self):
        self.lambda_variant = 0.0
        self.lambda_telomere = 0.0

    def solve(self, bg):

        copies_var = lambda edge: "copies_{0}_{1}".format(edge.edge_type, edge.id)
        likelihood_var = lambda edge: "llh_{0}_{1}".format(edge.edge_type, edge.id)

        cplex_variables = list()
        cplex_constraints = list()

        for edge in bg.edges:

            # Calculate objective coefficient based on regularization
            copies_var_obj = 0.0
            if edge.edge_type == 'variant' and not edge.unregularized:
                copies_var_obj = self.lambda_variant
            elif edge.edge_type == 'telomere' and not edge.unregularized:
                copies_var_obj = self.lambda_telomere

            # Calculate upper bound for enabled/disabled edges
            copies_var_ub = 100.0
            if edge.is_disabled:
                copies_var_ub = 0.0

            # Add regularized tumour coverage variables
            cplex_variables.append((copies_var(edge), copies_var_obj, 0.0, copies_var_ub))

            # No likelihood for zero length edges
            if edge.region_length <= 0.0:
                continue
            
            # Add total coverage likelihood variable and include in objective
            cplex_variables.append((likelihood_var(edge), 1.0, None, None))

            for line_idx, (slope, intercept) in enumerate(create_piecewise_linear(0.0, edge.readcount, edge.region_length)):

                # likelihood >= copies * slope + intercept
                cplex_constraints.append(([likelihood_var(edge), copies_var(edge)], [1.0, -slope], "G", intercept))

        for node_idx, node in enumerate(bg.nodes):

            node_edge_vars = list()
            node_edge_coeff = list()
            for edge_idx in node.get_in_edges():
                node_edge_vars.append(copies_var(bg.edges[edge_idx]))
                node_edge_coeff.append(1.0)
            for edge_idx in node.get_out_edges():
                node_edge_vars.append(copies_var(bg.edges[edge_idx]))
                node_edge_coeff.append(-1.0)

            # copies_in = copies_out
            cplex_constraints.append((node_edge_vars, node_edge_coeff, "E", 0.0))

        self.objective_value, self.solution_values = cplex_solve_min(cplex_variables, cplex_constraints)

        # Calculate expected counts
        self.expected_counts = list()
        for edge in bg.edges:
            copies = self.solution_values[copies_var(edge)]
            self.expected_counts.append(copies * edge.region_length)

        # Calculate tumour edge copies
        self.edge_copies = list()
        for edge in bg.edges:
            copies = self.solution_values[copies_var(edge)]
            self.edge_copies.append(copies)


SolutionInfo = namedtuple('SolutionInfo', ['objective_value',
                                           'edge_copies',
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

    return SolutionInfo(lp.objective_value, 
                        [copies for copies in lp.edge_copies],
                        [edge.edge_type for edge in bg.edges],
                        [edge.id for edge in bg.edges],
                        lambda_variant,
                        lambda_telomere)



