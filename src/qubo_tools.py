import numpy as np


def build_qubo_matrix(graph):
    # Number of nodes in the graph
    num_nodes = len(graph.nodes)

    # Initialize our QUBO matrix
    Q_matrix = np.zeros((num_nodes, num_nodes))

    # Update Q matrix for every edge in the graph
    for i, j in graph.edges:
        Q_matrix[i][i] += -1
        Q_matrix[j][j] += -1
        Q_matrix[i][j] += 2

    return Q_matrix


def brute_force(qubo_matrix):
    """Brute force solver for the QUBO problem

    Args:
        qubo_matrix: QUBO matrix representing the problem

    Returns:
        sort_solutions: Dictionary with the solutions in order
    """

    bitstrings = [np.binary_repr(i, len(qubo_matrix)) for i in range(2 ** len(qubo_matrix))]
    costs = []

    # this takes exponential time with the dimension of the QUBO

    for b in bitstrings:
        z = np.array(list(b), dtype=int)
        cost = z.T @ qubo_matrix @ z
        costs.append(cost)
        zipped = zip(bitstrings, costs)

    sort_solutions = sorted(zipped, key=lambda x: x[1])

    return sort_solutions
