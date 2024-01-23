import numpy as np
from pulser.devices import Chadoq2
from scipy.optimize import minimize
from scipy.spatial.distance import pdist, squareform


def evaluate_mapping(new_coords, *args):
    """Cost function to minimize.Function from the Pulser tutorial.

    Args:
        new_coords: the new coordinates for the encoding
        args: QUBO matrix and shape

    Returns:
        cost distance
    """
    Q, shape = args
    new_coords = np.reshape(new_coords, shape)
    new_Q = squareform(
        Chadoq2.interaction_coeff / pdist(new_coords) ** 6,
    )
    return np.linalg.norm(new_Q - Q)


def create_coordinates(
    qubo_matrix,
    optimizer="COBYLA",
    scale_factor=3,
):
    """Creation of coordinates for the register in the quantum device
       given a QUBO matrix

    Args:
        qubo_matrix: QUBO matrix representing the problem
        optimizer: Minimizer to use from scipy. Default COBYLA
        scale_factor: Factor to scale the QUBO matrix. Default 3

    Returns:
        coordinates: The coordinates for the register
    """

    size_matrix = len(qubo_matrix)

    shape = (size_matrix, 2)
    # costs = []
    np.random.seed(0)

    x0 = np.random.random(shape).flatten()

    res = minimize(
        evaluate_mapping,
        x0,
        args=(qubo_matrix * scale_factor, shape),
        method=optimizer,
        tol=1e-6,
        options={"maxiter": 200000},
    )

    coordinates = np.reshape(res.x, (size_matrix, 2))

    return coordinates
