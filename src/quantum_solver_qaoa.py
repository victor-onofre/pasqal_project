import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
from pulser import Pulse, Sequence
from pulser.devices import Chadoq2
from pulser_simulation import QutipEmulator
from scipy.optimize import minimize


def define_sequence(
    register,
    p_layers,
):
    """Definition of parametrized sequence to use for the QAOA.

    Args:
        register: The atomic register representing the problem in the quantum device
        p_layers: Number of layers for the QAOA

    Returns:
        seq: parametrized sequence
    """

    # Parametrized sequence
    seq = Sequence(register, Chadoq2)
    seq.declare_channel("ch0", "rydberg_global")

    beta_list = seq.declare_variable("beta_list", size=p_layers)
    gamma_list = seq.declare_variable("gamma_list", size=p_layers)

    if p_layers == 1:
        beta_list = [beta_list]
        gamma_list = [gamma_list]

    for beta, gamma in zip(beta_list, gamma_list):
        beta_pulse = Pulse.ConstantPulse(1000 * beta, 1.0, 0, 0)
        gamma_pulse = Pulse.ConstantPulse(1000 * gamma, 0, 1.0, 0)

        seq.add(beta_pulse, "ch0")
        seq.add(gamma_pulse, "ch0")

    seq.measure("ground-rydberg")

    return seq


# Building the quantum loop


def quantum_loop_qaoa(
    parameters,
    sequence,
    p_layers,
):
    """Simulation of the sequence in the quantum device

    Args:
        parameters: Parameters used in the sequence
        sequence: Definition of the sequence given a registers
        p_layers: Number of layers for the QAOA

    Returns:
        counts:Distribution of the samples
    """

    params = np.array(parameters)
    beta_params, gamma_params = np.reshape(params.astype(int), (2, p_layers))

    assigned_seq = sequence.build(beta_list=beta_params, gamma_list=gamma_params)

    sim = QutipEmulator.from_sequence(assigned_seq, sampling_rate=0.01)
    res = sim.run()

    # Samples from the state vector that results from running the simulation
    counts = res.sample_final_state(N_samples=6000)

    return counts


def plot_distribution(count_dist):
    count_dist = dict(sorted(count_dist.items(), key=lambda item: item[1], reverse=True))
    plt.figure(figsize=(6, 4))
    plt.xlabel("bitstings")
    plt.ylabel("counts")
    plt.bar(count_dist.keys(), count_dist.values(), width=0.5)
    plt.xticks(rotation="vertical")
    plt.show()


# Building the quantum loop


def get_cost(bitstring, Q):
    """Get the cost of the problem

    Args:
        Q: QUBO matrix
        bitstring: Solution in bitstring format

    Returns:
        cost
    """
    z = np.array(list(bitstring), dtype=int)
    cost = z.T @ Q @ z
    return cost


def get_cost_average(counts_distr, Q):
    """Get the average cost of the problem

    Args:
        counts_distr: Counts distribution from the sequence simulation
        Q: QUBO matrix

    Returns:
        Average cost
    """
    cost = sum(counts_distr[key] * get_cost(key, Q) for key in counts_distr)
    return cost / sum(counts_distr.values())  # Divide by total samples


def func_to_min_1(param, *args):
    """Function to minimize

    Args:
        param: Parameters beta and gamma for the sequence
        args: QUBO matrix, sequence and p layers

    Returns:
        Average cost
    """
    Q, sequence, p_l = args
    # sequence  = args[0][1]
    # p_l = args[0][2]
    C = quantum_loop_qaoa(param, sequence, p_l)
    cost = get_cost_average(C, Q)
    return cost


def QAOA_solver_for_max_cut(
    qubo_matrix,
    register,
    p_layers,
    optimizer_="COBYLA",
):
    """Function to minimize

    Args:
        qubo_matrix: QUBO matrix representing the problem
        register: The atomic register representing the problem in the quantum device
        p_layers: Number of layers for the QAOA
        optimizer_: Minimizer to use from scipy. Default COBYLA

    Returns:
        optimal_parameters
    """

    scores = []
    params = []
    testing = []
    sequence = define_sequence(register, p_layers)

    for repetition in range(10):
        testing.append(repetition)
        random_beta = np.random.uniform(1, 10, p_layers)
        random_gamma = np.random.uniform(1, 10, p_layers)

        try:
            res = minimize(
                func_to_min_1,
                args=(qubo_matrix, sequence, p_layers),
                x0=np.r_[random_beta, random_gamma],
                method=optimizer_,  # 'Nelder-Mead'
                tol=1e-5,
                options={"maxiter": 1000},
            )
            # print(res.fun)
            scores.append(res.fun)
            params.append(res.x)

            optimal_parameters = params[np.argmin(scores)]
        except:
            pass

    return optimal_parameters


def plot_solution_qaoa(graph,
                       Q_matrix,
                       optimal_parameters,
                       register,
                       p_layers,
                       plot_histogram= False):
    sequence = define_sequence(register, p_layers)
    optimial_count_dict = quantum_loop_qaoa(optimal_parameters, sequence, p_layers)
    best_solution = max(optimial_count_dict, key=optimial_count_dict.get)
    colors = [
        "b" if best_solution[node] == "0" else "r" for node in graph
    ]  # Define the colors of the nodes for the best solution
    best_cut = get_cost(best_solution, Q_matrix)
    print(f"Best solution: {best_solution} with {-best_cut} cuts")

    if plot_histogram:
        plot_distribution(optimial_count_dict)

    plt.figure(figsize=(3, 2))
    nx.draw(graph, with_labels=True, node_color=colors, font_weight="bold")
