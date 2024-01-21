import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
from pulser import Pulse, Sequence
from pulser.devices import Chadoq2
from pulser.waveforms import InterpolatedWaveform
from pulser_simulation import QutipEmulator
from scipy.optimize import minimize


def simple_adiabatic_sequence(
    device,
    register,
    time,
    Omega=3.271543,
    detuning=5,
) -> Sequence:
    """Creates the adiabatic sequence

    Args:
        device: physical device simulation
        Omega: Frecuency
        register: arrangement of atoms in the quantum processor
        time: time of the adiabatic process
        detuning: detuning use

    Returns:
        sequence
    """
    delta_0 = -detuning
    delta_f = -delta_0

    adiabatic_pulse = Pulse(
        InterpolatedWaveform(time, [1e-9, Omega, 1e-9]),
        InterpolatedWaveform(time, [delta_0, 0, delta_f]),
        0,
    )

    sequence = Sequence(register, device)
    sequence.declare_channel("ising", "rydberg_global")

    sequence.add(adiabatic_pulse, "ising")

    return sequence


# Building the quantum loop


def simple_quantum_loop(parameters, register):
    """Simulation of the sequence in the quantum device

    Args:
        parameters: Parameters used in the sequence
        register: arrangement of atoms in the quantum processor

    Returns:
        counts:Distribution of the samples
    """

    params = np.array(parameters)

    parameter_time, parameter_omega, parameter_detuning = np.reshape(params.astype(int), 3)
    seq = simple_adiabatic_sequence(
        Chadoq2,
        register,
        parameter_time,
        Omega=parameter_omega,
        detuning=parameter_detuning,
    )

    simul = QutipEmulator.from_sequence(seq, sampling_rate=0.1)
    res = simul.run()
    counts = res.sample_final_state(N_samples=5000)  # Sample from the state vector
    # print(counts)

    return counts


def plot_distribution(count_dist):
    count_dist = dict(sorted(count_dist.items(), key=lambda item: item[1], reverse=True))
    plt.figure(figsize=(6, 4))
    plt.xlabel("bitstings")
    plt.ylabel("counts")
    plt.bar(count_dist.keys(), count_dist.values(), width=0.5)
    plt.xticks(rotation="vertical")
    plt.show()


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


def func_to_min(param, *args):
    """Function to minimize

    Args:
        param: Parameters beta and gamma for the sequence
        args: QUBO matrix and register

    Returns:
        Average cost
    """
    Q, register = args
    quantum_counts = simple_quantum_loop(param, register)
    cost = get_cost_average(quantum_counts, Q)
    return cost


def VQAA(
    atomic_register,
    qubo_matrix,
    omega_range=(1, 5),
    detuning_range=(1, 5),
    time_range=(8, 25),
    minimizer_method="Nelder-Mead",
    repetitions=10,
) -> list:
    """Main function for the VQAA algorithm.

    Args:
        atomic_register: The atomic register representing the problem in the quantum device
        qubo_matrix: QUBO matrix representing the problem
        function_to_min: Function representing the problem to minimize
        omega_range: The range of frequencies to used for the optimizer parameters. Default (1,5)
        detuning_range: The range of detuning to used for the optimizer parameters. Default (1,5)
        time_range:Range of time evolution for QAA to used in optimizer parameters.Default (8,25)
        minimizer_method: Minimizer to use from scipy. Default Nelder-Mead
        repetitions: The number of times to repeat the optimization. Default(10)

    Returns:
        list: List of all final parameters.
    """
    scores = []
    params = []
    testing = []
    for repetition in range(repetitions):
        testing.append(repetition)
        random_omega = np.random.uniform(omega_range[0], omega_range[1])
        random_detuning = np.random.uniform(detuning_range[0], detuning_range[1])
        random_time = 1000 * np.random.randint(time_range[0], time_range[1])

        res = minimize(
            func_to_min,
            args=(qubo_matrix, atomic_register),
            x0=np.r_[random_time, random_omega, random_detuning],
            method=minimizer_method,
            tol=1e-5,
            options={"maxiter": 20},
        )

        # print(res.fun)
        scores.append(res.fun)
        params.append(res.x)

    optimal_parameters = params[np.argmin(scores)]

    return optimal_parameters


def plot_solution_vqaa(graph, Q_matrix, optimal_parameters, register):
    # sequence = define_sequence(register,p_layers )
    optimial_count_dict = simple_quantum_loop(optimal_parameters, register)
    best_solution = max(optimial_count_dict, key=optimial_count_dict.get)
    colors = [
        "b" if best_solution[node] == "0" else "r" for node in graph
    ]  # Define the colors of the nodes for the best solution
    best_cut = get_cost(best_solution, Q_matrix)
    print(f"Best solution: {best_solution} with {-best_cut} cuts")

    # plot_distribution(optimial_count_dict)
    plt.figure(figsize=(3, 2))
    nx.draw(graph, with_labels=True, node_color=colors, font_weight="bold")
