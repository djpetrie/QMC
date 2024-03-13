"""
This program uses a Variational Monte Carlo method to find the ground state energy of a helium atom as a quantum system.
Three subplots are displayed and updated at the end of each data collection period:
    Top Left: System energy level during the calibration period for the current run
    Bottom Left: Average system energy vs variational charge parameter values checked thus far (accumulated across runs)
    Right: Electron sample point distribution

Note: For ease of computation all values are considered to be in the system of atomic units.

While in principle one can determine the wave function and energy values in the course of solving Schrödinger equation
as the associated eigenstates and eigenvalues, for the many-body problem there is no closed form analytic solution.
Instead, a trial wave function is constructed as the product of two ground level hydrogen atom states, with nuclear
charge of 2. The Hamiltonian used accounts for the kinetic energy of the electrons, neglecting nuclear motion, and the
potential due to the electric field. The Hamiltonian and wave function are modified to admit a variational parameter q,
in place of the nuclear charge, to simulate charge shielding of the nucleus by the other electron.

The Monte Carlo method handles the evaluation of the multidimensional integrals required for the calculation of the
energy expectation value by sampling the configuration space. The sampling is done by constructing a Markov chain using
the Metropolis Hastings algorithm. The simulation is run multiple times and with differing values of the variational
parameter for effective nuclear charge, resulting in a curve of energy values. The minimum along this curve gives us the
energy value of the ground state at the corresponding value for effective nuclear charge.

"""

import qmcSimulation as qmcSim
import matplotlib.pyplot as plt
import numpy as np
import time
import sys

# Global Constants
DEFAULT_NUM_PAIRS = 75     # Default number of electron pairs
DEFAULT_RADIUS = 4          # Default radius for QmcSimulation, limits starting positions
DEFAULT_DELTA_R = 0.5       # Default upper limit on step size for electron position updates
DEFAULT_UPDATES = 2500      # Default number of updates to perform when collecting energy data for a given charge value
DEFAULT_CHARGE = 1.25       # Default starting value for effective charge parameter
DEFAULT_DELTA_Q = 0.05      # Default step size for effective charge parameter
DEFAULT_Q_STEPS = 15        # Default number of steps for effective charge parameter (number of charge data points)
DEFAULT_PLOTTING = True     # Default value for boolean flag that determines plotting behavior


def main():
    """
    Runs the qmcSimulation multiple times with varying levels of effective nuclear charge, as determined by the q_steps
    and delta_q settings, in order to collect data on the energy expectation values. Data from each run is plotted as
    the energy curve is constructed. Once all data has been collected across the space of effective charge values the
    resulting energy minimum is printed.
    :return: None
    """
    start = time.time()
    num_pairs = DEFAULT_NUM_PAIRS
    radius = DEFAULT_RADIUS
    delta_r = DEFAULT_DELTA_R
    updates = DEFAULT_UPDATES
    charge = DEFAULT_CHARGE
    delta_q = DEFAULT_DELTA_Q
    q_steps = DEFAULT_Q_STEPS
    plotting = DEFAULT_PLOTTING

    energy_data = np.zeros(q_steps, dtype=float)
    q_space = np.arange(charge, charge + q_steps * delta_q, delta_q)

    sim = qmcSim.Qmc(num_pairs, radius, delta_r, charge)

    for i in range(0, q_steps):
        if i != 0:
            sim.reset()
            plotting = False
        sim.set_charge(charge + i * delta_q)
        calibrate_system(sim, plotting)
        avg_sys_energy = 0.0
        for j in range(0, updates):
            sim.update()
            avg_sys_energy += sim.avg_energy()
        energy_data[i] = avg_sys_energy / updates
        sim.display.plot_system(sim)
        sim.display.plot_charge(q_space, energy_data, i, delta_q)
        sim.display.plot()
        print_progress(i, q_steps)

    poly = np.polyfit(q_space, energy_data, deg=2)
    curve_fit = np.poly1d(poly)
    sim.display.plot_charge_fit(q_space, curve_fit)
    sim.display.plot()

    end = time.time()
    t_total = end-start
    print_results(t_total, poly)
    plt.ioff()
    plt.show()


def calibrate_system(sim, plotting=True):
    """
    This method continuously updates a Qmc simulation until it approaches an energy equilibrium, where the distribution
    of electron sample points closely matches the probability distribution function given by the trial wave function.
    :param sim: The Qmc object (simulation) to calibrate
    :param plotting: A boolean flag to set plotting behavior. 'True' results in regular plot updates during calibration.
    :return: None
    """
    window_size = 25
    test_set = []
    buffer = np.zeros(window_size)
    equilibrium = False
    previous_avg = 100.0
    plt.ion()
    while not equilibrium:
        for i in range(0, window_size):
            current = sim.avg_energy()
            buffer[i] = current
            test_set.append(current)
            sim.update()
            if plotting and i % 10 == 0:
                sim.display.plot_system(sim)
                sim.display.plot()
        current_avg = np.average(buffer)
        if current_avg < previous_avg:
            previous_avg = current_avg
        else:
            equilibrium = True
        if plotting:
            sim.display.plot_calibration(test_set)
    sim.display.plot_calibration(test_set)


def print_progress(i, q_steps):
    """
    Prints a percentage value representing the current progress in calculating the grounds state energy value.
    :param i: The index of the current step
    :param q_steps: The total number of steps
    :return:
    """
    text = round(((i + 1) / q_steps) * 100, 1)
    if i == q_steps - 1:
        sys.stdout.write("\rRunning simulation: %d%%\n" % text)
    else:
        sys.stdout.write("\rRunning simulation: %d%%" % text)
    sys.stdout.flush()


def print_results(t_total, poly):
    """
    Prints the runtime of the simulation as well as the result for the ground state energy level.
    :param t_total: The total runtime of the simulation
    :param poly: A 3-vector containing the constant coefficients of the polynomial fit for the energy data
    :return: None
    """
    a, b, c = poly
    x_min = -b / (2 * a)
    y_min = a * x_min ** 2 + b * x_min + c
    print(x_min)
    print(y_min)
    time_str = str(t_total).rjust(6)
    min_energy_str = str(round(y_min, 5)).ljust(8)
    min_energy_str2 = str(round(y_min * 27.2114, 5)).ljust(8)
    charge_str = str(round(x_min, 3)).ljust(8)
    print("time: ", time_str, " s")
    print("Minimum Energy Value: ", min_energy_str, " Hartree")
    print("                      ", min_energy_str2, " EV")
    print("Effective charge:     ", charge_str, " e")


if __name__ == '__main__':
    main()
