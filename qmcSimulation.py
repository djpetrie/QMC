"""
@author David Petrie
"""

import random
import numpy as np
import qmcDisplay


class Electron:
    """
    A class representing an electron bound to a helium atom.
    """
    def __init__(self, r, norm):
        """
        Initializes the state for the electron
        :param r: A 3-vector (numpy array) representing the position of the electron in R3
        :param norm: The norm of the position vector (distance from the origin)
        """
        self.r = r
        self.norm = norm

    def get_disp(self):
        """
        returns the displacement from the origin
        :return: The displacement from the origin
        """
        return disp(self.r)


class Pair:
    """
    A class representing a pair of electrons bound to a helium atom
    """
    def __init__(self, e1, e2):
        """
        Initializes the state of the electron pair
        :param e1: The first electron to add to the pair
        :param e2: The second electron to add to the pair
        """
        self.e1 = e1
        self.e2 = e2

    def energy(self, charge):
        """
        Returns the system energy for the helium atom associated with the current electron pair and configuration.
        :param charge: The effective nuclear charge to use for the energy calculation
        :return: The system energy for the helium atom with the current electron pair and configuration.
        """
        q = charge
        _r1 = self.e1.norm
        _r2 = self.e2.norm
        _r12 = get_dist(self.e1.r, self.e2.r)
        f = -q ** 2 + (q - 2.0) / _r1 + (q - 2.0) / _r2 + 1.0 / _r12
        return f


class Qmc:
    """
    A simulation representing an ensemble helium atoms with bound electron pairs, in which each atom is simulated
    separately. The simulation's state is determined by the position vectors for both electrons in each of the electron
    pairs, and can be updated according to the probability density function given by the associated trial wave function.
    """
    def __init__(self, num_pairs, radius, step_size, charge):
        """
        Initializes the state of the QMC simulation
        :param num_pairs: Number of electron pairs
        :param radius: The maximum radius for initial electron positions with respect to the helium nucleus
        :param step_size: The maximum distance that an electron can be moved during an update
        :param charge: The effective nuclear charge seen by the electrons (as a result of shielding)
        """
        self.num_pairs = num_pairs
        self.pairs = []
        self.radius = radius
        self.step_size = step_size
        self.charge = charge
        self.x_data = np.zeros(2 * num_pairs, dtype=float)
        self.y_data = np.zeros(2 * num_pairs, dtype=float)
        self.z_data = np.zeros(2 * num_pairs, dtype=float)
        self.generate_pairs()
        self.display = qmcDisplay.QmcDisplay(self)

    def generate_pairs(self):
        """
        Populates the simulations list of electron pairs by generating new electrons with random positions.
        :return: None
        """
        for i in range(self.num_pairs):
            r1, r1_norm = self.random_vector()
            r2, r2_norm = self.random_vector()
            e1 = Electron(r1, r1_norm)
            e2 = Electron(r2, r2_norm)
            self.pairs.append(Pair(e1, e2))

    def get_data(self):
        """
        Updates and returns the simulations coordinate data for the all electron pairs.
        :return: A tuple of 3 numpy arrays containing the x, y, and z coordinate data respectively for all electrons in
        the simulation.
        """
        for i in range(self.num_pairs):
            self.x_data[2 * i], self.y_data[2 * i], self.z_data[2 * i] = self.pairs[i].e1.r
            self.x_data[2 * i + 1], self.y_data[2 * i + 1], self.z_data[2 * i + 1] = self.pairs[i].e2.r
        return self.x_data, self.y_data, self.z_data

    def reset(self):
        """
        Resets the simulation by assigning new random position vectors and norms to the electrons.
        :return:
        """
        for pair in self.pairs:
            pair.e1.r, pair.e1.norm = self.random_vector()
            pair.e1.r, pair.e2.norm = self.random_vector()

    def set_charge(self, charge):
        """
        Sets the effective nuclear charge seen by an electron (as a result of shielding).
        :param charge:
        :return:
        """
        self.charge = charge

    def random_vector(self):
        """
        Generates and returns a tuple containing a random position 3-vector and corresponding L2-norm with a maximum
        size given by the radius of the simulation.
        :return: A tuple containing a random position 3-vector and corresponding L2-norm within the simulation space.
        """
        radius = self.radius
        flag = True
        r = None
        norm = None
        while flag:
            r = np.random.uniform(-radius, radius, 3)
            norm = np.sqrt(r[0] ** 2 + r[1] ** 2 + r[2] ** 2)
            if norm <= radius:
                flag = False
        return r, norm

    def avg_energy(self):
        """
        Calculates and returns the average system energy across all helium atoms given their current configuration.
        :return: The average system energy across all helium atoms
        """
        avg_e = 0.0
        for pair in self.pairs:
            avg_e += pair.energy(self.charge)
        avg_e = avg_e / self.num_pairs
        return avg_e

    def update(self):
        """
        Updates the position of the electrons according to a Metropolis-Hastings algorithm. A new trial configuration is
        randomly generated for each pair, and tested against an acceptance criteria to determine if the move is allowed.
        The acceptance criteria biases the movements towards regions where the probability density from the given wave
        function is higher.
        :return: None
        """
        q = self.charge
        step = self.step_size
        for j in range(0, self.num_pairs):
            pair = self.pairs[j]
            r1 = pair.e1.r
            r2 = pair.e2.r
            t1 = r1 + np.random.uniform(-step, step, 3)
            t1_norm = np.sqrt(t1[0] ** 2 + t1[1] ** 2 + t1[2] ** 2)
            t2 = r2 + np.random.uniform(-step, step, 3)
            t2_norm = np.sqrt(t2[0] ** 2 + t2[1] ** 2 + t2[2] ** 2)
            if move_accept(pair.e1.norm, pair.e2.norm, t1_norm, t2_norm, q) > random.random():
                pair.e1.r = t1
                pair.e1.norm = t1_norm
                pair.e2.r = t2
                pair.e2.norm = t2_norm


def disp(r):
    """
    The displacement from the origin as given by taking the magnitude of the given position vector
    :param r: The position vector to find the magnitude of
    :return: The magnitude of the given position vector
    """
    return np.sqrt(r[0] ** 2 + r[1] ** 2 + r[2] ** 2)


def get_dist(r1, r2):
    """
    Returns the distance between a set of two given 3-vectors.
    :param r1: The first position vector of the set
    :param r2: The second position vector of the set
    :return: The magnitude of the difference between the two 3-vectors
    """
    return np.sqrt((r2[0] - r1[0]) ** 2 + (r2[1] - r1[1]) ** 2 + (r2[2] - r1[2]) ** 2)


def move_accept(r1, r2, t1, t2, charge):
    """
    Returns a value representing the ratio of the probability of the trial configuration over the initial configuration.
    This acts as an acceptance ratio for the Metropolis Hastings algorithm.
    The proper expression of this term is: |Psi'|^2 / |Psi|^2; where
        |Psi'|^2 is the squared modulus of the wave function evaluated at the trial configuration
        |Psi|^2 is the squared modulus of the wave function evaluated at the initial configuration
        Psi = e^(-q * (r1 + r2)); A trial wave function constructed as the product between two ground level hydrogen
            states, with a variational parameter q representing effective nuclear charge.


    :param r1: The first electron's distance from the origin
    :param r2: The second electron's distance from the origin
    :param t1: The distance from the origin for the trial position of the first electron
    :param t2: The distance from the origin for the trial position of the second electron
    :param charge: The effective nuclear charge felt by the electrons as a result of shielding
    :return: a value representing the ratio of the probability of the trial configuration over the initial configuration
    """
    # A simplified form of the expression for faster calculation
    return np.exp(2 * charge * (r1 + r2 - t1 - t2))
