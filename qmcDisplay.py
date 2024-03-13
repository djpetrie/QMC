import matplotlib.pyplot as plt
import numpy as np
from matplotlib.pyplot import pause


class QmcDisplay:
    def __init__(self, sim):
        self.fig = None
        self.ax_scatter3d = None
        self.image_scatter3d = None
        self.ax_calibration = None
        self.image_calibration = None
        self.ax_charge = None
        self.image_charge = None
        self.image = None
        self.setup(sim)

    def setup(self, sim):
        """
        Creates a new figure, configured for displaying three subplots to visualize different types of data gathered
        during simulation.
        :param sim: The simulation to pull data points from
        :return: The new figure
        """
        self.fig = plt.figure(figsize=(12, 8))
        self.fig.subplots_adjust(hspace=.35)
        self.ax_calibration = plt.subplot2grid((2, 3), (0, 0), colspan=1)
        self.ax_scatter3d = plt.subplot2grid((2, 3), (0, 1), colspan=2, rowspan=2, projection='3d')
        self.setup_scatter3d(sim)
        self.ax_charge = plt.subplot2grid((2, 3), (1, 0), colspan=1)
        self.ax_charge.set_xlim(0.0, 2.0)

        self.ax_calibration.set_xlabel('# of steps taken')
        self.ax_calibration.set_ylabel('Energy value (Hartree)')
        self.ax_calibration.set_title('Early energy change as the system settles')
        return self.fig

    def plot_calibration(self, data):
        """
        Plots the given data to the calibration subplot.
        :param data: The data set to plot (energy levels during the calibration period)
        :return: None
        """
        self.ax_calibration.cla()
        self.ax_calibration.set_xlabel('# of steps taken')
        self.ax_calibration.set_ylabel('Energy value (Hartree)')
        self.ax_calibration.set_title('Early energy change as the system settles')
        self.image_calibration = self.ax_calibration.plot(data)

    def plot_charge(self, q_space, energy_values, i, delta_q):
        """
        Plots the data for average energy level vs effective charge to the charge subplot.
        :param q_space: An array of evenly spaced points for the charge axis (horizontal axis, independent variable)
        :param energy_values: An array of data points taken as the average system energy for a given charge value
        :param i: The index marking the end of the data to plot
        :param delta_q: The step size between charge data points
        :return: None
        """
        self.ax_charge.cla()
        self.ax_charge.set_xlim(q_space[0], q_space[0] + delta_q * len(q_space))
        self.ax_charge.set_xlabel('Effective Nuclear Charge')
        self.ax_charge.set_xlabel('Energy minimum (Hartree)')
        self.ax_charge.set_title('Energy Configuration Curve')
        self.image_charge = self.ax_charge.scatter(q_space[0:i + 1], energy_values[0:i + 1], marker='.')

    def plot_charge_fit(self, q_space, curve):
        self.ax_charge.plot(q_space, curve(q_space))

    def setup_scatter3d(self, sim):
        """
        Sets the configurations needed for the 3D scatter plot of electron positions.
        :param sim: The simulation from which the pull data
        :return: None
        """
        radius = sim.radius
        self.ax_scatter3d.set_xlim(-radius, radius)
        self.ax_scatter3d.set_ylim(-radius, radius)
        self.ax_scatter3d.set_zlim(-radius, radius)
        self.ax_scatter3d.set_aspect('equal')
        self.ax_scatter3d.set_title('Electron positions')

    def plot_system(self, sim):
        """
        Plots the updated 3D scatter plot of the electron positions.
        :param sim: The simulation from which to pull the position data
        :return: None
        """
        if not (self.image_scatter3d is None):
            self.image_scatter3d.remove()
        x_data, y_data, z_data = sim.get_data()
        self.image_scatter3d = self.ax_scatter3d.scatter(x_data, y_data, z_data, marker='.', color='b')

    def plot(self):
        """
        Draws the current figure, updating the visuals on each subplot to the most recently constructed version.
        :return:
        """
        self.fig.canvas.draw_idle()
        pause(0.01)
