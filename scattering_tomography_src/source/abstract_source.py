from abc import ABC, abstractmethod
import numpy as np
import yaml
from mpl_toolkits.mplot3d import Axes3D
import mpl_toolkits.mplot3d as a3
import matplotlib.pyplot as plt
import utils.delete_files
import os


class AbstractSource(ABC):
    """
    Interface of the abstract source
    """
    @abstractmethod
    def __init__(self):
        """
        Initializes the source
        """
        # int: number of overall sources
        self._n_source_units = None
        # ndarray: _n_source_units x 3 (dims), the 3D coordinates of center of each source
        self._locations = None
        # ndarray: _n_source_units x 3 x 3, 3 orthogonal unit vectors which describe the spatial orientation of each source
        # the first column [ind, :, 0] describes the orientation in which the source opens with angle of theta, the second column is the orientation
        # in which the source opens with angle of delta, the third is the orientation in which the the source opens (i.e. the central beam of the cone)
        self._orientations = None
        # tuple: 2 x 1, the opening angles of the source (theta, delta), theta is the azimuthal (polar) angle, measured relatively to positive x axis, values are [0, pi/2]
        # delta is the latitude angle, measured relatively to a vector on the XY plane, values ar [0, pi/2]
        # it is related to phi (the zenith angle in spherical coordinates) by phi = pi/2 - delta
        self._focal_angles = None
        self._spectrum = None
        self._edge_energies = None

    def draw_3d(self, ax: Axes3D = None):
        """
        :param ax: handle of 3D axes if we wish to continue to use some figure
        :return: return handle of 3D axes in case we want to continue using the figure
        """

        if ax is None:
            fig = plt.figure()
            ax = Axes3D(fig)

        for source_ind in np.arange(self._n_source_units):
            self.draw_source(source_ind, ax)
        radius = np.linalg.norm(self._locations[0, :])
        ax.set_xlim3d(2*(-radius), 2*(radius))
        ax.set_ylim3d(2*(-radius), 2*(radius))
        ax.set_zlim3d(2*(-radius), 2*(radius))
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')

        return ax

    def draw_source(self, source_ind: int, ax: Axes3D):
        """
        calculate
        :param ax: handle of the 3D figure which we're plotting on
        :param source_ind:
        :return:
        """
        # create an elliptical cone concentric with Z axis, X axis is related to theta angle, Y axis is related to delta angle
        theta, delta = self._focal_angles
        radius = np.linalg.norm(self._locations[source_ind, :])
        source_loc = self._locations[source_ind, :]
        n_grid_points = 100
        a = 2 * radius * np.tan(theta)
        b = 2 * radius * np.tan(delta)
        x = np.linspace(-a, a, n_grid_points)
        y = np.linspace(-b, b, n_grid_points)
        xx, yy = np.meshgrid(x, y)
        zz = np.sqrt((2 * radius) ** 2 * ((np.square(xx) / np.square(a)) + (np.square(yy) / np.square(b))))
        zz[zz > 2*radius] = np.nan

        R = self._orientations[source_ind]  # the orientation matrix of each source (each vector corresponds to a column)
                                            # is also the rotation matrix from the canonical coordinate system to the
                                            # coordinate system that corresponds to the orientations of the cone

        # iterating through every point of the cone grid and rotate + translate it to the location of the source
        for ind_i in np.arange(n_grid_points):
            for ind_j in np.arange(n_grid_points):
                curr_vec = np.reshape((xx[ind_i, ind_j], yy[ind_i, ind_j], zz[ind_i, ind_j]), (3, 1))
                xx[ind_i, ind_j], yy[ind_i, ind_j], zz[ind_i, ind_j] = np.dot(R, curr_vec) + np.reshape(source_loc, (3, 1))

        ax.plot_wireframe(xx, yy, zz)
        return ax

    def read_spectrum(self, path="../../run_inputs/spectrum_const.txt"):
        self._spectrum = np.loadtxt(path)
        self._edge_energies = (0, len(self._spectrum))

    def create_120keV_spect(self, path="../"):
        file = path + '120keV_spectrum.yaml'
        with open(file, 'r') as f:
            cfg = yaml.load(f)
        self._spectrum = cfg['SPECTRUM']
        self._edge_energies = (0, len(self._spectrum))

    def create_120keV_spect_trimed(self, path="../"):
        file = path + '120keV_spectrum.yaml'
        with open(file, 'r') as f:
            cfg = yaml.load(f)
        spectrum = np.array(cfg['SPECTRUM'])
        spectrum[0:40] = 0
        spectrum = spectrum / np.sum(spectrum)
        self._spectrum = list(spectrum)
        self._edge_energies = (0, len(self._spectrum))

    def create_60keV_delta(self):
        self._spectrum = np.zeros((150,))
        self._spectrum[59] = 1.0
        self._edge_energies = (0, len(self._spectrum))

    def create_multi_delta(self, bins: list):
        self._spectrum = np.zeros((150,))
        self._spectrum[bins] = 1.0 / len(bins)
        self._edge_energies = (0, len(self._spectrum))

    ###########################
    # Writers
    ###########################

    def write_location_files(self, path, loc_file_name):
        file = os.path.join(path, loc_file_name)
        np.savetxt(file, self._locations, fmt='%.5e')
        return

    def write_orientation_file(self, path, orient_file_name):
        file = os.path.join(path, orient_file_name)
        orientations_2d_mat = self.get_orientations().reshape((self.get_orientations().shape[0], 9))
        # turn the orientation 3d mat into 2d mat nu_src x 9
        np.savetxt(file, orientations_2d_mat, fmt='%.5e')
        return

    def write_spectrum_file(self, path, spectrum_file_name):
        file = os.path.join(path, spectrum_file_name)
        np.savetxt(file, self._spectrum, fmt='%.5e')
        return

    def write_files(self, path="../run_inputs/",
                    loc_file="src_loc.txt",
                    orient_file="src_orient.txt",
                    spectrum_file="spectrum.txt"):

        self.write_location_files(path, loc_file)
        self.write_orientation_file(path, orient_file)
        self.write_spectrum_file(path, spectrum_file)

    ###########################
    # Getters
    ###########################

    def get_locations(self) -> np.ndarray:
        """
        API method
        :return: ndarray: _n_source_units x 3 (dims), the 3D coordinates of center of each source
        """
        return self._locations

    def get_orientations(self) -> np.ndarray:
        """
        API method
        :return: ndarray: _n_source_units x 3 x 3, 3 orthogonal unit vectors which describe the spatial orientation of each source
        the first column [ind, :, 1] describes the orientation in which the source opens with angle of theta, the second column is the orientation
        in which the source opens with angle of delta, the third is the orientation in which the the source opens (i.e. the central beam of the cone)

        """
        return self._orientations

    def get_spectrum(self) -> np.ndarray:
        """
        API method
        :return: TODO: document the exact structure
        """
        return self._spectrum

