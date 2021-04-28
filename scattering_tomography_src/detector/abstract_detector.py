from abc import ABC, abstractmethod
from mpl_toolkits.mplot3d import Axes3D
from utils.draw_cuboid import draw_cuboid
import matplotlib.pyplot as plt
import numpy as np
import utils.delete_files
import os


class AbstractDetector(ABC):
    """
    Interface of the abstract detector
    """
    @abstractmethod
    def __init__(self):
        """
        Initializes the detector
        """
        # int: overall number of detectors
        self._n_detector_units = None
        # ndarray: _n_detector_units x 3 (dims), the 3D coordinates of center of each detector
        self._locations = None
        # ndarray: _n_detector_units x 3 x 3, 3 orthogonal unit vectors which describe the spatial orientation of each detector
        # the first column [ind, :, 1] describes the orientation of the X direction face, the second column is the orientation
        # of the Y direction face, the third is the orientation in the Z direction face
        self._orientations = None
        # tuple: 3 x 1, the size of each detector unix in the (x,y,z) directions
        self._sizes = None

    def vertices(self, detector_ind: int):
        """"
        returns the  ndarray of size(3,8)
        of the vertices of a chest shaped detector given its sizes, center of mass coordinates and the orientation
        (i.e. the unit vector orthogonal to it's surface)
        *** IMPORTANT***
        We assume that we have 3 orthogonal orientations of every detector which means self._orientations is a self._n_detector_unitsXdimXdim ndarray
        and that each column is an orientation
        """
        det_location = self._locations[detector_ind, :]
        det_orientation_x = self._orientations[detector_ind, :, 0]
        det_orientation_y = self._orientations[detector_ind, :, 1]
        det_orientation_z = self._orientations[detector_ind, :, 2]
        Sx, Sy, Sz = np.array(self._sizes)/2

        vertices = np.zeros((3, 8))  # 8 corners in 3D
        vertices[:, 0] = det_location + Sx * det_orientation_x + Sy * det_orientation_y + Sz * det_orientation_z
        vertices[:, 1] = det_location + Sx * det_orientation_x + Sy * det_orientation_y - Sz * det_orientation_z
        vertices[:, 2] = det_location + Sx * det_orientation_x - Sy * det_orientation_y + Sz * det_orientation_z
        vertices[:, 3] = det_location + Sx * det_orientation_x - Sy * det_orientation_y - Sz * det_orientation_z
        vertices[:, 4] = det_location - Sx * det_orientation_x + Sy * det_orientation_y + Sz * det_orientation_z
        vertices[:, 5] = det_location - Sx * det_orientation_x + Sy * det_orientation_y - Sz * det_orientation_z
        vertices[:, 6] = det_location - Sx * det_orientation_x - Sy * det_orientation_y + Sz * det_orientation_z
        vertices[:, 7] = det_location - Sx * det_orientation_x - Sy * det_orientation_y - Sz * det_orientation_z

        return vertices

    def draw_3d(self, ax: Axes3D = None):
        """
        :param ax: handle of 3D axes if we wish to continue to use some figure
        :return: return handle of 3D axes in case we want to continue using the figure
        """

        if ax is None:
            fig = plt.figure()
            ax = Axes3D(fig)

        for det_ind in np.arange(self._n_detector_units):
            det_vert = self.vertices(det_ind)
            draw_cuboid(ax, det_vert)
        min_locations = np.min(self._locations, 0)
        max_locations = np.max(self._locations, 0)
        Sx, Sy, Sz = np.array(self._sizes) / 2
        ax.set_xlim3d(1.2*(min_locations[0]-Sx), 1.2*(max_locations[0]+Sx))
        ax.set_ylim3d(1.2*(min_locations[1]-Sy), 1.2*(max_locations[1]+Sy))
        ax.set_zlim3d(1.2*(min_locations[2]-Sz), 1.2*(max_locations[2]+Sz))
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')

        return ax

    ############################
    # Writers
    ############################
    def write_location_files(self, path, loc_file_name):
        file = os.path.join(path, loc_file_name )
        np.savetxt(file, self._locations, fmt='%.5e')
        return

    def write_orientation_file(self, path, orient_file_name):
        file = os.path.join(path, orient_file_name)
        orientations_2d_mat = self.get_orientations().reshape((self.get_orientations().shape[0], 9))
        # turn the orientation 3d mat into 2d mat nu_det x 9
        np.savetxt(file, orientations_2d_mat, fmt='%.5e')
        return

    def write_files(self, path="../run_inputs/",
                    orient_file_name='det_orient.txt',
                    loc_file_name='det_loc.txt'):
        self.write_orientation_file(path, orient_file_name)
        self.write_location_files(path, loc_file_name)

    ############################
    # Getters
    ############################
    def get_sizes(self):
        return self._sizes

    def get_locations(self) -> np.ndarray:
        """
        API method
        :return: ndarray. _n_detector_units x 3 (dims), the 3D coordinates of center of each detector
        """
        return self._locations

    def get_orientations(self) -> np.ndarray:
        """
        API method
        :return: ndarray. _n_detector_units x 3 x 3, 3 orthogonal unit vectors which describe the spatial orientation of each detector
        # the first column [ind, :, 1] describes the orientation of the X direction face, the second column is the orientation
        # of the Y direction face, the third is the orientation in the Z direction face

        """
        return self._orientations