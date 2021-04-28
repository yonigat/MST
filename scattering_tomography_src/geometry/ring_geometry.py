from geometry.abstract_geometry import AbstractGeometry
from detector.ring_detector import RingDetector
from source.ring_source import RingSource
from volume.abstract_volume import AbstractVolume
import numpy as np


class RingGeometry(AbstractGeometry):

    def __init__(self, detectors: 'RingDetector', sources: 'RingSource', volume: AbstractVolume):
        super().__init__()
        self._detectors = detectors
        self._sources = sources
        self._volume = volume

    def get_projection(self, source_idx: int):
        """
        return a vector of the columns and rows the cone of the source reaches
        :param source_idx: index of the source we want to get its projection
        :return: the upper and lower boundaries of the columns of the projection and the angle to the middle of the projection
        (polar angle [-pi, pi], relative ot postive side of X axis
        """
        num_cols = self._detectors.get_cols()
        angle_between_detectors = (2*np.pi) / num_cols  # angle between 2 detectors
        gamma = self._detectors.get_beta()  # angle of the width of the detector
        source_loc = self._sources._locations[source_idx, :]
        beta = self.set_angle_to_range(np.arctan2(source_loc[1], source_loc[0]))   # angle of the source
        angle_of_middle_det = self.set_angle_to_range(beta + np.pi)

        # vol_size = np.array(self._volume.get_grid().get_shape()) * np.array(self._volume.get_grid().get_voxel_size())
        # alpha = 0.9 * np.arctan2(vol_size[1]/2, np.linalg.norm(source_loc))
        """
        We now calculate the degree between the middle detector and the furthest detectors the source reaches.
        This is done according to cosine theroem, and can fit a case where the radius of the source (Rs) is smaller or larger than
        the detectors radius (Rd).
        In case the Rs>Rd, there is a critical angle for the opening of the cone that after that there is no intersection between the 
        edge of the cone and the detector array, so we use the critical angle in that case.
        """
        Rs = np.linalg.norm(self._sources.get_radius())  # radius of the source
        Rd = np.linalg.norm(self._detectors.get_radius())  # radius of the detectors
        alpha = self._sources.get_focal_angle()[0]  # polar angle of the source
        if Rd > Rs or alpha <= np.arcsin(Rd/Rs):  # check if were over the critical angle
            longest_beam_of_source = Rs*np.cos(alpha) + np.sqrt(Rd**2 - (Rs**2)*(np.sin(alpha)**2))
            beta_center = np.pi - np.arccos((longest_beam_of_source**2 - Rs**2 - Rd**2)/(-2*Rs*Rd))
        else:
            longest_beam_of_source = Rs * np.cos(np.arcsin(Rd/Rs)) + np.sqrt(Rd ** 2 - (Rs ** 2) * ((Rd/Rs) ** 2))
            beta_center = np.pi - np.arccos((longest_beam_of_source ** 2 - Rs ** 2 - Rd ** 2) / (-2 * Rs * Rd))

        upper_boundry = self.set_angle_to_range(angle_of_middle_det + beta_center)
        lower_boundry = self.set_angle_to_range(angle_of_middle_det - beta_center)

        # it's ok if the lower_boundry > upper_boundry, it means we need to rearrange the image we will
        # send as the detectors image

        # we now translate the angles to column index in the detector array

        upper_col = 0
        upper_col_angle = 0
        if upper_boundry >= 0:
            while upper_col_angle < upper_boundry:
                upper_col_angle += angle_between_detectors
                upper_col += 1
            if upper_col_angle - gamma > upper_boundry:
                upper_col -= 1
        elif upper_boundry < 0:
            upper_col = num_cols
            while upper_col_angle > upper_boundry:
                upper_col_angle -= angle_between_detectors
                upper_col -= 1
            if upper_col_angle + gamma < upper_boundry:
                upper_col += 1

        lower_col = 0
        lower_col_angle = 0
        if lower_boundry >= 0:
            while lower_col_angle < lower_boundry:
                lower_col_angle += angle_between_detectors
                lower_col += 1
            if lower_col_angle - gamma > lower_boundry:
                lower_col -= 1
        elif lower_boundry < 0:
            lower_col = num_cols
            while lower_col_angle > lower_boundry:
                lower_col_angle -= angle_between_detectors
                lower_col -= 1
            if lower_col_angle + gamma < lower_boundry:
                lower_col += 1

        lower_col = np.mod(lower_col, num_cols)

        cols = np.array([lower_col, upper_col], dtype=int)
        return cols, angle_of_middle_det

    def get_center_and_orientation(self, angle: float):
        """
        :param angle: the polar angle of the middle of the projections detector
        :return: 3 vectors: the center of the detector, unit vector that represents the (1,0) vector of the detector array
        and a unit vector of the (0,1) vector of the detector array (in the case of ring detector, always 'Z' direction
        """
        radius = self._detectors.get_radius()
        rot_mat = np.array([[np.cos(angle),  -np.sin(angle), 0],
                            [np.sin(angle), np.cos(angle), 0],
                            [0,              0,             1]])
        center_orig = np.array([radius, 0, 0])
        u_orig = np.array([0, -1, 0])
        v_orig = np.array([0, 0, 1])
        center = np.dot(rot_mat, center_orig)
        u = self._detectors.get_sizes()[1]*np.dot(rot_mat, u_orig)
        v = self._detectors.get_sizes()[2]*np.dot(rot_mat, v_orig)
        return center, u, v

    @staticmethod
    def set_angle_to_range(angle: float):
        """
        :param angle:
        :return: sets an angle (in radians) to the range of [pi, -pi)
        """
        new_angle = angle
        if new_angle > np.pi:
            while new_angle > np.pi:
                new_angle -= 2*np.pi
        elif new_angle <= -np.pi:
            while new_angle <= -np.pi:
                new_angle += 2*np.pi
        return new_angle
