from detector.abstract_detector import AbstractDetector
import numpy as np



class RingDetector(AbstractDetector):
    """
    Creates a ring detector
    """

    def __init__(self, radius: float, n_units_along_ring: int, n_rows: int, unit_size: tuple, plane_normal: str = 'Z'):
        """
        Creates the ring detector from the given specs
        :param radius: the ring radius
        :param plane_normal: the normal to the plane on which the ring lies, can be 'X', 'Y', 'Z'
        :param n_units_along_ring: number of detector elements along the ring
        :param n_rows: number of rows (concatenation of rings along the orthogonal axis to the radius)
        :param unit_size: tuple of size (x_size: float, y_size: float, z_size: float) of each detector unit
        """
        if plane_normal != 'Z':
            raise Exception("Currently plane normal supports only Z direction")

        super().__init__()
        self._n_detector_units = n_units_along_ring * n_rows
        self._sizes = unit_size
        self._n_rows = n_rows
        self._radius = radius
        self._beta = np.arctan(self._sizes[1]/(2*self._radius))  # angle between the line connecting the origin to the center
                                                             # of the detector to the line connecting the origin with edge
                                                             # of the detector, in the XY plane


        dims = 3
        alpha = 2*np.pi / n_units_along_ring  # angle between each two detectors
        trans_dict = {'X': np.array([0, radius, 0]), 'Y': np.array([radius, 0, 0]), 'Z': np.array([radius, 0, 0])}
        orientations = np.zeros((self._n_detector_units, dims, dims))  # each detector has 3 orthogonal orientations
        locations = np.zeros((self._n_detector_units, dims))  # each detectors location is determined by its center

        """
        creating positions and orientations of detectors
        we assume plane normal == 'Z'
        In future implementations if a different plane normal is needed - a simple linear transformation around the axis can be applied 
        """
        orig_orientation = np.eye(3)  # each column is an orientation of the chest in original position
        H = n_rows*unit_size[2]  # total length of each column along Z axis
        for col_ind in np.arange(n_units_along_ring):
            rot_mat = np.array([[np.cos(col_ind * alpha), -np.sin(col_ind * alpha), 0],
                               [np.sin(col_ind * alpha), np.cos(col_ind * alpha), 0],
                               [0,                         0,                      1]])
            for row_ind in np.arange(n_rows):
                locations[col_ind*n_rows + row_ind, :] = np.dot(rot_mat, trans_dict[plane_normal].T).T  # transform to correct location in XY plane
                orientations[col_ind*n_rows + row_ind, :, :] = np.dot(rot_mat, orig_orientation)  # translate the orientations (each column is an orientation now)
                locations[col_ind * n_rows + row_ind, :] += np.array([0, 0, H/2 - (unit_size[2]/2) * (2*row_ind + 1)])  # translate to correct location along Z axis

        # if needed we can now perform a transformation on the locations and orientations to fit a different normal
        # than Z
        orientations[np.abs(orientations) < 1e-10] = 0
        locations[np.abs(locations) < 1e-10] = 0
        self._locations = locations
        self._orientations = orientations

    def get_rows(self):
        return self._n_rows

    def get_cols(self):
        return self._n_detector_units / self._n_rows

    def get_radius(self):
        return self._radius

    def get_beta(self):
        return self._beta

