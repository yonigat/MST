from detector.abstract_detector import AbstractDetector
import numpy as np



class FlatPanelDetector(AbstractDetector):
    """
    Creates a ring detector
    """

    def __init__(self, radius: float, n_cols: int, n_rows: int, unit_size: tuple, plane_normal: str = 'X'):
        """
        Creates the ring detector from the given specs
        :param radius: the ring radius
        :param plane_normal: the normal to the plane on which the ring lies, can be 'X', 'Y', 'Z'
        :param n_units_along_ring: number of detector elements along the ring
        :param n_rows: number of rows (concatenation of rings along the orthogonal axis to the radius)
        :param unit_size: tuple of size (x_size: float, y_size: float, z_size: float) of each detector unit
        """
        if plane_normal != 'X':
            raise Exception("Currently plane normal supports only X direction")

        super().__init__()
        self._n_detector_units = n_cols * n_rows
        self._sizes = unit_size
        self._n_rows = n_rows
        self._radius = radius

        dims = 3

        trans_dict = {'X': np.array([0, radius, 0]), 'Y': np.array([radius, 0, 0]), 'Z': np.array([radius, 0, 0])}
        orientations = np.zeros((self._n_detector_units, dims, dims))  # each detector has 3 orthogonal orientations
        locations = np.zeros((self._n_detector_units, dims))  # each detectors location is determined by its center

        """
        creating positions and orientations of detectors
        we assume plane normal == 'Z'
        In future implementations if a different plane normal is needed - a simple linear transformation around the axis can be applied 
        """
        orig_orientation = np.eye(3)  # each column is an orientation of the chest in original position
        height = n_rows*unit_size[2]
        width = n_cols*unit_size[1]  # total length of each column along Z axis
        for col_ind in np.arange(n_cols):
            y_loc = -width/2 + unit_size[1]/2 + col_ind*unit_size[1]
            for row_ind in np.arange(n_rows):
                z_loc = -height/2 + unit_size[2]/2 + row_ind*unit_size[2]
                locations[col_ind * n_rows + row_ind, :] = np.array([radius, y_loc, z_loc])  # translate to correct location along Z axis
                orientations[col_ind * n_rows + row_ind, :, :] = np.eye(dims)
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


