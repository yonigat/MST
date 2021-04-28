from source.abstract_source import AbstractSource
import numpy as np


class RingSource(AbstractSource):
    """
    Creates a ring source
    """

    def __init__(self, radius: float, n_units_along_ring: int, unit_focal_angles: tuple, plane_normal: str = 'Z'):
        """
        Creates the ring source from the given specs
        :param radius: the ring radius
        :param plane_normal: the normal to the plane on which the ring lies, can be 'X', 'Y', 'Z'
        :param n_units_along_ring: number of source elements along the ring
        :param unit_focal_angles: tuple, the cone beam opening angles (theta, delta) where the opening is symmetric
                theta is the azimuthal (polar) angle, measured relatively to positive x axis, values are [0, pi]
                delta is the latitude angle, measured relatively to a vector on the XY plane, values ar [0, pi/2]
                it is related to phi (the zenith angle in spherical coordinates) by phi = pi/2 - delta
        """

        def totuple(a):
            try:
                return tuple(totuple(i) for i in a)
            except TypeError:
                return a

        if plane_normal != 'Z':
            raise Exception("Currently plane normal supports only Z direction")

        super().__init__()
        self._n_source_units = n_units_along_ring
        self._focal_angles = totuple((np.pi) / np.array(unit_focal_angles))
        self._radius = radius

        dims = 3
        alpha = 2 * np.pi / n_units_along_ring  # angle between each two detectors
        trans_dict = {'X': np.array([0, radius, 0]), 'Y': np.array([radius, 0, 0]), 'Z': np.array([-radius, 0, 0])}

        # for now we assume all the sources are evenly spaced around the middle of the ring (i.e. on the XY plane) and
        # are all oriented towards the origin
        orientations = np.zeros((self._n_source_units, dims, dims))  # each source has 3 orthogonal orientations
                                                                       # 1. The direction in which the cone opens with angle theta
                                                                       # 2. The direction in which the cone opens with angle delta
                                                                       # 3. The direction which the cone opnes
        locations = np.zeros((self._n_source_units, dims))

        """
                creating positions and orientations of detectors
                we assume plane normal == 'Z'
                In future implementations if a different plane normal is needed - a simple linear transformation around the axis can be applied 
                """

        orig_orientation = np.array(((0, 1, 0), (0, 0, 1), (1, 0, 0))).T  # each column is an orientation of the chest in original position
        for source_ind in np.arange(n_units_along_ring):
            rot_mat = np.array([[np.cos(source_ind * alpha),  -np.sin(source_ind * alpha), 0],
                                [np.sin(source_ind * alpha), np.cos(source_ind * alpha), 0],
                                [0,                           0,                          1]])
            locations[source_ind, :] = np.dot(rot_mat, trans_dict[plane_normal].T).T  # transform to correct location in XY plane
            orientations[source_ind, :] = np.dot(rot_mat, orig_orientation)  # translate the orientations (each column is an orientation)

        # if needed we can now perform a transformation on the locations and orientations to fit a different normal
        # than Z
        orientations[np.abs(orientations) < 1e-10] = 0
        locations[np.abs(locations) < 1e-10] = 0
        self._locations = locations
        self._orientations = orientations

    def get_radius(self):
        return self._radius

    def get_focal_angle(self):
        """
        :return: the focal angles of the sources as a ndarray
        """
        return np.array(self._focal_angles)


