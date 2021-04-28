from abc import ABC, abstractmethod
from source.abstract_source import AbstractSource
from detector.abstract_detector import AbstractDetector
from projection_code.projection_code import ProjectionCode
from geometry.ring_geometry import RingGeometry
import scipy.interpolate
import numpy as np
from numpy import ndarray
from volume.abstract_volume import AbstractVolume


class Abstract2DRingSolver(ABC):
    """
    Interface for linear tomography solver
    """

    def __init__(self, source: AbstractSource, detector: AbstractDetector, vol: AbstractVolume, code: ProjectionCode):
        """
        Initializes the linear tomography solver
        Creates the projection matrix A using the setup geometry
        :param source:
        :param detector:
        :param grid:

        """
        self._source = source
        self._detector = detector
        self._volume = vol
        self._geometry = RingGeometry(detector, source, vol)
        self._code = code._projection_code
        if len(code._projection_code.shape) == 1:
            self._code = code._projection_code.reshape((code._projection_code.size, 1))

    def code2astravec(self):
        """
        :return: an ndarray of size (L,12) where L is the overall number of sources used in the experiment
        """
        num_shots = self._code.shape[0]
        num_sources_in_shot = self._code.shape[1] if len(self._code.shape) > 1 else 1
        voxel_sizes = self._volume.get_grid().get_voxel_size()
        vec = np.zeros((num_shots*num_sources_in_shot, 6))
        for ind_shot in np.arange(num_shots):
            for ind_source_in_shot in np.arange(num_sources_in_shot):
                if self._code.shape[1] == 1:
                    curr_source = int(self._code[ind_shot])
                else:
                    curr_source = int(self._code[ind_shot, ind_source_in_shot])
                curr_borders, mid_angle = self._geometry.get_projection(curr_source)
                s = self._geometry._sources._locations[curr_source, :]
                d, u, v = self._geometry.get_center_and_orientation(mid_angle)
                vec[ind_source_in_shot + ind_shot * num_sources_in_shot, :2] = s[0:2] / voxel_sizes
                vec[ind_source_in_shot + ind_shot * num_sources_in_shot, 2:4] = d[0:2] / voxel_sizes
                vec[ind_source_in_shot + ind_shot * num_sources_in_shot, 4:] = -u[0:2] / voxel_sizes
        return vec

    def cut_frames(self, full_images: ndarray):
        """
        :param full_images: ndarray N x n x m: all images of all the detectors in every shot
        :return: the sub images of all the full images, each sub image is cut according to the geometry of the cone
        """
        num_shots, num_sources_in_shot = self._code.shape
        max_length = 0  # maximal size of all the subimages we cut from the full frame.
        borders = np.zeros((num_shots*num_sources_in_shot, 2), dtype=int)
        angles = np.zeros((num_shots*num_sources_in_shot,))

        for ind_shot in np.arange(num_shots):
            for ind_source_in_shot in np.arange(num_sources_in_shot):
                full_image = full_images[ind_shot, :, :]  # take the image relevant for the current shot
                curr_shot = ind_shot*num_sources_in_shot+ind_source_in_shot
                if self._code.shape[1] == 1:
                    curr_source = int(self._code[ind_shot])
                else:
                    curr_source = int(self._code[ind_shot, ind_source_in_shot])
                borders[curr_shot, :], angles[curr_shot] = self._geometry.get_projection(curr_source)
                curr_length = int(self.sub_image_size(borders[curr_shot, :], full_image.shape))
                if curr_length > max_length:
                    max_length = curr_length

        sub_images = np.zeros((full_images.shape[1], num_shots*num_sources_in_shot, max_length))
        for ind_shot in np.arange(num_shots):
            for ind_source_in_shot in np.arange(num_sources_in_shot):
                full_image = full_images[ind_shot, :, :]  # take the image relevant for the current shot
                curr_shot = ind_shot * num_sources_in_shot + ind_source_in_shot
                sub_image = self.get_sub_image(borders[curr_shot, :], full_image)
                sub_images[:, curr_shot, :sub_image.shape[1]] = sub_image

        return sub_images

    def project_ring_image_to_flat(self, I: ndarray):

        src_loc = self._geometry._sources.get_locations()[0,:]
        detector_unit_size = np.array(self._geometry._detectors.get_sizes()[1:])
        detector_size_in_mm =  detector_unit_size * I.T.shape
        cols = np.linspace(-(detector_size_in_mm[0]-detector_unit_size[0])/2, (detector_size_in_mm[0]-detector_unit_size[0])/2, num=I.shape[1],
                           endpoint=True)
        rows = np.linspace(-(detector_size_in_mm[1] - detector_unit_size[1]) / 2,
                           (detector_size_in_mm[1] - detector_unit_size[1]) / 2, num=I.shape[0],
                           endpoint=True)
        f = scipy.interpolate.interp2d(cols, rows, I, kind='linear')
        N = len(rows)*len(cols)
        detector_loc = np.zeros((N, 3))
        detector_loc[:,0] = self._geometry._detectors.get_radius()
        detector_loc[:,1] = np.repeat(cols, N//len(cols))
        detector_loc[:,2] = np.tile(rows, N//len(rows))

        projected_pixels = self.project_flat2cone(detector_loc, src_loc)

        ring_loc = self._geometry._sources.get_radius() * np.arctan2(projected_pixels[:, 1], projected_pixels[:, 0])
        I_new = np.zeros_like(I)
        for ind_row, row_loc in enumerate(rows):
            for ind_col, col_loc in enumerate(cols):
                eff_ind = ind_col * len(rows) + ind_row
                I_new[ind_row, ind_col] = f(ring_loc[eff_ind], projected_pixels[eff_ind,2])
        # I_new = I_new.reshape(I.shape)
        return I_new


    @staticmethod
    def project_flat2cone(C_prime: ndarray, S: ndarray):
        """
        This function returns the projection of a point on a cylinder (the detector array) from a line given by the location
        of the source and a point C_prime on a flat panel which is tangent to the cylinder.
        The projection was based on the calculation given in:
        https://johannesbuchner.github.io/intersection/intersection_line_cylinder.html
        :param S: size(3,) coordinates (in mm) of the source
        :param C_prime:  size(N,3) coordinates (in mm) of the point on the flat panel we wish to project
        :return: C size(N,3) the coordinates (in mm) of the projection of C_prime onto a cone.
        """
        R = C_prime[0,0]

        t = -S[0] + C_prime[:,0]
        k = (-S[1] + C_prime[:, 1]) / t
        l = (-S[2] + C_prime[:, 2]) / t

        t = (1 / (k ** 2 + 1)) * (- S[1] * k - S[0] + np.sqrt(
            R ** 2 * (k ** 2 + 1) - (k * S[0]) ** 2 + 2 * k * S[0] * S[1] - S[1] ** 2))
        line_param = np.array([t, k * t, l * t]).T

        C = S + line_param
        return C


    @staticmethod
    def sub_image_size(borders: ndarray, full_image_shape: tuple):
        """
        The method returns the size of the sub image according to the borders, the lower border can be larger than the higher,
        since the detector array is circular.

        :param borders:
        :param full_image_shape:
        :return:
        """
        Ny = full_image_shape[1]
        if np.any(borders < 0) or np.any(borders > Ny):
            raise Exception("Borders of the image do not fit the number of columns in the image")

        if borders[0] > borders[1]:
            length = Ny - borders[0] + borders[1]+1
        elif borders[0] < borders[1]:
            length = borders[1] - borders[0] + 1
        else:
            length = 1
        return length

    @staticmethod
    def get_sub_image(borders: ndarray, full_image: ndarray):
        """
        A method that extracs part of a full image, according to the borders given. It's possible that the lower boundry
        will be larger than the upper, this just means we make the image cyclic continues
        :param borders: 2x1 ndarray, contains the lower and upper bounds of the image we wish to extract
        :param full_image: N x M ndarray, the full image from which we wish to extract part of the image
        :return:
        """
        # TODO: should the image be given, or maybe part of the detectors class? for now we assume it's given.

        # check that the value of the borders are legal
        Nx, Ny = full_image.shape
        if np.any(borders < 0) or np.any(borders > Ny):
            raise Exception("Borders of the image do not fit the number of columns in the image")

        if borders[0] > borders[1]:
            image = np.zeros((Nx, Ny - borders[0] + borders[1]+1))
            image[:, :Ny-borders[0]] = full_image[:, borders[0]:]
            image[:, Ny-borders[0]:] = full_image[:, :borders[1]+1]

        elif borders[0] < borders[1]:
            image = full_image[:, borders[0]:(borders[1]+1)]
        else:
            image = full_image[:, borders[0]]
            image = image.reshape((image.shape[0], 1))
        # return np.fliplr(image)  # the flip is needed to match between how the image is seen on the detector array, and
                                # how it is stored in the image

        return np.flipud(image)

    @abstractmethod
    def inverse(self, projections, projections0):
        # TODO: the projections should match exactly the geometry setup - maybe it is better to form them under one object?
        """
        Gets the projections (measurments) that match the geometry and reconstructs a volume
        :param projections:
        :return: Reconstructed volume
        """
        pass
