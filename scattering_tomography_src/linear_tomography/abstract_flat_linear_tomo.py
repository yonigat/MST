from abc import ABC, abstractmethod
from source.abstract_source import AbstractSource
from detector.abstract_detector import AbstractDetector
from projection_code.projection_code import ProjectionCode
from volume.abstract_volume import AbstractVolume
import numpy as np


class AbstractFlatLinearTomoSolver(ABC):
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
        self._code = code
        if len(code._projection_code.shape) == 1:
            self._code = code._projection_code.reshape((code._projection_code.shape[0], 1))

    def code2astravec(self):
        """
        :return: an ndarray of size (L,12) where L is the overall number of sources used in the experiment
        """
        num_shots = self._code._projection_code.shape[0]
        num_sources_in_shot = self._code._projection_code.shape[1] if len(self._code._projection_code.shape) > 1 else 1
        voxel_sizes = self._volume.get_grid().get_voxel_size()
        vec = np.zeros((num_shots*num_sources_in_shot, 12))
        for ind_shot in np.arange(num_shots):
            for ind_source_in_shot in np.arange(num_sources_in_shot):
                if num_sources_in_shot == 1:
                    curr_source = int(self._code._projection_code[ind_shot])
                else:
                    curr_source = int(self._code._projection_code[ind_shot, ind_source_in_shot])
                s, d, u, v = self.get_center_and_orientation(ind_shot, curr_source)
                vec[ind_source_in_shot + ind_shot * num_sources_in_shot, :3] = s / voxel_sizes
                vec[ind_source_in_shot + ind_shot * num_sources_in_shot, 3:6] = d / voxel_sizes
                vec[ind_source_in_shot + ind_shot * num_sources_in_shot, 6:9] = -u / voxel_sizes
                vec[ind_source_in_shot + ind_shot * num_sources_in_shot, 9:] = v / voxel_sizes
        return vec

    def get_center_and_orientation(self, ind_shot: int, curr_source: int):
        sizes = self._detector.get_sizes()
        d_orig = np.mean(self._detector.get_locations(), axis = 0)
        u_orig = self._detector.get_orientations()[0, :, 1]*sizes[1]
        v_orig = self._detector.get_orientations()[0, :, 2]*sizes[2]
        s_orig = self._source.get_locations()[curr_source, :]

        alpha = -np.pi/90  # a difference of 2 degrees between shots, minus because in the simulation the vol rotates by 2 deg
                            # so the screen rotates by -2 deg.
        beta = alpha*ind_shot
        transform = np.array([[np.cos(beta) , -np.sin(beta), 0],
                               [np.sin(beta), np.cos(beta) , 0],
                               [0            , 0             , 1]])
        s = np.matmul(transform, s_orig)
        d = np.matmul(transform, d_orig)
        u = np.matmul(transform, u_orig)
        v = np.matmul(transform, v_orig)
        return s, d, u, v

    @abstractmethod
    def inverse(self, projections) -> AbstractVolume:
        # TODO: the projections should match exactly the geometry setup - maybe it is better to form them under one object?
        """
        Gets the projections (measurments) that match the geometry and reconstructs a volume
        :param projections:
        :return: Reconstructed volume
        """
        pass
