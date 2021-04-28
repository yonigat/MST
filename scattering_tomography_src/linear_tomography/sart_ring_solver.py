import astra
import numpy as np
from numpy import ndarray
from tqdm import tqdm

from linear_tomography.linear_reconstruct_util import cache_reconstruction
from projection_code.projection_code import ProjectionCode
from linear_tomography.abstract_ring_linear_tomography_solver import AbstractRingLinearTomoSolver
from source.abstract_source import AbstractSource
from detector.abstract_detector import AbstractDetector
from utils.identify_outliers import non_bone_voxel_filter
from utils.read_write_cache import get_experiment_path
from volume.grid import Grid
from volume.abstract_volume import AbstractVolume
import odl, proximal
import scipy.sparse


class SARTSolver(AbstractRingLinearTomoSolver):
    """
    Linear tomography using SART
    """

    def __init__(self, source: AbstractSource, detector: AbstractDetector, vol: AbstractVolume, code: ProjectionCode,
                 cfg: dict):
        """
        Initializes the linear tomography solver
        Creates the projection matrix A using the setup geometry
        :param source:
        :param detector:
        :param grid:
        """
        super().__init__(source, detector, vol, code)
        self._cfg = cfg

    @cache_reconstruction
    def inverse(self, projections: np.ndarray, projections0: np.ndarray, algorithm='SIRT3D_CUDA',
                return_projections=False, iterations: int=100, initializer=None,
                load_from_cache: bool = True, images_from_cache: bool = False, project=True) -> AbstractVolume:
        """
        Gets the projections (measurments) that match the geometry and reconstructs a volume
        :param projections: assuming N x n x m shape
        :param algorithm: 'CGLS3D_CUDA', 'FDK_CUDA', 'FP3D_CUDA', 'BP3D_CUDA', 'SIRT3D_CUDA'
        :return: Reconstructed volume attenuation coefficients in units [1/cm]
        """
        if len(projections.shape) < 3:
            projections = np.expand_dims(projections, axis=0)
            projections0 = np.expand_dims(projections0, axis=0)

        astra_vec = self.code2astravec()
        projections = self.cut_frames(projections)
        projections0 = self.cut_frames(projections0)

        detector_rows, num_projections, detector_cols = projections.shape

        for i in tqdm(np.arange(num_projections), desc='projectiong images onto flat panel'):
            projections[:, i, :] = self.project_ring_image_to_flat(projections[:, i, :])
            projections0[:, i, :] = self.project_ring_image_to_flat(projections0[:, i, :])

        BSZ = (3, 3)
        m, n = projections[:,0,:].shape
        # for i in range(projections.shape[1]):
            # averging every BSZ pixels and due to voxels shadowing more than one pixel
            # projections[m%BSZ[0]:, i, n%BSZ[1]:] = projections[m%BSZ[0]:, i, n%BSZ[1]:].reshape(m // BSZ[0], BSZ[0], n // BSZ[1], BSZ[1]).mean(axis=(1, 3)).repeat(BSZ[0], axis=0).repeat(BSZ[1], axis=1)
            # projections0[m%BSZ[0]:, i, n%BSZ[1]:] = projections0[m%BSZ[0]:, i, n%BSZ[1]:].reshape(m // BSZ[0], BSZ[0], n // BSZ[1], BSZ[1]).mean(axis=(1, 3)).repeat(BSZ[0], axis=0).repeat(BSZ[1], axis=1)
            # projections[:, i, :] = np.mean(projections[:, i, :], axis=0) * np.ones_like(projections[:, i, :])
            # projections0[:, i, :] = np.mean(projections0[:, i, :], axis=0) * np.ones_like(projections0[:, i, :])

        distance_source_origin = np.linalg.norm(self._source.get_locations()[0, :])  # [mm]
        distance_origin_detector = np.linalg.norm(self._detector.get_locations()[0, :2])  # [mm]
        detector_rows = projections.shape[0]  # Vertical size of detector [pixels].
        detector_cols = projections.shape[2]  # Horizontal size of detector [pixels].
        num_of_projections = self._source._n_source_units
        angles = np.linspace(np.pi/2, np.pi/2 - 2*np.pi, num=num_of_projections, endpoint=False)
        if self._code.size > num_of_projections:
            angles = np.append(angles, angles[:self._code.size-num_of_projections])
        grid = self._volume.get_grid()
        grid_shape = grid.get_shape()

        proj_geom = \
            astra.create_proj_geom('cone', self._detector.get_sizes()[1], self._detector.get_sizes()[2], detector_rows,
                                   detector_cols, angles,
                                   distance_source_origin, distance_origin_detector)

        grid = self._volume.get_grid()
        grid_shape = grid.get_shape()
        grid_size = grid.get_voxel_size()
        offset = self._volume.get_offset()

        # vol_geom = astra.create_vol_geom(grid_shape[1], grid_shape[0], grid_shape[2])
        vol_geom = astra.create_vol_geom(grid_shape[1], grid_shape[0], grid_shape[2],
                                         -0.5 * (grid_size[0] * grid_shape[0]),
                                         +0.5 * (grid_size[0] * grid_shape[0]),
                                         -0.5 * (grid_size[0] * grid_shape[1]),
                                         +0.5 * (grid_size[0] * grid_shape[1]),
                                         -0.5 * (grid_size[2] * grid_shape[2]),
                                         +0.5 * (grid_size[2] * grid_shape[2]))

        div = projections/projections0
        mask_0 = div == 0.0
        mask_nan = np.isnan(div)
        mask_inf = np.isinf(div)
        div[mask_0] = 1.0
        div[mask_nan] = 1.0
        div[mask_inf] = 1.0
        div[div>1.0] = 1.0
        log_div = -np.log(div)
        # for i in range(log_div.shape[1]):
        #     log_div[:, i, :] = np.mean(log_div[:, i, :], axis=0) * np.ones_like(log_div[:, i, :])
        sinogram_id = astra.data3d.create('-sino', proj_geom, log_div)
        projector_id = astra.create_projector('cuda3d', proj_geom, vol_geom)

        # Create a data object for the reconstruction
        rec_id = astra.data3d.create('-vol', vol_geom, data=initializer)
        # Set up the parameters for a reconstruction algorithm using GPU
        cfg = astra.astra_dict(algorithm)

        cfg['ReconstructionDataId'] = rec_id
        cfg['ProjectionDataId'] = sinogram_id
        cfg['option'] = {
            'MinConstraint': 0  # No negative values.
        }

        # Create the algorithm object from the configuration structure
        alg_id = astra.algorithm.create(cfg)

        # Run algorithm according to specified number of iterations to get initial results
        niter = iterations
        astra.algorithm.run(alg_id, niter)

        print('preforming recovery')
        recovery = astra.data3d.get(rec_id)

        astra.algorithm.delete(alg_id)
        astra.data3d.delete(rec_id)
        astra.data3d.delete(sinogram_id)
        astra.projector.delete(projector_id)

        # Transpose and flip to return in (x,y,z) vol. format
        ret_val = np.transpose(recovery, (2,1,0))
        # ret_val = np.flip(ret_val, 0)
        ret_val = ret_val[::-1, :, :]

        ret_val = ret_val*10 # from mm to cm
        # ret_val[ret_val < 0.1] = 0
        # ret_val[ret_val > 0.3] *= 1.4
        # ret_val*=1.2


        # a gamma correction for bone tissues
        # gamma_val_dict = {1: 0.785, 2: 0.6, 3: 0.6, 4: 0.49, 5: 0.48}
        # th_val_dict = {1: 0.31, 2: 0.2, 3: 0.2, 4: 0.18, 5: 0.18}
        # th_std_dict = {1: 3, 2: 1.7, 3: 0.8, 4: 1.5, 5: 2.5}
        # val_th = th_val_dict.get(self._cfg['NUM_PROJECTIONS_IN_SHOT'], 0.31)
        # std_th = th_std_dict.get(self._cfg['NUM_PROJECTIONS_IN_SHOT'], 0.31)
        # filtered = non_bone_voxel_filter(ret_val, th_val=val_th, std_th=std_th)
        #
        # pow_mat = np.power(ret_val, gamma_val_dict.get(self._cfg['NUM_PROJECTIONS_IN_SHOT'], 0.785))
        # ret_val[filtered > 0] = pow_mat[filtered > 0]
        return ret_val
