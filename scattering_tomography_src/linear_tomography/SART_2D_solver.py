import astra
import numpy as np
from numpy import ndarray
from tqdm import tqdm

from projection_code.projection_code import ProjectionCode
from linear_tomography.abstract_ring_2D_solver import Abstract2DRingSolver
from source.abstract_source import AbstractSource
from detector.abstract_detector import AbstractDetector
from volume.grid import Grid
from volume.abstract_volume import AbstractVolume
import odl, proximal
import scipy.sparse


class SART2DSolver(Abstract2DRingSolver):
    """
    Linear tomography using SART
    """

    def __init__(self, source: AbstractSource, detector: AbstractDetector, vol: AbstractVolume, code: ProjectionCode):
        """
        Initializes the linear tomography solver
        Creates the projection matrix A using the setup geometry
        :param source:
        :param detector:
        :param grid:
        """
        super().__init__(source, detector, vol, code)

    def inverse(self, projections: np.ndarray, projections0: np.ndarray, algorithm='SIRT3D_CUDA', return_projections=False, iterations: int=100, initializer=None) -> AbstractVolume:
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

        # for i in tqdm(np.arange(num_projections), desc='projectiong images onto flat panel'):
        #     projections[:, i, :] = self.project_ring_image_to_flat(projections[:, i, :])
        #     projections0[:, i, :] = self.project_ring_image_to_flat(projections0[:, i, :])

        proj_geom = astra.create_proj_geom('cone_vec',
                                           detector_rows,
                                           detector_cols,
                                           astra_vec)

        grid = self._volume.get_grid()
        grid_shape = grid.get_shape()
        grid_voxel_size = grid.get_voxel_size()
        offset = self._volume.get_offset()

        vol_geom = astra.create_vol_geom(grid_shape[0], grid_shape[1], grid_shape[2],
                                         offset[0] - 0.5 * (grid_shape[0]),
                                         offset[0] + 0.5 * (grid_shape[0]),
                                         offset[1] - 0.5 * (grid_shape[1]),
                                         offset[1] + 0.5 * (grid_shape[1]),
                                         offset[2] - 0.5 * (grid_shape[2]),
                                         offset[2] + 0.5 * (grid_shape[2]))

        div = projections/projections0
        mask_0 = div == 0.0
        mask_nan = np.isnan(div)
        mask_inf = np.isinf(div)
        div[mask_0] = 1.0
        div[mask_nan] = 1.0
        div[mask_inf] = 1.0
        log_div = -np.log(div)
        sinogram_id = astra.data2d.create('-sino', proj_geom, log_div)
        projector_id = astra.create_projector('cuda3d', proj_geom, vol_geom)

        # Create a data object for the reconstruction
        rec_id = astra.data3d.create('-vol', vol_geom, data=0.0)
        # Set up the parameters for a reconstruction algorithm using GPU
        # cfg = astra.astra_dict('CGLS3D_CUDA')
        cfg = astra.astra_dict(algorithm)
        # cfg = astra.astra_dict('FDK_CUDA')
        cfg['ReconstructionDataId'] = rec_id
        cfg['ProjectionDataId'] = sinogram_id
        cfg['option'] = {
            'MinConstraint': 0,  # No negative values.
            'VoxelSuperSampling': True
        }

        # Create the algorithm object from the configuration structure
        alg_id = astra.algorithm.create(cfg)

        # Run algorithm according to specified number of iterations to get initial results
        astra.algorithm.run(alg_id, iterations)

        print('preforming recovery')
        recovery = astra.data3d.get(rec_id)

        astra.algorithm.delete(alg_id)
        astra.data3d.delete(rec_id)
        astra.data3d.delete(sinogram_id)
        astra.projector.delete(projector_id)

        # Transpose and flip to return in (x,y,z) vol. format
        ret_val = np.transpose(recovery, (2,1,0))
        # ret_val = np.flip(ret_val, 0)
        ret_val = ret_val*10 # from mm to cm
        if return_projections:
            return ret_val, projections, projections0
        return ret_val
