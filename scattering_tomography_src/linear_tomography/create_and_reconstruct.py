import astra
import numpy as np

from detector.abstract_detector import AbstractDetector
from linear_tomography.abstract_ring_linear_tomography_solver import AbstractRingLinearTomoSolver
from projection_code.projection_code import ProjectionCode
from source.abstract_source import AbstractSource
from volume.abstract_volume import AbstractVolume


class SARTSolver(AbstractRingLinearTomoSolver):
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

    def create(self, vol_atten: np.ndarray):
        vol_atten /= 10  # convert to 1/mm
        distance_source_origin = np.linalg.norm(self._source.get_locations()[0, :])  # [mm]
        distance_origin_detector = np.linalg.norm(self._detector.get_locations()[0, :2])  # [mm]
        detector_rows = 80  # Vertical size of detector [pixels].
        detector_cols = 500  # Horizontal size of detector [pixels].
        num_of_projections = 1000  # self._source._n_source_units
        angles = np.linspace(0, 2 * np.pi, num=num_of_projections, endpoint=False)

        grid = self._volume.get_grid()
        grid_shape = grid.get_shape()
        grid_size = grid.get_voxel_size()

        # Create projections. With increasing angles, the projection are such that the
        # object is rotated clockwise. Slice zero is at the top of the object. The
        # projection from angle zero looks upwards from the bottom of the slice.
        vol_geom = astra.create_vol_geom(grid_shape[1], grid_shape[0], grid_shape[2],
                                         - (grid_size[0] * grid_shape[0]),
                                         + (grid_size[0] * grid_shape[0]),
                                         - (grid_size[0] * grid_shape[1]),
                                         + (grid_size[0] * grid_shape[1]),
                                         - (grid_size[2] * grid_shape[2]),
                                         + (grid_size[2] * grid_shape[2]))
        phantom_id = astra.data3d.create('-vol', vol_geom, data=np.transpose(vol_atten, (2, 1, 0)))
        proj_geom = \
            astra.create_proj_geom('cone', self._detector.get_sizes()[1], self._detector.get_sizes()[2], detector_rows,
                                   detector_cols, angles,
                                   distance_source_origin, distance_origin_detector)
        projections_id, projections = \
            astra.creators.create_sino3d_gpu(phantom_id, proj_geom, vol_geom)
        # projections /= np.max(projections)

        # Cleanup
        astra.data3d.delete(projections_id)
        astra.data3d.delete(phantom_id)

        return projections

    def reconstruct(self, projections: np.ndarray, algorithm='SIRT3D_CUDA', iterations: int = 100):
        distance_source_origin = np.linalg.norm(self._source.get_locations()[0, :])  # [mm]
        distance_origin_detector = np.linalg.norm(self._detector.get_locations()[0, :2])  # [mm]
        detector_rows = projections.shape[0]  # Vertical size of detector [pixels].
        detector_cols = projections.shape[2]  # Horizontal size of detector [pixels].
        num_of_projections = self._source._n_source_units
        angles = np.linspace(0, -2 * np.pi, num=num_of_projections, endpoint=False)

        grid = self._volume.get_grid()
        grid_shape = grid.get_shape()
        grid_size = grid.get_voxel_size()

        proj_geom = \
            astra.create_proj_geom('cone', 1, 1, detector_rows,
                                   detector_cols, angles,
                                   (distance_source_origin + distance_origin_detector) /
                                   self._detector.get_sizes()[1], 0)
        projections_id = astra.data3d.create('-sino', proj_geom, projections)
        vol_geom = astra.create_vol_geom(grid_shape[1], grid_shape[0], grid_shape[2])  # ,
                                         # -0.5 * (grid_size[0] * grid_shape[0]),
                                         # +0.5 * (grid_size[0] * grid_shape[0]),
                                         # -0.5 * (grid_size[1] * grid_shape[1]),
                                         # +0.5 * (grid_size[1] * grid_shape[1]),
                                         # -0.5 * (grid_size[2] * grid_shape[2]),
                                         # +0.5 * (grid_size[2] * grid_shape[2]))
        # Create reconstruction.
        # vol_geom = astra.creators.create_vol_geom(detector_cols, detector_cols,
        #                                           detector_rows)
        reconstruction_id = astra.data3d.create('-vol', vol_geom, data=0.1)
        alg_cfg = astra.astra_dict(algorithm)
        alg_cfg['ProjectionDataId'] = projections_id
        alg_cfg['ReconstructionDataId'] = reconstruction_id
        alg_cfg['option'] = {
            'MinConstraint': 0  # No negative values.
            # 'VoxelSuperSampling': True
        }
        algorithm_id = astra.algorithm.create(alg_cfg)
        astra.algorithm.run(algorithm_id, iterations=iterations)
        reconstruction = astra.data3d.get(reconstruction_id)

        # Cleanup.
        astra.algorithm.delete(algorithm_id)
        astra.data3d.delete(reconstruction_id)
        astra.data3d.delete(projections_id)

        ret_val = np.transpose(reconstruction, (1, 2, 0))
        # ret_val = np.flip(ret_val, 0)
        ret_val = ret_val * 10  # from mm to cm

        return ret_val

    def get_projections(self, projections: np.ndarray, projections0: np.ndarray):
        if len(projections.shape) < 3:
            projections = np.expand_dims(projections, axis=0)
            projections0 = np.expand_dims(projections0, axis=0)

        projections = self.cut_frames(projections)
        projections0 = self.cut_frames(projections0)

        # for i in tqdm(np.arange(projections.shape[1]), desc='projectiong images onto flat panel'):
        #     projections[:, i, :] = self.project_ring_image_to_flat(projections[:, i, :])
        #     projections0[:, i, :] = self.project_ring_image_to_flat(projections0[:, i, :])

        div = projections / projections0
        mask_0 = div == 0.0
        mask_nan = np.isnan(div)
        mask_inf = np.isinf(div)
        div[mask_0] = 1.0
        div[mask_nan] = 1.0
        div[mask_inf] = 1.0
        div[div > 1.0] = 1.0
        log_div = -np.log(div)

        return log_div

    def run_linear_tomography(self, projections: np.ndarray, projections0: np.ndarray, algorithm='SIRT3D_CUDA',
                              iterations: int = 1000):
        log_div = self.get_projections(projections, projections0)
        reconstruction = self.reconstruct(log_div, algorithm=algorithm, iterations=iterations)
        reconstruction[reconstruction < 0] = 0
        return reconstruction

    def inverse(self, projections: np.ndarray, projections0: np.ndarray, algorithm='SIRT3D_CUDA',
                return_projections=False, iterations: int = 100, initializer=None) -> AbstractVolume:
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

        vol_geom = astra.create_vol_geom(grid_shape[1], grid_shape[0], grid_shape[2],
                                         offset[0] - 0.5 * (grid_shape[0]),
                                         offset[0] + 0.5 * (grid_shape[0]),
                                         offset[0] - 0.5 * (grid_shape[1]),
                                         offset[0] + 0.5 * (grid_shape[1]),
                                         offset[2] - 0.5 * (grid_shape[2]),
                                         offset[2] + 0.5 * (grid_shape[2]))

        div = projections / projections0
        mask_0 = div == 0.0
        mask_nan = np.isnan(div)
        mask_inf = np.isinf(div)
        div[mask_0] = 1.0
        div[mask_nan] = 1.0
        div[mask_inf] = 1.0
        log_div = -np.log(div)
        sinogram_id = astra.data3d.create('-sino', proj_geom, log_div)
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
        ret_val = np.transpose(recovery, (2, 1, 0))
        # ret_val = np.flip(ret_val, 0)
        ret_val = ret_val * 10  # from mm to cm
        return ret_val
