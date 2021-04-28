import astra
import numpy as np
import odl
import proximal

from detector.abstract_detector import AbstractDetector
from linear_tomography.abstract_flat_linear_tomo import AbstractFlatLinearTomoSolver
from projection_code.projection_code import ProjectionCode
from source.abstract_source import AbstractSource
from volume.abstract_volume import AbstractVolume


class FlatPanelSolver(AbstractFlatLinearTomoSolver):
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
        distance_origin_detector = np.linalg.norm(self._detector.get_locations()[0, :1])  # [mm]
        detector_rows = 224  # Vertical size of detector [pixels].
        detector_cols = 272  # Horizontal size of detector [pixels].
        num_of_projections = 180  # self._source._n_source_units
        angles = np.linspace(np.pi/2, np.pi/2 + 2*np.pi, num=num_of_projections, endpoint=False)

        grid = self._volume.get_grid()
        grid_shape = grid.get_shape()
        grid_size = grid.get_voxel_size()

        # Create projections. With increasing angles, the projection are such that the
        # object is rotated clockwise. Slice zero is at the top of the object. The
        # projection from angle zero looks upwards from the bottom of the slice.
        vol_geom = astra.create_vol_geom(grid_shape[1], grid_shape[0], grid_shape[2],
                                         -0.5 * (grid_size[0] * grid_shape[0]),
                                         +0.5 * (grid_size[0] * grid_shape[0]),
                                         -0.5 * (grid_size[0] * grid_shape[1]),
                                         +0.5 * (grid_size[0] * grid_shape[1]),
                                         -0.5 * (grid_size[2] * grid_shape[2]),
                                         +0.5 * (grid_size[2] * grid_shape[2]))
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

    def reconstruct(self, projections: np.ndarray, algorithm='SIRT3D_CUDA', iterations: int = 50) -> AbstractVolume:
        distance_source_origin = np.linalg.norm(self._source.get_locations()[0, :])  # [mm]
        distance_origin_detector = np.linalg.norm(self._detector.get_locations()[0, :1])  # [mm]
        detector_rows = projections.shape[0]  # Vertical size of detector [pixels].
        detector_cols = projections.shape[2]  # Horizontal size of detector [pixels].
        num_of_projections = len(self._code.get_code())
        angles = np.linspace(np.pi / 2, np.pi / 2 + 2 * np.pi, num=num_of_projections, endpoint=False)

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
        sinogram_id = astra.data3d.create('-sino', proj_geom, projections)
        projector_id = astra.create_projector('cuda3d', proj_geom, vol_geom)

        # Create a data object for the reconstruction
        rec_id = astra.data3d.create('-vol', vol_geom, data=0)
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
        ret_val = np.transpose(recovery, (2, 1, 0))
        # ret_val = np.flip(ret_val, 0)
        # ret_val = ret_val[::-1, :, :]
        ret_val = ret_val * 10  # from mm to cm

        return ret_val

    def inverse(self, projections: np.ndarray, projections0: np.ndarray, algorithm='SIRT3D_CUDA', iterations: int = 50) -> AbstractVolume:
        """
        Gets the projections (measurments) that match the geometry and reconstructs a volume
        :param projections:
        :return: Reconstructed volume attenuation coefficients in units [1/cm]
        :param algorithm: 'CGLS3D_CUDA', 'FDK_CUDA', 'FP3D_CUDA', 'BP3D_CUDA', 'SIRT3D_CUDA'
        """
        if len(projections.shape) < 3:
            projections = np.expand_dims(projections, axis=0)
            projections0 = np.expand_dims(projections0, axis=0)
        projections = projections.transpose(1, 0, 2) # we assume projections is [N_SHOTS,N_ROWS,N_COLS] and ASTRA needs
                                                    # the projections in [N_ROWS, N_SHOTS, N_COLS] format
        projections0 = projections0.transpose(1, 0, 2)
        astra_vec = self.code2astravec()
        detector_rows, num_projections, detector_cols = projections.shape

        distance_source_origin = np.linalg.norm(self._source.get_locations()[0, :])  # [mm]
        distance_origin_detector = np.linalg.norm(self._detector.get_locations()[0, :1])  # [mm]
        detector_rows = projections.shape[0]  # Vertical size of detector [pixels].
        detector_cols = projections.shape[2]  # Horizontal size of detector [pixels].
        num_of_projections = len(self._code.get_code())
        angles = np.linspace(np.pi/2, np.pi/2 + 2*np.pi, num=num_of_projections, endpoint=False)

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
        sinogram_id = astra.data3d.create('-sino', proj_geom, log_div)
        projector_id = astra.create_projector('cuda3d', proj_geom, vol_geom)

        # Create a data object for the reconstruction
        rec_id = astra.data3d.create('-vol', vol_geom, data=0)
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
        # ret_val = ret_val[::-1, :, :]
        ret_val = ret_val*10 # from mm to cm
        return ret_val
