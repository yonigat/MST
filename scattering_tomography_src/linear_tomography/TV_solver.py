import astra
import numpy as np
from numpy import ndarray
from tqdm import tqdm

from projection_code.projection_code import ProjectionCode
from linear_tomography.abstract_ring_linear_tomography_solver import AbstractRingLinearTomoSolver
from source.abstract_source import AbstractSource
from detector.abstract_detector import AbstractDetector
from volume.grid import Grid
from volume.abstract_volume import AbstractVolume
import odl, proximal
import scipy.sparse


class TV_solver(AbstractRingLinearTomoSolver):
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
        distance_source_origin = np.linalg.norm(self._source.get_locations()[0, :])  # [mm]
        distance_origin_detector = np.linalg.norm(self._detector.get_locations()[0, :2])  # [mm]
        detector_pixel_size = self._detector.get_sizes()[1]  # [mm]
        detector_rows = projections.shape[0]  # Vertical size of detector [pixels].
        detector_cols = projections.shape[2]  # Horizontal size of detector [pixels].
        num_of_projections = self._source._n_source_units
        angles = np.linspace(0, -2 * np.pi, num=num_of_projections, endpoint=False)

        grid = self._volume.get_grid()
        grid_shape = grid.get_shape()
        grid_size = grid.get_voxel_size()

        proj_geom = \
            astra.create_proj_geom('cone', self._detector.get_sizes()[1], self._detector.get_sizes()[2], detector_rows,
                                   detector_cols, angles,
                                   distance_source_origin, distance_origin_detector)
        projections_id = astra.data3d.create('-sino', proj_geom, projections)
        vol_geom = astra.create_vol_geom(grid_shape[1], grid_shape[0], grid_shape[2],
                                         -0.5 * (grid_size[0] * grid_shape[0]),
                                         +0.5 * (grid_size[0] * grid_shape[0]),
                                         -0.5 * (grid_size[1] * grid_shape[1]),
                                         +0.5 * (grid_size[1] * grid_shape[1]),
                                         -0.5 * (grid_size[2] * grid_shape[2]),
                                         +0.5 * (grid_size[2] * grid_shape[2]))
        # Create reconstruction.
        # vol_geom = astra.creators.create_vol_geom(detector_cols, detector_cols,
        #                                           detector_rows)
        reconstruction_id = astra.data3d.create('-vol', vol_geom, data=0)
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

        proj_operator = astra.OpTomo(projections_id)
        """
        Tomography with TV regularization using the ProxImaL solver.
        Solves the optimization problem
            min_{0 <= x }  ||A(x) - g||_2^2 + alpha || |grad(x)| ||_1
        Where ``A`` is a parallel beam forward projector, ``grad`` the spatial
        gradient and ``g`` is given noisy data.
        """

        class Projection(odl.Operator):
            def __init__(self, OpTomo):
                self.OpTomo = OpTomo
                self.adjoint = OpTomo.adjoint()
                dom = odl.rn(OpTomo.shape[1])
                ran = odl.rn(OpTomo.shape[0])
                odl.Operator.__init__(self, dom, ran)

            def _call(self, x):
                return self.OpTomo(x)

            def adjoint(self, x):
                return self.adjoint(x)

        # Initialize the ray transform (forward projection).
        ray_trafo = Projection(proj_operator)

        # Convert ray transform to proximal language operator
        proximal_lang_ray_trafo = odl.as_proximal_lang_operator(ray_trafo)

        # Convert to array for ProxImaL
        rhs_arr = projections.ravel()

        # Set up optimization problem
        # Note that proximal is not aware of the underlying space and only works with
        # matrices. Hence the norm in proximal does not match the norm in the ODL space exactly.
        alpha = 1
        x = proximal.Variable(grid_shape[0] * grid_shape[1] * grid_shape[2])
        funcs = proximal.sum_squares(proximal_lang_ray_trafo(x) - rhs_arr) + proximal.nonneg(x) + \
                proximal.norm1(proximal.grad(proximal.reshape(x, (grid_shape[0], grid_shape[1], grid_shape[2]))))

        # Solve the problem using ProxImaL
        prob = proximal.Problem(funcs)
        prob.solve(verbose=True)
        recovery = x.value.reshape((grid_shape[20], grid_shape[0], grid_shape[1]))

        astra.algorithm.delete(algorithm_id)
        astra.data3d.delete(reconstruction_id)
        astra.data3d.delete(sinogram_id)
        astra.projector.delete(projections_id)

        # Transpose and flip to return in (x,y,z) vol. format
        # ret_val = np.transpose(recovery, (2,1,0))
        # ret_val = np.flip(ret_val, 0)
        ret_val = recovery*10 # from mm to cm
        if return_projections:
            return ret_val, projections, projections0
        return ret_val
