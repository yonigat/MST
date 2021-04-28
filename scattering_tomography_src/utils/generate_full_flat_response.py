import projections.projections as projection_class
import detector.ring_detector
import source.ring_source
import detector.flat_panel_detector
import volume.grid
import projection_code.projection_code
from nontrivial_projection_code import NonTrivialProjectionCode
from illum_code_optim.generate_src2pix_map import generate_illum_mat, generate_full_code
from illum_code_optim.illum_loss import optimize_illum
import numpy as np
import os
import g4driver.g4driver
import linear_tomography.sart_ring_solver
from linear_tomography.ct_to_vol_util import create_vol_simple_table, create_vol_from_ct_recovery, \
    construct_atten_coeff_vol
from optimization.ADAM import ADAM
from optimization.oracle import ring_oracle_1st_order_helper
from utils.params import read_configs
from volume.abstract_volume import load_volume, save_volume
from volume.xcat_volume_util import create_xcat
from utils.generate_run_name import generate_run_name
from utils.tf_logger import Logger
from utils.read_write_cache import get_experiment_path
from g4driver.g4driver_util import save_results_to_dir


def merge_images(cfg_file: str, single_flat_response, projections=None):
    class Shifter:
        def __init__(self, src: source.ring_source.RingSource, detectors: detector.ring_detector.RingDetector):
            self._src = src
            self._detector = detectors
            self._orig_src = 19
        def mid_detec_ind(self, src):
            src_loc = self._src.get_locations()[src, 0:2]  # only x and y are needed
            src_angle = self.set_angle_to_range(np.arctan2(src_loc[1], src_loc[0]))   # angle of the source
            num_cols = self._detector.get_cols()
            angle_between_detectors = (2 * np.pi) / num_cols  # angle between 2 detectors
            gamma = self._detector.get_beta()
            angle_of_middle_det = self.set_angle_to_range(src_angle + np.pi)

            mid_col_ind = 0
            mid_col_ind_angle = 0
            if angle_of_middle_det >= 0:
                while mid_col_ind_angle < angle_of_middle_det:
                    mid_col_ind_angle += angle_between_detectors
                    mid_col_ind += 1
                if mid_col_ind_angle - gamma > angle_of_middle_det:
                    mid_col_ind -= 1
            elif angle_of_middle_det < 0:
                mid_col_ind = int(num_cols)
                while mid_col_ind_angle > angle_of_middle_det:
                    mid_col_ind_angle -= angle_between_detectors
                    mid_col_ind -= 1
                if mid_col_ind_angle + gamma < angle_of_middle_det:
                    mid_col_ind += 1
            return mid_col_ind

        @staticmethod
        def set_angle_to_range(angle: float):
            """
            :param angle:
            :return: sets an angle (in radians) to the range of [pi, -pi)
            """
            new_angle = angle
            if new_angle > np.pi:
                while new_angle > np.pi:
                    new_angle -= 2 * np.pi
            elif new_angle <= -np.pi:
                while new_angle <= -np.pi:
                    new_angle += 2 * np.pi
            return new_angle

        def __call__(self, image: np.ndarray, curr_src):
            mid_det_col_ind = self.mid_detec_ind(curr_src)
            mid_det_col_orig = self.mid_detec_ind(self._orig_src)
            new_image = np.zeros_like(image)
            for ind in np.arange(image.shape[0]):
                new_image[ind] = np.roll(image[ind], mid_det_col_ind- mid_det_col_orig, axis=0)
            return new_image

    # multi_GT_dir = 'bla'
    cfg = read_configs(cfg_file)  # create a dictionary of all the paramters, it will also be sent to the G4 sim
    single_GT_path = get_experiment_path(cfg, GT_data=True)

    # parameters

    # sources
    radius_src = cfg['SOURCE_TO_CENTER']
    num_sources = cfg['NUM_OF_SOURCES']
    angle_theta = cfg['MAX_THETA']
    angle_delta = cfg["MAX_PHI"]

    # detectors
    radius_det = cfg["CENTER_TO_DET"]
    num_det_cols = cfg["NUM_DETECTOR_COLS"]
    num_det_rows = cfg["NUM_DETECTOR_ROWS"]
    size_det_X = cfg["DETECTOR_X"] * 2
    size_det_Y = cfg["DETECTOR_Y"] * 2
    size_det_Z = cfg["DETECTOR_Z"] * 2
    num_scorers = cfg["NUM_OF_SCORERS"]

    # grid
    num_voxels_X = cfg["NUM_OF_VOXELS_X"]
    num_voxels_Y = cfg["NUM_OF_VOXELS_Y"]
    num_voxels_Z = cfg["NUM_OF_Z_SLICES"]
    size_voxel_X = cfg["VOXEL_HALF_X"] * 2
    size_voxel_Y = cfg["VOXEL_HALF_Y"] * 2
    size_voxel_Z = cfg["VOXEL_HALF_Z"] * 2

    # volume
    active_elements = cfg['ACTIVE_ELEMENTS']
    vol_offset_X = cfg["PHANTOM_OFFSET_X"]
    vol_offset_Y = cfg["PHANTOM_OFFSET_Y"]
    vol_offset_Z = cfg["PHANTOM_OFFSET_Z"]
    volume_offset = np.array([vol_offset_X, vol_offset_Y, vol_offset_Z])

    # projection code
    num_shots = cfg["NUM_OF_SHOTS_GT"]
    num_GT_photons = cfg["NUM_OF_PHOTONS_GT"]
    num_overall_photons = num_GT_photons * num_shots

    num_projection_in_shot = cfg["NUM_PROJECTIONS_IN_SHOT"]
    energy = 50  # Kev
    src = source.ring_source.RingSource(radius_src, num_sources, (angle_theta, angle_delta))
    detect = detector.ring_detector.RingDetector(radius_det, num_det_cols, num_det_rows,
                                                 (size_det_X, size_det_Y, size_det_Z))
    grid = volume.grid.Grid((num_voxels_X, num_voxels_Y, num_voxels_Z,), (size_voxel_X, size_voxel_Y, size_voxel_Z))
    vol = create_xcat(cfg, 'xcat_util', 'stomach', grid, offset=volume_offset, num_of_materials=cfg["NUM_OF_MATERIALS"],
                      load_from_cache=True)
    # vol.blur_vol()
    # vol.AGWN(sigma=0.03)
    vol.reduce_materials(cfg["NUM_OF_MATERIALS"])
    # load spectrum array
    # src.create_120keV_spect()
    src.create_120keV_spect_trimed()
    # src.create_60keV_delta()

    code = projection_code.projection_code.ProjectionCode(num_sources, num_shots, num_projection_in_shot)
    code.create_lin_space_code()
    code.rand_perm_existing_code(code.get_code())

    projections = projection_class.Projections(active_scorer=2)  # We want to compare using direct transmission only

    G4drive = g4driver.g4driver.G4Driver(cfg, detect, src, projections)

    print('Running GT with phantom')
    (I_GT, _), loaded_from_cache = G4drive.runG4(vol, code, grad=False, GT=True, build_phantom=True)
    # I_GT *= 1e3  # 1000 - turn from MeV to KeV, num_GT_photons for normalization (assumig single source)
    # I_GT /= num_sources
    projections.set_I_GT(I_GT, built_phantom=True)
    print('Running GT without phantom')

    I0_GT = np.zeros_like(I_GT)
    shifted_image = Shifter(src, detect)
    for shot_ind, shot_code in enumerate(code.get_code()):
        I0_curr = I0_GT[shot_ind, :, :, :]
        for curr_source in shot_code:
            I0_curr += shifted_image(single_flat_response, curr_source)
        I0_GT[shot_ind] = I0_curr / num_overall_photons
    projections.set_I_GT(I0_GT, built_phantom=False)

    I0_multi_meta_data = ((single_GT_path, vol, code), {'grad': False, 'GT': True, 'build_phantom': False})
    save_results_to_dir(single_GT_path, I0_multi_meta_data, cfg, (I0_GT, _), 'I0_GT')

    class Shifter:
        def __init__(self, src: source.ring_source.RingSource, detectors: detector.ring_detector.RingDetector):
            self._src = src
            self._detect = detectors

        def __call__(self, image: np.ndarray, curr_src):
            mid_det_col_ind = self.mid_detec_ind(curr_src)
            new_image = np.zeros_like(image)
            for ind in np.arange(image.shape[0]):
                new_image[ind] = np.roll(image[ind], mid_det_col_ind, axis=1)
            return new_image

        def mid_detec_ind(self, src):
            src_loc = self._src.get_locations()[src, 0:2]  # only x and y are needed
            src_angle = self.set_angle_to_range(np.arctan2(src_loc[1], src_loc[0]))   # angle of the source
            num_cols = self._detector.get_cols()
            angle_between_detectors = (2 * np.pi) / num_cols  # angle between 2 detectors
            gamma = self._detector.get_beta()
            angle_of_middle_det = self.set_angle_to_range(src_angle + np.pi)

            mid_col_ind = 0
            mid_col_ind_angle = 0
            if angle_of_middle_det >= 0:
                while mid_col_ind_angle < angle_of_middle_det:
                    mid_col_ind_angle += angle_between_detectors
                    mid_col_ind += 1
                if mid_col_ind_angle - gamma > angle_of_middle_det:
                    mid_col_ind -= 1
            elif angle_of_middle_det < 0:
                mid_col_ind = num_cols
                while mid_col_ind_angle > angle_of_middle_det:
                    mid_col_ind_angle -= angle_between_detectors
                    mid_col_ind -= 1
                if mid_col_ind_angle + gamma < angle_of_middle_det:
                    mid_col_ind += 1

        @staticmethod
        def set_angle_to_range(angle: float):
            """
            :param angle:
            :return: sets an angle (in radians) to the range of [pi, -pi)
            """
            new_angle = angle
            if new_angle > np.pi:
                while new_angle > np.pi:
                    new_angle -= 2 * np.pi
            elif new_angle <= -np.pi:
                while new_angle <= -np.pi:
                    new_angle += 2 * np.pi
            return new_angle

cfg_file = '../tests/GT_generate_new_stomach_3D_cfg.yaml'
cfg = read_configs(cfg_file)  # create a dictionary of all the paramters, it will also be sent to the G4 sim
single_GT_path = get_experiment_path(cfg, GT_data=True)
for i in range(5):
    file = f'run_0outputEnergyDep_{i}.csv'
    projection = np.loadtxt(os.path.join(single_GT_path, file), delimiter=',', skiprows=1)
    if i == 0:
        I0 = np.zeros((5, projection.shape[0], projection.shape[1]))
    I0[i] = projection
merge_images(cfg_file, I0)
