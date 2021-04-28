import sys
from copy import deepcopy

sys.path.append("/home/yonatangat/Projects/scattering_tomo/scattering_tomography_src")
import projections.projections
import detector.ring_detector
import source.ring_source
import detector.flat_panel_detector
import volume.grid
import projection_code.projection_code
import numpy as np
import os
import g4driver.g4driver
import linear_tomography.sart_ring_solver
import linear_tomography.TV_solver
import linear_tomography.create_and_reconstruct
from linear_tomography.ct_to_vol_util import create_vol_simple_table, create_vol_from_ct_recovery, \
    construct_atten_coeff_vol
from optimization.ADAM import ADAM
from optimization.oracle import ring_oracle_1st_order_helper, constrain_func_helper, adjust_func_helper
from utils.params import read_configs
from volume.abstract_volume import load_volume, save_volume
from volume.xcat_volume_util import create_xcat
from volume.toy_volume import ToyVolume
from utils.generate_run_name import generate_run_name
from utils.tf_logger import Logger
from utils.read_write_cache import get_experiment_path

import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.use('Agg')

cfg = read_configs('full_optimization_cfg.yaml')  # create a dictionary of all the paramters, it will also be sent to the G4 sim

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
size_det_X = cfg["DETECTOR_X"]*2
size_det_Y = cfg["DETECTOR_Y"]*2
size_det_Z = cfg["DETECTOR_Z"]*2
num_scorers = cfg["NUM_OF_SCORERS"]

# grid
num_voxels_X = cfg["NUM_OF_VOXELS_X"]
num_voxels_Y = cfg["NUM_OF_VOXELS_Y"]
num_voxels_Z = cfg["NUM_OF_Z_SLICES"]
size_voxel_X = cfg["VOXEL_HALF_X"]*2
size_voxel_Y = cfg["VOXEL_HALF_Y"]*2
size_voxel_Z = cfg["VOXEL_HALF_Z"]*2

# volume
active_elements = cfg['ACTIVE_ELEMENTS']
vol_offset_X = cfg["PHANTOM_OFFSET_X"]
vol_offset_Y = cfg["PHANTOM_OFFSET_Y"]
vol_offset_Z = cfg["PHANTOM_OFFSET_Z"]
volume_offset = np.array([vol_offset_X, vol_offset_Y, vol_offset_Z])

# projection code
num_shots_opt = cfg["NUM_OF_SHOTS_OPTIMIZATION"]
num_shots = cfg["NUM_OF_SHOTS_GT"]

num_projection_in_shot = cfg["NUM_PROJECTIONS_IN_SHOT"]
energy = 50  # Kev


src = source.ring_source.RingSource(radius_src, num_sources, (angle_theta, angle_delta))
detect = detector.ring_detector.RingDetector(radius_det, num_det_cols, num_det_rows, (size_det_X, size_det_Y, size_det_Z))
grid = volume.grid.Grid((num_voxels_X, num_voxels_Y, num_voxels_Z,), (size_voxel_X, size_voxel_Y, size_voxel_Z))
vol = create_xcat(cfg, 'xcat_util', cfg.get('ORGAN', 'knee'), grid, offset=volume_offset, num_of_materials=cfg["NUM_OF_MATERIALS"], load_from_cache=True)
vol.AGWN(sigma=0.03)
vol.reduce_materials(cfg["NUM_OF_MATERIALS"])

# load spectrum array
src.create_120keV_spect()

# GT_atten_coeff = construct_atten_coeff_vol(vol, src._spectrum, os.path.join(cfg['MAIN_FILE_DIR'], cfg['MASS_ATTEN_DIR']))


code = projection_code.projection_code.ProjectionCode(num_sources, num_shots, num_projection_in_shot)
code.create_lin_space_code()
code.rand_perm_existing_code(code.get_code())

projections = projections.projections.Projections(active_scorer=2)  # We want to compare using direct transmission only

G4drive = g4driver.g4driver.G4Driver(cfg, detect, src, projections)

print('Running GT with phantom')
(I_GT, _), loaded_from_cache = G4drive.runG4(vol, code, grad=False, GT=True, build_phantom=True)
projections.set_I_GT(I_GT, built_phantom=True)
print('Running GT without phantom')
(I0_GT, _), loaded_from_cache0 = G4drive.runG4(vol, code, grad=False, GT=True, build_phantom=False)
projections.set_I_GT(I0_GT, built_phantom=False)

projections.calc_direct_ray_mask(eps=0.25)


if cfg['LOAD_INIT']:
    initial_volume = load_volume(cfg)
elif cfg['LOAD_CHECKPOINT']:
    initial_volume = load_volume(cfg, name='volume_checkpoint')
else:
    print('Running linear tomography')
    linear_solver = linear_tomography.sart_ring_solver.SARTSolver(src, detect, vol, code, cfg)
    reconstruct = linear_solver.inverse(np.transpose(I_GT[:, 2, :, :], (0, 2, 1)),
                                        np.transpose(I0_GT[:, 2, :, :], (0, 2, 1)),
                                        algorithm='SIRT3D_CUDA', iterations=100, initializer=0.5, load_from_cache=True,
                                        images_from_cache=loaded_from_cache and loaded_from_cache0)

    print('Volume initialization')
    initial_volume = create_vol_from_ct_recovery(cfg, reconstruct, grid, src._spectrum, offset=volume_offset,
                                                     load_from_cache=False, arch_mat=True)

    print('Saving volume')
    save_volume(get_experiment_path(cfg), cfg, initial_volume)

print('Creating logger')
run_name = generate_run_name(cfg)
slices2plot = 0
logger = Logger(os.path.join(get_experiment_path(cfg, data=False), 'logs'), run_name)
# logger.log_multiple_figures('GT', 0, I_GT[:, 2, :, :], title='GT with phantom', colorbar=True, dim_permute=(0, 2, 1))
# logger.log_multiple_figures('GT_vol', 0, vol.get_density()[:,:,slices2plot], title='density of GT volume', colorbar=True, dim_permute=(2, 0, 1))
# dens_per_el = np.zeros((len(cfg['ACTIVE_ELEMENTS']), num_voxels_X, num_voxels_Y))
# for index, el in enumerate([0, [5, 6, 7], [13, 18]]):
#     if type(el) is list:
#         dens_per_el[index, :, :] = vol.get_density()[:, :,slices2plot]*np.sum(vol._fractional_mass[el, :, :, slices2plot], axis=0)
#     else:
#         dens_per_el[index, :, :] = vol.get_density()[:, :, slices2plot] * vol._fractional_mass[el, :, :, slices2plot]
# logger.log_multiple_figures('GT dens per element', 0, dens_per_el, colorbar=True)
x0 = initial_volume.get_density_vec()

print('Optimization process begins')
oracle = ring_oracle_1st_order_helper(cfg, G4drive, initial_volume, vol, code, projections, logger)
dummy_vol = deepcopy(initial_volume)  # we want to know what are the bone voxels
constrain_func = adjust_func_helper(dummy_vol, cfg)
vol_optimization = ADAM(oracle, x0, nIter=cfg["OPTIMIZATION_ITER"], stepSize=cfg.get('STEP_SIZE',[0.005, 0.02, 0.005, 0.005]), beta1=0.7,
                        beta2=0.99, constrain_func=constrain_func, cfg=cfg, logger=logger, constrain=True)
#[0.003, 0.012, 0.003, 0.003]

final_vol = volume.toy_volume.ToyVolume(grid, active_elements, offset=volume_offset, arch_mat=False)
final_vol.load_density_vec(vol_optimization[:, -1], num_materials=cfg["NUM_OF_MATERIALS"])
I_final, _ = G4drive.runG4(final_vol, code, grad=False, GT=False, build_phantom=True)
logger.log_multiple_figures('I_final', 0, I_final[:, 2, :, :], title='I_final', colorbar=True, dim_permute=(0, 2, 1))
save_volume(get_experiment_path(cfg), cfg, final_vol, name='final_vol')
