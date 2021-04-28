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
mpl.use('TkAgg')

cfg = read_configs('GT_generate_new_flat_cfg.yaml')  # create a dictionary of all the paramters, it will also be sent to the G4 sim

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
vol = ToyVolume(grid, active_elements, arch_mat=True)

radius = 30


def in_sphere(loc: tuple, center: tuple, radius: float = 2):
    # radius units in number of voxels
    norm = np.linalg.norm(np.array(loc)-np.array(center))
    return norm < radius

def in_cube(loc: tuple, center: tuple, radius: float = 2):
    # radius units in number of voxels
    is_in = True
    for i in range(len(loc)):
        is_in = is_in and center[i] - radius <= loc[i] <= center[i] + radius
    return is_in


bone_center = (15, 15, 15)
bone_radius = 10
soft_center = (40,45,30)
soft_radius = 12

for ix in np.arange(num_voxels_X):
    for iy in np.arange(num_voxels_Y):
        for iz in np.arange(num_voxels_Z):
            loc = (ix, iy, iz)
            vol.add_voxel_element(loc, 'O', 0.0012) # air voxels
            if in_sphere(loc, bone_center, bone_radius):
                vol.add_voxel_element(loc, 'O', 0.4)
                vol.add_voxel_element(loc, 'Ca', 0.2)
            if in_cube(loc, soft_center, soft_radius):
                vol.add_voxel_element(loc, 'H', 0.1)
                vol.add_voxel_element(loc, 'O', 0.3)
            # if (ix - num_voxels_X/2)**2 + (iy - num_voxels_Y/2)**2 < radius**2:
            #     vol.add_voxel_element((ix, iy, iz), 'O', 0.8)
# vol.blur_vol()

code = projection_code.projection_code.ProjectionCode(num_sources, num_shots, num_projection_in_shot)
code.create_lin_space_code()
code.rand_perm_existing_code(code.get_code())

# load spectrum array
# src.create_120keV_spect()
src.create_60keV_delta()
# src.create_multi_delta([15])
# src.create_120keV_spect_trimed()
GT_atten_coeff = construct_atten_coeff_vol(vol, src._spectrum, os.path.join(cfg['MAIN_FILE_DIR'], cfg['MASS_ATTEN_DIR']))
vol.reduce_materials(cfg["NUM_OF_MATERIALS"])


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
                                        algorithm='SIRT3D_CUDA', iterations=2500, initializer=0.0, load_from_cache=True,
                                        images_from_cache=loaded_from_cache and loaded_from_cache0)

    print('Volume initialization')
    initial_volume = create_vol_from_ct_recovery(cfg, reconstruct, grid, src._spectrum, offset=volume_offset,
                                                     load_from_cache=False, arch_mat=True)

    print('Saving volume')
    save_volume(get_experiment_path(cfg), cfg, initial_volume)




print('Creating logger')
run_name = generate_run_name(cfg)

logger = Logger(os.path.join(get_experiment_path(cfg, data=False), 'logs'), run_name)
# logger.log_multiple_figures('GT', 0, I_GT[:, 2, :, :], title='GT with phantom', colorbar=True, dim_permute=(0, 2, 1))
logger.log_multiple_figures('GT_vol', 0, vol.get_density()[:, :, cfg['SLICES_TO_PLOT']], title='density of GT volume', colorbar=True, dim_permute=(2, 0, 1))
dens_per_el = np.zeros((len(cfg['ACTIVE_ELEMENTS']), num_voxels_X, num_voxels_Y, len(cfg['SLICES_TO_PLOT'])))
for index, el in enumerate([0, 1, 2]):
    if type(el) is list:
        dens_per_el[index, :, :] = vol.get_density()[:, :, cfg['SLICES_TO_PLOT']]*np.sum(vol._fractional_mass[el, :, :, cfg['SLICES_TO_PLOT']], axis=0)
    else:
        dens_per_el[index, :, :] = vol.get_density()[:, :, cfg['SLICES_TO_PLOT']] * vol._fractional_mass[el, :, :, cfg['SLICES_TO_PLOT']].transpose(1,2,0)
logger.log_multiple_figures('GT dens per element', 0, dens_per_el, colorbar=True, dim_permute=(2, 0, 1))
x0 = initial_volume.get_density_vec()

print('Optimization process begins')
oracle = ring_oracle_1st_order_helper(cfg, G4drive, initial_volume, vol, code, projections, logger)
dummy_vol = deepcopy(initial_volume)  # we want to know what are the bone voxels
constrain_func = adjust_func_helper(dummy_vol, cfg)
vol_optimization = ADAM(oracle, x0, nIter=cfg["OPTIMIZATION_ITER"], stepSize=0.005, beta1=0.9, beta2=0.99, constrain_func=constrain_func, cfg=cfg, logger=logger, constrain=False)

final_vol = volume.toy_volume.ToyVolume(grid, active_elements, offset=volume_offset, arch_mat=False)
final_vol.load_density_vec(vol_optimization[:, -1], num_materials=cfg["NUM_OF_MATERIALS"])
I_final, _ = G4drive.runG4(final_vol, code, grad=False, GT=False, build_phantom=True)
logger.log_multiple_figures('I_final', 0, I_final[:, 2, :, :], title='I_final', colorbar=True, dim_permute=(0, 2, 1))
save_volume(get_experiment_path(cfg), cfg, final_vol, name='final_vol')
