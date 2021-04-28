from copy import deepcopy

import projections.projections
import detector.ring_detector
import source.ring_source
import detector.flat_panel_detector
import volume.grid
import projection_code.projection_code
import numpy as np
import pickle
import os
import g4driver.g4driver
import linear_tomography.sart_ring_solver
import linear_tomography.create_and_reconstruct
from linear_tomography.ct_to_vol_util import create_vol_simple_table, create_vol_from_ct_recovery, \
    construct_atten_coeff_vol
from utils.plotting_util import calc_delta, calc_epsilon
from utils.generate_run_name import generate_run_name
from utils.identify_outliers import non_bone_voxel_filter, cond_median_filter
from utils.params import read_configs
from utils.tf_logger import Logger
from volume.abstract_volume import load_volume, save_volume
from volume.toy_volume import ToyVolume
from volume.xcat_volume_util import create_xcat
from utils.read_write_cache import get_experiment_path

import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.use('TkAgg')

cfg = read_configs(
    'projections_compare.yaml')  # create a dictionary of all the paramters, it will also be sent to the G4 sim

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
num_shots_opt = cfg["NUM_OF_SHOTS_OPTIMIZATION"]
num_shots = cfg["NUM_OF_SHOTS_GT"]

num_projection_in_shot = cfg["NUM_PROJECTIONS_IN_SHOT"]
energy = 50  # Kev

src = source.ring_source.RingSource(radius_src, num_sources, (angle_theta, angle_delta))
detect = detector.ring_detector.RingDetector(radius_det, num_det_cols, num_det_rows,
                                             (size_det_X, size_det_Y, size_det_Z))
grid = volume.grid.Grid((num_voxels_X, num_voxels_Y, num_voxels_Z,), (size_voxel_X, size_voxel_Y, size_voxel_Z))
vol = ToyVolume(grid, active_elements, arch_mat=True)

vol = create_xcat(cfg, 'xcat_util', 'knee', grid, offset=volume_offset, num_of_materials=cfg["NUM_OF_MATERIALS"],
                  load_from_cache=True)
# vol.blur_vol()
vol.AGWN(sigma=0.03)
vol.reduce_materials(cfg["NUM_OF_MATERIALS"])
# load spectrum array
# src.create_120keV_spect()
# src.create_120keV_spect_trimed()
# src.create_60keV_delta()
src.create_multi_delta([60])
GT_atten_file = os.path.join(get_experiment_path(cfg, GT_data=True), 'GT_vol_atten')

if os.path.exists(GT_atten_file) and True:
    with open(GT_atten_file, 'rb') as f:
        GT_atten_coeff = pickle.load(f)
else:
    GT_atten_coeff = construct_atten_coeff_vol(vol, src._spectrum,
                                               os.path.join(cfg['MAIN_FILE_DIR'], cfg['MASS_ATTEN_DIR']))
    with open(GT_atten_file, 'wb') as f:
        pickle.dump(GT_atten_coeff, f, protocol=4)


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

print('Running linear tomography')
linear_solver = linear_tomography.sart_ring_solver.SARTSolver(src, detect, vol, code, cfg)


reconstruct_with_projection = linear_solver.inverse(np.transpose(I_GT[:, 2, :, :], (0, 2, 1)),
                                    np.transpose(I0_GT[:, 2, :, :], (0, 2, 1)),
                                    algorithm='SIRT3D_CUDA', iterations=100000, initializer=0.5, load_from_cache=True,
                                    images_from_cache=loaded_from_cache and loaded_from_cache0)

reconstruct_without_projection = linear_solver.inverse(np.transpose(I_GT[:, 2, :, :], (0, 2, 1)),
                                    np.transpose(I0_GT[:, 2, :, :], (0, 2, 1)),
                                    algorithm='SIRT3D_CUDA', iterations=15000, initializer=0.0, load_from_cache=False,
                                    images_from_cache=loaded_from_cache and loaded_from_cache0, project=False)

plt.figure()
plt.subplot(1,3,1)
plt.imshow(reconstruct_without_projection[:,:,0])
plt.title('without projection')
plt.subplot(1,3,2)
plt.imshow(reconstruct_with_projection[:,:,0])
plt.title('with projection')
plt.subplot(1,3,3)
plt.imshow(reconstruct_with_projection[:,:,0]-reconstruct_without_projection[:,:,0])
plt.title('difference')
plt.show()
a=1

