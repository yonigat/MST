import copy
import sys


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

from utils.params import read_configs

from volume.xcat_volume_util import create_xcat


import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.use('TkAgg')

cfg = read_configs('../tests/GT_generate_new_3D_cfg.yaml')  # create a dictionary of all the paramters, it will also be sent to the G4 sim

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
# src.create_120keV_spect()
src.create_multi_delta([60])
# src.create_120keV_spect_trimed()
# src.create_60keV_delta()

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

images2show = np.linspace(0, len(code.get_code()), 90, dtype=int, endpoint=False)
first_in_code = code.get_code()[:,0] if len(code.get_code().shape) ==2 else code.get_code()
unpremed = np.argsort(first_in_code.flatten())

for ind, im_ind in enumerate(unpremed[images2show]):
    plt.figure()
    im2show = copy.copy(I_GT[im_ind, 2, :, :].T)
    norm_val = np.max(I_GT[im_ind, 2, :, :].T)
    im2show = im2show / norm_val
    im2show = np.power(im2show, 0.5)
    plt.imshow(im2show, cmap='bone_r')
    plt.axis('off')
    if not os.path.isdir(f'/home/yonatangat/Figures/{num_projection_in_shot}'):
        os.mkdir(f'/home/yonatangat/Figures/{num_projection_in_shot}')
    plt.savefig(f'/home/yonatangat/Figures/{num_projection_in_shot}/{ind}.png', dpi=500)
# plt.show()
# a=1