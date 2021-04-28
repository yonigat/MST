import projections.projections
import detector.ring_detector
import source.ring_source
import detector.flat_panel_detector
import volume.grid
from projection_code.projection_code import ProjectionCode
from nontrivial_projection_code import NonTrivialProjectionCode
import numpy as np
import g4driver.g4driver
from illum_code_optim.generate_src2pix_map import generate_illum_mat, generate_full_code
import illum_code_optim.illum_loss
from utils.params import read_configs
from volume.xcat_volume_util import create_xcat

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
vol = create_xcat(cfg, 'xcat_util', 'knee', grid, offset=volume_offset, num_of_materials=cfg["NUM_OF_MATERIALS"], load_from_cache=True)
vol.blur_vol()
vol.reduce_materials(cfg["NUM_OF_MATERIALS"])

# code = ProjectionCode(num_sources, num_shots, num_projection_in_shot)
# code.create_lin_space_code()
# code.rand_perm_existing_code(code.get_code())

opt_code = NonTrivialProjectionCode(num_sources, num_shots, num_projection_in_shot)
W = np.array([1.0, 0.3, 0.6, 0.3]).reshape(1,4)
full_code = generate_full_code(W, num_shots)
opt_code.load_full_multiplex_code(full_code, cfg)


# load spectrum array
src.create_120keV_spect()

projections = projections.projections.Projections(active_scorer=2)  # We want to compare using direct transmission only

G4drive = g4driver.g4driver.G4Driver(cfg, detect, src, projections)

print('Running GT with phantom')
I_GT, _ = G4drive.runG4(vol, opt_code, grad=False, GT=True, build_phantom=True)
projections.set_I_GT(I_GT, built_phantom=True)
print('Running GT without phantom')
I0_GT, _ = G4drive.runG4(vol, opt_code, grad=False, GT=True, build_phantom=False)
projections.set_I_GT(I0_GT, built_phantom=False)

# illum_mat = generate_illum_mat(src, detect)
# S = np.zeros((num_sources,))
# # for idx in range(len(S)):
# #     S[idx] = random.choice([0.2, 0.8])
# S[code.get_code()[0,:]] = 1.0
# codes, loss = illum_code_optim.illum_loss.optimize_illum(S, illum_mat, C=9.0, eta=1e-3, lam=5)
# opt_code = NonTrivialProjectionCode(num_sources, num_shots, num_projection_in_shot)
# full_code = generate_full_code(codes[-1], num_shots)
# opt_code.load_full_multiplex_code(full_code, cfg)
# G4drive = g4driver.g4driver.G4Driver(cfg, detect, src, projections)
#
# print('Running GT with phantom')
# I_GT, _ = G4drive.runG4(vol, opt_code, grad=False, GT=True, build_phantom=True)
# projections.set_I_GT(I_GT, built_phantom=True)
# print('Running GT without phantom')
# I0_GT, _ = G4drive.runG4(vol, opt_code, grad=False, GT=True, build_phantom=False)
# projections.set_I_GT(I0_GT, built_phantom=False)
# a=1
