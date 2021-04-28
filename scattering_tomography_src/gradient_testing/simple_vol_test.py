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
from linear_tomography.ct_to_vol_util import create_vol_simple_table
from optimization.ADAM import ADAM
from optimization.vanilla_SGD import VSGD
from optimization.oracle import ring_oracle_1st_order_helper
from utils.params import read_configs
from utils.plotting_util import calc_epsilon, calc_delta
from volume.abstract_volume import load_volume, save_volume
from volume.xcat_volume_util import create_xcat
from utils.generate_run_name import generate_run_name
from utils.tf_logger import Logger
from utils.read_write_cache import get_experiment_path
import matplotlib
import matplotlib.pyplot as plt

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
vol = volume.toy_volume.ToyVolume(grid, active_elements, offset=volume_offset)
for i in range(num_voxels_X):
    for j in range(num_voxels_Y):
        for k in range(num_voxels_Z):
            vol.add_voxel_element((i, j, k), 'H', 0.4)
            vol.add_voxel_element((i, j, k), 'O', 0.25)
            vol.add_voxel_element((i, j, k), 'Ca', 0.1)

code = projection_code.projection_code.ProjectionCode(num_sources, num_shots, num_projection_in_shot)
code.create_lin_space_code()
code.rand_perm_existing_code(code.get_code())

# load spectrum array
src.create_120keV_spect()

projections = projections.projections.Projections(active_scorer=2)  # We want to compare using direct transmission only

G4drive = g4driver.g4driver.G4Driver(cfg, detect, src, projections)

print('Running GT with phantom')
I_GT, _ = G4drive.runG4(vol, code, grad=False, GT=True, build_phantom=True)
projections.set_I_GT(I_GT, built_phantom=True)
print('Running GT without phantom')
I0_GT, _ = G4drive.runG4(vol, code, grad=False, GT=True, build_phantom=False)
projections.set_I_GT(I0_GT, built_phantom=False)

# projections.calc_cone_only_mask()

if cfg['LOAD_INIT']:
    initial_volume = load_volume(cfg)
elif cfg['LOAD_CHECKPOINT']:
    initial_volume = load_volume(cfg, name='volume_checkpoint')
else:
    initial_volume = volume.toy_volume.ToyVolume(grid, active_elements)
    # noise_level = np.reshape(np.abs(np.random.randn(vol._fractional_mass.size)*0.3 + 1.0), vol._fractional_mass.shape)
    for ix in range(num_voxels_X):
        for iy in range(num_voxels_Y):
            for iz in range(num_voxels_Z):
                for el in active_elements:
                    el_vol_ind = vol._active_elements[el][0]
                    el_dens = vol.get_density()[ix, iy, iz] * vol._fractional_mass[el_vol_ind,ix,iy, iz]
                    el_dens = el_dens# * 2.5 # * np.abs(np.random.randn()*0.3 + 1.0)
                    initial_volume.add_voxel_element((ix, iy, iz), el, el_dens) #* noise_level[el_vol_ind,ix,iy, iz])

print('Creating logger')
run_name = generate_run_name(cfg)
logger = Logger(os.path.join(get_experiment_path(cfg, data=False), 'logs'), run_name)
logger.log_multiple_figures('GT', 0, I_GT[:, 2, :, :], title='GT with phantom', colorbar=True, dim_permute=(0, 2, 1))
logger.log_multiple_figures('GT_vol', 0, vol.get_density(), title='density of GT volume', colorbar=True, dim_permute=(2, 0, 1))

x0 = initial_volume.get_density_vec()

print('Optimization process begins')
oracle = ring_oracle_1st_order_helper(cfg, G4drive, initial_volume, vol, code, projections, logger)
vol_optimization = ADAM(oracle, x0, nIter=cfg["OPTIMIZATION_ITER"], stepSize=0.3)

final_vol = volume.toy_volume.ToyVolume(grid, active_elements, offset=volume_offset)
final_vol.load_density_vec(vol_optimization[:, -1], num_materials=cfg["NUM_OF_MATERIALS"])
I_final, _ = G4drive.runG4(final_vol, code, grad=False, GT=False, build_phantom=True)
I0_final, _ = G4drive.runG4(final_vol, code, grad=False, GT=False, build_phantom=False)

# logging data for final volume
sampled_shots = range(cfg["NUM_OF_SHOTS_OPTIMIZATION"])
loss = projections.calc_loss(angles=sampled_shots, I=I_final, I0=I0_final)
eps = calc_epsilon(vol, final_vol)
delta = calc_delta(vol, final_vol)
logger.log_scalar('loss', loss, projections.itr)
logger.log_scalar('epsilon', eps, projections.itr)
logger.log_scalar('delta', delta, projections.itr)

if loss <= logger._optimal_loss:
    save_volume(get_experiment_path(cfg), cfg, vol, name='volume_checkpoint')

logger.log_multiple_figures('I_final', 0, I_final[:, 2, :, :], title='I_final', colorbar=True, dim_permute=(0, 2, 1))

save_volume(get_experiment_path(cfg), cfg, final_vol, name='final_vol')




