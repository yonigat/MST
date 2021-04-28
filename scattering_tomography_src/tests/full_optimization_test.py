import utils.params
import utils.write_mac_files
import detector.flat_panel_detector
import projections.projections
# import detector.ring_detector
import source.ring_source
import volume.toy_volume
import volume.grid
import projection_code.projection_code
import numpy as np
from optimization.oracle import oracle_1st_order_helper
from optimization.ADAM import ADAM
import linear_tomography.flat_linear_tomo_solver
import matplotlib.pyplot as plt
import g4driver.g4driver



cfg = utils.params.read_configs('test_optimization_config.yaml')  # create a dictionary of all the paramters, it will also be sent to the G4 sim


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

src = source.ring_source.RingSource(radius_src, num_sources, (angle_theta, angle_delta))
detect = detector.flat_panel_detector.FlatPanelDetector(radius_det, num_det_cols, num_det_rows, (size_det_X, size_det_Y, size_det_Z))
grid = volume.grid.Grid((num_voxels_X, num_voxels_Y, num_voxels_Z,), (size_voxel_X, size_voxel_Y, size_voxel_Z))
vol = volume.toy_volume.ToyVolume(grid, active_elements, offset=volume_offset)
code = projection_code.projection_code.ProjectionCode(num_sources, num_shots, num_projection_in_shot)

code.create_one_source_code()

# load materials
# materials = np.loadtxt(cfg["INPUT_DIR"]+'materials.txt')
# vol._materials = materials
# vol._num_materials = materials.shape[1]
#
# # load id and density
# for indZ in np.arange(num_voxels_Z):
#     # vol._id[:, :, indZ] = np.zeros((1, 1))
#     # vol._density[:, :, indZ] = np.ones((1, 1)) * 0.01
#     vol._id[:, :, indZ] = np.loadtxt(cfg["INPUT_DIR"]+'Y/id'+str(indZ)+'.dat').T
#     vol._density[:, :, indZ] = np.loadtxt(cfg["INPUT_DIR"]+'Y/dens' + str(indZ) + '.dat').T

vol.create_simple_phantom(cfg["NUM_OF_MATERIALS"])

# load spectrum array
src.create_120keV_spect()

projections = projections.projections.Projections(num_photons_forward=cfg["NUM_OF_PHOTONS_FORWARD"])

G4drive = g4driver.g4driver.G4Driver(cfg, detect, src, projections)

I_GT = G4drive.runG4(vol, code, grad=False, GT=True, build_phantom=True)
projections.set_I_GT(I_GT, built_phantom=True)

optimization_code = projection_code.projection_code.ProjectionCode(num_sources, num_shots_opt, num_projection_in_shot)
optimization_code.create_one_source_code()

# create initialization
x0 = vol.get_density_vec()
x0 = 1.5 * x0

oracle = oracle_1st_order_helper(cfg, G4drive, vol, optimization_code, projections)
final_vol = ADAM(oracle, x0, nIter=cfg["OPTIMIZATION_ITER"])

plt.figure()
plt.plot(final_vol[0, :], label=r"$voxel_{0}$ $el_{0}$")
plt.plot(final_vol[1, :], label=r"$voxel_{1}$ $el_{0}$")
plt.xlabel("iteration")
plt.ylabel(r"voxel density [$\frac{g}{cm^{3}}$]")
plt.title(r"$Element_{0}$")
plt.legend(loc='best')

plt.figure()
plt.plot(final_vol[2, :], label=r"$voxel_{0}$ $el_{1}$")
plt.plot(final_vol[3, :], label=r"$voxel_{1}$ $el_{1}$")
plt.xlabel("iteration")
plt.ylabel(r"voxel density [$\frac{g}{cm^{3}}$]")
plt.title(r"$Element_{1}$")
plt.legend(loc='best')

plt.figure()
plt.plot(projections._loss)
plt.xlabel(r"iteration")
plt.ylabel(r"loss")
plt.title(r"Loss vs. Iterations")

plt.show()