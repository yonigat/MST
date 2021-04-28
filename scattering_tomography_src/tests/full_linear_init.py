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
from optimization.oracle import ring_oracle_1st_order_helper, constrain_func_helper
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


cfg = read_configs('GT_generate_new_cfg.yaml')  # create a dictionary of all the paramters, it will also be sent to the G4 sim

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
# vol.blur_vol()
vol.reduce_materials(cfg["NUM_OF_MATERIALS"])
# load spectrum array
src.create_120keV_spect()

GT_atten_coeff = construct_atten_coeff_vol(vol, src._spectrum, os.path.join(cfg['MAIN_FILE_DIR'], cfg['MASS_ATTEN_DIR']), elements=cfg['ACTIVE_ELEMENTS'])


code = projection_code.projection_code.ProjectionCode(num_sources, num_shots, num_projection_in_shot)
code.create_lin_space_code()
code.rand_perm_existing_code(code.get_code())

projections = projections.projections.Projections(active_scorer=2)  # We want to compare using direct transmission only

G4drive = g4driver.g4driver.G4Driver(cfg, detect, src, projections)

print('Running GT with phantom')
I_GT, _ = G4drive.runG4(vol, code, grad=False, GT=True, build_phantom=True)
projections.set_I_GT(I_GT, built_phantom=True)
print('Running GT without phantom')
I0_GT, _ = G4drive.runG4(vol, code, grad=False, GT=True, build_phantom=False)
projections.set_I_GT(I0_GT, built_phantom=False)

projections.calc_direct_ray_mask(eps=0.25)

cfg['NUM_OF_MATERIALS'] = 5

if cfg['LOAD_INIT']:
    initial_volume = load_volume(cfg)
elif cfg['LOAD_CHECKPOINT']:
    initial_volume = load_volume(cfg, name='volume_checkpoint')
else:
    print('Running linear tomography')
    linear_solver = linear_tomography.sart_ring_solver.SARTSolver(src, detect, vol, code)
    # ann_proj = linear_solver.create(GT_atten_coeff)
    # reconstruct = linear_solver.run_linear_tomography(np.transpose(I_GT[:, 2, :, :], (0, 2, 1)), np.transpose(I0_GT[:, 2, :, :], (0, 2, 1)),
    #                                                   algorithm='FDK_CUDA', iterations=200)

    reconstruct = linear_solver.inverse(np.transpose(I_GT[:, 0, :, :], (0, 2, 1)), np.transpose(I0_GT[:, 0, :, :], (0, 2, 1)),
                                        algorithm='SIRT3D_CUDA', iterations=10000, initializer=10)
    # reconstruct = GT_atten_coeff
    print('Volume initialization')
    initial_volume = create_vol_from_ct_recovery(cfg, reconstruct, grid, src._spectrum, offset=volume_offset,
                                                     load_from_cache=False, arch_mat=True)
    GT = np.hstack(GT_atten_coeff[:, :, 0])
    rec = np.hstack(reconstruct[:, :, 0])
    err = (rec - GT)

    fig = plt.figure()
    ax = fig.add_subplot(111)

    ax.scatter(GT, rec)
    ax.plot(GT, GT, 'r--')
    ax.set_xlabel('Ground Truth')
    ax.set_ylabel('Linear Initialization')

    plt.title('Attenuation coefficients scatter plot')

    fig = plt.figure()

    clim = [0, 0]
    clim[0] = np.minimum(np.min(GT_atten_coeff[:, :, 0]), np.min(reconstruct[:, :, 0]))
    clim[1] = np.maximum(np.max(GT_atten_coeff[:, :, 0]), np.max(reconstruct[:, :, 0]))

    plt.subplot(1, 2, 1)
    im = plt.imshow(GT_atten_coeff[:, :, 0], clim=clim)
    plt.colorbar()
    plt.title('Ground truth atten. coeff. [1/cm]')
    plt.subplot(1, 2, 2)
    im2 = plt.imshow(reconstruct[:, :, 0], clim=clim)
    plt.colorbar()
    plt.title('Linear tomography atten. coeff. [1/cm]')

    # plt.subplot(1, 3, 3)
    # im2 = plt.imshow((GT_atten_coeff[:, :, 0] - recovery[:, :, 0]))
    #
    # plt.title('GT-initialization [1/cm]')
    # plt.colorbar()

    gt_frac_mass = vol._fractional_mass[[0, 5, 6, 7, 13, 18], :, :, :]
    init_frac_mass = initial_volume._fractional_mass
    fig = plt.figure()
    elements = {0: 'H', 1: 'C', 2: 'N', 3: 'O', 4: 'P', 5: 'Ca'}
    plt.suptitle('Comparing elements density in ground truth and initialization')
    for el in np.arange(3):
        if el == 0:
            el_GT_frac = gt_frac_mass[el, :, :, 0] * vol.get_density()[:, :, 0]
            el_init_frac = init_frac_mass[el, :, :, 0] * initial_volume.get_density()[:, :, 0]
            el_name = 'H'
        elif el == 1:
            el_GT_frac = np.sum(gt_frac_mass[[1,2,3], :, :, 0], axis=0) * vol.get_density()[:, :, 0]
            el_init_frac = np.sum(init_frac_mass[[1,2,3], :, :, 0], axis=0) * initial_volume.get_density()[:, :, 0]
            el_name = 'O'
        elif el == 2:
            el_GT_frac = np.sum(gt_frac_mass[[4, 5], :, :, 0], axis=0) * vol.get_density()[:, :, 0]
            el_init_frac = np.sum(init_frac_mass[[4, 5], :, :, 0], axis=0) * initial_volume.get_density()[:, :, 0]
            el_name = 'Ca'
        clim = [0, 0]
        clim[0] = np.minimum(np.min(el_GT_frac), np.min(el_init_frac))
        clim[1] = np.maximum(np.max(el_GT_frac), np.max(el_init_frac))

        plt.subplot(3, 2, 2 * el+1)
        plt.imshow(el_GT_frac, clim=clim)
        plt.colorbar()
        plt.axis('off')
        plt.title(f'Ground truth {elements[el]} dens [g/cm^3]')

        plt.subplot(3, 2, 2 + 2 * el)
        plt.imshow(el_init_frac, clim=clim)
        plt.colorbar()
        plt.axis('off')
        plt.title(f'Initialization {elements[el]} dens [g/cm^3]')

        # plt.subplot(len(elements), 3, 3 + 3 * el)
        # plt.imshow((el_GT_frac - el_init_frac))
        # plt.axis('off')
        # plt.title(f'GT-initialization for {elements[el]} [g/cm^3]')
        # plt.colorbar()

    plt.show()
    a=1

