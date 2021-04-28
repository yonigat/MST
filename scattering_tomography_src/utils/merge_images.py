from nontrivial_projection_code.nontrivial_projection_code import NonTrivialProjectionCode


def merge_images(cfg_file: str, multi_GT_dir: str, num_sources_multiplexed: int, non_triv: bool = False):
    import projections.projections as projection_class
    import detector.ring_detector
    import source.ring_source
    import detector.flat_panel_detector
    import volume.grid
    import projection_code.projection_code
    # from nontrivial_projection_code import NonTrivialProjectionCode
    # from nontrivial_projection_code.nontrivial_projection_code import NonTrivialProjectionCode
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

    # multi_GT_dir = 'bla'
    cfg = read_configs(cfg_file)  # create a dictionary of all the paramters, it will also be sent to the G4 sim
    single_GT_path = get_experiment_path(cfg, GT_data=True)
    multi_GT_path = os.path.join(single_GT_path, os.pardir, multi_GT_dir)

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
    num_GT_photons = cfg["NUM_OF_PHOTONS_GT"]
    num_overall_photons = num_GT_photons * num_shots

    num_projection_in_shot = cfg["NUM_PROJECTIONS_IN_SHOT"]


    src = source.ring_source.RingSource(radius_src, num_sources, (angle_theta, angle_delta))
    detect = detector.ring_detector.RingDetector(radius_det, num_det_cols, num_det_rows,
                                                 (size_det_X, size_det_Y, size_det_Z))
    grid = volume.grid.Grid((num_voxels_X, num_voxels_Y, num_voxels_Z,), (size_voxel_X, size_voxel_Y, size_voxel_Z))
    vol = create_xcat(cfg, 'xcat_util', cfg.get('ORGAN', 'knee'), grid, offset=volume_offset, num_of_materials=cfg["NUM_OF_MATERIALS"],
                      load_from_cache=True)

    vol.AGWN(sigma=0.03)
    vol.reduce_materials(cfg["NUM_OF_MATERIALS"])

    # load spectrum array
    src.create_120keV_spect()

    code = projection_code.projection_code.ProjectionCode(num_sources, num_shots, num_projection_in_shot)
    code.create_lin_space_code()
    code.rand_perm_existing_code(code.get_code())

    projections = projection_class.Projections(active_scorer=2)  # We want to compare using direct transmission only

    G4drive = g4driver.g4driver.G4Driver(cfg, detect, src, projections)

    print('Running GT with phantom')
    (I_GT, _), loaded_from_cache = G4drive.runG4(vol, code, grad=False, GT=True, build_phantom=True)

    projections.set_I_GT(I_GT, built_phantom=True)
    print('Running GT without phantom')
    (I0_GT, _), loaded_from_cache = G4drive.runG4(vol, code, grad=False, GT=True, build_phantom=False)

    projections.set_I_GT(I0_GT, built_phantom=False)
    projections.calc_direct_ray_mask(eps=0.25)
    multiplexed_code = projection_code.projection_code.ProjectionCode(num_sources,
                                                                      np.int(np.ceil(num_sources / num_sources_multiplexed)),
                                                                      num_sources_multiplexed)
    multiplexed_code.create_lin_space_code()
    multiplexed_code.rand_perm_existing_code(multiplexed_code.get_code())

    if non_triv:
        illum_mat = generate_illum_mat(src, detect)
        S = np.zeros((num_sources,))

        S[multiplexed_code.get_code()[0, :]] = 1.0

        codes, loss = optimize_illum(S, illum_mat, C=num_sources_multiplexed, eta=1e-3, lam=20000.0)
        multiplexed_code = NonTrivialProjectionCode(num_sources, num_shots // num_sources_multiplexed,
                                                    num_projection_in_shot)
        full_code = generate_full_code(codes[-1], num_shots // num_sources_multiplexed)
        multiplexed_code.load_full_multiplex_code(full_code, cfg)

    num_multi_shots = int(np.ceil(num_sources / num_sources_multiplexed))
    I_multi_GT = np.zeros((num_multi_shots, I_GT.shape[1], I_GT.shape[2], I_GT.shape[3]))
    I0_multi_GT = np.zeros((num_multi_shots, I_GT.shape[1], I_GT.shape[2], I_GT.shape[3]))

    multi_cfg = cfg
    multi_cfg["NUM_OF_SHOTS_GT"] = num_multi_shots
    multi_cfg["NUM_PROJECTIONS_IN_SHOT"] = multiplexed_code._num_projection_in_shot

    for shot_ind, shot_code in enumerate(multiplexed_code.get_code()):
        sorter = np.argsort(code.get_code(), axis=0).reshape((-1,))
        idx_in_single_code = sorter[np.searchsorted(code.get_code().reshape((-1,)), shot_code, sorter=sorter)]
        intens = np.ones((len(idx_in_single_code), 1, 1, 1))
        if non_triv:
            intens = multiplexed_code._intens_code[shot_ind, :].reshape((len(idx_in_single_code), 1, 1, 1))
        I_multi_GT[shot_ind, :, :, :] = np.sum(I_GT[idx_in_single_code, :, :, :] * intens,
                                               axis=0)  # num_sources_multiplexed
        I0_multi_GT[shot_ind, :, :, :] = np.sum(I0_GT[idx_in_single_code, :, :, :] * intens,
                                                axis=0)  # num_sources_multiplexed

    multi_projections = projection_class.Projections(
        active_scorer=2)  # We want to compare using direct transmission only

    multi_G4drive = g4driver.g4driver.G4Driver(multi_cfg, detect, src, multi_projections)

    multi_projections.set_I_GT(I_multi_GT, built_phantom=True)
    multi_projections.set_I_GT(I0_multi_GT, built_phantom=False)

    I_multi_meta_data = ((multi_GT_dir, vol, multiplexed_code), {'grad': False, 'GT': True, 'build_phantom': True})
    save_results_to_dir(multi_GT_path, I_multi_meta_data, multi_cfg, (I_multi_GT, _), 'I_GT')
    I0_multi_meta_data = ((multi_GT_dir, vol, multiplexed_code), {'grad': False, 'GT': True, 'build_phantom': False})
    save_results_to_dir(multi_GT_path, I0_multi_meta_data, multi_cfg, (I0_multi_GT, _), 'I0_GT')
    if non_triv:
        np.savetxt(os.path.join(multi_GT_path, 'full_opt_code.txt'), full_code, fmt='%1.5f')