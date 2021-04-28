def generate_illumination_matrix(cfg_file: str):
    import projections.projections as projection_class
    import detector.ring_detector
    import source.ring_source
    import detector.flat_panel_detector
    import volume.grid
    import projection_code.projection_code
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
    vol = create_xcat(cfg, 'xcat_util', 'knee', grid, offset=volume_offset, num_of_materials=cfg["NUM_OF_MATERIALS"],
                      load_from_cache=True)
    vol.blur_vol()
    GT_atten_coeff = construct_atten_coeff_vol(vol)
    vol.reduce_materials(cfg["NUM_OF_MATERIALS"])

    code = projection_code.projection_code.ProjectionCode(num_sources, num_shots, num_projection_in_shot)
    code.create_lin_space_code()
    code.rand_perm_existing_code(code.get_code())

    # load spectrum array
    src.create_120keV_spect()

    projections = projection_class.Projections(active_scorer=2)  # We want to compare using direct transmission only

    G4drive = g4driver.g4driver.G4Driver(cfg, detect, src, projections)

    print('Running GT with phantom')
    I_GT, _ = G4drive.runG4(vol, code, grad=False, GT=True, build_phantom=True)
    projections.set_I_GT(I_GT, built_phantom=True)
    print('Running GT without phantom')
    I0_GT, _ = G4drive.runG4(vol, code, grad=False, GT=True, build_phantom=False)
    projections.set_I_GT(I0_GT, built_phantom=False)

    illum_mat = np.zeros((I0_GT.shape[0], I0_GT.shape[2]*I0_GT.shape[3]))
    for sample_ind, sample in enumerate(I0_GT):
        illum_mat[code._projection_code[sample_ind][0], :] = np.hstack(sample[0])
    return illum_mat

