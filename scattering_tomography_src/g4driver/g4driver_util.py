import numpy as np
import os
import pickle
from utils.read_write_cache import save_data, get_experiment_path, read_data

def create_g4_params_dict(run_config: dict, grad=False, n_photons=10000, build_phantom=True, n_elements=6, GT=False):
    """
    Creates the parameters dict that is used by the G4 simulator.
    This method takes the run configuration and modifies it to G4 params_dict format.
    :param: run_config: The configuration of the entire run
    :param: grad: bool. whether to run with/without gradient
    :param: n_photons. Number of photons to generate
    :return: dict with G4 config
    """
    params_dict = dict()
    params_dict["i_MULTI"] = run_config['MULTI_THREADED']
    params_dict["i_NUM_OF_SCORERS"] = run_config['NUM_OF_SCORERS']
    params_dict["i_NUM_OF_BASE_SCORERS"] = 5
    params_dict["i_SINGLE_SCATTERING"] = run_config['SINGLE_SCATTERING']
    params_dict["i_NUM_OF_SOURCES"] = run_config['NUM_OF_SOURCES']
    params_dict["i_NUM_OF_SPECTRUM_BINS"] = run_config['NUM_OF_SPECTRUM_BIN']
    params_dict["i_NUM_OF_ELEMENTS"] = n_elements
    params_dict["i_NUM_OF_THREADS"] = run_config['NUM_OF_THREADS']
    params_dict["i_NUM_OF_BASE_MATERIALS"] = run_config['NUM_OF_MATERIALS']
    params_dict["i_NUM_OF_PHOTONS"] = n_photons
    params_dict["i_ROTATE_PHANTOM"] = run_config['ROTATE_PHANTOM']

    # verbose
    params_dict["i_VERBOSE_SCORING"] = 0
    params_dict["i_VERBOSE_PHYSICS_LIST"] = 0
    params_dict["i_VERBOSE_ANALYSIS_MANAGER"] = 0
    params_dict["i_VERBOSE_VIS"] = 0
    params_dict["i_VERBOSE_RUN"] = 1
    params_dict["i_VERBOSE_EVENT"] = 0
    params_dict["i_VERBOSE_TRACK"] = 0
    params_dict["i_VERBOSE_G4NAVIGATOR"] = 0
    params_dict["i_BIASING"] = run_config['BIASING']
    params_dict["i_SPLITTING_FACTOR_NP"] = 70
    params_dict["i_SPLITTING_FACTOR_NS"] = 70
    params_dict["f_DETECTOR_SPECIAL_CUTS"] = 0.0
    params_dict["f_PHANTOM_PRODUCTION_CUTS"] = 0.0
    params_dict["i_KILL_ELECTRONS"] = 1
    params_dict["i_RECORD_HIST"] = 0
    params_dict["i_PRINT_ELEMENTS_XS"] = 0

    # particle gun
    params_dict["f_PARTICLE_ENERGY"] = 60
    params_dict["f_MIN_THETA"] = run_config['MIN_THETA']
    params_dict["f_MAX_THETA"] = np.pi / run_config['MAX_THETA']
    params_dict["f_MIN_PHI"] = run_config['MIN_PHI']
    params_dict["f_MAX_PHI"] = np.pi / run_config['MAX_PHI']
    params_dict["i_CONE_BEAM"] = run_config['CONE_BEAM']
    params_dict["f_FILL_FACTOR"] = 0.7
    params_dict["f_THETA_CUT"] = 6
    params_dict["i_ALT_SOURCES"] = 1

    # projection code
    params_dict["i_NUM_SHOTS_GT"] = run_config['NUM_OF_SHOTS_GT']
    params_dict["i_NUM_SHOTS"] = run_config['NUM_OF_SHOTS_GT'] if GT else run_config['NUM_OF_SHOTS_OPTIMIZATION']
    params_dict["i_NUM_SHOTS_GT"] = run_config['NUM_OF_SHOTS_GT']
    params_dict["i_NUM_PROJECTION_IN_SHOT"] = run_config['NUM_PROJECTIONS_IN_SHOT']

    # geometry
    params_dict["i_BUILD_DETECTORS"] = 1
    params_dict["i_CALC_GRADIENT"] = 1 if grad else 0
    params_dict["i_BUILD_PHANTOM"] = 1 if build_phantom else 0
    params_dict["i_XCAT"] = 0
    params_dict["i_CT_PHANTOM"] = 1

    params_dict["i_USE_DICOM_INIT"] = 0
    params_dict["i_ROTATE_SEQUENTIALLY"] = 1 if GT else 0
    params_dict["f_WORLD_XY"] = run_config['WORLD_XY']
    params_dict["f_WORLD_Z"] = run_config['WORLD_Z']
    params_dict["f_DETECTOR_X"] = run_config['DETECTOR_X']
    params_dict["f_DETECTOR_Y"] = run_config['DETECTOR_Y']
    params_dict["f_DETECTOR_Z"] = run_config['DETECTOR_Z']
    params_dict["i_NUM_OF_DET_COLS"] = run_config['NUM_DETECTOR_COLS']
    params_dict["i_NUM_OF_DET_ROWS"] = run_config['NUM_DETECTOR_ROWS']
    params_dict["f_CENTER_TO_DET"] = run_config['CENTER_TO_DET']
    params_dict["f_SOURCE_TO_CENTER"] = run_config['SOURCE_TO_CENTER']
    params_dict["f_RADIUS"] = run_config['RADIUS']
    params_dict["f_OFFSET_U"] = run_config['OFFSET_U']
    params_dict["f_OFFSET_V"] = run_config['OFFSET_V']
    params_dict["f_SHIFT"] = run_config['SHIFT']
    params_dict["i_GT_MODE"] = 1
    params_dict["i_LIV_MODE"] = 0
    params_dict["f_PHANTOM_OFFSET_X"] = run_config['PHANTOM_OFFSET_X']
    params_dict["f_PHANTOM_OFFSET_Y"] = run_config['PHANTOM_OFFSET_Y']
    params_dict["f_PHANTOM_OFFSET_Z"] = run_config['PHANTOM_OFFSET_Z']

    # XCAT
    params_dict["i_NUM_OF_Z_SLICES"] = run_config['NUM_OF_Z_SLICES']
    params_dict["i_NUM_OF_VOXELS_X"] = run_config['NUM_OF_VOXELS_X']
    params_dict["i_NUM_OF_VOXELS_Y"] = run_config['NUM_OF_VOXELS_Y']
    params_dict["i_NUM_OF_VOXELS"] = run_config['NUM_OF_VOXELS']
    params_dict["i_NUM_OF_PIXELS_SLICE"] = run_config['NUM_OF_VOXELS_X'] * run_config['NUM_OF_VOXELS_Y']
    params_dict["f_VOXEL_HALF_X"] = run_config['VOXEL_HALF_X']
    params_dict["f_VOXEL_HALF_Y"] = run_config['VOXEL_HALF_Y']
    params_dict["f_VOXEL_HALF_Z"] = run_config['VOXEL_HALF_Z']

    # files
    params_dict["i_WRITE_TO_FILES"] = run_config['WRITE_TO_FILES']

    input_path = os.path.join(run_config['MAIN_FILE_DIR'], run_config['INPUT_DIR'])
    output_path = os.path.join(run_config['MAIN_FILE_DIR'], run_config['OUTPUT_DIR'])
    grad_path = os.path.join(run_config['MAIN_FILE_DIR'], run_config['GRADIENT_DIR'])
    error_path = os.path.join(run_config['MAIN_FILE_DIR'], run_config['ERROR_DIR'])

    params_dict["s_FILE_PROJECTION_CODE"] = os.path.join(input_path, run_config['FILE_PROJECTION_CODE'])
    params_dict["s_FILE_INTENS_CODE"] = os.path.join(input_path, run_config['FILE_INTENS_CODE'])
    params_dict["s_FILE_FIXED_SOURCE"] = os.path.join(input_path, run_config['FILE_FIXED_SOURCE'])
    params_dict["s_FILE_SOURCE_TO_DET"] = os.path.join(run_config['MAIN_FILE_DIR'], run_config['FILE_SOURCE_TO_DET'])
    params_dict["s_FILE_SOURCE_POS"] = os.path.join(run_config['MAIN_FILE_DIR'], run_config['FILE_SOURCE_POS'])
    params_dict["s_FILE_DET_POS"] = os.path.join(run_config['MAIN_FILE_DIR'], run_config['FILE_DET_POS'])
    params_dict["s_INPUT_DIR"] = input_path
    params_dict["s_OUTPUT_DIR"] = output_path
    params_dict["s_GRADIENT_DIR"] = grad_path
    params_dict["s_ERROR_DIR"] = error_path
    params_dict["s_FILE_SPECTRUM"] = os.path.join(input_path, run_config['FILE_SPECTRUM'])
    params_dict["s_FILE_MATERIALS"] = os.path.join(input_path, run_config['FILE_MATERIALS'])
    params_dict["s_ELEMENTS_FILE"] = os.path.join(input_path, run_config['FILE_ELEMENTS'])
    params_dict["s_A_FILE"] = os.path.join(input_path, run_config['FILE_A'])
    params_dict["s_Z_FILE"] = os.path.join(input_path, run_config['FILE_Z'])
    params_dict["s_FILE_MATERIALS_DICOM_BASIC"] = os.path.join(input_path, run_config['FILE_MATERIALS_DICOM_BASIC'])
    params_dict["s_FILE_VOXEL_TO_MATERIALS"] = os.path.join(input_path, run_config['FILE_VOXEL_TO_MATERIALS'])
    params_dict["s_FILE_VOXEL_TO_MATERIALS_ID"] = os.path.join(input_path, run_config['FILE_VOXEL_TO_MATERIALS_ID'])
    params_dict["s_FILE_VOXEL_TO_MATERIALS_TEST"] = os.path.join(input_path, run_config['FILE_VOXEL_TO_MATERIALS_TEST'])
    params_dict["s_FILE_VOXEL_TO_MATERIALS_Y"] = os.path.join(input_path, run_config['FILE_VOXEL_TO_MATERIALS_Y'])
    params_dict["s_FILE_VOXEL_TO_MATERIALS_DENS"] = os.path.join(input_path, run_config['FILE_VOXEL_TO_MATERIALS_DENS'])
    params_dict["s_FILE_XCAT_ID_TO_COMP"] = os.path.join(input_path, run_config['FILE_XCAT_ID_TO_COMP'])
    params_dict["s_FILE_XCAT_SLICE_PREFIX"] = os.path.join(input_path, run_config['FILE_XCAT_SLICE_PREFIX'])
    params_dict["s_FILE_SOURCES"] = os.path.join(input_path, run_config['FILE_SOURCES'])
    params_dict["s_FILE_SOURCES_ORIENTATION"] = os.path.join(input_path, run_config['FILE_SOURCES_ORIENTATION'])
    params_dict["s_FILE_DETECTORS"] = os.path.join(input_path, run_config['FILE_DETECTORS'])
    params_dict["s_FILE_DETECTORS_ORIENTATION"] = os.path.join(input_path, run_config['FILE_DETECTORS_ORIENTATION'])
    params_dict["s_FILE_PHANTOM_OFFSET_ORIENTATION"] = os.path.join(input_path, run_config['FILE_PHANTOM_OFFSET_ORIENTATION'])

    return params_dict


def write_g4_params_dict(param_dict: dict, path: str = '/home/yonatangat/Geant4-workspace/geant4.10.05/examples/basic/ring_sim_debug/'):
    file = os.path.join(path, 'dict.txt')
    text = open(file, 'w')
    for key in param_dict:
        text.write("%s %s\n" % (key, param_dict[key]))


def create_inputs_dir(path: str="../../run_inputs"):
    """

    :param path:
    :return:
    """
    if not os.path.isdir(path):
        os.makedirs(path)


def cache_g4(func):
    def wrapper_func(*args, **kwargs):
        # First check if we need to cache/ check for previous cache:
        if not kwargs['GT']:
            I_res = func(*args, **kwargs)
            return I_res
        else:
            loaded_from_cache = False
            # Check if there is previous caching results
            phantom_state = 'I_GT' if kwargs['build_phantom'] else 'I0_GT'

            # args[0] is the self (G4driver)
            run_cfg = args[0].get_cfg()
            exp_dir_path = get_experiment_path(run_cfg, GT_data=True)
            meta_data = [args, kwargs]

            if not (os.path.isdir(exp_dir_path) and
                    os.path.exists(os.path.join(exp_dir_path, 'run_cfg_' + phantom_state)) and
                    os.path.exists(os.path.join(exp_dir_path, 'meta_data_' + phantom_state)) and
                    os.path.exists(os.path.join(exp_dir_path, phantom_state))):
                # Need to cache
                I_res = func(*args, **kwargs)
                save_results_to_dir(exp_dir_path, meta_data, run_cfg, I_res, phantom_state)
            else:
                # Old results exists, compare fun configurations, and arguments
                if check_valid_results(exp_dir_path, meta_data, run_cfg, phantom_state):
                    I_res = read_data(exp_dir_path, phantom_state)
                    loaded_from_cache = True
                else:
                    # Need to cache
                    I_res = func(*args, **kwargs)
                    save_results_to_dir(exp_dir_path, meta_data, run_cfg, I_res, phantom_state)

            return I_res, loaded_from_cache
    return wrapper_func


def save_results_to_dir(path, meta_data, run_cfg, I_res, phantom_state):
    if not os.path.isdir(path):
        os.makedirs(path)
    save_data(path, meta_data, 'meta_data_' + str(phantom_state))
    save_data(path, run_cfg, 'run_cfg_' + str(phantom_state))
    save_data(path, I_res, str(phantom_state))


def read_previous_results(path,phantom_state):
    if not os.path.isdir(path):
        os.makedirs(path)
    prev_meta_data = read_data(path, 'meta_data_' + str(phantom_state))
    prev_cfg = read_data(path, 'run_cfg_' + str(phantom_state))

    return prev_cfg, prev_meta_data


def check_valid_results(path, meta_data, run_cfg, phantom_state):
    fields_to_compare = ['NUM_OF_SCORERS', 'SINGLE_SCATTERING', 'NUM_OF_SOURCES', 'NUM_OF_SPECTRUM_BIN', 'NUM_OF_PHOTONS_GT', 'ROTATE_PHANTOM',
                         'MIN_THETA', 'MAX_THETA', 'MIN_PHI', 'MAX_PHI', 'CONE_BEAM', 'NUM_OF_SHOTS_GT', 'NUM_PROJECTIONS_IN_SHOT', 'WORLD_XY','WORLD_Z', 'DETECTOR_X',
                         'DETECTOR_Y', 'DETECTOR_Z' ,'NUM_DETECTOR_COLS' ,'NUM_DETECTOR_ROWS' ,'CENTER_TO_DET' , 'SOURCE_TO_CENTER','OFFSET_U' ,'OFFSET_V' ,'SHIFT',
                         'PHANTOM_OFFSET_X','PHANTOM_OFFSET_Y','PHANTOM_OFFSET_Z','NUM_OF_Z_SLICES','NUM_OF_VOXELS_X','NUM_OF_VOXELS_Y','VOXEL_HALF_X','VOXEL_HALF_Y','VOXEL_HALF_Z','NUM_OF_MATERIALS']

    prev_cfg, prev_meta_data = read_previous_results(path, phantom_state)

    # Compare the configs:
    equal_cfgs = True
    unequal_fields = {}
    for field in fields_to_compare:
        if run_cfg[field] != prev_cfg[field]:
            unequal_fields[field] = (run_cfg[field], prev_cfg[field])
        equal_cfgs = equal_cfgs and run_cfg[field] == prev_cfg[field]

    equal_args_meta_data = meta_data[0][1:] == prev_meta_data[0][1:]
    equal_kwargs_meta_data = meta_data[1] == prev_meta_data[1]

    if len(unequal_fields) > 0:
        for field, values in unequal_fields.items():
            print(f'{field} - new: {values[0]}, old: {values[1]}')
    if not equal_args_meta_data:
        print('***problem in args metadata')
    if not equal_kwargs_meta_data:
        print('***problem in kwargs metadata')
    return equal_cfgs and equal_args_meta_data and equal_kwargs_meta_data


