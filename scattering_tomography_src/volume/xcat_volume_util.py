from utils.read_write_cache import read_data, save_data, get_experiment_path
import os
from volume.xcat_volume import XCATVolume
import numpy as np


def cache_xcat(func):
    def wrapper_func(*args, **kwargs):
        # First check if we need to cache/ check for previous cache:
        exp_dir_path = get_experiment_path(args[0], GT_data=True)
        meta_data = [args, kwargs]
        if not kwargs['load_from_cache']:
            res = func(*args, **kwargs)
            save_results_to_dir(exp_dir_path, meta_data, res)
        else:
            if not (os.path.isdir(exp_dir_path) and
                    os.path.exists(os.path.join(exp_dir_path, 'meta_data_xcat')) and
                    os.path.exists(os.path.join(exp_dir_path, 'xcat'))):
                # Need to cache
                res = func(*args, **kwargs)
                save_results_to_dir(exp_dir_path, meta_data, res)
            else:
                # Old results exists, compare arguments
                if check_valid_results(exp_dir_path, meta_data):
                    res = read_data(exp_dir_path, 'xcat')
                else:
                    # Need to cache
                    res = func(*args, **kwargs)
                    save_results_to_dir(exp_dir_path, meta_data, res)

        return res
    return wrapper_func


@cache_xcat
def create_xcat(cfg, path, organ, grid, offset=None, orientations=None, num_of_materials=10, load_from_cache=False):
    full_xcat_path = os.path.join(cfg['MAIN_FILE_DIR'], path)
    xcat_vol = XCATVolume(full_xcat_path, organ, grid, offset, orientations)
    return xcat_vol


def save_results_to_dir(path, meta_data, res):
    if not os.path.isdir(path):
        os.makedirs(path)
    save_data(path, meta_data, 'meta_data_xcat')
    save_data(path, res, 'xcat')


def read_previous_results(path):
    if not os.path.isdir(path):
        os.makedirs(path)
    prev_meta_data = read_data(path, 'meta_data_xcat')
    return prev_meta_data


def check_valid_results(path, meta_data):
    fields_to_compare = ['PHANTOM_OFFSET_X', 'PHANTOM_OFFSET_Y', 'PHANTOM_OFFSET_Z', 'NUM_OF_Z_SLICES',
                         'NUM_OF_VOXELS_X', 'NUM_OF_VOXELS_Y', 'VOXEL_HALF_X', 'VOXEL_HALF_Y',
                         'VOXEL_HALF_Z']

    prev_meta_data = read_previous_results(path)
    unequal_fields = {}
    equal_args_meta_data = True
    for ind in np.arange(len(meta_data[0])):
        if type(meta_data[0][ind]) is dict:
            for field in fields_to_compare:
                if meta_data[0][0][field] != prev_meta_data[0][0][field]:
                    unequal_fields[field] = (meta_data[0][0][field], prev_meta_data[0][0][field])
                equal_args_meta_data = equal_args_meta_data and meta_data[0][0][field] == prev_meta_data[0][0][field]
        else:
            equal_args_meta_data = equal_args_meta_data and meta_data[0][ind] == prev_meta_data[0][ind]

    prev_meta_data[1].pop('load_from_cache', None)
    meta_data_load_cache = meta_data[1].pop('load_from_cache', None)

    equal_kwargs_meta_data = True
    for key in prev_meta_data[1]:
        if type(prev_meta_data[1][key]) == np.ndarray:
            if not np.array_equal(meta_data[1][key], prev_meta_data[1][key]):
                unequal_fields[key] = (meta_data[1][key], prev_meta_data[1][key])
            eq = np.array_equal(prev_meta_data[1][key], meta_data[1][key])
        else:
            if meta_data[1][key] != prev_meta_data[1][key]:
                unequal_fields[key] = (meta_data[1][key], prev_meta_data[1][key])
            eq = prev_meta_data[1][key] == meta_data[1][key]

        equal_kwargs_meta_data = equal_kwargs_meta_data and eq
    meta_data[1]['load_from_cache'] = meta_data_load_cache
    if len(unequal_fields) > 0:
        for field, values in unequal_fields.items():
            print(f'{field} - new: {values[0]}, old: {values[1]}')
    return equal_args_meta_data and equal_kwargs_meta_data


