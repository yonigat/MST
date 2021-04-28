import os

import numpy as np

from utils.read_write_cache import read_data, save_data, get_experiment_path


def cache_reconstruction(func):
    def wrapper_func(*args, **kwargs):
        # First check if we need to cache/ check for previous cache:
        exp_dir_path = get_experiment_path(args[0]._cfg, GT_data=True)
        meta_data = [args, kwargs]
        if not kwargs['load_from_cache']:
            res = func(*args, **kwargs)
            save_results_to_dir(exp_dir_path, meta_data, res)
        else:
            if not (os.path.isdir(exp_dir_path) and
                    os.path.exists(os.path.join(exp_dir_path, 'meta_data_reconstruction')) and
                    os.path.exists(os.path.join(exp_dir_path, 'linear_reconstruction'))):
                # Need to cache
                res = func(*args, **kwargs)
                save_results_to_dir(exp_dir_path, meta_data, res)
            else:
                # Old results exists, compare arguments
                if check_valid_results(exp_dir_path, meta_data):
                    res = read_data(exp_dir_path, 'linear_reconstruction')
                else:
                    # Need to cache
                    res = func(*args, **kwargs)
                    save_results_to_dir(exp_dir_path, meta_data, res)

        return res

    return wrapper_func


def save_results_to_dir(path, meta_data, res):
    if not os.path.isdir(path):
        os.makedirs(path)
    save_data(path, meta_data, 'meta_data_reconstruction')
    save_data(path, res, 'linear_reconstruction')


def read_previous_results(path):
    if not os.path.isdir(path):
        os.makedirs(path)
    prev_meta_data = read_data(path, 'meta_data_reconstruction')
    return prev_meta_data


def check_valid_results(path, meta_data):
    fields_to_compare = ['algorithm', 'iterations', 'initializer']
    prev_meta_data = read_previous_results(path)
    unequal_fields = {}

    equal_kwargs_meta_data = True
    for key in prev_meta_data[1]:
        if key in fields_to_compare:
            if type(prev_meta_data[1][key]) == np.ndarray:
                if not np.array_equal(meta_data[1][key], prev_meta_data[1][key]):
                    unequal_fields[key] = (meta_data[1][key], prev_meta_data[1][key])
                eq = np.array_equal(prev_meta_data[1][key], meta_data[1][key])
            else:
                if meta_data[1][key] != prev_meta_data[1][key]:
                    unequal_fields[key] = (meta_data[1][key], prev_meta_data[1][key])
                eq = prev_meta_data[1][key] == meta_data[1][key]
            equal_kwargs_meta_data = equal_kwargs_meta_data and eq

    if len(unequal_fields) > 0:
        for field, values in unequal_fields.items():
            print(f'{field} - new: {values[0]}, old: {values[1]}')

    return meta_data[1].get('images_from_cache', False) and equal_kwargs_meta_data
