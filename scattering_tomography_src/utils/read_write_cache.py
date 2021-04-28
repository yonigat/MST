import pickle
import os


def save_data(path, obj, file_name):
    with open(os.path.join(path, file_name), 'wb') as f:
        pickle.dump(obj, f, protocol=4)


def read_data(path, file_name):
    with open(os.path.join(path, file_name), 'rb') as f:
        obj = pickle.load(f)

    return obj


def get_experiment_path(cfg: dict, data=True, GT_data = False):
    if GT_data and cfg.get('GT_DATA_DIR') is not None:
        path = os.path.join(cfg['MAIN_FILE_DIR'], 'data', cfg['GT_DATA_DIR'])
    else:
        path = os.path.join(cfg['MAIN_FILE_DIR'], 'data' if data else 'experiments', cfg['DATA_DIR'] if data else cfg['EXP_DIR'])
    return path
