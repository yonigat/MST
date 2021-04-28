import yaml
import os
from utils.read_write_cache import get_experiment_path


def read_configs(file_path: str):
    """

    :param file_path: config file path relative to the main "scattering_tomography_src" dir
    :return:
    """
    with open(file_path, 'r') as f:
        cfg = yaml.load(f)
    return cfg


def write_configs(cfg: dict):
    """

    :param file_path: full path
    :return:
    """
    exp_dir_path = get_experiment_path(cfg)
    if not os.path.isdir(exp_dir_path):
        os.makedirs(exp_dir_path)
    with open(os.path.join(exp_dir_path, 'curr_run_cfg.yml'), 'w') as f:
        yaml.dump(cfg, f)
