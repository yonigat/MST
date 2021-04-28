import numpy as np
from volume.abstract_volume import AbstractVolume


def calc_error_tmp(gt_vol: AbstractVolume, vol: AbstractVolume):
    assert gt_vol.get_grid() == vol.get_grid()
    assert gt_vol.get_active_elements() == vol.get_active_elements()
    error = np.linalg.norm(gt_vol.get_density_vec() - vol.get_density_vec())
    return error
