import numpy as np


def calc_cone_angle(src_radius: float, voxel_size: tuple, vol_dim: tuple):
    R = src_radius
    r_x = vol_dim[0] * voxel_size[0] / 2
    r_y = vol_dim[1] * voxel_size[1] / 2
    alpha = np.arctan(r_y/(R - r_x))

    return np.pi/alpha
