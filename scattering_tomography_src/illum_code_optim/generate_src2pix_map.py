from source.abstract_source import AbstractSource
from detector.ring_detector import RingDetector
from projection_code.projection_code import ProjectionCode
import numpy as np

def generate_illum_mat(srcs: AbstractSource, ring_detec: RingDetector):
    num_sources = srcs._n_source_units
    num_rows, num_cols = ring_detec.get_rows(), ring_detec.get_cols()
    W = np.zeros((num_rows*num_cols, num_sources))

    for src in range(num_sources):
        W[:, src] = find_pix_in_cone(srcs, ring_detec, src).reshape((-1,), order='F')

    return W


def generate_src2pix_map(srcs: AbstractSource, ring_detec: RingDetector, proj_code: ProjectionCode):

    code = proj_code.get_code()
    num_shots, num_projection_in_shot = code.shape
    num_rows, num_cols = ring_detec.get_rows(), ring_detec.get_cols()
    I = np.zeros((num_shots, num_rows, num_cols))

    for shot_idx, shot_code in enumerate(code):
        for src_in_shot in shot_code:
            I[shot_idx] += find_pix_in_cone(srcs, ring_detec, src_in_shot)

    return I.reshape((num_shots, num_rows*num_cols), order='F')


def find_pix_in_cone(srcs: AbstractSource, ring_detec: RingDetector, src_in_shot: int):
    num_rows, num_cols = ring_detec.get_rows(), ring_detec.get_cols()
    I = np.zeros((num_rows*num_cols,))
    src_loc = -srcs.get_locations()[src_in_shot, :]
    det_loc = src_loc.T + ring_detec.get_locations()
    src_det_dot_prod = np.matmul(det_loc, src_loc.T)
    angle = np.arccos(src_det_dot_prod / (np.linalg.norm(det_loc, axis=1)*np.linalg.norm(src_loc)))
    I[angle<=srcs._focal_angles[0]] = 1
    I = I.reshape((num_rows, num_cols), order='F')
    return I

def generate_full_code(single_shot_code: np.ndarray, num_shots: int):
    full_code = np.zeros((num_shots, single_shot_code.shape[0]))
    for idx in range(num_shots):
        full_code[idx, :] = np.roll(single_shot_code, idx)
    return full_code
