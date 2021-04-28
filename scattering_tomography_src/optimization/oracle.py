from copy import deepcopy

from g4driver.g4driver import G4Driver
from volume.abstract_volume import AbstractVolume
from projection_code.projection_code import ProjectionCode
from projections.projections import Projections
import numpy as np
from utils.tf_logger import Logger
from utils.plotting_util import calc_delta, calc_epsilon
from volume.abstract_volume import save_volume
from utils.read_write_cache import get_experiment_path


def oracle_1st_order(x: np.ndarray,
                     iter: int,
                     cfg: dict,
                     G4drive: G4Driver,
                     vol: AbstractVolume,
                     gt_vol: AbstractVolume,
                     code: ProjectionCode,
                     projections: Projections,
                     logger: Logger):
    # 1. load the given volume
    vol.load_density_vec(x, cfg["NUM_OF_MATERIALS"])

    # 2. run forward and calc error
    I, angles = G4drive.runG4(vol, code, grad=False, GT=False)

    loss = projections.calc_loss(angles=angles, I=I)

    # 3. run inverse
    grad, P = G4drive.runG4(vol, code, grad=True, GT=False, I_previous=I, angles=angles)
    grad = np.sum(-grad, axis=0)  # sum grad over all shots

    # error = calc_error_tmp(gt_vol, vol)
    logger.log_scalar('loss', loss, projections.itr)
    # logger.log_scalar('error', error, projections.itr)
    # all_density = vol.get_density()
    # logger.log_images('slices', [all_density[:, :, 0], all_density[:, :, 1]], projections.itr)

    grad = np.hstack(grad)
    grad[vol.get_air_voxels_mask().astype(bool)] = 0
    return


def oracle_1st_order_helper(cfg: dict, G4drive: G4Driver, vol: AbstractVolume, gt_vol: AbstractVolume,
                            code: ProjectionCode, projections: Projections, logger: Logger):
    return lambda x, y: oracle_1st_order(x, y, cfg, G4drive, vol, gt_vol, code, projections, logger)


def ring_oracle_1st_order(x: np.ndarray,
                          iter: int,
                          cfg: dict,
                          G4drive: G4Driver,
                          vol: AbstractVolume,
                          gt_vol: AbstractVolume,
                          full_code: ProjectionCode,
                          projections: Projections,
                          logger: Logger):
    # 1. sample random projection code
    code = ProjectionCode(cfg['NUM_OF_SOURCES'], cfg['NUM_OF_SHOTS_OPTIMIZATION'], cfg["NUM_PROJECTIONS_IN_SHOT"])
    sampled_shots = code.subsample_projection_code(full_code.get_code(), projections.itr)

    # 2. load the given volume
    vol.load_density_vec(x, cfg["NUM_OF_MATERIALS"])

    # 3. run forward and calc error
    I, _ = G4drive.runG4(vol, code, grad=False, GT=False)
    # I0, _ = G4drive.runG4(vol, code, grad=False, GT=False, build_phantom=False)

    # 4. run inverse
    grad, P = G4drive.runG4(vol, code, grad=True, GT=False, I_previous=I, angles=sampled_shots)
    grad = grad / (np.expand_dims(P, -1) + 1)
    grad = np.sum(grad, axis=0) / cfg['NUM_OF_SHOTS_OPTIMIZATION']  # sum grad over all shots

    vol_shape = vol._air_voxels.shape

    P_to_log = np.zeros((P.shape[0], vol_shape[0], vol_shape[1], len(cfg['SLICES_TO_PLOT'])))
    for i in np.arange(P.shape[0]):
        P_to_log[i, :, :, :] = np.reshape(P[i, :], (vol_shape[0], vol_shape[1], vol_shape[2]), order='F')[:, :, cfg['SLICES_TO_PLOT']]
    dens_per_el = np.zeros((grad.shape[1], vol_shape[0], vol_shape[1], len(cfg['SLICES_TO_PLOT'])))
    for el in range(grad.shape[1]):
        dens_per_el[el, :, :, :] = vol.get_density()[:, :, cfg['SLICES_TO_PLOT']] * vol._fractional_mass[el, :, :, cfg['SLICES_TO_PLOT']].transpose(1, 2, 0)

    # sampled_shots = code.subsample_projection_code(full_code.get_code(), 0)
    # I, _ = G4drive.runG4(vol, code, grad=False, GT=False)
    error = projections.calc_error(I=I, angles=sampled_shots,
                                   n_photons_GT=cfg['NUM_OF_PHOTONS_GT'] * cfg['NUM_OF_SOURCES'],
                                   n_photons_fward=cfg['NUM_OF_PHOTONS_FORWARD'] * cfg[
                                       'NUM_OF_SOURCES'])

    loss = projections.calc_loss(angles=sampled_shots, I=I,
                                 n_photons_fward=
                                 cfg['NUM_PROJECTIONS_IN_SHOT'] * cfg['NUM_OF_SHOTS_GT'] * cfg[
                                     'NUM_OF_PHOTONS_FORWARD'])
    eps = calc_epsilon(gt_vol, vol)
    delta = calc_delta(gt_vol, vol)

    logger.log_scalar('loss', loss, projections.itr)
    logger.log_scalar('epsilon', eps, projections.itr)
    logger.log_scalar('delta', delta, projections.itr)

    eps_el = calc_epsilon(gt_vol, vol, elements=True)
    delta_el = calc_delta(gt_vol, vol, elements=True)
    for el in vol.get_active_elements():
        logger.log_scalar(f'{el} epsilon', eps_el[el], projections.itr)
        logger.log_scalar(f'{el} delta', delta_el[el], projections.itr)

    if loss <= logger._optimal_loss:
        save_volume(get_experiment_path(cfg), cfg, vol, name='volume_checkpoint')
        logger._optimal_loss = loss

    logger.log_multiple_figures('vol X face', projections.itr, vol.get_density()[:, :, cfg['SLICES_TO_PLOT']], colorbar=True, dim_permute=(2,0,1))
    logger.log_multiple_figures('dens per elenment', projections.itr, dens_per_el, colorbar=True, dim_permute=(2,0,1))
    # logger.log_multiple_figures('Error term', projections.itr, error, colorbar=True)
    # logger.log_multiple_figures('gradient per element', projections.itr, grad_to_log, colorbar=True)
    logger.log_multiple_figures('P', projections.itr, P_to_log, colorbar=True, dim_permute=(2,0,1))
    return np.reshape(grad, (-1,), order='F')


def ring_oracle_1st_order_helper(cfg: dict, G4drive: G4Driver, vol: AbstractVolume, gt_vol: AbstractVolume,
                                 code: ProjectionCode, projections: Projections, logger: Logger):
    return lambda x, y: ring_oracle_1st_order(x, y, cfg, G4drive, vol, gt_vol, code, projections, logger)


def constrain_func(x: np.ndarray, vol: AbstractVolume, cfg: dict, kernel_size: int = 7, th_factor: float = 3.0):
    vol.load_density_vec(x, cfg["NUM_OF_MATERIALS"])
    vol.filter_outliers(kernel_size=kernel_size, th_factor=th_factor)
    return vol.get_density_vec()


def adjust_func(x: np.ndarray, vol: AbstractVolume, cfg: dict):
    vol.load_density_vec(x, cfg["NUM_OF_MATERIALS"])
    vol._density[vol.get_bone_voxels() == 0] += vol._density[vol.get_bone_voxels() == 0] * 0.003 * (cfg['NUM_PROJECTIONS_IN_SHOT'])

    return vol.get_density_vec()


def adjust_func_helper(vol: AbstractVolume, cfg: dict):
    return lambda x: adjust_func(x, vol, cfg)


def constrain_func_helper(vol: AbstractVolume, cfg: dict, kernel_size: int = 7, th_factor: float = 3.0):
    return lambda x: constrain_func(x, vol, cfg, kernel_size=kernel_size, th_factor=th_factor)
