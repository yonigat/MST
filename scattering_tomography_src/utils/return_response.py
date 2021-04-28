import glob
import numpy as np


def returnResponse(path: str, norm_factor: float = 1.0, scorer_num: int = 2):
    """
    reads and loads into ndarray files of the run from G4
    :param path: string of the path to the directory
    :param norm_factor: normalizing factor for the readings
    :param scorer_num:
    :return:
    """
    projections_pattern = path + '/run*_' + str(scorer_num) + '.csv'
    projection_files = glob.glob(projections_pattern)

    for ind, file in enumerate(projection_files):
        projection = np.loadtxt(file, delimiter=',', skiprows=1).T
        spl = file.split('run_')
        file_ind = int(spl[2][:-21])
        if ind == 0:
            projections = np.zeros((projection.shape[0], projection.shape[1], len(projection_files)))
        projections[:, :, file_ind] = np.fliplr(projection) / norm_factor
    return projections


def returnGradResponse(path: str, norm_factor: float = 1.0):
    """
    reads and loads into ndarray files of the run from G4
    :param path: string of the path to the directory
    :param norm_factor: normalizing factor for the readings
    :param scorer_num:
    :return:
    """
    projections_pattern = path + '*run_gradient.csv'
    projection_files = glob.glob(projections_pattern)

    for ind, file in enumerate(projection_files):
        grad = np.loadtxt(file, delimiter=',', )
        spl = file.split('run_gradient')
        file_ind = int(spl[0][20:])
        if ind == 0:
            Total_grad = np.zeros((grad.shape[0], grad.shape[1], len(projection_files)))
        Total_grad[:, :, file_ind] = np.fliplr(grad) / norm_factor
    return Total_grad

