import os
import numpy as np
from utils.delete_files import clear_dir


def read_angles(path="../../run_outputs/angles"):
    angles = []
    angle_files = os.listdir(path)
    angle_files.sort()
    for file in angle_files:
        angle = np.loadtxt(path + file, delimiter=',', dtype='int')
        angles.append(angle.reshape((angle.size,))[0])
    clear_dir(path)
    return angles


def write_angles(angles: list, path="../../run_inputs/angles"):
    clear_dir(path)
    for ind, angle in enumerate(angles):
        file = os.path.join(path, str(ind) + '.csv')
        np.savetxt(file, np.array(angle, ndmin=1), delimiter=" ")
