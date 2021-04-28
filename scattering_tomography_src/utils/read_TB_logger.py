import os

from tensorflow.tensorboard.backend.event_processing import event_accumulator
import numpy as np


def load_run(path: str):
    event_acc = event_accumulator.EventAccumulator(path)
    event_acc.Reload()
    data = {}

    for tag in sorted(event_acc.Tags()["scalars"]):
        x, y = [], []

        for scalar_event in event_acc.Scalars(tag):
            x.append(scalar_event.step)
            y.append(scalar_event.value)

        data[tag] = (np.asarray(x), np.asarray(y))
    return data


def read_all_exp_loggers(dir_path: str):
    if not os.path.exists(dir_path):
        raise Exception('Path given does not exist')
    exp_dirs = os.listdir(dir_path)

    data = {}

    for curr_dir in exp_dirs:
        data[curr_dir] = load_run(os.path.join(dir_path, curr_dir))

    return data