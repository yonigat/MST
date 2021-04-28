"""Simple example on how to log scalars and images to tensorboard without tensor ops.
License: BSD License 2.0
"""
import time
import os

import PIL
from pygments.formatters import img

__author__ = "Michael Gygli"

import tensorflow as tf
from io import StringIO, BytesIO
from utils.gen_image import gen_image_buffer
from io import StringIO, BytesIO
import matplotlib.pyplot as plt
import numpy as np
from tensorboard.backend.event_processing import event_accumulator
import PIL.Image
from torchvision.transforms import ToTensor


class Logger(object):
    """Logging in tensorboard without tensorflow ops."""

    def __init__(self, log_dir, run_name = None):
        """Creates a summary writer logging to log_dir."""
        if run_name is None:
            log_dir = os.path.join(log_dir, time.strftime("%d_%B_%Y_%H_%M_%S"))
        else:
            log_dir = os.path.join(log_dir, run_name + time.strftime("_%d_%B_%Y_%H_%M_%S"))
        self._log_dir = log_dir

        self.writer = tf.summary.FileWriter(log_dir)

        self._optimal_loss = np.inf

    def log_scalar(self, tag, value, step):
        """Log a scalar variable.
        Parameter
        ----------
        tag : basestring
            Name of the scalar
        value
        step : int
            training iteration
        """
        summary = tf.Summary(value=[tf.Summary.Value(tag=tag,
                                                     simple_value=value)])
        self.writer.add_summary(summary, step)
        self.writer.flush()

    def log_multiple_figures(self,
                             tag: str,
                             step: int,
                             data: np.ndarray,
                             title: str = None,
                             xlable: str = None,
                             ylable: str = None,
                             colorbar: bool = False,
                             dim_permute=(0, 1, 2)):
        """
        loas a 3D matrix as pyplot figures, we iterate over the first dimension (after changing order if we need)
        :param tag: the tag of the data in the tensorboard
        :param step:
        :param data: 3D ndarray (can be more than 3D but the effective dimension has to be 3)
        :param title:
        :param xlable:
        :param ylable:
        :param colorbar:
        :param dim_permute: with this we can control which two dimensions we plot and on which we iterate
        :return:
        """

        # data = np.squeeze(data)
        while len(data.shape) < 3:
            data = np.expand_dims(data, axis=0)
        extra_dim = 1 if len(data.shape) == 3 else data.shape[0]
        extra_dim2el = {0: 'H', 1: 'O', 2: 'Ca'}
        for i in range(extra_dim):
            if extra_dim > 1:
                data_curr = np.transpose(data[i], dim_permute)
            else:
                data_curr = np.transpose(data, dim_permute)
            for ind, curr_image in enumerate(data_curr):
                buff = gen_image_buffer(curr_image, title, xlable, ylable, colorbar)
                el_name = extra_dim2el[i] if extra_dim > 1 else i
                full_tag = f'{tag}/{i}/{ind}'
                self.log_figure(full_tag, buff, step)

    def log_figure(self, tag, plot_buf, step):
        """
        loads a pyplot figure into tensorboard
        :param tag:
        :param plot_buf: a buffer saving a pyplot figure
        :param step:
        :return:
        """
        image = PIL.Image.open(plot_buf)
        image = ToTensor()(image)
        # image = image.permute((1, 2, 0))  # change order to [k, w, h, c] (k num of images)

        img_sum = tf.Summary.Image(encoded_image_string=plot_buf.getvalue(),
                                   height=image.shape[1],
                                   width=image.shape[2])
        summary = tf.Summary(value=[tf.Summary.Value(tag=tag, image=img_sum)])
        self.writer.add_summary(summary, step)
        self.writer.flush()


    def log_images(self, tag, images, step):
        """Logs a list of images."""

        im_summaries = []
        for nr, img in enumerate(images):
            # Write the image to a string
            s = BytesIO()

            plt.imsave(s, img, format='png', vmin=0, vmax=1)

            # Create an Image object
            img_sum = tf.Summary.Image(encoded_image_string=s.getvalue(),
                                       height=img.shape[0],
                                       width=img.shape[1])
            # Create a Summary value
            im_summaries.append(tf.Summary.Value(tag='%s/%d' % (tag, nr),
                                                 image=img_sum))

        # Create and write Summary
        summary = tf.Summary(value=im_summaries)
        self.writer.add_summary(summary, step)
        self.writer.flush()

    def log_histogram(self, tag, values, step, bins=1000):
        """Logs the histogram of a list/vector of values."""
        # Convert to a numpy array
        values = np.array(values)

        # Create histogram using numpy
        counts, bin_edges = np.histogram(values, bins=bins)

        # Fill fields of histogram proto
        hist = tf.HistogramProto()
        hist.min = float(np.min(values))
        hist.max = float(np.max(values))
        hist.num = int(np.prod(values.shape))
        hist.sum = float(np.sum(values))
        hist.sum_squares = float(np.sum(values ** 2))

        # Requires equal number as bins, where the first goes from -DBL_MAX to bin_edges[1]
        # See https://github.com/tensorflow/tensorflow/blob/master/tensorflow/core/framework/summary.proto#L30
        # Thus, we drop the start of the first bin
        bin_edges = bin_edges[1:]

        # Add bin edges and counts
        for edge in bin_edges:
            hist.bucket_limit.append(edge)
        for c in counts:
            hist.bucket.append(c)

        # Create and write Summary
        summary = tf.Summary(value=[tf.Summary.Value(tag=tag, histo=hist)])
        self.writer.add_summary(summary, step)
        self.writer.flush()

    def load_run(self, path: str = None):
        if path is None:
            path = self._log_dir
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
