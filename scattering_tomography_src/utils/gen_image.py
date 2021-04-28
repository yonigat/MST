import matplotlib
import numpy as np
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import io


def gen_image_buffer(data, title=None, xlable=None, ylable=None, colorbar=False):
    """Create a pyplot plot and save to buffer."""
    fig = plt.figure()
    if len(np.squeeze(data).shape) == 1:
        plt.plot(np.squeeze(data))
    else:
        plt.imshow(data)
        if colorbar:
            plt.colorbar()
    plt.xlabel(xlable)
    plt.ylabel(ylable)
    plt.title(title)

    buf = io.BytesIO()
    plt.savefig(buf, format='jpeg', dpi=200)
    buf.seek(0)
    plt.close(fig)
    return buf
