from mpl_toolkits.mplot3d import Axes3D
import mpl_toolkits.mplot3d as a3
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import numpy as np


def draw_cuboid(ax: Axes3D, vertices: np.ndarray, color: tuple = None, density: float = np.inf, beta: float = 1., draw_edges: bool=True):
    """
    :param ax: 3D axes handle where we plot the polygons of the detector
    :param vertices: an 3x8 ndarray of the coordinates of the corners of the detector
    :param color: a tuple of 2 strings, first is the color of the voxel, second is the color of the edge of the voxel
    :param beta: parameter to control the voxels opacity rate
    :param density: the density of the voxel
    :return:
    """
    edges = 1.0 if draw_edges else 0.0
    faces = np.array(((0, 2, 6, 4),
                      (1, 3, 7, 5),
                      (0, 1, 3, 2),
                      (4, 5, 7, 6),
                      (0, 1, 5, 4),
                      (2, 3, 7, 6)))  # vertices of every face of the chest, the order of vertices in each face is crucial
    for face_ind in np.arange(faces.shape[0]):
        curr_face = faces[face_ind]
        face_coord = vertices[:, curr_face].T
        chest = a3.art3d.Poly3DCollection([face_coord])
        if color is None:
            chest.set_color('g')
            chest.set_edgecolor((0, 0, 1, edges))

        else:
            chest.set_color(tuple(color))
            # chest.set_edgecolor('b')
            chest.set_edgecolor((0, 0, 1, edges))
        chest.set_alpha(np.tanh(beta*density))
        ax.add_collection3d(chest)
    return ax
