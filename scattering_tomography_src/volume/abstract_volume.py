from abc import ABC, abstractmethod

from utils.read_write_cache import get_experiment_path, read_data, save_data
from volume.grid import Grid
from utils.draw_cuboid import draw_cuboid
import numpy as np
from numpy import ndarray
import scipy.linalg as linal
from mpl_toolkits.mplot3d import Axes3D
from sklearn.cluster import KMeans
import matplotlib.pyplot as plt
from periodictable import elements
import utils.delete_files
import operator
import os
from tqdm import tqdm
from scipy.ndimage import gaussian_filter
from utils.identify_outliers import conditional_median_filter

class AbstractVolume(ABC):
    """
    Defines the basic interface for volume/ phantom
    """
    def __init__(self, grid: Grid, active_elements: list, offset: ndarray = None, orientations: ndarray = None, arch_mat=False):
        # Grid: A grid object, containing number of voxels in each dimension and 3d size of each voxel
        self._grid = grid
        # list: a list containing the symbols of the elements that will be used all along the simulation, written into a ditc
        # each key is a a string with elements symbol, the value is a tuple containing the index of the element in the fractional mass array
        # and en element class of this element, containing the physical data about it.
        self._arch_mat = arch_mat
        self._P2Ca_ratio = 0.45
        self._active_elements = {}
        self.create_elements_dict(active_elements)

        # Creates the arrays: _density, _fractional_mass
        self._fractional_mass_fractional_mass = None
        self._density = None
        self._id = None # ndarray: num_voxel_x X num_voxel_y X num_voxel_z need to consist to the G4 input array
        self.create_density_and_fractional_mass_and_id()

        self._Z = None
        self._A = None
        self._num_materials = None
        self._materials = None
        if orientations is None:
            self._orientations = np.eye(3)
        else:
            self._orientations = orientations
        # ndarray: 3 (dims) x 3 (dims), 3 orthogonal orientations of all the voxels, first column responds to the orientation of
        # X face, second column to the orientation of the Y face, third column to the orientation of the Z face

        if offset is None:
            self._offset = np.zeros((3,))
        else:
            self._offset = offset
         # ndarray: 3 (dims) X 1, offset from the axes origin where the middle of the grid is located

        self._air_voxels = np.zeros(grid.get_shape())
        self._bone_voxels = np.zeros(grid.get_shape())

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            equal = (self._grid == other._grid) and \
                    (self._active_elements == other._active_elements) and \
                    (np.allclose(self.get_density_vec(), other.get_density_vec(), atol=0.01)) and \
                    (np.array_equal(self._orientations, other._orientations)) and \
                    (np.array_equal(self._offset, other._offset))

            return equal
        else:
            return False

    @abstractmethod
    def create_density_and_fractional_mass_and_id(self):
        """
        self._fractional_mass: ndarray: num_active_elements X num_voxels_x X num_voxels_y  X num_voxels_z, the fractional mass of every active
        element in every voxel. the some of every [ind1, ind2, ind3, :] has to be 1.

        self._density: ndarray: num_voxels_x X num_voxels_y  X num_voxels_z, the overall density of all the active elements in the volume
        :return:
        """
        raise NotImplemented


    def draw_3d(self, ax: Axes3D = None, beta_color: float = 1.0, beta_density: float = 1.0, draw_edges: bool=True):
        """
        Draws the volume in 3D
        :param beta:
        :param ax: handle of Axes3D figure to draw the slice in
        :return:
        """
        if ax is None:
            fig = plt.figure()
            ax = Axes3D(fig)
        num_materials = len(self._active_elements.keys())
        color_mat = linal.orth(np.random.rand(max(num_materials, 3), min(num_materials, 3)))  # an orthogonal set of 3 vectors in the size of number of elements
                                                                      # s.t we can project each voxel and get RGB colors
        if num_materials > 3:
            color_mat = color_mat.T
        length = self.get_grid().get_shape()[2]
        for ind in np.arange(length):
            self.draw_slice('Z', ind, ax, color_mat, beta_color=beta_color, beta_density=beta_density, draw_edges=draw_edges)

        Sx, Sy, Sz = np.array(self.get_grid().get_voxel_size()) / 2
        Nx, Ny, Nz = np.array(self.get_grid().get_shape())

        ax.set_xlim3d(1.2*(-(Sx * Nx) + self._offset[0]), 1.2*(+(Sx * Nx) + self._offset[0]))
        ax.set_ylim3d(1.2*(-(Sy * Ny) + self._offset[1]), 1.2*(+(Sy * Ny) + self._offset[1]))
        ax.set_zlim3d(1.2*(-(Sz * Nz) + self._offset[2]), 1.2*(+(Sz * Nz) + self._offset[2]))
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')

        return ax

    def draw_slice_2d(self, plane_normal: str = 'Z', slice_num: int = 0, color_mat: ndarray = None, beta_color: float = 1.0, beta_dens: float = 1.0):
        """
        :param plane_normal:  The normal to the slice plane, can be either 'X', 'Y', 'Z' (ignoring sign)
        :param slice_num: Slice number along the normal axis
        :param color_mat: n x 3 ndarray of 3 non-trivial orthogonal vectors which describe the basic color components
        :return:
        """
        fig_2d = plt.figure()

        if color_mat is None:
            num_materials = len(self._active_elements.keys())
            color_mat = linal.orth(np.random.rand(max(num_materials, 3), min(num_materials, 3)))  # an orthogonal set of 3 vectors in the size of number of elements
            # s.t we can project each voxel and get RGB colors
            if num_materials > 3:
                color_mat = color_mat.T

        axis_to_draw_dict = {'X': np.array((1, 2)), 'Y': np.array((0, 2)), 'Z': np.array((0, 1))}
        num_voxels_arr = self.get_grid().get_shape()  # number of voxels in each dimension
        iter_ind = axis_to_draw_dict[plane_normal]
        length, width = np.array(num_voxels_arr)[iter_ind]  # the dimensions of the slice we show
        slice_image = np.zeros([length, width, 4])  # each pixel in this image is an RGBA coordinate

        for ind_dim0 in np.arange(length):
            for ind_dim1 in np.arange(width):
                curr_composition = self.get_curr_composition(iter_ind, np.array((ind_dim0, ind_dim1)), slice_num)
                slice_image[ind_dim0, ind_dim1, 0:3] = np.reshape(np.abs(np.tanh(beta_color * np.dot(color_mat, curr_composition))), (1, 1, 3))
                slice_image[ind_dim0, ind_dim1, 3] = np.tanh(self.get_curr_density(iter_ind, np.array((ind_dim0, ind_dim1)), slice_num) * beta_dens)
                plt.title('Slice number: '+ str(slice_num)+', Plane normal: '+plane_normal)
                plt.tick_params(axis='both', which='both', bottom=False, top=False, labelbottom=False, right=False,
                                left=False, labelleft=False)
                plt.imshow(slice_image)

        return

    def draw_slice(self, plane_normal: str, slice_num: int, ax: Axes3D = None, color_mat: ndarray = None, beta_color: float = 1.0, beta_density: float = 1.0, draw_edges: bool=True):
        """
        Draws a plane slice of the volume
        :param beta:
        :param color_mat:
        :param ax: handle of Axes3D figure to draw the slice in
        :param plane_normal: The normal to the slice plane, can be either 'X', 'Y', 'Z' (ignoring sign)
        :param slice_num: Slice number along the normal axis
        :return ax: handle of Axes3D with the slice drawn in
        """

        # TODO: add some protection from user input, check if the slice really exists?

        if ax is None:
            fig = plt.figure()
            ax = Axes3D(fig)

        if color_mat is None:
            num_materials = len(self._active_elements.keys())
            color_mat = linal.orth(np.random.rand(max(num_materials, 3), min(num_materials, 3)))  # an orthogonal set of 3 vectors in the size of number of elements
                                                                                                  # s.t we can project each voxel and get RGB colors
            if num_materials > 3:  # necessary to set right dimensions
                color_mat = color_mat.T

        axis_to_draw_dict = {'X': np.array((1, 2)), 'Y': np.array((0, 2)), 'Z': np.array((0, 1))}

        num_voxels_arr = self.get_grid().get_shape()  # number of voxels in each dimension
        iter_ind = axis_to_draw_dict[plane_normal]
        iter_dim1, iter_dim2 = np.array(num_voxels_arr)[iter_ind]   # the size of dimensions we need to iterate over

        for ind_dim0 in np.arange(iter_dim1):
            for ind_dim1 in np.arange(iter_dim2):
                vertices = self.vertices(iter_ind, np.array((ind_dim0, ind_dim1)), slice_num)
                curr_density = self.get_curr_density(iter_ind, np.array((ind_dim0, ind_dim1)), slice_num) * beta_density
                curr_composition = self.get_curr_composition(iter_ind, np.array((ind_dim0, ind_dim1)), slice_num)
                color = np.abs(np.tanh(beta_color * np.dot(color_mat, curr_composition)))
                ax = draw_cuboid(ax, vertices, density=curr_density, color=color, draw_edges=draw_edges)

        edges = self.get_edge_points(iter_ind, slice_num)
        min_x, min_y, min_z = np.min(edges, 1)
        max_x, max_y, max_z = np.max(edges, 1)
        ax.set_xlim3d(1.2*min_x, 1.2*max_x)
        ax.set_ylim3d(1.2*min_y, 1.2*max_y)
        ax.set_zlim3d(1.2*min_z, 1.2*max_z)
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')

        return ax

    def vertices(self, iter_ind: ndarray, ind_dim: ndarray, slice_num: int):
        """
        :param iter_ind: the two axis we are iterating over, they are the axis orthogonal to the plane normal
        :param ind_dim: the index of voxel we are currently drawing in the slice
        :param slice_num: the index of the slice in the plane normal axis
        :return: vertices: an ndarray of 3 x 8 of the 3D coordinates of the 8 vertices of the voxel
        """
        size_arr = np.array(self.get_grid().get_voxel_size()) / 2  # half sizes of voxels relative to the main X,Y,Z axis
        shape_arr = np.array(self.get_grid().get_shape())
        center_of_slice = size_arr * (2 * (slice_num - np.floor(shape_arr/2)))  # since one of the coordinates are constant, we create this vector
                                                        # and then change in the relevant coordinates for each voxel
        center_of_slice[iter_ind[0]] = size_arr[iter_ind[0]] * (2 * (ind_dim[0] - np.floor(shape_arr[iter_ind[0]]/2)))
        center_of_slice[iter_ind[1]] = size_arr[iter_ind[1]] * (2 * (ind_dim[1] - np.floor(shape_arr[iter_ind[1]]/2)))
        center_of_slice = center_of_slice.T + self._offset  # add offset of volume if exists
        orientation_x = self._orientations[:, 0]
        orientation_y = self._orientations[:, 1]
        orientation_z = self._orientations[:, 2]

        vertices = np.zeros((3, 8))  # 8 corners in 3D
        vertices[:, 0] = center_of_slice + size_arr[0] * orientation_x + size_arr[1] * orientation_y + size_arr[2] * orientation_z
        vertices[:, 1] = center_of_slice + size_arr[0] * orientation_x + size_arr[1] * orientation_y - size_arr[2] * orientation_z
        vertices[:, 2] = center_of_slice + size_arr[0] * orientation_x - size_arr[1] * orientation_y + size_arr[2] * orientation_z
        vertices[:, 3] = center_of_slice + size_arr[0] * orientation_x - size_arr[1] * orientation_y - size_arr[2] * orientation_z
        vertices[:, 4] = center_of_slice - size_arr[0] * orientation_x + size_arr[1] * orientation_y + size_arr[2] * orientation_z
        vertices[:, 5] = center_of_slice - size_arr[0] * orientation_x + size_arr[1] * orientation_y - size_arr[2] * orientation_z
        vertices[:, 6] = center_of_slice - size_arr[0] * orientation_x - size_arr[1] * orientation_y + size_arr[2] * orientation_z
        vertices[:, 7] = center_of_slice - size_arr[0] * orientation_x - size_arr[1] * orientation_y - size_arr[2] * orientation_z

        return vertices

    def get_edge_points(self, iter_ind: ndarray, slice_num: int):
        """
        :param iter_ind: ndarray of shape (2,), holds the indices of the two variables of the slice currently drawn
        :param slice_num: int, the index of the slice in the plane normal variable
        :return: a (3,8) ndarray of the coordinates of the current slice
        """
        edges = np.zeros((3, 8))
        size_arr = np.array(self.get_grid().get_voxel_size()) / 2  # half sizes of voxels relative to the main X,Y,Z axis
        shape_arr = np.array(self.get_grid().get_shape())
        step_size = np.copy(size_arr)
        step_size[iter_ind[0]] *= shape_arr[iter_ind[0]]
        step_size[iter_ind[1]] *= shape_arr[iter_ind[1]]

        center_of_slice = size_arr * (2 * (slice_num - np.floor(shape_arr / 2)))  # since one of the coordinates are constant, we create this vector
                                                                                # and then change in the relevant coordinates for each voxel
        center_of_slice[iter_ind[0]] = 0
        center_of_slice[iter_ind[1]] = 0
        center_of_slice = center_of_slice.T + self._offset  # add offset of volume if exists
        orientation_x = self._orientations[:, 0]
        orientation_y = self._orientations[:, 1]
        orientation_z = self._orientations[:, 2]

        edges[:, 0] = center_of_slice + step_size[0] * orientation_x + step_size[1] * orientation_y + step_size[2] * orientation_z
        edges[:, 1] = center_of_slice + step_size[0] * orientation_x + step_size[1] * orientation_y - step_size[2] * orientation_z
        edges[:, 2] = center_of_slice + step_size[0] * orientation_x - step_size[1] * orientation_y + step_size[2] * orientation_z
        edges[:, 3] = center_of_slice + step_size[0] * orientation_x - step_size[1] * orientation_y - step_size[2] * orientation_z
        edges[:, 4] = center_of_slice - step_size[0] * orientation_x + step_size[1] * orientation_y + step_size[2] * orientation_z
        edges[:, 5] = center_of_slice - step_size[0] * orientation_x + step_size[1] * orientation_y - step_size[2] * orientation_z
        edges[:, 6] = center_of_slice - step_size[0] * orientation_x - step_size[1] * orientation_y + step_size[2] * orientation_z
        edges[:, 7] = center_of_slice - step_size[0] * orientation_x - step_size[1] * orientation_y - step_size[2] * orientation_z

        return edges

    def reduce_materials(self, num_materials: int, seed: int = 42):
        print("**********Reduce Materials**********")
        all_densities_before_quantization = self.get_density_vec()
        self._num_materials = num_materials
        if num_materials <= 0:
            raise Exception("The number of materials given is illegal")
        if self._fractional_mass is None:
            raise Exception("No mass was defined to the volume")
        num_elements, shape0, shape1, shape2 = self._fractional_mass.shape
        # unique_materials = np.unique(self._fractional_mass.reshape(num_elements, shape0*shape1*shape2).round(decimals=10), axis=1)
        reduced_elements = KMeans(n_clusters=num_materials, random_state=seed).fit(self._fractional_mass.reshape(num_elements, shape0*shape1*shape2).T)
        self._id = reduced_elements.labels_.reshape((shape0, shape1, shape2))
        self._materials = reduced_elements.cluster_centers_
        mask = np.abs(self._materials) < 1e-3
        self._materials[mask] = 0.0
        self._materials = self._materials / np.expand_dims(np.sum(self._materials, axis=1), axis=1)
        # In case all fractional masses were 0 materials will return nan:
        self._materials[np.isnan(self._materials)] = 1.0 / len(self._active_elements)
        for ind0 in tqdm(np.arange(shape0)):
            for ind1 in np.arange(shape1):
                for ind2 in np.arange(shape2):
                    self._fractional_mass[:, ind0, ind1, ind2] = self._materials[self._id[ind0, ind1, ind2], :]
        all_densities_after_quantization = self.get_density_vec()
        # assert np.allclose(all_densities_before_quantization, all_densities_after_quantization, atol=0.05), 'Reduce materials changed the volume too much!'

    def get_curr_density(self, iter_ind: ndarray, ind_dim: ndarray, slice_num: int):
        """
        :param iter_ind: the two axis we are iterating over, they are the axis orthogonal to the plane normal
        :param ind_dim: the index of voxel we are currently drawing in the slice
        :param slice_num: the index of the slice in the plane normal axis
        :return: vertices: an ndarray of 3 x 8 of the 3D coordinates of the 8 vertices of the voxel
        """
        curr_voxel = slice_num * np.ones([3, ], dtype=np.int8)
        curr_voxel[iter_ind[0]] = int(ind_dim[0])
        curr_voxel[iter_ind[1]] = int(ind_dim[1])
        return self._density[curr_voxel[0]][curr_voxel[1]][curr_voxel[2]]

    def get_curr_composition(self, iter_ind: ndarray, ind_dim: ndarray, slice_num: int):
        """
        :param iter_ind: ndarray of shape (2,), holds the indices of the two variables of the slice currently drawn ('x' =1, 'y'=2, 'z'=3)
        :param ind_dim:  ndarray of shape (2,), the index of the voxel we wish to find its fractional mass
        :param slice_num:
        :return: returns an nd_array of the the fractional mass of all active elements in a given voxel used in the simulation.
        """
        curr_voxel = slice_num * np.ones([3, ], dtype=np.int8)
        curr_voxel[iter_ind[0]] = int(ind_dim[0])
        curr_voxel[iter_ind[1]] = int(ind_dim[1])
        return self._fractional_mass[:, curr_voxel[0], curr_voxel[1], curr_voxel[2]]

    def add_voxel_element(self, loc: tuple, element: str, density: float):
        """
        Adds the input element with its input density to the location in the volume
        :param loc: tuple with the (x: int, y: int, z: int) location in the volume of the voxel to update
        :param element: Atomic number of the element to add
        :param density: The density of the element to add
        :return:
        """
        if density < 0.0:
            raise ValueError("Density should be a non-negative number")
            return
        loc_x, loc_y, loc_z = loc
        el_ind, el = self._active_elements[element]  # divide the value of the dict to the index in which we need to updated the matrices, and the element class itself
        new_dens = self._density[loc] + density  # calculate new density of the the voxel
        if new_dens == 0.0:
            return
        new_farc_mass = np.zeros((self._fractional_mass.shape[0]))
        new_farc_mass[el_ind] = density
        overall_mass = self._density[loc_x, loc_y, loc_z] * self._fractional_mass[:, loc_x, loc_y, loc_z] + new_farc_mass  # create new mass vector of the voxel
        overall_mass = overall_mass / np.sum(overall_mass)  # normalize the masses so it fractional mass
        self._fractional_mass[:, loc_x, loc_y, loc_z] = overall_mass
        self._density[loc] = new_dens

    def create_elements_dict(self, element_symbol: list):
        """
        :param element_symbol: a list with the symbols
        :return:
        create a dictionary of the active elements in the simulation, the keys are the symbol of the elements given by user,
        the values are tuples, the first element is the index this element appear in the list, the second is an
        element class of the "periodictable" package containing all the information about the element, to learn more
        about it please read the official documentation:
        http://www.reflectometry.org/danse/docs/elements/index.html
        """

        for ind, el in enumerate(element_symbol):
            self._active_elements[el] = (ind, elements.symbol(el))

    def create_A_Z_list(self):
        self._A = []
        self._Z = []
        O_item = None
        if self._arch_mat:
            O_item = self._active_elements.get('O')
            el_list_to_change = ['C', 'N']

        for val in sorted(self._active_elements.items(), key=operator.itemgetter(1)):
            el = val[1][1]
            el_name = val[0]
            if O_item is not None and el_name in el_list_to_change:
                el = O_item[1]
            self._A.append(el.mass)
            self._Z.append(el.number)

    def elements_mass_ratio(self):
        total_mass = 0
        mass_ratio = np.zeros((len(self._active_elements), 1))
        for el_ind, element in self._active_elements.values():
            total_mass += element.mass
            mass_ratio[el_ind] = element.mass
        return mass_ratio / total_mass

    def set_zero_density(self):
        shape = self._density.shape
        self._density = np.zeros(shape)
        shape = self._fractional_mass.shape
        self._fractional_mass = np.zeros(shape)

    def get_density_vec(self):

        x0 = self._fractional_mass.transpose((1,2,3,0)).reshape((-1,1), order='F').reshape((-1, len(self.get_active_elements())),order='F') * \
             self._density.reshape((self.get_grid().get_num_voxels(), 1), order='F') # set order to be x->y->z
            # creates a num_voxels x num_elements matrix of the density of each element in each voxel
        x0 = x0.reshape((x0.size,), order='F') # do a row stach of the matrix
        return x0

    def get_air_voxels_mask(self):
        air_voxels = np.repeat(np.expand_dims(self._air_voxels, axis=0), len(self.get_active_elements()))
        mask_air = air_voxels.reshape((-1, 1), order='F').reshape((-1, len(self.get_active_elements())))
        mask_air = mask_air.reshape((mask_air.size,))
        return mask_air

    def load_density_vec(self, density_vec: np.ndarray, num_materials: int):
        num_elements = len(self.get_active_elements())
        density_vec = np.reshape(density_vec,
                                 (self.get_grid().get_shape()[0],
                                  self.get_grid().get_shape()[1],
                                  self.get_grid().get_shape()[2],
                                  num_elements),
                                 order="F")
        self.set_zero_density()
        Ca_item = self._active_elements.get("Ca")
        for i in np.arange(num_elements):
            curr_element = self.get_i_element(i)
            for j in np.arange(self.get_grid().get_shape()[0]):
                for k in np.arange(self.get_grid().get_shape()[1]):
                    for l in np.arange(self.get_grid().get_shape()[2]):

                        factor = 1.0
                        if self._arch_mat and curr_element == "P" and Ca_item is not None:
                            i = Ca_item[0]
                            factor = self._P2Ca_ratio
                        if density_vec[j, k, l, i] == 0:
                            continue
                        self.add_voxel_element((j, k, l), curr_element, factor * density_vec[j, k, l, i])

    def get_i_element(self, i: int):
        for key, val in self._active_elements.items():
            if val[0] == i:
                return key

    def blur_vol(self, sigma=1):
        self._density = gaussian_filter(self._density,
                                        sigma=sigma)
        el_num, voxels_x, voxels_y, voxels_z = self._fractional_mass.shape
        for el in np.arange(el_num):
            self._fractional_mass[el, :, :, :] = gaussian_filter(self._fractional_mass[el, :, :, :],
                                                                 sigma=sigma)
        for ix in np.arange(voxels_x):
            for iy in np.arange(voxels_y):
                for iz in np.arange(voxels_z):
                    self._fractional_mass[:, ix, iy, iz] /= np.max(self._fractional_mass[:, ix, iy, iz])

    def AGWN(self, sigma=1, seed=42):
        np.random.seed(seed)
        el_num, voxels_x, voxels_y, voxels_z = self._fractional_mass.shape
        dens_sig = sigma * np.std(self._density)
        # Generate Gaussian noise
        gauss_dens = np.random.normal(0, dens_sig, self._density.size)
        gauss_dens = gauss_dens.reshape(self._density.shape[0], self._density.shape[1], self._density.shape[2])
        gauss_frac = np.random.normal(0, sigma, self._fractional_mass.size)
        gauss_frac = gauss_frac.reshape(el_num, voxels_x, voxels_y, voxels_z)

        # Add noise if the density in that voxel isn't zero
        self._density = np.maximum(np.where(self._density > 0, self._density + gauss_dens, self._density), 0)
        # self._fractional_mass = np.where(self._fractional_mass > 0, self._fractional_mass + gauss_frac, self._fractional_mass)

        # Renormalize the fractional mass
        for ix in np.arange(voxels_x):
            for iy in np.arange(voxels_y):
                for iz in np.arange(voxels_z):
                    self._fractional_mass[:, ix, iy, iz] /= np.max(self._fractional_mass[:, ix, iy, iz])

    def filter_outliers(self, kernel_size: int = 3, th_factor: float = 3.0):
        for el_ind, _ in self._active_elements.values():
            for iz in np.arange(self.get_grid().get_shape()[2]):
                new_frac_mass = conditional_median_filter(self._fractional_mass[el_ind, :, :, iz] * self._density[:, :, iz],
                                                          kernel_size=kernel_size, th_factor=th_factor)
                self._density[:, :, iz] -= self._fractional_mass[el_ind, :, :, iz] * self._density[:, :, iz]
                self._density[:, :, iz] += new_frac_mass


    ######################################
    # Writers
    ######################################
    def write_location_orientation_file(self, path, phantom_location_orientation_name):
        # writes in a single file the offset of the volume (int the first row)
        file = os.path.join(path, phantom_location_orientation_name)
        np.savetxt(file, np.concatenate((np.reshape(self._offset, (1, 3)), self._orientations), axis=0), fmt='%.5e')
        return

    def write_density_files (self, path):
        utils.delete_files.clear_dir(path)
        num_files = self.get_grid().get_shape()[2]  # number of Z slices
        for ind in np.arange(num_files):
            file = os.path.join(path, "dens" + str(ind) + ".txt")
            np.savetxt(file, self._density[:, :, ind].T/2, fmt='%.5e')
        return

    def write_id_files (self, path="../run_inputs/Y/id/"):
        utils.delete_files.clear_dir(path)
        num_files = self.get_grid().get_shape()[2]  # number of Z slices
        for ind in np.arange(num_files):
            file = os.path.join(path, "id" + str(ind) + ".txt")
            np.savetxt(file, self._id[:, :, ind].T, fmt='%d')
        return

    def write_material_file (self, path, mat_file_name):
        if self._materials is None:
            raise Exception("Materials matrix does not exist, create one using reduce_materials method")
        file = os.path.join(path, mat_file_name)
        np.savetxt(file, self._materials, fmt='%.5e')
        return

    def write_elements_file(self, path="../run_inputs/"):
        if self._active_elements is None:
            raise Exception("Active elements list does not exist, create one using reduce_materials method")
        file = path + "elements.txt"
        elements = []
        for val in sorted(self._active_elements.items(), key=operator.itemgetter(1)):
            if self._arch_mat:
                if self._active_elements.get('O') is not None and (val[0] == 'C' or val[0] == 'N'):
                    elements.append('O')
                    continue
                if self._active_elements.get('Ca') is not None and val[0] == 'P':
                    elements.append('Ca')
                    continue
            elements.append(val[0])
        with open(file, 'w+') as f:
            f.write(' '.join(elements))

    def write_A_Z_file(self, path="../run_inputs/"):
        if self._Z is None or self._A is None:
            raise Exception("Z list does not exist, create one using reduce_materials method")
        A_file = path + "A.txt"
        Z_file = path + "Z.txt"
        f = open(Z_file, 'w+')
        for i in self._Z:
            f.write(str(i) + " ")
        f = open(A_file, 'w+')
        for i in self._A:
            f.write(str(i) + " ")

    def add_bone_voxel(self, ix, iy, iz):
        self._bone_voxels[ix, iy, iz] = 1

    def add_air_voxel(self, ix, iy, iz):
        self._air_voxels[ix, iy, iz] = 1

    def write_files(self, path="../run_inputs/",
                    phantom_loc_ori_file="phantom_location_orientation.txt",
                    Y_dir="Y/",
                    mat_file="materials.txt",
                    num_materials=10):
        self.reduce_materials(num_materials)
        self.write_location_orientation_file(path, phantom_loc_ori_file)
        dens_path = os.path.join(path, Y_dir, "dens")
        id_path = os.path.join(path, Y_dir, "id")
        self.write_density_files(dens_path)
        self.write_id_files(id_path)
        self.write_material_file(path, mat_file)

    ######################################
    # Getters
    ######################################
    def get_grid(self) -> Grid:
        """
        :return: Grid object
        """
        return self._grid

    def get_active_elements(self) -> dict:
        """
        API method
        :return: dictionary - the keys are the symbol of the elements given by user,
        the values are tuples, the first element is the index this element appear in the list, the second is an
        element class of the "periodictable" package containing all the information about the element
        """
        return self._active_elements

    def get_materials(self) -> np.ndarray:
        """
        API method.
        TODO: document the exact structure of the returned array. This is an API method.
        :return:
        """
        return self._materials

    def get_ids(self) -> np.ndarray:
        """
        API method.
        Returns the voxels ids array. ndarray of type int.
        The shape of the array is: num_voxel_z X num_voxel_y X num_voxel_x.
        """
        return self._id

    def get_density(self) -> np.ndarray:
        """
        API method.
        Returns the voxels density array. ndarray of type float.
        The shape of the array is: num_voxel_z X num_voxel_y X num_voxel_x.
        """
        return self._density

    def get_offset(self):
        return self._offset

    def get_orientations(self) -> np.ndarray:
        """
        API method
        :return: ndarray: 3 (dims) x 3 (dims), 3 orthogonal orientations of all the voxels, first column responds to the orientation of
        # X face, second column to the orientation of the Y face, third column to the orientation of the Z face
        """
        return self._orientations

    def get_air_voxels(self):
        return self._air_voxels

    def get_bone_voxels(self):
        return self._bone_voxels

    ######################################
    # Loader
    ######################################


def load_volume(cfg: dict, name: str ='init_vol'):
    if name == 'init_vol':
        exp_dir_path = get_experiment_path(cfg)
    else:
        exp_dir_path = get_experiment_path(cfg)

    if not (os.path.isdir(exp_dir_path) and
        os.path.exists(os.path.join(exp_dir_path, 'run_cfg_'+name)) and
        os.path.exists(os.path.join(exp_dir_path, name))):

        raise Exception('Missing files to load the volume for')
    elif check_valid_results(exp_dir_path, cfg):
        vol = read_data(exp_dir_path, name)
    else:
        raise Exception('Saved volume does not match current volume')
    return vol


def check_valid_results(path, run_cfg):
    fields_to_compare = ['PHANTOM_OFFSET_X','PHANTOM_OFFSET_Y','PHANTOM_OFFSET_Z','NUM_OF_Z_SLICES','NUM_OF_VOXELS_X',
                         'NUM_OF_VOXELS_Y','VOXEL_HALF_X','VOXEL_HALF_Y','VOXEL_HALF_Z','NUM_OF_MATERIALS']

    prev_cfg = read_previous_results(path)

    # Compare the configs:
    equal_cfgs = True
    for field in fields_to_compare:
        equal_cfgs = equal_cfgs and run_cfg[field] == prev_cfg[field]

    return equal_cfgs


def read_previous_results(path):
    if not os.path.isdir(path):
        os.makedirs(path)
    prev_cfg = read_data(path, 'run_cfg_init_vol')
    return prev_cfg


def save_volume(path, cfg, res, name='init_vol'):
    if not os.path.isdir(path):
        os.makedirs(path)
    save_data(path, cfg, 'run_cfg_' + name)
    save_data(path, res, name)
