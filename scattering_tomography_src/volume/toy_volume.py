from volume.abstract_volume import AbstractVolume
from volume.grid import Grid
import numpy as np
from numpy import ndarray

import matplotlib.pyplot as plt


class ToyVolume(AbstractVolume):
    """
    Creates a toy volume
    """

    def create_density_and_fractional_mass_and_id(self):
        """
        self._fractional_mass: ndarray: num_voxels_x X num_voxels_y  X num_voxels_z X num_active_elements, the fractional mass of every active
        element in every voxel. the some of every [ind1, ind2, ind3, :] has to be 1.

        self._density: ndarray: num_voxels_x X num_voxels_y  X num_voxels_z, the overall density of all the active elements in the volume
        :return:
        """
        self._density = np.zeros(self.get_grid().get_shape())
        self._fractional_mass = np.zeros((len(self._active_elements),) + self.get_grid().get_shape())
        self._id = np.zeros(self.get_grid().get_shape(), dtype='int64') # ndarray: num_voxel_z X num_voxel_y X num_voxel_x need to consist to the G4 input array

    def load_initialization(self, linear_tomo: ndarray):
        """
        takes an ndarray of betas from the linear tomography and translates it into fractional mass of materials of all
        active elements.
        :param linear_tomo:
        :return:
        """
        grid_shape = self.get_grid().get_shape()
        mass_ratio = self.elements_mass_ratio()
        O_item = self._active_elements.get("O")
        for idx_x in np.arange(grid_shape[0]):
            for idx_y in np.arange(grid_shape[1]):
                for idx_z in np.arange(grid_shape[2]):
                    if linear_tomo[idx_x, idx_y, idx_z] == 0:
                        continue
                    for el_ind, el in self._active_elements.values():
                        if self._arch_mat and el.symbol is 'P' and O_item is not None:
                            self.add_voxel_element((idx_x, idx_y, idx_z), el.symbol,
                                                   self._P2Ca_ratio * self.mass_ratio[O_item[0]] * linear_tomo[idx_x, idx_y, idx_z])
                        else:
                            self.add_voxel_element((idx_x, idx_y, idx_z), el.symbol, mass_ratio[el_ind]*linear_tomo[idx_x, idx_y, idx_z])

    def create_simple_phantom(self, num_materials: int, seed: int = 42):
        # self._materials = np.array([[0.5,  0.5]])
        np.random.set_state(seed)
        self._materials = np.random.rand(num_materials, len(self._active_elements.keys()))
        self._materials = self._materials / np.expand_dims(np.sum(self._materials, axis=1), axis=1)

        for ix in np.arange(self.get_grid().get_shape()[0]):
            for iy in np.arange(self.get_grid().get_shape()[1]):
                for iz in np.arange(self.get_grid().get_shape()[2]):
                    mat_ind = np.random.randint(0, num_materials)
                    # density = np.random.rand()
                    density = 0.5
                    for key, val in self._active_elements.items():
                        self.add_voxel_element((ix, iy, iz), key, density*self._materials[mat_ind, val[0]])

    def set_vol_from_arrays(self, dens: np.ndarray, fracs: np.ndarray, num_materials: int):
        """
        A method that sets the internal vol arrays manually by overiding the original arrays
        :param dens: Nx * Ny * Nz
        :param fracs: N_el * Nx * Ny * Nz
        :param ids: Nx * Ny * Nz
        :param num_materials: number of materials to reduce the vol
        :return:
        """
        assert np.shape(dens) == np.shape(self._density)
        assert np.shape(fracs) == np.shape(self._fractional_mass)
        self._density = dens
        self._fractional_mass = fracs
