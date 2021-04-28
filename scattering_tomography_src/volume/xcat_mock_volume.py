from tqdm import tqdm

from volume.abstract_volume import AbstractVolume
from volume.grid import Grid
import numpy as np
import pandas as pd
import os
# from tqdm import tqdm


class XCATMockVolume(AbstractVolume):
    """
    Creates a volume object from XCAT generated files
    """
    def __init__(self,
                 path: str,
                 grid: Grid,):
        """

        :param path: path to the originam XCAT files. Should include Materials.csv, id_to_comp.txt, and a phantom dir with data files
        :param organ: Can be either 'hand' or 'knee'
        :param grid:
        :param offset:
        :param orientations:
        """

        self._xcat_path = path
        df = pd.read_csv(os.path.join(self._xcat_path, 'Materials.csv'))
        df = df.rename(columns={df.columns[0]: 'Element'})
        df.fillna(0, inplace=True)
        human_elements = list(df['Element'].values[:-1])
        self._xcat_materials_df = df

        with open(os.path.join(self._xcat_path, 'id_to_comp.txt'), 'rb') as f:
            id_to_comp = np.genfromtxt(f, delimiter=' ')

        self._xcat_intensity_to_mat_id = dict(zip(id_to_comp[:, 1], id_to_comp[:, 2]))

        super().__init__(grid, human_elements)

    def create_density_and_fractional_mass_and_id(self):
        """
        self._fractional_mass: ndarray: num_voxels_x X num_voxels_y  X num_voxels_z X num_active_elements, the fractional mass of every active
        element in every voxel. the some of every [ind1, ind2, ind3, :] has to be 1.

        self._density: ndarray: num_voxels_x X num_voxels_y  X num_voxels_z, the overall density of all the active elements in the volume
        :return:
        """

        def extract_num(elem):
            """
            Helper method for sorting
            :param elem:
            :return:
            """
            return int(elem.split('_')[-1].split('.')[0])

        # Extract all XCAT data files
        files = [f for f in os.listdir(os.path.join(self._xcat_path,'phantom')) if
                     os.path.isfile(os.path.join(self._xcat_path,'phantom',f))]

        # Sort files by name
        files.sort(key=extract_num)

        grid_shape = self._grid.get_shape()

        self._density = np.zeros(self.get_grid().get_shape())
        self._fractional_mass = np.zeros((len(self._active_elements),) + self.get_grid().get_shape())
        self._material_per_voxel = np.chararray(self.get_grid().get_shape(), itemsize=15)
        # Iterate over the sub vol. and create the density and fractional mass array
        print('Building XCAT...')
        for z_ind in tqdm(range(grid_shape[2])):
            for x_ind in range(grid_shape[0]):
                for y_ind in range(grid_shape[1]):
                    id_ind = self._xcat_materials_df.shape[1]-2
                    # intensity = sub_vol[x_ind, y_ind, z_ind]
                    # material_id = self._xcat_intensity_to_mat_id[intensity]
                    material = self._xcat_materials_df.iloc[:, int(id_ind)]

                    # Fill density
                    self._density[x_ind, y_ind, z_ind] = material.iloc[-1]
                    self._material_per_voxel[x_ind, y_ind, z_ind] = self._xcat_materials_df.columns[int(id_ind)+1]
                    # Fill fractional mass
                    fracs = np.array(material.iloc[:-1])
                    self._fractional_mass[:, x_ind, y_ind, z_ind] = fracs

    def _get_xcat_shifting(self, grid: Grid):
        """
        Helper method to get the correct shifting used for each human body part
        :return:
        """
        half_grid_shape = np.array(grid.get_shape()) // 2
        if self._organ == 'hand':
            return 118 + (50 - half_grid_shape[0]), 0 + (50 - half_grid_shape[1]), 330 + (40 - half_grid_shape[2])
        if self._organ == 'knee':
            return 150 + (35 - half_grid_shape[0]), 90 + (33 - half_grid_shape[1]), 215 + (15 - half_grid_shape[2])
