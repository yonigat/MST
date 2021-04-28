import numpy as np
import pandas as pd
import os
from utils.read_write_cache import read_data, save_data, get_experiment_path
from volume.grid import Grid
from volume.abstract_volume import AbstractVolume
from volume.toy_volume import ToyVolume
from tqdm import tqdm

def cache_vol(func):
    def wrapper_func(*args, **kwargs):
        # First check if we need to cache/ check for previous cache:
        exp_dir_path = get_experiment_path(args[0])
        meta_data = [args, kwargs]
        if not kwargs['load_from_cache']:
            res = func(*args, **kwargs)
            save_results_to_dir(exp_dir_path, meta_data, res)
        else:
            if not (os.path.isdir(exp_dir_path) and
                    os.path.exists(os.path.join(exp_dir_path, 'meta_data_init_vol')) and
                    os.path.exists(os.path.join(exp_dir_path, 'init_vol'))):
                # Need to cache
                res = func(*args, **kwargs)
                save_results_to_dir(exp_dir_path, meta_data, res)
            else:
                # Old results exists, compare arguments
                if check_valid_results(exp_dir_path, meta_data):
                    res = read_data(exp_dir_path, 'init_vol')
                else:
                    # Need to cache
                    res = func(*args, **kwargs)
                    save_results_to_dir(exp_dir_path, meta_data, res)

        return res
    return wrapper_func


def save_results_to_dir(path, meta_data, res):
    if not os.path.isdir(path):
        os.makedirs(path)
    save_data(path, meta_data, 'meta_data_init_vol')
    save_data(path, res, 'init_vol')


def read_previous_results(path):
    if not os.path.isdir(path):
        os.makedirs(path)
    prev_meta_data = read_data(path, 'meta_data_init_vol')
    return prev_meta_data


def check_valid_results(path, meta_data):
    prev_meta_data = read_previous_results(path)
    for i in np.arange(len(meta_data[0])):
        equal_args_meta_data = np.all(meta_data[0][i] == prev_meta_data[0][i])
        if not equal_args_meta_data:
            break
    prev_meta_data[1].pop('load_from_cache', None)
    meta_data_load_cache = meta_data[1].pop('load_from_cache', None)

    equal_kwargs_meta_data = True
    for key in prev_meta_data[1]:
        if type(prev_meta_data[1][key]) == np.ndarray:
            eq = np.array_equal(prev_meta_data[1][key], meta_data[1][key])
        else:
            eq = prev_meta_data[1][key] == meta_data[1][key]

        equal_kwargs_meta_data = equal_kwargs_meta_data and eq
    meta_data[1]['load_from_cache'] = meta_data_load_cache
    return equal_args_meta_data and equal_kwargs_meta_data


def get_data_arrays():
    # https://www.nde-ed.org/EducationResources/CommunityCollege/Radiography/Physics/attenuationCoef.html
    # https://physics.nist.gov/PhysRefData/XrayMassCoef/ComTab/water.html
    # First column are energy values in Mev. Second column is mass attenuation coefficient in cm^2/g values.
    water_atten_coef = np.array([[1.00000E-03, 4.078E+03],
                        [1.50000E-03,  1.376E+03],
                        [2.00000E-03,  6.173E+02],
                        [3.00000E-03,  1.929E+02],
                        [4.00000E-03,  8.278E+01],
                        [5.00000E-03,  4.258E+01],
                        [6.00000E-03,  2.464E+01],
                        [8.00000E-03,  1.037E+01],
                        [1.00000E-02,  5.329E+00],
                        [1.50000E-02,  1.673E+00],
                        [2.00000E-02,  8.096E-01],
                        [3.00000E-02,  3.756E-01],
                        [4.00000E-02,  2.683E-01],
                        [5.00000E-02,  2.269E-01],
                        [6.00000E-02,  2.059E-01],
                        [8.00000E-02,  1.837E-01],
                        [1.00000E-01,  1.707E-01],
                        [1.50000E-01,  1.505E-01],
                        [2.00000E-01,  1.370E-01],
                        [3.00000E-01,  1.186E-01],
                        [4.00000E-01,  1.061E-01],
                        [5.00000E-01,  9.687E-02],
                        [6.00000E-01,  8.956E-02],
                        [8.00000E-01,  7.865E-02],
                        [1.00000E+00,  7.072E-02],
                        [1.25000E+00,  6.323E-02],
                        [1.50000E+00,  5.754E-02],
                        [2.00000E+00,  4.942E-02],
                        [3.00000E+00,  3.969E-02],
                        [4.00000E+00,  3.403E-02],
                        [5.00000E+00,  3.031E-02],
                        [6.00000E+00,  2.770E-02],
                        [8.00000E+00,  2.429E-02],
                        [1.00000E+01,  2.219E-02],
                        [1.50000E+01,  1.941E-02],
                        [2.00000E+01,  1.813E-02]])

    # Conversion between the materials in the Materials.csv file to densities
    dict_mat_to_dens ={
        'Air': [0.0, 0.207],
        'Lung Inhale': [0.207, 0.919],
        'Adipose tissue': [0.919, 0.979],
        'Breat Mammary': [0.979, 1.004],
        'Water': [1.004, 1.015],
        'Red marrow': [1.015, 1.039],
        'Muscle': [1.039, 1.109],
        'Liver': [1.109, 1.2],
        # 'Trabecular': [1.113, 1.496],
        # 'dry rib water': [1.496, 1.654],
        'dry rib water': [1.2, 5.0]
    }

    # Conversion between CT numbers to density
    ct_to_dens = np.array([[-5000,	0.0],
                   [-1000,	0.01],
                   [-82,	0.95],
                   [0,	    1.000],
                   [53,     1.05],
                   [68, 	1.06],
                   [1380,	1.55],
                   [66000,  7.8]])

    return water_atten_coef, dict_mat_to_dens, ct_to_dens

def get_data_array_new():
    atten2tissue = {
        'knee': {
            'Air': [0.0, 0.1],
            'Adipose tissue': [0.1, 0.32],
            'Muscle': [0.32, 0.45],
            'dry rib water': [0.45, 1.5]
        }
    }

    return atten2tissue


def read_materials_mass_atten_tables(dir: str, el_list: list):
    mass_atten_dict = {}
    
    for el in el_list:
        file = os.path.join(dir,f'{el}.csv')
        mass_atten_dict[el] = np.loadtxt(file, delimiter=',')

    return mass_atten_dict


def mass_atten_from_spectrum(atten_coeff_table: np.ndarray, spectrum: np.ndarray):
    spectrum = np.array(spectrum)
    atten_coeff_per_bin = np.interp(np.arange(0.0, float(len(spectrum)))/1000, atten_coeff_table[:, 0], atten_coeff_table[:, 1])
    atten_coeff_tot = np.dot(atten_coeff_per_bin, np.array(spectrum))
    return atten_coeff_tot


def mass_atten_coeff_for_different_elements(elments_list):
    """
    all attenuation coefficient are for 50 KeV, and have units of cm^2/g.
    (50keV because when using the 120keV spectrum source, the effective atten. coeff. fits ~ the 50keV values.
    :param elments_list:
    :return:
    """
    mass_atten_coeff_all = {
        'H': 3.355e-1,
        'He': 1.703e-1,
        'Li': 1.488e-1,
        'Be': 1.554e-1,
        'B': 1.665e-1,
        'C': 1.871e-1,
        'N': 1.98e-1,
        'O': 2.132e-1,
        'F': 2.214e-1,
        'Ne': 2.579e-1,
        'Na': 2.804e-1,
        'Mg': 3.292e-1,
        'Al': 3.681e-1,
        'P': 4.916e-1,
        'S': 5.849e-1,
        'Cl': 6.483e-1,
        'Ar': 7.012e-1,
        'K': 8.679e-1,
        'Ca': 1.019,
        'Sc': 1.087,
        'Ti': 1.213,
        'V': 1.347,
        'Cr': 1.55,
        'Mn': 1.714,
        'Fe': 1.958,
        'I': 1.232e1,
        'Pb': 8.041
    }


    mass_atten_coeff = {}
    for element in elments_list:
        mass_atten_coeff[element] = mass_atten_coeff_all[element]
    return mass_atten_coeff

def get_xcat_materials_df(path='/home/yonatangat/Projects/scattering_tomo/xcat_util'):
    df = pd.read_csv(os.path.join(path, 'Materials.csv'))
    df = df.rename(columns={df.columns[0]: 'Element'})
    df.fillna(0, inplace=True)
    human_elements = list(df['Element'].values[:-1])
    return df, human_elements


@cache_vol
def create_vol_from_ct_recovery(cfg: dict,
                                recovery: np.ndarray,
                                grid: Grid,
                                spectrum: list,
                                offset: np.ndarray = None,
                                orientations: np.ndarray = None,
                                load_from_cache = False,
                                arch_mat = False):
    """

    :param recovery: Nx * Ny * Nz
    :param grid:
    :param energy: KeV
    :param num_of_materials
    :param offset:
    :param orientations:
    :return: ToyVolume object
    """
    df_materials, human_el = get_xcat_materials_df()


    # find the index of the active elements from the full list of elements
    d = {k: v for v, k in enumerate(human_el)}
    el_ind = []
    for k in cfg['ACTIVE_ELEMENTS']:
        if arch_mat:
            if k == 'O':
                el_ind.append(d['C'])
                el_ind.append(d['N'])
            # if k == 'Ca':
            #     el_ind.append(d['P'])
        el_ind.append(d[k])

    vol = ToyVolume(grid, cfg['ACTIVE_ELEMENTS'], offset, orientations, arch_mat=arch_mat)
    vol._material_per_voxel = np.chararray(grid.get_shape(), itemsize=15)
    water_atten_coef, dict_mat_to_dens, ct_to_dens = get_data_arrays()

    # Attenuation coef to CT numbers
    atten_coeff_per_bin = np.interp(np.arange(1.0, float(len(spectrum))+1)/1000, water_atten_coef[:, 0], water_atten_coef[:, 1])
    water_atten = np.dot(atten_coeff_per_bin, np.array(spectrum))
    recovery_ct = (recovery - water_atten) / water_atten * 1e3
    vol._ct_num = recovery_ct
    # CT numbers to density
    recovery_dens = np.interp(recovery_ct, ct_to_dens[:, 0], ct_to_dens[:, 1])
    vol._recover_dens= recovery_dens
    id = np.zeros_like(recovery_dens)
    # Density to materials
    fractional_mass = np.zeros((len(cfg['ACTIVE_ELEMENTS']),) + grid.get_shape())
    for ix in tqdm(range(grid.get_shape()[0])):
        for iy in range(grid.get_shape()[1]):
            for iz in range(grid.get_shape()[2]):
                mat_name = [key for key in dict_mat_to_dens if
                            (dict_mat_to_dens[key][0] <= recovery_dens[ix, iy, iz] < dict_mat_to_dens[key][1])]
                ind = ix + iy * grid.get_shape()[0] + iz * grid.get_shape()[0] * grid.get_shape()[1]
                # TODO: delete later, this is temp workaround to assign all materials in table to some material
                if len(mat_name) == 0:
                    mat_name.append('dry rib water')
                mat_name = mat_name[0]
                vol._material_per_voxel[ix, iy, iz] = mat_name
                id[ix, iy, iz] = df_materials.columns.get_loc(mat_name)
                material = df_materials.loc[el_ind, mat_name]
                fracs = np.array(material)

                if arch_mat and mat_name!='Air':
                    material_new = []
                    # material_new = [0.04, 0.632, 0.328] if mat_name=='dry rib water' else [0.112, 0.888, 0]
                    if 'H' in cfg['ACTIVE_ELEMENTS']:
                        material_new.append(fracs[0])
                    if 'O' in cfg['ACTIVE_ELEMENTS']:
                        material_new.append(np.sum(fracs[1:4]))
                    if 'P' in cfg['ACTIVE_ELEMENTS']:
                        material_new.append(np.sum(fracs[4]))
                    if 'Ca' in cfg['ACTIVE_ELEMENTS']:
                        material_new.append(np.sum(fracs[5]))
                    fracs = np.array(material_new)
                elif arch_mat and mat_name=='Air':
                    fracs = np.array([0, 1.0, 0, 0])

                if mat_name == 'Air':
                    vol.add_air_voxel(ix, iy, iz)
                elif mat_name == 'dry rib water':
                    vol.add_bone_voxel(ix, iy, iz)
                # Fill fractional mass

                # Ensure fracs sum up to 1:
                fracs = fracs * (1/np.sum(fracs))
                fractional_mass[:, ix, iy, iz] = fracs

    vol.set_vol_from_arrays(recovery_dens, fractional_mass, cfg["NUM_OF_MATERIALS"])

    return vol

@cache_vol
def create_vol_simple_table(cfg: dict,
                                recovery: np.ndarray,
                                grid: Grid,
                                organ: str,
                                offset: np.ndarray = None,
                                orientations: np.ndarray = None,
                                load_from_cache = False,
                                arch_mat: bool = False):
    """

    :param recovery: Nx * Ny * Nz
    :param grid:
    :param energy: KeV
    :param num_of_materials
    :param offset:
    :param orientations:
    :return: ToyVolume object
    """
    df_materials, human_el = get_xcat_materials_df()

    # find the index of the active elements from the full list of elements
    d = {k: v for v, k in enumerate(human_el)}
    el_ind = []
    for k in cfg['ACTIVE_ELEMENTS']:
        el_ind.append(d[k])

    density = np.zeros(grid.get_shape())
    vol = ToyVolume(grid, cfg['ACTIVE_ELEMENTS'], offset, orientations, arch_mat=arch_mat)
    vol._material_per_voxel = np.chararray(grid.get_shape(), itemsize=15)
    atten_coeff_to_vol = get_data_array_new()
    atten_coeff_to_vol = atten_coeff_to_vol[organ]


    # CT numbers to density

    # vol._recover_dens= recovery_dens
    # Density to materials
    fractional_mass = np.zeros((len(cfg['ACTIVE_ELEMENTS']),) + grid.get_shape())
    for ix in tqdm(range(grid.get_shape()[0])):
        for iy in range(grid.get_shape()[1]):
            for iz in range(grid.get_shape()[2]):
                mat_name = [key for key in atten_coeff_to_vol if
                            (atten_coeff_to_vol[key][0] <= recovery[ix, iy, iz] < atten_coeff_to_vol[key][1])]
                ind = ix + iy * grid.get_shape()[0] + iz * grid.get_shape()[0] * grid.get_shape()[1]
                ##TODO: delete later, this is temp workaround to assign all materials in table to some material
                if len(mat_name) == 0:
                    mat_name.append('dry rib water')
                mat_name = mat_name[0]
                vol._material_per_voxel[ix, iy, iz] = mat_name
                material = df_materials.loc[el_ind, mat_name]
                density[ix, iy, iz] = float(df_materials.loc[df_materials.shape[0]-1, mat_name])
                if mat_name == 'Air':
                    vol.add_air_voxel(ix,iy,iz)
                elif mat_name == 'dry rib water':
                    vol.add_bone_voxel(ix, iy, iz)
                # Fill fractional mass
                fracs = np.array(material)
                # Ensure fracs sum up to 1:
                fracs = fracs * (1/np.sum(fracs))
                fractional_mass[:, ix, iy, iz] = fracs

    vol.set_vol_from_arrays(density, fractional_mass, cfg["NUM_OF_MATERIALS"])

    return vol


def construct_atten_coeff_vol(vol: AbstractVolume, spectrum: np.ndarray, dir: str, elements: list=None, x_idx=None, y_idx=None, z_idx=None):
    if elements is None:
        elements = list(vol.get_active_elements().keys())
    elements_atten_coeff = read_materials_mass_atten_tables(dir, elements)
    grid = vol.get_grid()
    Nx, Ny, Nz = grid.get_shape()
    Nx = np.arange(Nx)
    Ny = np.arange(Ny)
    Nz = np.arange(Nz)
    if x_idx is not None:
        Nx = x_idx
    if y_idx is not None:
        Ny = y_idx
    if z_idx is not None:
        Nz = z_idx
    atten_coeff_vol = np.zeros((len(Nx), len(Ny), len(Nz)))

    for ind_x, ix in enumerate(Nx):
        for ind_y, iy in enumerate(Ny):
            for ind_z, iz in enumerate(Nz):
                for element in elements:
                    atten = mass_atten_from_spectrum(elements_atten_coeff[element], spectrum)
                    atten_coeff_vol[ind_x, ind_y, ind_z] += atten * vol.get_density()[ix, iy, iz] * \
                                                  vol._fractional_mass[vol.get_active_elements()[element][0], ix, iy, iz]

    return atten_coeff_vol

