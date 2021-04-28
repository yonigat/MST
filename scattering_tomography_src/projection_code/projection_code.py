import random
import warnings
import numpy as np
import os


class ProjectionCode:

    def __init__(self, num_sources, num_shots, num_projections_in_shot):
        self._num_sources = num_sources
        self._num_shots = num_shots
        self._num_projection_in_shot = num_projections_in_shot
        self._projection_code = None

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            equal = (self._num_sources == other._num_sources) and \
                    (self._num_shots == other._num_shots) and \
                    (self._num_projection_in_shot == other._num_projection_in_shot) and \
                    (np.array_equal(self._projection_code, other._projection_code))

            return equal
        else:
            return False

    def create_random_code(self, use_seed: bool = True, seed: int = 42):
        if use_seed:
            np.random.seed(seed)
        arr = np.random.permutation(np.arange(self._num_sources))
        arr = arr[:self._num_shots * self._num_projection_in_shot].reshape((self._num_shots, self._num_projection_in_shot))
        self._projection_code = np.random.permutation(arr)
        return self._projection_code

    def create_serial_code(self):
        self._projection_code = np.arange(self._num_shots * self._num_projection_in_shot).reshape((self._num_shots, self._num_projection_in_shot))
        return self._projection_code

    def create_one_source_code(self):
        self._projection_code = np.zeros((self._num_shots, self._num_projection_in_shot))
        return self._projection_code

    def sample_projection_code(self, full_code):
        sampled_shots = np.random.permutation(full_code.shape[0])[:self._num_shots]
        self._projection_code = full_code[sampled_shots, :]
        return sampled_shots

    def subsample_projection_code(self, full_code: np.ndarray, ind: int):
        """
        cyclic sampling of rows from a big array
        :param full_code:
        :param ind: the row index from where we start sampling (can be larger than the number of rows, then we use mod)
        :return:
        """
        full_num_shots = full_code.shape[0]

        sampled_code = full_code[np.arange(ind*self._num_shots, ind*self._num_shots+self._num_shots) % full_num_shots, :]
        self._projection_code = sampled_code
        return np.arange(ind*self._num_shots, ind*self._num_shots+self._num_shots) % full_num_shots

    def random_subsample_projection_code(self, full_code: np.ndarray):
        """
        cyclic sampling of rows from a big array
        :param full_code:
        :param ind: the row index from where we start sampling (can be larger than the number of rows, then we use mod)
        :return:
        """
        full_num_shots = full_code.shape[0]
        sampled_angles = random.sample(range(1, full_num_shots), self._num_shots)
        sampled_code = full_code[sampled_angles, :]
        self._projection_code = sampled_code
        return sampled_angles

    def rand_perm_existing_code(self, existing_code: np.ndarray, use_seed: bool = True, seed: int = 42):
        if use_seed:
            np.random.seed(seed)
        rand_perm = np.random.permutation(np.arange(existing_code.shape[0]))
        self._projection_code = existing_code[rand_perm, :]
        return existing_code[rand_perm, :]

    def create_lin_space_code(self):
        """
        creates a projection code that spreads the sources in each shot linearly
        for example if there are 16 sources and we do 4 shots, the sources in the 1st shot will be 0, 4, 8, 12
        the second shot will be with sources 1, 5, 9, 13 etc.
        :return:
        """
        if self._num_sources % self._num_shots is not 0:
            warnings.warn("Warning........... The number of shots should be a divisor of the number of sources, using same sources twice.")
            sources2use = np.append(np.arange(self._num_sources), np.arange(self._num_shots*self._num_projection_in_shot - self._num_sources))
            self._projection_code = np.reshape(sources2use,
                                               (self._num_projection_in_shot, self._num_shots)).T
        else:
            self._projection_code = np.reshape(np.arange(self._num_sources),
                                               (self._num_projection_in_shot, self._num_shots)).T
        return self._projection_code

    def write_code_file(self, path, code_file_name):
        file = os.path.join(path, code_file_name)
        if self._projection_code.shape[0] == 1:
            np.savetxt(file, self._projection_code, fmt='%1i')
        else:
            np.savetxt(file, self._projection_code, fmt='%1i')
        return

    def write_intens_file(self, path, code_file_name, num_photons):
        file = os.path.join(path, code_file_name)
        intens_mat = np.ones_like(self._projection_code)*num_photons
        if hasattr(self, '_intens_code'):
            intens_mat = self._intens_code*num_photons
        if intens_mat.shape[0] == 1:
            np.savetxt(file, intens_mat, fmt='%1i')
        else:
            np.savetxt(file, intens_mat, fmt='%1i')
        return

    def write_files(self,  num_photons: int, path="../run_inputs/", code_file_name=["projection_code.txt", "intensity_code.txt"]):
        self.write_code_file(path, code_file_name[0])
        self.write_intens_file(path, code_file_name[1], num_photons)

    def get_code(self):
        return self._projection_code
