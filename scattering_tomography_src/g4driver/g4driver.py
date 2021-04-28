from detector.abstract_detector import AbstractDetector
from g4driver.g4driver_util import create_g4_params_dict, write_g4_params_dict, create_inputs_dir, cache_g4
from projection_code.projection_code import ProjectionCode
from source.abstract_source import AbstractSource
from utils.read_write_angles import read_angles, write_angles
from utils.write_mac_files import write_mac_file
from volume.abstract_volume import AbstractVolume
from projections.projections import Projections
import utils.delete_files
import numpy as np
import importlib
import operator
import os
import subprocess


class G4Driver:

    def __init__(self, cfg: dict, detector: AbstractDetector, source: AbstractSource, projections: Projections):
        """

        :param cfg:
        :param
        :param detector:
        :param source:
        """
        self._cfg = cfg
        self._source = source
        self._detector = detector
        self._projection = projections
        self._empty_error = np.zeros(
            (self._cfg["NUM_OF_SHOTS_OPTIMIZATION"], self._cfg["NUM_DETECTOR_ROWS"], self._cfg["NUM_DETECTOR_COLS"]))
        # we create this zero error array to send it in every forward call (where it isn't needed).

        # Create the inputs dir
        create_inputs_dir(os.path.join(self._cfg["MAIN_FILE_DIR"], self._cfg["INPUT_DIR"]))
        self._source.write_files(os.path.join(self._cfg["MAIN_FILE_DIR"], self._cfg["INPUT_DIR"]),
                                 self._cfg['FILE_SOURCES'],
                                 self._cfg['FILE_SOURCES_ORIENTATION'],
                                 self._cfg['FILE_SPECTRUM'])
        self._detector.write_files(os.path.join(self._cfg["MAIN_FILE_DIR"], self._cfg["INPUT_DIR"]),
                                   self._cfg['FILE_DETECTORS_ORIENTATION'],
                                   self._cfg['FILE_DETECTORS'])

    @cache_g4
    def runG4(self, volume: AbstractVolume, projection_code: ProjectionCode, grad=False, build_phantom=True, GT=False,
              I_previous: np.ndarray = None, angles: list = None):

        full_path_input_dir = os.path.join(self._cfg["MAIN_FILE_DIR"], self._cfg["INPUT_DIR"])
        full_path_output_dir = os.path.join(self._cfg["MAIN_FILE_DIR"], self._cfg["OUTPUT_DIR"])
        full_path_grad_dir = os.path.join(self._cfg["MAIN_FILE_DIR"], self._cfg["GRADIENT_DIR"])

        full_path_angles_in = os.path.join(full_path_input_dir, self._cfg["ANGLES_DIR"])
        full_path_angles_out = os.path.join(full_path_output_dir, self._cfg["ANGLES_DIR"])

        error_path = os.path.join(self._cfg['MAIN_FILE_DIR'], self._cfg['ERROR_DIR'])
        g4_dir = os.path.join(self._cfg['MAIN_FILE_DIR'], self._cfg['G4_DIR'])

        # Delete outputs before run
        self.clear_all_outputs(full_path_output_dir, full_path_grad_dir)

        # num projections in shot:
        num_proj_in_shot = self._cfg["NUM_PROJECTIONS_IN_SHOT"]
        num_shots = self._cfg["NUM_OF_SHOTS_GT"]

        if GT:
            n_photons = self._cfg['NUM_OF_PHOTONS_GT']
            error = self._empty_error

        else:
            if grad is True:
                n_photons = self._cfg['NUM_OF_PHOTONS_INVERSE']
                error = self._projection.calc_error(I=I_previous, angles=angles,
                                                    n_photons_GT=self._cfg['NUM_OF_PHOTONS_GT'] * self._cfg[
                                                        'NUM_OF_SOURCES'],
                                                    n_photons_fward=self._cfg['NUM_OF_PHOTONS_FORWARD'] * self._cfg[
                                                        'NUM_OF_SOURCES'])
                write_angles(angles, full_path_angles_in)
            else:
                n_photons = self._cfg['NUM_OF_PHOTONS_FORWARD']
                error = self._empty_error

        overall_photons = num_shots * num_proj_in_shot * n_photons
        params_dict = create_g4_params_dict(self._cfg, grad=grad, n_photons=n_photons, build_phantom=build_phantom,
                                            GT=GT, n_elements=len(volume.get_active_elements()))
        write_g4_params_dict(params_dict, full_path_input_dir)
        self._create_g4_active_elements_dict(volume.get_active_elements(),
                                             write_files=True,
                                             path=full_path_input_dir,
                                             arch_mat=volume._arch_mat)

        self._create_g4_col_loc_ori(volume, write_files=True, path=full_path_input_dir)

        volume.write_files(full_path_input_dir,
                           self._cfg["FILE_PHANTOM_OFFSET_ORIENTATION"],
                           self._cfg["FILE_VOXEL_TO_MATERIALS_Y"],
                           self._cfg["FILE_MATERIALS"],
                           self._cfg['NUM_OF_MATERIALS'])

        self.write_error(error, path=error_path)
        projection_code.write_files(n_photons, full_path_input_dir,
                                    [self._cfg['FILE_PROJECTION_CODE'], self._cfg['FILE_INTENS_CODE']])
        num_shots = self._cfg['NUM_OF_SHOTS_GT'] if GT else self._cfg['NUM_OF_SHOTS_OPTIMIZATION']
        write_mac_file(n_photons, num_shots, path=g4_dir)

        direction = "inverse" if grad else "forward"
        print("***G4 " + direction + " run")

        p = subprocess.Popen(['./exampleB1'], cwd=g4_dir)
        p.wait()

        if grad is True:
            grad, P = self.read_grad(path=full_path_grad_dir)
            if grad is None:  # sometimes G4 has some seg error?
                grad = np.zeros((self._cfg["NUM_OF_SHOTS_OPTIMIZATION"], self._cfg["NUM_OF_VOXELS"],
                                 len(self._cfg['ACTIVE_ELEMENTS'])))
                P = np.zeros((self._cfg["NUM_OF_SHOTS_OPTIMIZATION"], self._cfg["NUM_OF_VOXELS"]))

            air_voxels = volume.get_air_voxels()
            bone_voxels = volume.get_bone_voxels()

            if np.any(air_voxels):
                grad[:, np.reshape(air_voxels == 1, (-1,), order='F'), :] = 0

            if volume._arch_mat:
                Ca_item = volume._active_elements.get('Ca')
                P_item = volume._active_elements.get('P')
                H_item = volume._active_elements.get('H')
                O_item = volume._active_elements.get('O')
                if P_item is not None and Ca_item is not None:
                    O2H_factor = 0.66
                    O_H_grad_sum = (grad[:, np.reshape(bone_voxels == 1, (-1,), order='F'), P_item[0]]
                                    + grad[:, np.reshape(bone_voxels == 1, (-1,), order='F'), Ca_item[0]])
                    grad[:, np.reshape(bone_voxels == 1, (-1,), order='F'), Ca_item[0]] = O_H_grad_sum
                    grad[:, np.reshape(bone_voxels == 1, (-1,), order='F'), P_item[0]] = 0
                # if Ca_item is not None and O_item is not None:
                #     O2Ca_factor = 0.66
                #     O_Ca_grad_sum = (grad[:, np.reshape(bone_voxels == 1, (-1,), order='F'), O_item[0]]
                #                      + grad[:, np.reshape(bone_voxels == 1, (-1,), order='F'), Ca_item[0]])
                #     grad[:, np.reshape(bone_voxels == 1, (-1,), order='F'), O_item[0]] = O2Ca_factor * O_Ca_grad_sum
                #     grad[:, np.reshape(bone_voxels == 1, (-1,), order='F'), Ca_item[0]] = (1 - O2H_factor) * O_Ca_grad_sum

            if np.any(bone_voxels):
                el2null = []
                # grad[:, np.reshape(bone_voxels == 1, (-1,), order='F'), volume._active_elements['Ca'][0]] /= 10
                # for el in ['H', 'C', 'N']:
                #     if volume._active_elements.get(el) is not None:
                #         grad[:, np.reshape(bone_voxels == 1, (-1,), order='F'), volume._active_elements[el][0]] = 0

            # if np.any(np.logical_not(bone_voxels)):
            #     for el in ['P', 'Ca']:
            #         if volume._active_elements.get(el) is not None:
            #             grad[:, np.reshape(bone_voxels == 0, (-1,), order='F'), volume._active_elements[el][0]] = 0
            #            # el2null.append(volume._active_elements[el][0])
            #     # grad[:, np.reshape(np.logical_not(bone_voxels) == 1, (-1,), order='F'), el2null] = 0

            # TODO: which elements to null if bone voxel?
            # before returning - normalize the values by the number of photons:
            # grad = grad / (n_photons * num_proj_in_shot)  # / overall_photons
            # P = P / (n_photons * num_proj_in_shot)  # / overall_photons
            return grad, P
        else:
            I = self.read_images(path=full_path_output_dir)
            angles = read_angles(full_path_angles_out)
            # I *= 1e3  # translate from MeV
            I = I / overall_photons  # / (n_photons * num_proj_in_shot)
            # for i in range(I.shape[0]):
            #     for j in range(I.shape[1]):
            #         I[i, j, :, :] = np.mean(I[i, j, :, :], axis=1, keepdims=True) * np.ones_like(I[i, j, :, :])
            return I, angles

    @staticmethod
    def _create_g4_col_loc_ori(volume: AbstractVolume, write_files=False, path="../run_inputs/") -> np.ndarray:
        """
        Concatenates the location and orientation arrays to a format that is required by G4.
        :param volume:
        :return: np.ndarray. A 4x3 ndarray, the describes the offset and spatial orientation of the volume
                the first row is the offset coordinates, the rows 1-3 is the spatial orientation matrix, each column represents
                orientation in X, Y and Z axis.
        """
        offset = volume.get_offset()
        orientations = volume.get_orientations()
        vol_loc_orient = np.concatenate((np.reshape(offset, (1, 3)), orientations), axis=0)
        if write_files is True:
            file = path + "phantom_location_orientation.txt"
            np.savetxt(file, vol_loc_orient, fmt='%.5e')
        return vol_loc_orient

    @staticmethod
    def _create_g4_active_elements_dict(vol_elements: dict, write_files=False, path="../run_inputs/",
                                        arch_mat=False) -> dict:
        """
        This method is in charge of changing the format between the active elements dict that is returned from AbstractVolume
        to the format required by G4.
        :param: vol_elements - dictionary.
        :return: elements: Dictionary of the available elements. The materials in the volume consists of these elements.
        This dict has the following keys: 'symbol' - list with the elements symbol.
                                          'Z' - list with the Z number.
                                          'A' - list with the molar mass.
        All lists are in the same length which is the number of available elements.
        """
        symbol = []
        A = []
        Z = []
        el_file = path + "elements.txt"
        A_file = path + "A.txt"
        Z_file = path + "Z.txt"

        if arch_mat:
            O_item = vol_elements.get('O')
            Ca_item = vol_elements.get('Ca')
            O_list_to_change = ['C', 'N']
            Ca_list_to_change = ['P']

        for val in sorted(vol_elements.items(), key=operator.itemgetter(1)):
            el_name = val[0]
            symbol.append(val[0])
            el = val[1][1]
            if arch_mat:
                if O_item is not None and el_name in O_list_to_change:
                    el = O_item[1]
                # if Ca_item is not None and el_name in Ca_list_to_change:
                #     el = Ca_item[1]
                symbol.append(el.symbol)
            else:
                symbol.append(val[0])

            A.append(el.mass)
            Z.append(el.number)

        elements = {'symbol': symbol, 'A': A, 'Z': Z}

        if write_files is True:
            with open(el_file, 'w') as f:
                f.write(' '.join(symbol))
            with open(Z_file, 'w') as f:
                for i in Z:
                    f.write(str(i) + " ")
            with open(A_file, 'w') as f:
                for i in A:
                    f.write(str(i) + " ")

        return elements

    @staticmethod
    def write_error(error_term, path: str = "../run_outputs/"):
        utils.delete_files.clear_dir(path)
        num_files = error_term.shape[0]
        if len(error_term.shape) == 2 or error_term.shape[1] == 1 or error_term.shape[2] == 1:
            for ind in np.arange(num_files):
                file = path + "error_" + str(ind) + ".txt"
                np.savetxt(file, error_term[ind], fmt='%.5e', delimiter=" ")
            return
        else:
            error_term = np.squeeze(error_term)
            while len(error_term.shape) <= 2:
                error_term = np.expand_dims(error_term, axis=0)
            num_files = error_term.shape[0]
            for ind in np.arange(num_files):
                file = path + "error_" + str(ind) + ".txt"
                np.savetxt(file, error_term[ind, :, :], fmt='%.5e', delimiter=" ")
        return

    def read_images(self, path: str = "../run_outputs/"):
        path = os.path.join(path, 'I')
        num_scorers = self._cfg["NUM_OF_SCORERS"]
        projection_files = os.listdir(path)
        for ind, file in enumerate(projection_files):
            projection = np.loadtxt(os.path.join(path, file), delimiter=',', skiprows=1)
            while len(projection.shape) < 2:
                projection = np.expand_dims(projection, axis=0)
                projection = projection.T
            spl = file.split('outputEnergyDep_')
            file_ind = int(spl[0][4:])
            scorer_ind = int(spl[1][:-4])
            if ind == 0:
                projections = np.zeros(
                    (int(len(projection_files) / num_scorers), num_scorers, projection.shape[0], projection.shape[1]))
            projections[file_ind, scorer_ind, :, :] = projection
        return projections

    def read_grad(self, path: str = "../run_outputs_grad/"):
        grad_dir = os.path.join(path, 'grad')
        P_dir = os.path.join(path, 'P')
        grad_files = os.listdir(grad_dir)
        P_files = os.listdir(P_dir)
        if grad_files == []:
            print("***seg error accured!")
            grads = None
            Ps = None
            return grads, Ps
        for ind, file in enumerate(grad_files):
            grad = np.loadtxt(os.path.join(grad_dir, file), delimiter=',')
            if len(grad.shape) < 2:
                grad = grad.reshape((grad.size, 1))
            spl = file.split('run_gradient')
            file_ind = int(spl[0])
            if ind == 0:
                grads = np.zeros((len(grad_files), grad.shape[0], grad.shape[1]))
            grads[file_ind, :, :] = grad
        for ind, file in enumerate(P_files):
            P = np.loadtxt(os.path.join(P_dir, file), delimiter=',')
            if len(P.shape) < 2:
                P = P.reshape((P.size, 1))
            spl = file.split('run_P')
            file_ind = int(spl[0])
            if ind == 0:
                Ps = np.zeros(
                    (len(P_files), P.shape[0]))
            Ps[file_ind, :] = np.squeeze(P)
        return grads, Ps

    @staticmethod
    def clear_all_outputs(out_dir, grad_dir):
        # TODO: delete hard coded I,grad,P
        utils.delete_files.clear_dir(os.path.join(out_dir, 'I'))
        utils.delete_files.clear_dir(os.path.join(grad_dir, 'grad'))
        utils.delete_files.clear_dir(os.path.join(grad_dir, 'P'))
        utils.delete_files.clear_dir(os.path.join(out_dir, 'angles'))

    def get_cfg(self):
        return self._cfg
