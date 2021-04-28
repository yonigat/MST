import os
import pickle
import random

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import yaml

from linear_tomography.ct_to_vol_util import construct_atten_coeff_vol

mpl.use('TkAgg')
from utils.read_write_cache import read_data
from volume.abstract_volume import load_volume
num_srcs = [1,5]
main_dir = '/home/yonatangat/Projects/scattering_tomo/data'
mass_atten_dir = '../elmental_media_info'
# GT_data_dir = 'knee_1_M1E7_NOISY_60'
# dirs_list = [f'knee_{i}_M1E7_NOISY_60' for i in num_srcs]
GT_data_dir = 'stomach_1_M1E8_1D'
dirs_list = [f'stomach_{i}_M1E8_1D' for i in num_srcs]
GT_vol = read_data(os.path.join(main_dir, GT_data_dir), 'GT_vol')
GT_atten_file = os.path.join(main_dir, GT_data_dir, 'GT_vol_atten')

if os.path.exists(GT_atten_file):
    with open(GT_atten_file, 'rb') as f:
        GT_atten_coeff_full = pickle.load(f)
else:
    GT_atten_coeff_full = None

init_atten_vols = []
init_vols = []
final_vols = []
for dir in dirs_list:
    init_atten_vols.append(read_data(os.path.join(main_dir, dir), 'lin_tomo'))
    init_vols.append(read_data(os.path.join(main_dir, dir), 'init_vol'))
    # final_vols.append(read_data(os.path.join(main_dir, dir), 'final_vol'))
    final_vols.append(read_data(os.path.join(main_dir, dir), 'volume_checkpoint'))

slices = [0]
partial = True
full_body = False
if partial:
    for ind, vols in enumerate(zip(init_vols, final_vols, init_atten_vols)):
        for slice in slices:
            GT_dens = GT_vol.get_density()[:, :, slice]
            init_dens = vols[0].get_density()[:, :, slice]
            final_dens = vols[1].get_density()[:, :, slice]
            file = '/home/yonatangat/Projects/scattering_tomo/scattering_tomography_src/120keV_spectrum.yaml'
            with open(file, 'r') as f:
                cfg = yaml.load(f)
            spectrum = np.array(cfg['SPECTRUM'])
            # # spectrum[0:40] = 0
            spectrum = spectrum / np.sum(spectrum)
            spectrum = list(spectrum)
            # spectrum = np.zeros((150,))
            # spectrum[60] = 1
            # spectrum = list(spectrum)

            if GT_atten_coeff_full is None:
                GT_atten_coeff = construct_atten_coeff_vol(GT_vol, spectrum,
                                                           os.path.join(main_dir, mass_atten_dir), z_idx=[slice])[:, :, 0]
            else:
                GT_atten_coeff = GT_atten_coeff_full[:,:,slice]
            init_atten_coeff = vols[2][:, :, slice]
            # init_atten_coeff = construct_atten_coeff_vol(vols[0], spectrum,
            #                                            os.path.join(main_dir, mass_atten_dir), z_idx=[slice])[:, :, 0]
            final_atten_coeff = construct_atten_coeff_vol(vols[1], spectrum,
                                                       os.path.join(main_dir, mass_atten_dir), z_idx=[slice])[:, :, 0]


            init_sum = np.sum(np.abs(init_atten_coeff))
            final_sum = np.sum(np.abs(final_atten_coeff))
            gt_sum = np.sum(np.abs(GT_atten_coeff))
            init_delta = (init_sum - gt_sum) / gt_sum
            final_delta = (final_sum - gt_sum) / gt_sum

            init_DiffNorm = np.sum(np.abs(GT_atten_coeff - init_atten_coeff))
            final_DiffNorm = np.sum(np.abs(GT_atten_coeff - final_atten_coeff))
            init_epsilonError = init_DiffNorm / gt_sum
            final_epsilonError = final_DiffNorm / gt_sum

            print(f'***Num srcs: {ind}. \n delta_init: {init_delta} \n delta_final: {final_delta}\n eps_init: {init_epsilonError}\n eps_final: {final_epsilonError}')

            # GT_el_dens = np.hstack(GT_vol.get_density() * GT_vol._fractional_mass[ind_el, :, :, :])
            # init_el_dens = np.hstack(vols[0].get_density() * vols[0]._fractional_mass[ind_el, :, :, :])
            # final_el_dens = np.hstack(vols[1].get_density() * vols[1]._fractional_mass[ind_el, :, :, :])

            plt.figure()

            clim = [0, 0]
            clim[0] = np.min([np.min(GT_atten_coeff), np.min(init_atten_coeff), np.min(final_atten_coeff)])
            clim[1] = np.max([np.max(GT_atten_coeff), np.max(init_atten_coeff), np.max(final_atten_coeff)])
            plt.subplot(1, 3, 1)
            im_GT = plt.imshow(GT_atten_coeff, clim=clim, cmap='bone')
            # plt.colorbar(im_GT,
            plt.axis('off')
            # f_ax0.set_title(f'Ground truth density [g/cm^3]')

            plt.subplot(1, 3, 2)
            im_GT = plt.imshow(init_atten_coeff, clim=clim, cmap='bone')
            # plt.colorbar(im_GT,
            plt.axis('off')

            plt.subplot(1, 3, 3)
            im_GT = plt.imshow(final_atten_coeff, clim=clim, cmap='bone')
            # plt.colorbar()
            plt.axis('off')

            plt.figure()

            clim = [0, 0]
            clim[0] = np.min([np.min(GT_atten_coeff), np.min(init_atten_coeff), np.min(final_atten_coeff)])
            clim[1] = np.max([np.max(GT_atten_coeff), np.max(init_atten_coeff), np.max(final_atten_coeff)])
            plt.subplot(1, 3, 1)
            im_GT = plt.imshow(GT_atten_coeff, clim=clim, cmap='bone')
            # plt.colorbar(im_GT,
            plt.axis('off')
            # f_ax0.set_title(f'Ground truth density [g/cm^3]')

            plt.subplot(1, 3, 2)
            im_GT = plt.imshow(init_atten_coeff, clim=clim, cmap='bone')
            # plt.colorbar(im_GT,
            plt.axis('off')

            plt.subplot(1, 3, 3)
            im_GT = plt.imshow(final_atten_coeff, clim=clim, cmap='bone')
            plt.colorbar()
            plt.axis('off')

            fig = plt.figure()
            gs3 = fig.add_gridspec(6, 6)

            f_ax0 = fig.add_subplot(gs3[0:2, 0:2])
            f_ax1 = fig.add_subplot(gs3[0:2, 2:4])
            f_ax2 = fig.add_subplot(gs3[0:2, 4:])
            f_ax3 = fig.add_subplot(gs3[2:, 0:3])
            f_ax4 = fig.add_subplot(gs3[2:, 3:])

            # Densities:
            fig.suptitle(f'Slice {slice} - Number of multiplexed sources {num_srcs[ind]}' if num_srcs[ind]>0 else 'Single source')

            clim = [0, 0]
            clim[0] = np.min([np.min(GT_dens), np.min(init_dens), np.min(final_dens)])
            clim[1] = np.max([np.max(GT_dens), np.max(init_dens), np.max(final_dens)])

            im_GT = f_ax0.imshow(GT_dens, clim=clim)
            plt.colorbar(im_GT, ax=f_ax0)
            f_ax0.axis('off')
            f_ax0.set_title(f'Ground truth density [g/cm^3]')

            im_init = f_ax1.imshow(init_dens, clim=clim)
            plt.colorbar(im_init, ax=f_ax1)
            f_ax1.axis('off')
            f_ax1.set_title(f'Linear CT Initialization density [g/cm^3]')

            im_final = f_ax2.imshow(final_dens, clim=clim)
            plt.colorbar(im_final, ax=f_ax2)
            f_ax2.axis('off')
            f_ax2.set_title(f'final volume density [g/cm^3]')

            # scatter plot - linear init and final vol
            max_y = np.max([np.max(init_dens), np.max(final_dens), np.max(GT_dens)])

            scatter_init = f_ax3.scatter(np.hstack(GT_dens), np.hstack(init_dens))
            f_ax3.set_xlabel('Ground Truth')
            f_ax3.set_ylabel('Linear Initialization')
            max_val = np.max([np.max(GT_dens), np.max(init_dens)])
            f_ax3.plot(np.linspace(0, max_val, GT_vol.get_density().size), np.linspace(0, max_val, GT_vol.get_density().size), 'r--')
            f_ax3.set_xlim(left=-0.02, right=max_y * 1.05)
            f_ax3.set_ylim(bottom=-0.02, top=max_y * 1.05)

            scatter_final = f_ax4.scatter(np.hstack(GT_dens), np.hstack(final_dens))
            f_ax4.set_xlabel('Ground Truth')
            f_ax4.set_ylabel('final volume')
            max_val = np.max([np.max(GT_dens), np.max(final_dens)])
            f_ax4.plot(np.linspace(0, max_val, GT_vol.get_density().size), np.linspace(0, max_val, GT_vol.get_density().size), 'r--')
            f_ax4.set_xlim(left=-0.02, right=max_y * 1.05)
            f_ax4.set_ylim(bottom=-0.02, top=max_y * 1.05)

            fig = plt.figure()
            gs3 = fig.add_gridspec(6, 6)

            f_ax0 = fig.add_subplot(gs3[0:2, 0:2])
            f_ax1 = fig.add_subplot(gs3[0:2, 2:4])
            f_ax2 = fig.add_subplot(gs3[0:2, 4:])
            f_ax3 = fig.add_subplot(gs3[2:, 0:3])
            f_ax4 = fig.add_subplot(gs3[2:, 3:])

            # Attenuation coefficients:
            fig.suptitle(f'Slice {slice} - Number of multiplexed sources {num_srcs[ind]}' if num_srcs[ind] > 0 else 'Single source')

            clim = [0, 0]
            clim[0] = np.min([np.min(init_atten_coeff), np.min(final_atten_coeff)])
            clim[1] = np.max([np.max(init_atten_coeff), np.max(final_atten_coeff)])

            im_init = f_ax0.imshow(init_atten_coeff, clim=clim)
            plt.colorbar(im_init, ax=f_ax0)
            f_ax0.axis('off')
            f_ax0.set_title(f'Initial mass atten. [1/cm]')

            im_final = f_ax1.imshow(final_atten_coeff, clim=clim)
            plt.colorbar(im_final, ax=f_ax1)
            f_ax1.axis('off')
            f_ax1.set_title(f'Final mass atten. [1/cm]')

            im_diff = f_ax2.imshow(final_atten_coeff - init_atten_coeff)
            plt.colorbar(im_diff, ax=f_ax2)
            f_ax2.axis('off')
            f_ax2.set_title(f'Difference between final and initial [1/cm]')

            # scatter plot - linear init and final vol
            max_y = np.max([np.max(init_atten_coeff), np.max(final_atten_coeff), np.max(GT_atten_coeff)])

            scatter_init = f_ax3.scatter(np.hstack(GT_atten_coeff), np.hstack(init_atten_coeff))
            f_ax3.set_xlabel('Ground Truth')
            f_ax3.set_ylabel('Linear Initialization')
            max_val = np.max([np.max(GT_atten_coeff), np.max(init_atten_coeff)])
            f_ax3.plot(np.hstack(GT_atten_coeff), np.hstack(GT_atten_coeff), 'r--')
            # f_ax3.plot(np.linspace(0, max_val, GT_vol.get_density().size), np.linspace(0, max_val, GT_vol.get_density().size),
            #            'r--')
            f_ax3.set_xlim(left=-0.02, right=max_y * 1.05)
            f_ax3.set_ylim(bottom=-0.02, top=max_y * 1.05)

            scatter_final = f_ax4.scatter(np.hstack(GT_atten_coeff), np.hstack(final_atten_coeff))
            f_ax4.set_xlabel('Ground Truth')
            f_ax4.set_ylabel('final volume')
            max_val = np.max([np.max(GT_atten_coeff), np.max(final_atten_coeff)])
            f_ax4.plot(np.hstack(GT_atten_coeff), np.hstack(GT_atten_coeff), 'r--')
            # f_ax4.plot(np.linspace(0, max_val, GT_vol.get_density().size), np.linspace(0, max_val, GT_vol.get_density().size),
            #            'r--')
            f_ax4.set_xlim(left=-0.02, right=max_y * 1.05)
            f_ax4.set_ylim(bottom=-0.02, top=max_y * 1.05)

            # plt.show()
            elements = {}
            # elements = {0: ('H', [0], [0]), 1: ('O', [5, 6, 7], [1]), 2: ('Ca', [13, 18], [2]), 6: ['density']}
            elements = {0: ('H', [0], [0]), 1: ('O', [5, 6, 7], [1]), 2: ('Ca', [13, 18], [2, 3]), 6: ['density']}
            # elements = {0: ('H', [0], [0]), 1: ('O', [1], [1]), 2: ('Ca', [2], [2]), 6: ['density']}
            # for el_ind, el in enumerate(vols[0]._active_elements):
            #     if el == 'P':
            #         continue
            #     if el == 'Ca':
            #         el_ind = 2
            #     gt_frac_mass = np.sum(GT_vol._fractional_mass[elements[el_ind][1], :, :, slice], axis=0) * GT_dens
            #     initial_frac_mass = np.sum(vols[0]._fractional_mass[elements[el_ind][2], :, :, slice], axis=0) * init_dens
            #     final_frac_mass = np.sum(vols[1]._fractional_mass[elements[el_ind][2], :, :, slice], axis=0) * final_dens
            #
            #     fig = plt.figure()
            #     gs3 = fig.add_gridspec(6, 6)
            #
            #     fig.suptitle(f'Slice {slice} - Number of multiplexed sources {num_srcs[ind]}' if num_srcs[ind] > 0 else 'Single source')
            #
            #     f_ax0 = fig.add_subplot(gs3[0:2, 0:2])
            #     f_ax1 = fig.add_subplot(gs3[0:2, 2:4])
            #     f_ax2 = fig.add_subplot(gs3[0:2, 4:])
            #     f_ax3 = fig.add_subplot(gs3[2:, 0:3])
            #     f_ax4 = fig.add_subplot(gs3[2:, 3:])
            #
            #     clim = [0, 0]
            #     clim[0] = np.min([np.min(gt_frac_mass), np.min(initial_frac_mass), np.min(final_frac_mass)])
            #     clim[1] = np.max([np.max(gt_frac_mass), np.max(initial_frac_mass), np.max(final_frac_mass)])
            #
            #     im_GT = f_ax0.imshow(gt_frac_mass, clim=clim)
            #     plt.colorbar(im_GT, ax=f_ax0)
            #     f_ax0.axis('off')
            #     f_ax0.set_title(f'Ground truth {el} dens [g/cm^3]')
            #
            #     im_init = f_ax1.imshow(initial_frac_mass, clim=clim)
            #     plt.colorbar(im_init, ax=f_ax1)
            #     f_ax1.axis('off')
            #     f_ax1.set_title(f'Linear CT Initialization {el} dens [g/cm^3]')
            #
            #     im_final = f_ax2.imshow(final_frac_mass, clim=clim)
            #     plt.colorbar(im_final, ax=f_ax2)
            #     f_ax2.axis('off')
            #     f_ax2.set_title(f'final volume {el} dens [g/cm^3]')
            #
            #     # scatter plot - linear init and final vol
            #     max_y = np.max([np.max(initial_frac_mass), np.max(final_frac_mass), np.max(gt_frac_mass)])
            #
            #     scatter_init = f_ax3.scatter(np.hstack(gt_frac_mass), np.hstack(initial_frac_mass))
            #     f_ax3.set_xlabel('Ground Truth')
            #     f_ax3.set_ylabel('Linear Initialization')
            #     max_val = np.max([np.max(gt_frac_mass), np.max(initial_frac_mass)])
            #     f_ax3.plot(np.linspace(0, max_val, GT_vol.get_density().size), np.linspace(0, max_val, GT_vol.get_density().size), 'r--')
            #     f_ax3.set_xlim(left=-0.02, right=max_y * 1.05)
            #     f_ax3.set_ylim(bottom=-0.02, top=max_y * 1.05)
            #
            #     scatter_final = f_ax4.scatter(np.hstack(gt_frac_mass), np.hstack(final_frac_mass))
            #     f_ax4.set_xlabel('Ground Truth')
            #     f_ax4.set_ylabel('final volume')
            #     max_val = np.max([np.max(gt_frac_mass), np.max(final_frac_mass)])
            #     f_ax4.plot(np.linspace(0, max_val, GT_vol.get_density().size), np.linspace(0, max_val, GT_vol.get_density().size), 'r--')
            #     f_ax4.set_xlim(left=-0.02, right=max_y * 1.05)
            #     f_ax4.set_ylim(bottom=-0.02, top=max_y * 1.05)

if full_body:
    points2plot = random.sample(range(GT_atten_coeff_full.size), int(GT_atten_coeff_full.size*0.05))
    for ind, vols in enumerate(zip(init_vols, final_vols)):
        GT_dens = GT_vol.get_density()
        init_dens = vols[0].get_density()
        final_dens = vols[1].get_density()

        file = '/home/yonatangat/Projects/scattering_tomo/scattering_tomography_src/120keV_spectrum.yaml'
        with open(file, 'r') as f:
            cfg = yaml.load(f)
        spectrum = np.array(cfg['SPECTRUM'])
        # spectrum[0:40] = 0
        spectrum = spectrum / np.sum(spectrum)
        spectrum = list(spectrum)
        # spectrum = np.zeros((150,))
        # spectrum[59] = 1
        # spectrum = list(spectrum)
        if GT_atten_coeff_full is None:
            GT_atten_coeff = construct_atten_coeff_vol(GT_vol, spectrum,
                                                       os.path.join(main_dir, mass_atten_dir))
        else:
            GT_atten_coeff = GT_atten_coeff_full
        init_atten_coeff = construct_atten_coeff_vol(vols[0], spectrum,
                                                     os.path.join(main_dir, mass_atten_dir))
        final_atten_coeff = construct_atten_coeff_vol(vols[1], spectrum,
                                                      os.path.join(main_dir, mass_atten_dir))

        plt.figure()
        plt.suptitle(f'Full Volume mass attenuation - '+ f'number of multiplexed sources {num_srcs[ind]}' if num_srcs[ind] > 0 else 'Single source' )

        gt_frac_mass = GT_atten_coeff.flatten()[points2plot]
        initial_frac_mass = init_atten_coeff.flatten()[points2plot]
        final_frac_mass = final_atten_coeff.flatten()[points2plot]

        max_y = np.max([np.max(initial_frac_mass), np.max(final_frac_mass), np.max(gt_frac_mass)])

        plt.subplot(1, 2, 1)

        scatter_init = plt.scatter(gt_frac_mass, initial_frac_mass)
        plt.xlabel('Ground Truth')
        plt.ylabel('Linear Initialization')
        max_val = np.max([np.max(gt_frac_mass), np.max(initial_frac_mass)])
        plt.plot(np.linspace(0, max_val, GT_vol.get_density().size),
                 np.linspace(0, max_val, GT_vol.get_density().size), 'r--')
        plt.xlim(left=-0.02, right=max_y * 1.05)
        plt.ylim(bottom=-0.02, top=max_y * 1.05)

        plt.subplot(1, 2, 2)

        scatter_final = plt.scatter(gt_frac_mass, final_frac_mass)
        plt.xlabel('Ground Truth')
        plt.ylabel('final volume')
        max_val = np.max([np.max(gt_frac_mass), np.max(final_frac_mass)])
        plt.plot(np.linspace(0, max_val, GT_vol.get_density().size),
                 np.linspace(0, max_val, GT_vol.get_density().size), 'r--')
        plt.xlim(left=-0.02, right=max_y * 1.05)
        plt.ylim(bottom=-0.02, top=max_y * 1.05)

        plt.figure()
        plt.suptitle(f'Full Volume density Scatter Plot - '+ f'number of multiplexed sources {num_srcs[ind]}' if num_srcs[ind] > 0 else 'Single source')

        gt_frac_mass = GT_dens.flatten()[points2plot]
        initial_frac_mass = init_dens.flatten()[points2plot]
        final_frac_mass = final_dens.flatten()[points2plot]

        max_y = np.max([np.max(initial_frac_mass), np.max(final_frac_mass), np.max(gt_frac_mass)])

        plt.subplot(1, 2, 1)

        scatter_init = plt.scatter(gt_frac_mass, initial_frac_mass)
        plt.xlabel('Ground Truth')
        plt.ylabel('Linear Initialization')
        max_val = np.max([np.max(gt_frac_mass), np.max(initial_frac_mass)])
        plt.plot(np.linspace(0, max_val, GT_vol.get_density().size),
                 np.linspace(0, max_val, GT_vol.get_density().size), 'r--')
        plt.xlim(left=-0.02, right=max_y * 1.05)
        plt.ylim(bottom=-0.02, top=max_y * 1.05)

        plt.subplot(1, 2, 2)

        scatter_final = plt.scatter(gt_frac_mass, final_frac_mass)
        plt.xlabel('Ground Truth')
        plt.ylabel('final volume')
        max_val = np.max([np.max(gt_frac_mass), np.max(final_frac_mass)])
        plt.plot(np.linspace(0, max_val, GT_vol.get_density().size),
                 np.linspace(0, max_val, GT_vol.get_density().size), 'r--')
        plt.xlim(left=-0.02, right=max_y * 1.05)
        plt.ylim(bottom=-0.02, top=max_y * 1.05)

        elements = {0: ('H', [0], [0]), 1: ('O', [5, 6, 7], [1]), 2: ('Ca', [13, 18], [2, 3])}
        for el_ind, el in enumerate(vols[0]._active_elements):
            if el == 'P':
                continue
            if el == 'Ca':
                el_ind = 2
            gt_frac_mass = (np.sum(GT_vol._fractional_mass[elements[el_ind][1], :, :, :], axis=0) * GT_dens).flatten()[points2plot]
            initial_frac_mass = (np.sum(vols[0]._fractional_mass[elements[el_ind][2], :, :, :], axis=0) * init_dens).flatten()[points2plot]
            final_frac_mass = (np.sum(vols[1]._fractional_mass[elements[el_ind][2], :, :, ], axis=0) * final_dens).flatten()[points2plot]

            max_y = np.max([np.max(initial_frac_mass), np.max(final_frac_mass), np.max(gt_frac_mass)])

            fig = plt.figure()

            plt.suptitle(f'Full Volume {el} Scatter Plot')

            plt.subplot(1,2,1)

            scatter_init = plt.scatter(gt_frac_mass, initial_frac_mass)
            plt.xlabel('Ground Truth')
            plt.ylabel('Linear Initialization')
            max_val = np.max([np.max(gt_frac_mass), np.max(initial_frac_mass)])
            plt.plot(np.linspace(0, max_val, GT_vol.get_density().size),
                       np.linspace(0, max_val, GT_vol.get_density().size), 'r--')
            plt.xlim(left=-0.02, right=max_y * 1.05)
            plt.ylim(bottom=-0.02, top=max_y * 1.05)

            plt.subplot(1,2,2)

            scatter_final = plt.scatter(gt_frac_mass, final_frac_mass)
            plt.xlabel('Ground Truth')
            plt.ylabel('final volume')
            max_val = np.max([np.max(gt_frac_mass), np.max(final_frac_mass)])
            plt.plot(np.linspace(0, max_val, GT_vol.get_density().size),
                       np.linspace(0, max_val, GT_vol.get_density().size), 'r--')
            plt.xlim(left=-0.02, right=max_y * 1.05)
            plt.ylim(bottom=-0.02, top=max_y * 1.05)

plt.show()

a=1