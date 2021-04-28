import pickle
from copy import deepcopy
import matplotlib.pyplot as plt
import numpy as np
from copy import copy
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.io import loadmat
import csv

avogadro = 6.02214076e23



def plot_slice(ax, gt_vol, initial_vol, final_vol, vmin=0, vmax=1.7, slice=40):

    im0 = ax[0].imshow(gt_vol._density[:, :, slice], vmin=vmin, vmax=vmax)
    im1 = ax[1].imshow(initial_vol._density[:, :, slice], vmin=vmin, vmax=vmax)
    im2 = ax[2].imshow(final_vol._density[:, :, slice], vmin=vmin, vmax=vmax)
    im = (im0, im1, im2)
    for i in range(len(ax)):
        divider = make_axes_locatable(ax[i])
        cax = divider.append_axes("right", size="5%", pad=0.05)
        plt.colorbar(im[i], cax=cax)

def plot_slice_element(ax, vol, el, slice=40):
    vol_densities = vol._density * vol._fractional_mass
    el_vol_density = vol_densities[final_vol._active_elements[el][0], :, :, slice]

    im = ax.imshow(el_vol_density)

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(im, cax=cax)


def plot_scatter_el(gt_vol, final_vol, el, ax, skip=100):
    final_vol_densities = final_vol._density * final_vol._fractional_mass
    el_final_vol_density = final_vol_densities[final_vol._active_elements[el][0], :, :, :]

    total_densities_GT = gt_vol._density * gt_vol._fractional_mass
    el_gt_density = total_densities_GT[gt_vol._active_elements[el][0], :, :, :]

    air_vol = final_vol._air_voxels
    air_vol_mask = air_vol.astype(bool)

    el_no_air = el_final_vol_density[np.logical_not(air_vol_mask)]
    el_no_air_GT = el_gt_density[np.logical_not(air_vol_mask)]

    ax.scatter(el_no_air[::skip], el_no_air_GT[::skip], facecolors='none', edgecolors='b')
    max_val = np.max(np.hstack((el_no_air, el_no_air_GT)))
    ax.set_xlim(0, max_val)
    ax.set_ylim(0, max_val)
    ax.plot([0, max_val], [0, max_val])
    ax.set_aspect('equal', adjustable='box')


def plot_scatter_density(gt_vol, final_vol, ax, skip=100, color='b', marker='o', facecolors='b'):
    final_vol_densities = final_vol._density
    total_densities_GT = gt_vol._density

    air_vol = final_vol._air_voxels
    air_vol_mask = air_vol.astype(bool)

    dens_no_air = final_vol_densities[np.logical_not(air_vol_mask)]
    dens_no_air_GT = total_densities_GT[np.logical_not(air_vol_mask)]

    handle = ax.scatter(dens_no_air[::skip], dens_no_air_GT[::skip], color=color, marker=marker, facecolors=facecolors)
    max_val = np.max(np.hstack((dens_no_air, dens_no_air_GT)))
    ax.set_xlim(0, max_val)
    ax.set_ylim(0, max_val)
    ax.plot([0, max_val], [0, max_val])
    ax.set_xlim(0, max_val)
    ax.set_ylim(0, max_val)
    ax.set_aspect('equal', 'box')
    return handle

def plot_scatter_from_vol(gt_vals, vol_vals, air_voxels, ax, skip=100, color='b', marker='o', facecolors='b'):

    air_vol_mask = air_voxels.astype(bool)

    vol_vals_no_air = vol_vals[np.logical_not(air_vol_mask)]
    gt_vals_no_air_GT = gt_vals[np.logical_not(air_vol_mask)]

    handle = ax.scatter(vol_vals_no_air[::skip], gt_vals_no_air_GT[::skip], color=color, marker=marker, facecolors=facecolors)
    max_val = np.percentile(np.hstack((vol_vals.flatten(), gt_vals.flatten())), 99.99)


    ax.set_xlim(0, max_val)
    ax.set_ylim(0, max_val)
    ax.plot([0, max_val], [0, max_val])
    ax.set_xlim(0, max_val)
    ax.set_ylim(0, max_val)
    ax.set_aspect('equal', 'box')

    return handle

def analyze_precision_recall(gt_vol, vol):
    air_vol = vol._air_voxels
    air_vol_mask = air_vol.astype(bool)

    active_elements_in_vol = vol._active_elements
    for el in active_elements_in_vol:

        el_gt_frac = gt_vol._fractional_mass[gt_vol._active_elements[el][0], :, :, :][np.logical_not(air_vol_mask)]
        el_gt_mask = el_gt_frac.astype(bool)

        # Get initial mask:
        el_vol_frac = vol._fractional_mass[active_elements_in_vol[el][0], :, :, :][np.logical_not(air_vol_mask)]
        el_vol_mask = el_vol_frac.astype(bool)

        in_gt_and_in_vol = np.logical_and(el_gt_mask, el_vol_mask)
        recall_rate = np.sum(in_gt_and_in_vol) / np.sum(el_gt_mask)
        precision = np.sum(in_gt_and_in_vol) / np.sum(el_vol_mask)

        print('el {}, recall: {}, precision: {}'.format(el, recall_rate, precision))


def calc_epsilon(gt_vol, vol, elements=False):

    gt_dens = gt_vol._density
    vol_dens = vol._density
    if elements:
        elements_eps = {}
        for ind, el in enumerate(vol.get_active_elements()):
            # elements_dict = {0: ('H', [0], [0]), 1: ('O', [1], [1]), 2: ('Ca', [2], [2])}
            elements_dict = {0: ('H', [0], [0]), 1: ('O', [5, 6, 7], [1]), 2: ('P', [13], [2]), 3: ('Ca', [18], [3])}
            # elements_dict = {0: ('H', [0], [0]), 1: ('O', [5, 6, 7], [1]), 2: ('Ca', [13, 18], [2]), }
            gt_dens = np.sum(gt_vol._fractional_mass[elements_dict[ind][1], :, :, :], axis=0) * gt_vol.get_density()
            vol_dens = vol._fractional_mass[ind, :, :, :] * vol.get_density()
            estimateDiffNorm = np.sum(np.abs(gt_dens - vol_dens))
            GTNorm = np.sum(np.abs(gt_dens))
            epsilonError = estimateDiffNorm / GTNorm
            elements_eps[el] = epsilonError
        return elements_eps
    estimateDiffNorm = np.sum(np.abs(gt_dens - vol_dens))
    GTNorm = np.sum(np.abs(gt_dens))
    epsilonError = estimateDiffNorm / GTNorm
    print('episilon_error: {}'.format(epsilonError))
    return epsilonError

def calc_delta(gt_vol, vol, elements=False):
    gt_dens = gt_vol._density
    vol_dens = vol._density
    if elements:
        elements_delta = {}
        for ind, el in enumerate(vol.get_active_elements()):
            # elements_dict = {0: ('H', [0], [0]), 1: ('O', [1], [1]), 2: ('Ca', [2], [2])}
            elements_dict = {0: ('H', [0], [0]), 1: ('O', [5, 6, 7], [1]), 2: ('P', [13], [2]), 3: ('Ca', [18], [3])}
            # elements_dict = {0: ('H', [0], [0]), 1: ('O', [5, 6, 7], [1]), 2: ('Ca', [13, 18], [2]) }
            gt_dens = np.sum(gt_vol._fractional_mass[elements_dict[ind][1], :, :, :], axis=0) * gt_vol.get_density()
            vol_dens = vol._fractional_mass[ind, :, :, :] * vol.get_density()
            estimateNorm = np.sum(np.abs(vol_dens))
            GTNorm = np.sum(np.abs(gt_dens))
            deltaMass = (estimateNorm - GTNorm) / GTNorm
            elements_delta[el] = deltaMass
        return elements_delta
    estimateNorm = np.sum(np.abs(vol_dens))
    GTNorm = np.sum(np.abs(gt_dens))
    deltaMass = (estimateNorm - GTNorm) / GTNorm
    print('delta_mass: {}'.format(deltaMass))
    return deltaMass

def calc_contrast(vol_before, vol_after, vol_gt, slice, a_val=1.55, b_val=1.05):
    mask_a = vol_gt._density[:, :, slice] == a_val
    mask_b = vol_gt._density[:, :, slice] == b_val

    patch_a_before = vol_before._density[mask_a, slice]
    patch_b_before = vol_before._density[mask_b, slice]

    patch_a_after = vol_after._density[mask_a, slice]
    patch_b_after = vol_after._density[mask_b, slice]

    std_a_before = np.std(patch_a_before)
    std_a_after = np.std(patch_a_after)

    std_b_before = np.std(patch_b_before)
    std_b_after = np.std(patch_b_after)

    i_a_before = np.mean(patch_a_before)
    i_a_after = np.mean(patch_a_after)
    i_b_before = np.mean(patch_b_before)
    i_b_after = np.mean(patch_b_after)

    contrast_before = (i_a_before - i_b_before) / i_b_before
    contrast_after = (i_a_after - i_b_after) / i_b_after

    noise_before = std_b_before / i_b_before
    noise_after = std_b_after / i_b_after

    CNR_before = contrast_before / noise_before
    CNR_after = contrast_after / noise_after

    print('CNR_before: {}, CNR_after: {}'.format(CNR_before, CNR_after))


def correct_gamma(image, gamma=0.5):
    im = image.astype('float') / 255.0
    im_gamma = np.power(im, gamma)
    return im_gamma


def calc_mac_sigma_o_tilde_gt(gt_vol, sigma_type, elements_data):
    vol_densities = gt_vol._density * gt_vol._fractional_mass
    Oxygen_vol_density = vol_densities[gt_vol._active_elements['O'][0], :, :, :]
    Carbon_vol_density = vol_densities[gt_vol._active_elements['C'][0], :, :, :]
    Nitrogen_vol_density = vol_densities[gt_vol._active_elements['N'][0], :, :, :]
    if sigma_type == 'scatter':
        sigma_s_Oxygen = elements_data['O']['compt'] + elements_data['O']['rayl']
        # sigma_s_Carbon = elements_data['C']['compt'] + elements_data['C']['rayl']
        # sigma_s_Nytrogen = elements_data['N']['compt'] + elements_data['N']['rayl']
    elif sigma_type == 'total':
        sigma_s_Oxygen = elements_data['O']['compt'] + elements_data['O']['rayl'] + elements_data['O']['phot']
        # sigma_s_Carbon = elements_data['C']['compt'] + elements_data['C']['rayl'] + elements_data['C']['phot']
        # sigma_s_Nytrogen = elements_data['N']['compt'] + elements_data['N']['rayl'] + elements_data['N']['phot']

    total_sigma = (avogadro / elements_data['O']['A']) * Oxygen_vol_density * sigma_s_Oxygen
                  # (avogadro / elements_data['C']['A']) * Carbon_vol_density * sigma_s_Carbon + \
                  # (avogadro / elements_data['N']['A']) * Nitrogen_vol_density * sigma_s_Nytrogen

    return total_sigma


def calc_mac_sigma_Ca_tilde_gt(gt_vol, sigma_type, elements_data):
    vol_densities = gt_vol._density * gt_vol._fractional_mass
    Calcium_vol_density = vol_densities[gt_vol._active_elements['Ca'][0], :, :, :]
    Phosphorus_vol_density = vol_densities[gt_vol._active_elements['P'][0], :, :, :]
    if sigma_type == 'scatter':
        sigma_s_Calcium = elements_data['Ca']['compt'] + elements_data['Ca']['rayl']
        sigma_s_Phosphorus = elements_data['P']['compt'] + elements_data['P']['rayl']
    elif sigma_type == 'total':
        sigma_s_Calcium = elements_data['Ca']['compt'] + elements_data['Ca']['rayl'] + elements_data['Ca']['phot']
        sigma_s_Phosphorus = elements_data['P']['compt'] + elements_data['P']['rayl'] + elements_data['P']['phot']

    total_sigma = (avogadro / elements_data['Ca']['A']) * Calcium_vol_density * sigma_s_Calcium + \
                  (avogadro / elements_data['P']['A']) * Phosphorus_vol_density * sigma_s_Phosphorus
    return total_sigma


def calc_mac_sigma_H_tilde_gt(gt_vol, sigma_type, elements_data):
    vol_densities = gt_vol._density * gt_vol._fractional_mass
    Hydrogen_vol_density = vol_densities[gt_vol._active_elements['H'][0], :, :, :]
    if sigma_type == 'scatter':
        sigma_s_Hydrogen = elements_data['H']['compt'] + elements_data['H']['rayl']
    elif sigma_type == 'total':
        sigma_s_Hydrogen = elements_data['H']['compt'] + elements_data['H']['rayl'] + elements_data['H']['phot']

    total_sigma = (avogadro / elements_data['H']['A']) * Hydrogen_vol_density * sigma_s_Hydrogen

    return total_sigma


# #####################################################################################################################
def calc_mac_sigma_o_tilde(vol, sigma_type, elements_data):
    vol_densities = vol._density * vol._fractional_mass
    Oxygen_vol_density = vol_densities[vol._active_elements['O'][0], :, :, :]
    if sigma_type == 'scatter':
        sigma_s_Oxygen = elements_data['O']['compt'] + elements_data['O']['rayl']
    elif sigma_type == 'total':
        sigma_s_Oxygen = elements_data['O']['compt'] + elements_data['O']['rayl'] + elements_data['O']['phot']

    total_sigma = (avogadro / elements_data['O']['A']) * Oxygen_vol_density * sigma_s_Oxygen
    return total_sigma


def calc_mac_sigma_Ca_tilde(vol, sigma_type, elements_data):
    vol_densities = vol._density * vol._fractional_mass
    Calcium_vol_density = vol_densities[vol._active_elements['Ca'][0], :, :, :]
    if sigma_type == 'scatter':
        sigma_s_Calcium = elements_data['Ca']['compt'] + elements_data['Ca']['rayl']
    elif sigma_type == 'total':
        sigma_s_Calcium = elements_data['Ca']['compt'] + elements_data['Ca']['rayl'] + elements_data['Ca']['phot']

    total_sigma = (avogadro / elements_data['Ca']['A']) * Calcium_vol_density * sigma_s_Calcium
    return total_sigma


def calc_mac_sigma_H_tilde(vol, sigma_type, elements_data):
    vol_densities = vol._density * vol._fractional_mass
    Hydrogen_vol_density = vol_densities[vol._active_elements['H'][0], :, :, :]
    if sigma_type == 'scatter':
        sigma_s_Hydrogen = elements_data['H']['compt'] + elements_data['H']['rayl']
    elif sigma_type == 'total':
        sigma_s_Hydrogen = elements_data['H']['compt'] + elements_data['H']['rayl'] + elements_data['H']['phot']

    total_sigma = (avogadro / elements_data['H']['A']) * Hydrogen_vol_density * sigma_s_Hydrogen

    return total_sigma


def calc_mac_sigma(vol, sigma_type, elements_data):
    vol_densities = vol._density * vol._fractional_mass

    total_sigma = 0
    for el in elements_data:
        if el not in vol._active_elements:
            continue
        el_vol_density = vol_densities[vol._active_elements[el][0], :, :, :]
        sigma_s = elements_data[el][sigma_type]
        total_sigma = total_sigma + (avogadro / elements_data[el]['A']) * el_vol_density * sigma_s
    return total_sigma

def calc_mac_sigma_scatter(vol, elements_data):
    vol_densities = vol._density * vol._fractional_mass

    total_sigma = 0
    for el in elements_data:
        if el not in vol._active_elements:
            continue
        el_vol_density = vol_densities[vol._active_elements[el][0], :, :, :]
        sigma_s = elements_data[el]['rayl'] + elements_data[el]['compt']
        total_sigma = total_sigma + (avogadro / elements_data[el]['A']) * el_vol_density * sigma_s
    return total_sigma

def calc_mac_sigma_total(vol, elements_data):
    vol_densities = vol._density * vol._fractional_mass

    total_sigma = 0
    for el in elements_data:
        if el not in vol._active_elements:
            continue
        el_vol_density = vol_densities[vol._active_elements[el][0], :, :, :]
        sigma_s = elements_data[el]['rayl'] + elements_data[el]['compt'] + elements_data[el]['phot']
        total_sigma = total_sigma + (avogadro / elements_data[el]['A']) * el_vol_density * sigma_s
    return total_sigma


