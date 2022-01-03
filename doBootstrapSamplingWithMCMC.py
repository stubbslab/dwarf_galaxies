import numpy as np
import math
from readInPotentialInterpolator import readInPotentialInterpolator
from doBootstrapSampling import doBootstrapOnMCMC
import time

if __name__ == '__main__':
    master_start = time.time()
    galaxy = 'fornax'
    halos_to_compute = ['nfw_obl','nfw_pro','burkert_obl','burkert_pro']
    #halos_to_compute = ['nfw_obl', 'burkert_obl']
    halo_types = ['nfw','nfw','burkert','burkert']
    #halo_types = ['nfw', 'burkert']
    disk_types = ['sech_disk','sech_disk','sech_disk','sech_disk']
    #disk_types = ['sech_disk','sech_disk']
    halos_prolate_or_oblate = ['either', 'either','either', 'either']
    #halos_prolate_or_oblate = ['either', 'either']
    specified_param_ranges_set = [{'el':[0.1, 1.0]}, {'el':[1.0, 10.0]},
                                  {'el':[0.1, 1.0]}, {'el':[1.0, 10.0]}
                                 ]
    #specified_param_ranges_set = [{'el':[0.1, 1.0]}, {'el':[1.0, 10.0]}]

    #specified_param_ranges_set = [
    #                              {'el':[0.1, 1.0]}, {'el':[1.0, 10.0]}
    #                             ]
    #halos_to_compute = ['burkert_pro']
    nfw_interp = readInPotentialInterpolator ('nfw')
    burkert_interp = readInPotentialInterpolator ('burkert')
    sech_disk_interp = readInPotentialInterpolator('sech_disk')
    halo_interp_array = [nfw_interp,nfw_interp, burkert_interp,burkert_interp]
    #halo_interp_array = [nfw_interp,nfw_interp]
    disk_interp_array = [sech_disk_interp for interp in halo_interp_array]

    param_sets_to_measure = [['h_x_center', 'h_z_center', 'el', 'phi', 'sigsqr_RR', 'omega_phi'], ['h_x_center', 'h_z_center', 'el', 'phi', 'sigsqr_RR', 'omega_phi'],
                             ['h_x_center', 'h_z_center', 'el', 'phi', 'sigsqr_RR', 'omega_phi'], ['h_x_center', 'h_z_center', 'el', 'phi', 'sigsqr_RR', 'omega_phi'] ]
    #param_sets_to_measure = [['h_x_center', 'h_z_center', 'el', 'phi', 'sigsqr_RR', 'omega_phi'], ['h_x_center', 'h_z_center', 'el', 'phi', 'sigsqr_RR', 'omega_phi']]

    extra_iterator_save_number = 1

    w_replace = 1
    n_MCMC_iterations = 3000
    n_interim_MCMC_iters_to_save = 500
    n_MCMC_chains = 2
    n_randomizations = 1
    target_n_sky_pix = 500 #500
    mask_R_bins = 100
    mask_z_bins = 100
    max_sky_angle_deg = None
    n_los_bins = 21 #21
    include_morph_prob = 1
    include_kinem_prob = 1
    smallest_max_val_to_be_on_curve = n_MCMC_iterations / 1000
    extra_save_str = 'nfw_burkert_elliptical_prolate_oblate_'
    rs_M_box_size = [500.0, 1.0 * 10.0 ** 9.0]
    pop_independent_params = ['M', 'rs', 'halo_center', 'el', 'phi', 'theta', 'eps', 'Rd', 'd_x_center', 'd_z_center', 'a', 'b']

    rs_M_box_size = [500.0, 1.0 * 10.0 ** 9.0]
    compute_with_disks = [0, 0, 0, 0]
    #compute_with_disks = [0, 0]
    print ('about to start bootstrapping')
    rand_bootstrap_of_data_from_MCMC = doBootstrapOnMCMC(galaxy, halos_to_compute, halo_types, disk_types, halo_interp_array, disk_interp_array, compute_with_disks, param_sets_to_measure,
                                                         w_replace = w_replace, n_MCMC_iterations = n_MCMC_iterations, n_MCMC_chains = n_MCMC_chains,
                                                         specified_param_ranges_set = specified_param_ranges_set,
                                                         n_randomizations = n_randomizations, n_interim_iters_to_save = n_interim_MCMC_iters_to_save,
                                                         target_n_sky_pix = target_n_sky_pix, mask_R_bins = mask_R_bins, mask_z_bins = mask_z_bins, max_sky_angle_deg = max_sky_angle_deg, n_los_bins = n_los_bins,
                                                         prolate_or_oblate_set = halos_prolate_or_oblate,
                                                         extra_iterator_save_number = extra_iterator_save_number, extra_save_str = extra_save_str,
                                                         smallest_max_val_to_be_on_curve = smallest_max_val_to_be_on_curve, include_morph_prob = include_morph_prob, include_kinem_prob = include_kinem_prob ,
                                                         curve_param_box_size = rs_M_box_size, pop_independent_params = pop_independent_params)

    print ('results: ')
    print (rand_bootstrap_of_data_from_MCMC)
    master_end = time.time()
    print ('Ended bootstrapping.  Took:' + str(master_end - master_start) + 's')
