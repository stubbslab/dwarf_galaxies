import numpy as np
import math
from readInPotentialInterpolator import readInPotentialInterpolator
from doSphericalBootstrapSampling import doBootstrapOnMCMC
import time

if __name__ == '__main__':
    master_start = time.time()
    galaxy = 'fornax'
    halos_to_compute = ['nfw_iso','nfw_aniso','burkert_iso','burkert_aniso']
    halos_to_compute = ['nfw_iso', 'burkert_iso']
    halo_types = ['nfw','nfw','burkert','burkert']
    halo_types = ['nfw', 'burkert']
    halos_iso_aniso = ['iso','aniso','iso','aniso']
    halo_iso_aniso = ['iso', 'iso']
    fixed_params_set = [{'sigsqr_rr_inf_to_0_rat':1.0, 'r_sigsqr_rr0':1000.0, 'alpha_sigsqr_rr':1.0, 'r_beta0':1000.0, 'alpha_beta':1.0, 'gamma_for_beta_inf':0.0},
                        {'sigsqr_rr_inf_to_0_rat':1.0, 'r_sigsqr_rr0':1000.0, 'alpha_sigsqr_rr':1.0, 'alpha_beta':1.0},
                        {'sigsqr_rr_inf_to_0_rat':1.0, 'r_sigsqr_rr0':1000.0, 'alpha_sigsqr_rr':1.0, 'r_beta0':1000.0, 'alpha_beta':1.0, 'gamma_for_beta_inf':0.0},
                        {'sigsqr_rr_inf_to_0_rat':1.0, 'r_sigsqr_rr0':1000.0, 'alpha_sigsqr_rr':1.0, 'alpha_beta':1.0}
                       ]
    fixed_params_set = [{'sigsqr_rr_inf_to_0_rat':1.0, 'r_sigsqr_rr0':1000.0, 'alpha_sigsqr_rr':1.0, 'r_beta0':1000.0, 'alpha_beta':1.0, 'gamma_for_beta_inf':0.0},
                        {'sigsqr_rr_inf_to_0_rat':1.0, 'r_sigsqr_rr0':1000.0, 'alpha_sigsqr_rr':1.0, 'r_beta0':1000.0, 'alpha_beta':1.0, 'gamma_for_beta_inf':0.0},
                     ]
    param_sets_to_measure = [['h_x_center', 'h_z_center', 'sigsqr_rr_0'], ['h_x_center', 'h_z_center', 'sigsqr_rr_0', 'gamma_for_beta_inf', 'r_beta0'],
                         ['h_x_center', 'h_z_center', 'sigsqr_rr_0'], ['h_x_center', 'h_z_center', 'sigsqr_rr_0', 'gamma_for_beta_inf', 'r_beta0'] ]
    param_sets_to_measure = [ ['h_x_center', 'h_z_center', 'sigsqr_rr_0'], ['h_x_center', 'h_z_center', 'sigsqr_rr_0']]


    extra_iterator_save_number = 2

    w_replace = 1
    n_MCMC_iterations = 4000
    n_interim_MCMC_iters_to_save = 1000
    n_MCMC_chains = 2
    n_randomizations = 1
    target_n_sky_pix = 500 #500
    mask_R_bins = 100
    mask_z_bins = 100
    max_sky_angle_deg = None
    n_los_bins = 21 #21
    include_morph_prob = 1
    include_kinem_prob = 1
    smallest_max_val_to_be_on_curve = n_MCMC_iterations / 500
    extra_save_str = 'nfw_burkert_spherical_iso_aniso_test_'
    rs_M_box_size = [500.0, 1.0 * 10.0 ** 9.0]
    pop_independent_params = ['M', 'rs', 'halo_center']
    n_ignore_frac = 0.2

    print ('about to start bootstrapping')

    rand_bootstrap_of_data_from_MCMC = doBootstrapOnMCMC(galaxy, halos_to_compute, halo_types, param_sets_to_measure,
                                                                  w_replace = w_replace, n_MCMC_iterations = n_MCMC_iterations, n_MCMC_chains = n_MCMC_chains,
                                                                  fixed_params_set = fixed_params_set, n_ignore_frac = n_ignore_frac,
                                                                  n_randomizations = n_randomizations, n_interim_iters_to_save  = n_interim_MCMC_iters_to_save,
                                                                  target_n_sky_pix = target_n_sky_pix, mask_R_bins = mask_R_bins, mask_z_bins = mask_z_bins, max_sky_angle_deg = max_sky_angle_deg, n_los_bins = n_los_bins,
                                                                  extra_iterator_save_number = extra_iterator_save_number, extra_save_str = extra_save_str,
                                                                  smallest_max_val_to_be_on_curve = smallest_max_val_to_be_on_curve, include_morph_prob = include_morph_prob, include_kinem_prob = include_kinem_prob ,
                                                                  curve_param_box_size = rs_M_box_size, pop_independent_params = pop_independent_params)
    print ('results: ')
    print (rand_bootstrap_of_data_from_MCMC)
    master_end = time.time()
    print ('Ended bootstrapping.  Took:' + str(master_end - master_start) + 's')
