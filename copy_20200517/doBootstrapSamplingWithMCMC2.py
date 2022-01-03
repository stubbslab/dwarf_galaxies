import numpy as np
import math
from readInPotentialInterpolator import readInPotentialInterpolator 
from doBootstrapSampling import doBootstrapOnMCMC
import time

if __name__ == '__main__':
    master_start = time.time() 
    galaxy = 'fornax'
    halos_to_compute = ['nfw_obl','nfw_pro','cored_obl','cored_pro','burkert_obl','burkert_pro']
    halo_types = ['nfw','nfw','cored','cored','burkert','burkert']
    halos_prolate_or_oblate = ['oblate','prolate','oblate','prolate','oblate','prolate']
    specified_param_ranges_set = [{'el':[0.1, 1.0]}, {'el':[1.0, 10.0]},
                                  {'el':[0.1, 1.0]}, {'el':[1.0, 10.0]},
                                  {'el':[0.1, 1.0]}, {'el':[1.0, 10.0]}
                                 ]
    #halos_prolate_or_oblate = ['oblate','prolate']
    #halos_to_compute = ['burkert_obl','burkert_pro']
    #halo_types = ['burkert','burkert']
    #specified_param_ranges_set = [
    #                              {'el':[0.1, 1.0]}, {'el':[1.0, 10.0]}
    #                             ]
    #halos_to_compute = ['burkert_pro']
    nfw_interp = readInPotentialInterpolator ('nfw')
    cored_interp = readInPotentialInterpolator ('cored')
    burkert_interp = readInPotentialInterpolator ('burkert')
    sech_disk_interp = readInPotentialInterpolator('sech_disk')
    halo_interp_array = [nfw_interp,nfw_interp,cored_interp,cored_interp,burkert_interp,burkert_interp]
    #halo_interp_array = [burkert_interp,burkert_interp]
    disk_interp_array = [sech_disk_interp for interp in halo_interp_array]
    w_replace = 1
    n_MCMC_iterations = 1000
    n_MCMC_chains = 4
    n_randomizations = 2
    outer_step = 2.0
    disk_MCMC_multiplier = 1
    extra_iterator_save_number = 400
    smallest_max_val_to_be_on_curve = n_MCMC_iterations / 500 
    extra_save_str = 'allConfigs_wRot_Replace2'

    rs_M_box_size = [500.0, 1.0 * 10.0 ** 9.0]
    compute_with_disks = [0, 0, 0, 0, 0, 0]
    #compute_with_disks = [0, 0] 
    print 'about to start bootstrapping'  
    rand_bootstrap_of_data_from_MCMC = doBootstrapOnMCMC(galaxy, halos_to_compute, halo_interp_array, disk_interp_array, compute_with_disks, 
                                                         w_replace = w_replace, n_MCMC_iterations = n_MCMC_iterations, n_MCMC_chains = n_MCMC_chains, 
                                                         specified_param_ranges_set = specified_param_ranges_set, 
                                                         n_randomizations = n_randomizations, outer_step = outer_step,
                                                         disk_MCMC_multiplier = disk_MCMC_multiplier, prolate_or_oblate_set = halos_prolate_or_oblate, 
                                                         extra_iterator_save_number = extra_iterator_save_number, extra_save_str = extra_save_str,
                                                         smallest_max_val_to_be_on_curve = smallest_max_val_to_be_on_curve,
                                                         curve_param_box_size = rs_M_box_size)
    print 'results: '
    print rand_bootstrap_of_data_from_MCMC
    master_end = time.time() 
    print 'Ended bootstrapping.  Took:' + str(master_end - master_start) + 's' 
    
