

#Do a sequence of the MCMC Series algorithms, with varying run conditions.
#You can change various parameters (present set listed in commented line right above parameter_serires assignment), and have program execute in sequence
#This script can be run from command line without having to be in python environment (bash$ python doMCMCSeriesAllPops.py)
#This method tries to be as efficient as possible, reading in the various potentialFunctionArrays only once, as that can be a time consuming part of the process.

#from PotentialFunctionArray import PotentialFunctionArray
from PotentialArchive import PotentialArchive
from ComputationalArchive import ComputationalArchive
from SurfaceBrightnessProfile import SurfaceBrightnessProfile
from DwarfGalaxyParametersStorer import DwarfGalaxyParametersStorer
from AstronomicalParameterArchive import AstronomicalParameterArchive
from DwarfGalDataArchive import DwarfGalDataArchive
from ObservedGalaxyStarData import ObservedGalaxyStarData
from readInPotentialInterpolator import readInPotentialInterpolator
from compRZLos import compRZLos
from GalaxyMask import GalaxyMask
import math
import numpy as np

if __name__ == '__main__':
    #parameters at present are: [el, lam, galaxy, number of iterations, population selection method, start parameter index]
    print ('Starting computation of likelihood series. ')
    astro_archive = AstronomicalParameterArchive ()
    compute_params_archive = ComputationalArchive()
    dSph_archive = DwarfGalDataArchive()
    r_half_light = 791.0
    M_star = 10.0 ** 7.39

    saverun = 1
    save_array = []
    results_number = '0'
    gamma = astro_archive.getGamma()

    #                     halo       disk          gal       pops               pop_select      step,  WD, apply_observation_mask,  M                rs                      halo_sym_axis                  halo_center                        Rd                     el     eps         disk_sym_axis                disk_center                        lam

    halo_index = 0
    gal_index = halo_index + 1
    pops_index = gal_index + 1
    pop_select_index = pops_index + 1
    step_index = pop_select_index + 1 #step for determining density of profile sampling
    apply_observation_mask_index = withDisk_index + 1
    M_index = apply_observation_mask_index + 1
    rs_index = M_index + 1
    halo_center_index = halo_sym_axis_index + 1
    sigsqr_rr_0_index = halo_center_index + 1
    sigsqr_rr_inf_index = sigsqr_rr_0_index + 1
    r_sigsqr_rr0_index = sigsqr_rr_inf_index + 1
    alpha_sigsqr_rr_index = r_sigsqr_rr0_index + 1
    gamma_for_beta_inf_index = alpha_sigsqr_rr_index + 1
    r_beta0_index = gamma_for_beta_inf_index + 1
    alpha_beta_index = r_beta0_index + 1




    #                     halo       gal       pops               pop_select      step,  apply_observation_mask,  M        rs       halo_center      sigsqr_rr_0 sigsqr_rr_inf  r_sigsqr_rr0  alpha_sigsqrr_rr gamma_for_beta_inf  r_beta0  alpha_beta
    #Gaussian-on-one-variable peaks
    parameter_series = [ ['nfw',     'fornax', ['MR','IM','MP'],  'metal_rigid',  2.0,   0,                       9.0e+09, 3.8e+03, [0.0, 0.0, 0.0], 100.0,      100.0,         1.0,          2.0,             1.0,                1.0,     2.0,],
                         ['nfw',     'fornax', ['MR','IM','MP'],  'metal_rigid',  2.0,   0,                       9.0e+09, 3.8e+03, [0.0, 0.0, 0.0], 100.0,      100.0,         1.0,          2.0,             1.0,                1.0,     2.0,],
                         ['cored',   'fornax', ['MR','IM','MP'],  'metal_rigid',  2.0,   0,                       2.0e+09, 3.0e+03, [0.0, 0.0, 0.0], 100.0,      100.0,         1.0,          2.0,             1.0,                1.0,     2.0,],
                         ['cored',   'fornax', ['MR','IM','MP'],  'metal_rigid',  2.0,   0,                       2.0e+09, 3.0e+03, [0.0, 0.0, 0.0], 100.0,      100.0,         1.0,          2.0,             1.0,                1.0,     2.0,],
                         ['burkert', 'fornax', ['MR','IM','MP'],  'metal_rigid',  2.0,   0,                       1.2e+09, 5.6e+03, [0.0, 0.0, 0.0], 100.0,      100.0,         1.0,          2.0,             1.0,                1.0,     2.0,],
                         ['burkert', 'fornax', ['MR','IM','MP'],  'metal_rigid',  2.0,   0,                       2.0e+09, 3.0e+03, [0.0, 0.0, 0.0], 100.0,      100.0,         1.0,          2.0,             1.0,                1.0,     2.0,],
                         ['nfw',     'fornax', ['MR','IM','MP'],  'metal_rigid',  2.0,   0,                       8.0e+09, 3.8e+03, [0.0, 0.0, 0.0], 100.0,      100.0,         1.0,          2.0,             1.0,                1.0,     2.0,],
                         ['nfw',     'fornax', ['MR','IM','MP'],  'metal_rigid',  2.0,   0,                       8.0e+09, 5.2e+03, [0.0, 0.0, 0.0], 100.0,      100.0,         1.0,          2.0,             1.0,                1.0,     2.0,],
                         ['burkert', 'fornax', ['MR','IM','MP'],  'metal_rigid',  2.0,   0,                       3.5e+09, 4.9e+03, [0.0, 0.0, 0.0], 100.0,      100.0,         1.0,          2.0,             1.0,                1.0,     2.0,],
                         ['burkert', 'fornax', ['MR','IM','MP'],  'metal_rigid',  2.0,   0,                       2.5e+09, 5.7e+03, [0.0, 0.0, 0.0], 100.0,      100.0,         1.0,          2.0,             1.0,                1.0,     2.0,],
                         ['cored',   'fornax', ['MR','IM','MP'],  'metal_rigid',  2.0,   0,                       2.0e+09, 3.0e+03, [0.0, 0.0, 0.0], 100.0,      100.0,         1.0,          2.0,             1.0,                1.0,     2.0,],
                         ['cored',   'fornax', ['MR','IM','MP'],  'metal_rigid',  2.0,   0,                       2.0e+09, 3.0e+03, [0.0, 0.0, 0.0], 100.0,      100.0,         1.0,          2.0,             1.0,                1.0,     2.0,],
                         ]


    pot_archive = PotentialArchive()
    unique_halo_types = []
    unique_disk_types = []
    #Create all of the potential files used in the MCMC algorithms.
    # That way, we don't need to reread in a potential file every time it is used.
    for parameter_set in parameter_series:
        if parameter_set[halo_index] not in unique_halo_types:
            print ('assigning halo type ' + str(parameter_set[halo_index]) + ' to unique_halo_types')
            unique_halo_types = unique_halo_types + [parameter_set[halo_index]]
        if parameter_set[disk_index] not in unique_disk_types:
            print ('assigning disk type ' + str(parameter_set[disk_index]) + ' to unique_disk_types')
            unique_disk_types = unique_disk_types + [parameter_set[disk_index]]
    halo_funct_array = {}
    disk_funct_array = {}
    for unique_halo_type in unique_halo_types:
        print ('Loading potential interpolator for halo type ' + unique_halo_type)
        halo_funct_array[unique_halo_type] = readInPotentialInterpolator(unique_halo_type)
    for unique_disk_type in unique_disk_types:
        print ('Loading potential interpolator for disk type ' + unique_disk_type)
        disk_funct_array[unique_disk_type] = readInPotentialInterpolator(unique_disk_type)

    #arcmin_limits_z =50.0
    likelihood_array = []
    save_array = []
    #Now actually compute the likelihood for each of the specified point.
    for parameter_set in parameter_series:
        print ('Computing likelihood for parameter_set: ' + str(parameter_set))
        parameter_storers = []
        surface_profiles = []

        halo_type = parameter_set[halo_index]
        galaxy = parameter_set[gal_index]
        pops  = parameter_set[pops_index]
        if len(pops) == 0:
            pops = dSph_archive.getPopulations(galaxy)
        pop_selection_method = parameter_set[pop_select_index]
        outer_step = parameter_set[step_index]
        apply_observation_mask = parameter_set[apply_observation_mask_index]
        M = parameter_set[M_index]
        rs = parameter_set[rs_index]
        dist = dSph_archive.getDistanceFromSun([galaxy,'dummy_pop_variable'])
        halo_center = parameter_set[halo_center_index]
        sigsqr_rr_0 = parameter_set[sigsqr_rr_0_index]
        sigsqr_rr_inf = parameter_set[sigsqr_rr_inf_index]
        r_sigsqr_rr0 = parameter_set[r_sigsqr_rr0_index]
        alpha_sigsqr_rr = parameter_set[alpha_sigsqr_rr_index]
        gamma_for_beta_inf = parameter_set[gamma_for_beta_inf_index]
        r_beta0 = parameter_set[r_beta0_index]
        alpha_beta = parameter_set[alpha_beta_index]

        #print 'max, min R in arcminutes are: ' + str(max(R_outer)) + ', ' + str(min(R_outer))
        #print 'max, min z in arcminutes are: ' + str(max(z_outer)) + ', ' + str(min(z_outer))
        storer = 0
        log_likelihood = 0.0

        arcmin_limits_R,arcmin_limits_z = dSph_archive.getObservationBounds([galaxy,'dummy_var'], return_unit = 'arcmin')
        #arcmin_limits_R = [-60.0,60.0]
        #arcmin_limits_z = [-60.0,60.0]
        arcmin_limits_los = [-50.0,50.0]
        #R,z,los_bins = compRZLos(zeta, zeta * outer_step, outer_step, arcmin_limits_R, arcmin_limits_z, arcmin_limits_los, rs, dist)
        R,z,los_bins = compRZLos(1.0, 1.0 * outer_step, outer_step, arcmin_limits_R, arcmin_limits_z, arcmin_limits_los, rs, dist)
        Rmesh_degrees, zmesh_degrees = np.meshgrid(R * (rs / dist) * (180.0 / math.pi), z * (rs / dist) * (180.0 / math.pi))
        if apply_observation_mask:
            observation_mask = GalaxyMask(Rmesh_degrees, zmesh_degrees, galaxy, mask_types = ['n_vel_meas']).final_mask
        else:
            observation_mask = np.zeros(np.shape(Rmesh_degrees)) + 1.0
        for pop in pops:
            population = [galaxy,pop]
            star_data = ObservedGalaxyStarData(population, pop_selection_method = pop_selection_method)
            #print 'for population ' + pop + ', len(star_data.corrRa) ' + str(len(star_data.corrRa))
            storer =  DwarfGalaxyParametersStorer(population, withDisk,
                                                  M = M, rs = rs, dist = dist, halo_center = halo_center,
                                                  dispersion_rr_params = [sigsqr_rr_0, sigsqr_rr_inf, r_sigsqr_rr0, alpha_sigsqr_rr], dispersion_beta_params = [gamma_for_beta_inf, r_beta0, alpha_beta ]
                                                  stellar_data = star_data, )
            #print 'storer.sigsqr = ' + str(storer.sigsqr) + ' for population ' + pop
            print ('storer.printContents() = ')
            storer.printContents()

            parameter_storers = parameter_storers + [storer]
            def __init__(self, x, z, los_bins, parameter_storer, halo_type, disk_file=None, halo_file=None, C=1.0, observation_mask = None, masking_fields = [], population = None)
            surface_profile = SurfaceBrightnessProfile(R, z, los_bins, storer, halo_type, observation_mask = observation_mask )
            print ('surface_profile.onSkyInterpolator([[0.0, 0.1, 0.0, 0.1], [0.0, 0.0, 0.1, 0.1]] ) = ' +str(surface_profile.onSkyInterpolator( np.dstack([[0.0, 0.1, 0.0, 0.1], [0.0, 0.0, 0.1, 0.1]]) [0]  )))
            #if pop == 'MR':
            #        #print 'pop == ' + pop
            #        #print 'surface_profile.Rlikeaxis  = '
            #        #print surface_profile.Rlikeaxis
            #        #print 'surface_profile.zlikeaxis  = '
            #        #print surface_profile.zlikeaxis
            #        brightness_mesh = surface_profile.surfaceBrightness
            #        #print 'surface_profile.surfaceBrightness = '
            #        #print surface_profile.surfaceBrightness
            surface_profiles = surface_profiles + [surface_profile]
            log_likelihood = log_likelihood + surface_profile.sumLogSurfaceBrightness(star_data.proj_x, star_data.proj_y)

            #plt.scatter(np.array(star_data.corrRa) * 60.0 ,np.array(star_data.corrDec) * 60.0)
        pop_string = ' '.join(pops)
        halo_center_string = ''
        if halo_center is 'none':
            halo_center_string = 'none'
        else:
            halo_center_string = ' '.join([str(elem) for elem in halo_center])
        #                     halo       gal       pops               pop_select      step, apply_observation_mask,  dist,  omegaphi, sigsqr_rr_0, sigsqr_rr_inf, r_sigsqr_rr0, alpha_sigsqr_rr, gamma_for_beta_inf, r_beta0, alpha_beta, c, M  rs halo_center
        save_array = save_array + [[halo_type,galaxy, '/'.join(pops), pop_selection_method, outer_step, withDisk, apply_observation_mask]
                                    + [storer.dist, storer.omegaphi, storer.sigsqr_rr_0, storer.sigsqr_rr_inf, storer.r_sigsqr_rr0, storer.alpha_sigsqr_rr, storer.gamma_for_beta_inf, storer.r_beta0, storer.alpha_beta, storer.c]
                                    + [storer.M, storer.rs]
                                    + storer.halo_center
                                    + [log_likelihood]
                                   ]

                                    #+ parameter_set[0:halo_sym_axis_index] + [halo_sym_string]
                                    #+ parameter_set[halo_sym_axis_index+1:halo_center_index] + [halo_center_string]
                                    #+ parameter_set[halo_center_index+1:disk_sym_axis_index] + [disk_sym_string]
                                    #+ parameter_set[disk_sym_axis_index+1:disk_center_index] + [disk_center_string]
                                    #+ parameter_set[disk_center_index+1:pops_index] + [pop_string] + parameter_set[pops_index+1:] + [log_likelihood]]
        likelihood_array = likelihood_array + [log_likelihood]


    #print likelihood_array
    #save_array = [parameter_series[i] + [likelihood_array[i]] for i in range(len(likelihood_array)) ]
    #print save_array
    #print np.array(save_array)

    if saverun:
        save_dir = compute_params_archive.getSpotProbabilityDir()
        print ('save_dir = ' + save_dir )
        file_name = 'test_spherical_probabilities_number_' + str(results_number) + '.csv'
        print ('Saving series to ' + save_dir + file_name)
        header = 'halo, gal, pops, pop_select, step, apply_observation_mask, dist, omegaphi, sigsqr_rr_0, sigsqr_rr_inf, r_sigsqr_rr0, alpha_sigsqr_rr, gamma_for_beta_inf, r_beta0, alpha_beta, c, M, rs, halo_center, logLikelihood'
        np.savetxt(save_dir + file_name, np.array(save_array), delimiter = ',',fmt = '%10s',header = header)
