#Write MCMC type procedure for comparing stellar distributions
#Due to computational constraints, assumes a fixed value of
#  halo ellipticity and disk height-to-radius scale.
#Once I get this working, perhaps we can generalize

from ObservedGalaxyStarData import ObservedGalaxyStarData
from PotentialFunctionArray import PotentialFunctionArray
#from generateGalaxyDensityProfile import GalaxyDensityProfile
import SphericalSharedObjectHolder as soh
from PotentialArchive import PotentialArchive
import SurfaceBrightnessProfile as sbp
from DwarfGalaxyParametersStorer import DwarfGalaxyParametersStorer
from DwarfGalDataArchive import DwarfGalDataArchive
#from DwarfGalaxyParametersStorer import getMCMCInfo
from AstronomicalParameterArchive import AstronomicalParameterArchive
from ComputationalArchive import ComputationalArchive
from VisitedDsphMCMCPoint import VisitedDsphMCMCPoint
from StartMCMCParameterStorer import StartMCMCParameterStorer
from logList import logList
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import copy
import random
import math
import numpy as np
import time
from compRZLos import compRZLos
from distributeRandomStars import distributeRandomStars
from ArtificialParamsStorer import ArtificialParamsStorer
from ArtificialStarPositions import ArtificialStarPositions
from GalaxyMask import GalaxyMask
from readInPotentialInterpolator import readInPotentialInterpolator
import cantrips as cant
import os

def runMCMCForDwarfGalaxyProfilesAllPops(galaxy,
                                         disk_funct = 0, halo_funct = 0, withDisk=0, nIterations = 2, saverun = 0, show_best_fit = 0, pop_selection_method = 'none', start_param_index = 0,
                                         outer_step = 1.0, use_artificial_data = 0, generation_params_storer_index = 'None',generation_disk_funct = 0, generation_halo_funct = 0,
                                         mask_R_bins = 200, mask_z_bins = 200, max_sky_angle_deg = None, target_n_sky_pixels = 1250, n_los_bins = 21,
                                         pop_independent_params = ['M', 'rs', 'el','halo_center', 'halo_sym_axis', 'eps', 'Rd', 'lam','disk_center', 'disk_sym_axis'],
                                         max_R_bins = 200, max_z_bins = 200, apply_observation_mask = 1, h_type = '', d_type = '', halo_edge_on = 0, disk_edge_on = 0, apply_rotation = 1,
                                         fixed_disk_angle = None, fixed_halo_angle = None, fixed_params = {}, apply_crowding_mask = 0, extra_save_str = '', star_data = 'real',
                                         start_param_draw_type = 'random', arcmin_limits_los = [-50.0, 50.0], extra_dir = '', specified_param_ranges = {}, prolate_or_oblate = 'either',
                                         include_morph_prob = 1, include_kinem_prob = 1, n_iters_to_save = 200, file_suffix = '.csv',
                                         big_jump_multiplier = 10, big_jump_frequency = 100):

    master_start = time.time()
    astro_archive = AstronomicalParameterArchive ()
    computational_params_archive = ComputationalArchive()
    start_param_storer = StartMCMCParameterStorer(draw = start_param_draw_type )
    dSph_archive = DwarfGalDataArchive(pop_selection_method = pop_selection_method)
    pops = dSph_archive.getPopulations(galaxy)

    if apply_rotation:
        vphi = 'none'
    else:
        vphi = 0.0

    start_param_storer = StartMCMCParameterStorer( draw = start_param_draw_type )
    start_params = [start_param_storer.getStartParameters(start_param_index) for pop in pops]
    print ('start_params = ' + str(start_params))
    if start_param_draw_type == 'best_fit_random':
        n_kinem_params = 6
        new_start_params = [{} for pop in pops]
        for param_key in start_params[0].keys():
            for i in range(len(start_params)):
                if param_key in pop_independent_params:
                    new_start_params[i][param_key] = start_params[i][param_key]
                else:
                    new_start_params[i][param_key] = start_params[i][param_key][i]
        start_params = new_start_params
    #allow user to look for only prolate or oblate halo, fixing start params accordingly
   # (it is up to the user to specify 'el' param range to keep it within the appropriate range)


    #allow user to look for only prolate or oblate halo, fixing start params accordingly
   # (it is up to the user to specify 'el' param range to keep it within the appropriate range)

    if prolate_or_oblate.lower() in ['p','prolate', 'prol']:
        for i in range(len(start_params)):
            if start_params[i]['el'] < 1.0: start_params[i]['el'] = 1.0 / start_params[i]['el']
    elif prolate_or_oblate.lower() in ['o','oblate', 'obl']:
        for i in range(len(start_params)):
            if start_params[i]['el'] > 1.0: start_params[i]['el'] = 1.0 / start_params[i]['el']
    #print ('start_params = ' + str(start_params))
    print ('start_params = ' + str(start_params))

    for i in range(len(pops)):
        for param in specified_param_ranges.keys():
            if start_params[i][param] < specified_param_ranges[param][0]: start_params[i][param] = specified_param_ranges[param][0]
            elif start_params[i][param] > specified_param_ranges[param][1]: start_params[i][param] = specified_param_ranges[param][0]
        for param in fixed_params:
            start_params[i][param] = fixed_params[param]

    deg_to_rad = astro_archive.getDegToRad()
    deg_limits_R,deg_limits_z = dSph_archive.getObservationBounds([galaxy,'dummy_pop_var'],return_unit = 'degrees')

    if star_data is 'real':
        star_data = []
        for i in range(len(pops)):
            pop = pops[i]
            star_data = star_data + [ObservedGalaxyStarData([galaxy,pop], pop_selection_method = pop_selection_method)]

    popGalPairs = []
    start_dsph_parameter_storers  = []

    for i in range(len(pops)):
        pop = pops[i]
        popGalPairs = popGalPairs + [[galaxy,pop]]


        #We presently don't want to have the mass parameter varied
        #start_dsph_parameter_storers = start_dsph_parameter_storers + [
        #    DwarfGalaxyParametersStorer([galaxy,pop], withDisk=withDisk, el=start_params['el'], lam = start_params['lam'],
        #                                M = start_params['M'], rs = start_params['rs'], phi = start_params['phi'], theta = start_params['theta'],
        #                                Rd = start_params['Rd'], eps = start_params['eps'], a = start_params['a'], b= start_params['b'],
        #                                sigsqr = star_data[i].sigSqr)
        #     ]
        #So we do it without the mass
        start_dsph_parameter_storers = start_dsph_parameter_storers + [
            DwarfGalaxyParametersStorer([galaxy,pop], withDisk=withDisk, el=start_params[i]['el'], lam = start_params[i]['lam'],
                                        halo_sym_axis = start_params[i]['halo_sym_axis'], halo_center = start_params[i]['halo_center'], halo_edge_on = halo_edge_on,
                                        M = start_params[i]['M'], rs = start_params[i]['rs'], Rd = start_params[i]['Rd'], eps = start_params[i]['eps'],
                                        disk_sym_axis = start_params[i]['disk_sym_axis'], disk_center = start_params[i]['disk_center'], disk_edge_on = disk_edge_on, fixed_disk_angle = fixed_disk_angle, fixed_halo_angle = fixed_halo_angle,
                                        sigsqr_RR = start_params[i]['sigsqr_RR'], omega_phi = start_params[i]['omega_phi'], specified_param_ranges = specified_param_ranges, stellar_data = star_data[i])
             ]

    if use_artificial_data:
        if generation_params_storer_index == 'None':
            print ('generation_params_storer was not entered.  Using real data.')
        else:
            art_params_storer = ArtificialParamsStorer()
            generation_params_storer = art_params_storer.getStorer(generation_params_storer_index)
            if not generation_disk_funct and generation_params_storer.lam == lam:
                generation_disk_funct = disk_funct
            if not generation_halo_funct and generation_params_storer.el == el:
                generation_halo_funct = halo_funct

            for i in range(len(pops)):
                single_generation_params_storer = generation_params_storer
                single_generation_params_storer.sigsqr = start_dsph_parameter_storers[i].sigsqr
                art_stellar_positions = ArtificialStarPositions(len(star_data[i].proj_x), single_generation_params_storer,
                                                                                       max_R_bins, max_z_bins, outer_step,
                                                                                       arcmin_limits_R, arcmin_limits_z, arcmin_limits_los,
                                                                                       generation_disk_funct, generation_halo_funct, withDisk )
                star_data[i].proj_x, star_data[i].proj_y = art_stellar_positions.getStarPositions()
                star_data[i].proj_x = np.array(star_data[i].proj_x)
                star_data[i].proj_y = np.array(star_data[i].proj_y)
                print ('For population ' + pops[i] + ', we have ' + str(len(star_data[i].proj_x)) + ' artificial stars. ')

    potential_archive = PotentialArchive()
    if not disk_funct:
        print ('Creating the disk_funct from file: ' + disk_file)
        disk_funct = PotentialFunctionArray(disk_file)
    if not halo_funct:
        print ('Creating the halo_funct from file: ' + halo_file)
        halo_funct = PotentialFunctionArray(halo_file)
    #print ('Halo and disk interpolating functions created. ')

    c = start_dsph_parameter_storers[0].c
    dist = start_dsph_parameter_storers[0].dist
    rs = start_dsph_parameter_storers[0].rs
    Rd = start_dsph_parameter_storers[0].Rd
    gamma = astro_archive.getGamma()
    (param_indeces_to_vary, param_functions) =  start_dsph_parameter_storers[0].getMCMCInfo()
    param_indeces_to_fix = {}
    for param in fixed_params.keys():
        index = start_dsph_parameter_storers[0].param_order_array.index(param)
        param_indeces_to_fix[index] = fixed_params[param]
        print ('[index, param] = ' + str([index, param] ))

    pop_independent_param_indeces = [start_dsph_parameter_storers[0].param_order_array.index(param) for param in pop_independent_params]

    #print ('Generating shared objects (including sky mask)...')
    if (max_sky_angle_deg) is None: max_sky_angle_deg = np.max(np.abs([deg_limits_R ,deg_limits_z]))
    mask_R = np.linspace(-max_sky_angle_deg, max_sky_angle_deg, mask_R_bins)
    mask_z = mask_R[:]
    #v_los_to_calc = np.linspace(0.0, max_v_los, n_los_vel_bins)
    #shared_object_holder = ssoh.SharedObjects(start_dsph_parameter_storers[0].dist, mask_R, mask_z, np.linspace(-1.0, 1.0, n_los_bins), target_n_sky_pixels, max_sky_angle_deg * deg_to_rad, v_los_to_calc, compute_mask = apply_observation_mask)
    shared_object_holder = soh.SharedObjects(start_dsph_parameter_storers[0].dist, mask_R, mask_z, np.linspace(-1.0, 1.0, n_los_bins), target_n_sky_pixels, max_sky_angle_deg * deg_to_rad, compute_mask = apply_observation_mask)


    current_parameter_storers = start_dsph_parameter_storers
    current_brightness_profiles = [sbp.SurfaceBrightnessProfile(storer, shared_object_holder, h_type, d_type, disk_interpolating_function=disk_funct, halo_interpolating_function=halo_funct, include_morph_prob = include_morph_prob, include_kinem_prob = include_kinem_prob ) for storer in current_parameter_storers ]

    #current_observation_mask = dSph_archive.getObservationMask(popGalPairs[0],R * rs / dist * 1.0 / deg_to_rad, z * rs / dist * 1.0 / deg_to_rad )
    #Rmesh_degrees, zmesh_degrees = np.meshgrid(R * (rs / dist) * (180.0 / math.pi), z * (rs / dist) * (180.0 / math.pi))
    #universal_observation_mask = GalaxyMask(Rmesh_degrees, zmesh_degrees, galaxy, mask_types = ['n_vel_meas'])
    #for storer in current_parameter_storers:
    #    new_profile = SurfaceBrightnessProfile(R, z, los_bins, gamma, storer, disk_interpolating_function = disk_funct,
    #                                           halo_interpolating_function = halo_funct, observation_mask = universal_observation_mask.final_mask)
    #    if apply_crowding_mask:
    #        specific_crowding_mask = GalaxyMask(Rmesh_degrees, zmesh_degrees, 'fornax', mask_types = ['crowding'], probability_profile = new_profile)
    #        print ('crowding_mask min = ' + str(np.min(specific_crowding_mask.final_mask) ) )
    #        print ('crowding_mask max = ' + str(np.max(specific_crowding_mask.final_mask) ))
    #        new_profile.normalizeSurfaceProfile(specific_crowding_mask.final_mask, renormalizing = 1)
    #    current_brightness_profiles = current_brightness_profiles + [new_profile]

    current_log_likelihood = 0.0
    for i in range(len(current_brightness_profiles)):
        #print 'current_likelihood for pop ' + popGalPairs[i][1] + ' ' + str(current_brightness_profiles[i].sumSurfaceBrightness(star_data[i].proj_x,star_data[i].proj_y))

        current_log_likelihood = current_log_likelihood + current_brightness_profiles[i].getStarProbabilityValues(star_data[i].proj_x, star_data[i].proj_y, star_data[i].corr_Vhel)

    #Note that the parameters we vary as part of our MCMC algorithm of interest in varying MCMC
    visited_points = [VisitedDsphMCMCPoint(current_parameter_storers, current_log_likelihood,1)]

    if saverun:
        save_dir = computational_params_archive.getMCMCOutDir() + extra_dir
        fixed_param_str = '_'.join(fixed_params.keys())
        file_name_prefix = 'MCMC_out_' + extra_save_str + h_type + '_start_' + str(start_param_index) + '_' + galaxy + '_simul_' + 'PSM_' + str(pop_selection_method) + 'mask_' + str(apply_crowding_mask) + '_N_'
        if use_artificial_data and generation_params_storer_index != 'None':
            file_name_prefix = 'artificial_params_' + str(generation_params_storer_index) + '_' + file_name_prefix
        file_name =  file_name_prefix + '0' + file_suffix
        #header = 'dist, c ,M, el, rs , phi , theta , h_xhat, h_yhat, h_zhat, h_x_center, h_y_center, h_z_center, sigsqr_RR, omega_phi, dispersion_b, lam, Rd, eps, a, b, d_xhat, d_yhat, d_zhat, d_x_center, d_y_center, d_z_center, nVisits, logLikelihood'
        header = 'gal, halo, disk, pops, pop_select, mask?, n_mask_bins_R, n_mask_bins_z, sky_angle, target_n_sky_pix, n_los_bins, dist, c, M, rs, el, hsym_x, hsym_y, hsym_z, phi, theta, hc_x, hc_y, hc_z, ' + ','.join(cant.flattenListOfLists([[disp_str + '_' + pop for disp_str in ['sigsqr_RR', 'omega_phi', 'dispersion_b']] for pop in pops])) + ', lam, Rd, eps, a, b, d_xhat, d_yhat, d_zhat, d_x_center, d_y_center, d_z_center, nVisits, logLikelihood'

    #Begin MCMC loop
    for i in range(nIterations):
        #print ('[target_n_sky_pixels, n_los_bins] = ' + str([target_n_sky_pixels, n_los_bins] ))
        start = time.time()
        jump_multiplier = 1
        if i % big_jump_frequency == big_jump_frequency - 1:
            jump_multiplier = big_jump_multiplier
        current_parameter_sets = []
        test_parameter_sets = []
        for storer in current_parameter_storers:
            current_parameter_sets = current_parameter_sets + [storer.getParameterSet()]
            test_parameter_sets = test_parameter_sets + [storer.getParameterSet()[:]]
        #rot_index = storer.whereEqualsomega_phistorer.'omega_phi'
        rot_index = 2

        print ('Here!!')
        print ('test_parameter_sets = ' + str(test_parameter_sets))
        print ('param_indeces_to_vary = ' + str(param_indeces_to_vary))
        print ('param_functions = ' + str(param_functions))
        #Note that the indeces we vary are underlying properties of potential, and so must be the same for each population
        for j in range(len(param_indeces_to_vary)):
            index = param_indeces_to_vary[j]
            varying_function = param_functions[j]
            #If parameter is in list of parameters to keep fixed, don't move it
            if index in param_indeces_to_fix.keys():
                new_parameter_value = param_indeces_to_fix[index]
                new_parameter_values = [new_parameter_value for pop in pops]
            elif index in pop_independent_param_indeces:
                print ('[index, varying_function, j] = ' + str([index, varying_function, j]))
                new_parameter_value = varying_function(test_parameter_sets[0][index], jump_multiplier)
                #print ('Assigning new value  ' + str(new_parameter_value ) + ' for shared param index ' + str(index) + ' corresponding to variable ' + start_dsph_parameter_storers[0].param_order_array[index])
                new_parameter_values = [new_parameter_value for pop in pops]
            #Otherwise, vary it
            else:
                new_parameter_values = [varying_function(test_parameter_set[index], jump_multiplier) for test_parameter_set in test_parameter_sets]


            for k in range(len(pops)):
                test_parameter_sets[k][index] = new_parameter_values[k]


        #test_parameter_set = [ random.gauss(current_parameter_set[j],gauss_sampling_width[j]) if j in param_indeces_to_vary else current_parameter_set[j] for j in range(len(current_parameter_set)) ]

        test_parameter_storers = []
        test_brightness_profiles = []
        test_log_likelihood = 0.0
        for j in range(len(pops)):
            pop = pops[j]
            test_param_set = test_parameter_sets[j]
            test_param_storer = DwarfGalaxyParametersStorer([galaxy,pop], withDisk, *test_param_set, halo_edge_on = halo_edge_on, disk_edge_on = disk_edge_on, fixed_disk_angle = fixed_disk_angle, fixed_halo_angle = fixed_halo_angle, specified_param_ranges = specified_param_ranges, stellar_data = star_data[j])
            test_parameter_storers = test_parameter_storers + [test_param_storer]
            test_brightness_profile = sbp.SurfaceBrightnessProfile(test_param_storer, shared_object_holder, h_type, d_type, disk_interpolating_function = disk_funct, halo_interpolating_function = halo_funct, include_morph_prob = include_morph_prob, include_kinem_prob = include_kinem_prob)
            test_brightness_profiles = test_brightness_profiles + [test_brightness_profile]
            prob_of_pop = test_brightness_profile.getStarProbabilityValues(star_data[j].proj_x, star_data[j].proj_y, star_data[j].corr_Vhel)
            test_log_likelihood = test_log_likelihood + prob_of_pop

        log_rel_likelihood = test_log_likelihood - current_log_likelihood
        #print 'The log_likelihood of the new configuration is: ' + str(test_log_likelihood)
        #print 'The log_relative likelihood for the new configuration is: ' + str(log_rel_likelihood)

        if log_rel_likelihood > math.log(random.random()):
            #print 'Moving to new point. '
            visited_points = visited_points + [VisitedDsphMCMCPoint(test_parameter_storers, test_log_likelihood, 1)]
            current_brightness_profiles = test_brightness_profiles
            current_parameter_storers = test_parameter_storers
            current_log_likelihood = test_log_likelihood
        else:
            #print 'Not moving to another point.  '
            visited_points[-1].n_visits = visited_points[-1].n_visits + 1

        if saverun and (i+1) % n_iters_to_save == n_iters_to_save - 1:
            if os.path.isfile(save_dir + file_name):
                os.remove(save_dir + file_name)
            file_name = file_name_prefix + str(i) + file_suffix
            #header = 'dist, M, omega_phi, sigsqr ,el, rs , phi , theta , h_xhat, h_yhat, h_zhat, h_x_center, h_y_center, h_z_center, c, lam, Rd, eps, a, b, d_xhat, d_yhat, d_zhat, d_x_center, d_y_center, d_z_center, nVisits, logLikelihood'
            print ('Saving run so far to file: ' + save_dir + file_name)
            parameter_result_array = []
            for j in range(len(visited_points)):
                point = visited_points[j]
                storers = point.parameter_storers

                #For those parameters that are actually arrays (eg, the halo and disk symmetry axes),
                # we need to flatten them for writing to a file.
                # save_params =               'gal,  halo,      disk,          pops,          pop_select,           mask?,                  n_mask_bins_R, n_mask_bins_z, sky_angle,        target_n_sky_pix,    n_los_bins,   dist,            c,            M,            rs,     el,              hsym_x, hsym_y, hsym_z,       phi, theta,                        hc_x, hc_y, hc_z, ' + ','.join(cant.flattenListOfLists([[disp_str + '_' + pop for disp_str in ['sigsqr_RR', 'omega_phi', 'dispersion_b']] for pop in pops])) + ', nVisits, logLikelihood'
                flattened_params_to_write = [galaxy, h_type, d_type] + ['/'.join(pops), pop_selection_method, apply_observation_mask, mask_R_bins,   mask_z_bins,   max_sky_angle_deg, target_n_sky_pixels, n_los_bins, storers[0].dist, storers[0].c, storers[0].M, storers[0].rs, storers[0].el] + storers[0].halo_sym_axis + [storers[0].phi , storers[0].theta] + storers[0].halo_center + cant.flattenListOfLists([[storer.sigsqr_RR, storer.omega_phi, storer.dispersion_b] for storer in storers]) + [storers[0].lam, storers[0].Rd, storers[0].eps, storers[0].a, storers[0].b] +  storers[0].disk_sym_axis + storers[0].disk_center
                flattened_params_to_write = [str(save_term) for save_term in flattened_params_to_write ]
                #for param in point.parameter_storer.getParameterSet():
                #    #print 'type(param) = ' + str(type(param))
                #    #if np.size(np.array(param)) == 1: flattened_params_to_write = flattened_params_to_write + [param]
                #    #else:  flattened_params_to_write = flattened_params_to_write + param

                parameter_result_array = parameter_result_array + [flattened_params_to_write + [str(point.n_visits)] + [str(point.log_likelihood)] ]

            parameter_result_array = np.array(parameter_result_array)
            np.savetxt(save_dir + file_name, parameter_result_array, delimiter=", ",header = header, fmt = '%s')

        end = time.time()
        if i % 50 == 1-1:
            print ('Took ' + str(end - start) + ' seconds for iteration ' + str(i) + '. ')

    parameter_result_array = []
    max_visits_index = 0
    max_visits = 0
    for j in range(len(visited_points)):
        point = visited_points[j]
        storers = point.parameter_storers

        #For those parameters that are actually arrays (eg, the halo and disk symmetry axes),
        # we need to flatten them for writing to a file.
        # save_params =    'halo, disk, gal, pops, pop_select, step, WD, obs_mask, dist, omega_phi, sigsqr, c, M, rs, hsym_x, hsym_y, hsym_z, phi, theta, hc_x, hc_y, hc_z, Rd, el, eps, hsym_x, hsym_y, hsym_z, a, b, hc_x, hc_y, hc_z, lam, nVisits, logLikelihood'
        flattened_params_to_write = [galaxy, h_type, d_type] + ['/'.join(pops), pop_selection_method, apply_observation_mask, mask_R_bins, mask_z_bins, max_sky_angle_deg, target_n_sky_pixels, n_los_bins, storers[0].dist, storers[0].c, storers[0].M, storers[0].rs, storers[0].el] + storers[0].halo_sym_axis + [storers[0].phi , storers[0].theta] + storers[0].halo_center + cant.flattenListOfLists([[storer.sigsqr_RR, storer.omega_phi, storer.dispersion_b] for storer in storers]) + [storers[0].lam,         storers[0].Rd, storers[0].eps, storers[0].a, storers[0].b] +  storers[0].disk_sym_axis + storers[0].disk_center
        flattened_params_to_write = [str(save_term) for save_term in flattened_params_to_write ]
        #for param in point.parameter_storer.getParameterSet():
        #    #print 'type(param) = ' + str(type(param))
        #    #if np.size(np.array(param)) == 1: flattened_params_to_write = flattened_params_to_write + [param]
        #    #else:  flattened_params_to_write = flattened_params_to_write + param

        parameter_result_array = parameter_result_array + [flattened_params_to_write + [str(point.n_visits)] + [str(point.log_likelihood)] ]
        if point.n_visits > max_visits:
            max_visits = point.n_visits
            max_visits_index = j

    parameter_result_array = np.array(parameter_result_array)
    #print 'max_visits = ' + str(max_visits)
    #print 'max_visits index = ' + str(max_visits_index)
    max_point = visited_points[max_visits_index]
    max_parameter_storers = max_point.parameter_storers
    max_visit_parameters = max_parameter_storers[0].getParameterSet()

    print ('max_visit_parameters = ')
    print (max_visit_parameters)

    if show_best_fit:

        best_fit_brightness_profile = sbp.SurfaceBrightnessProfile(R_best, z_best, los_bins, gamma, max_parameter_storer,
                                                               disk_interpolating_function = disk_funct,
                                                               halo_interpolating_function = halo_funct,
                                                               observation_mask = universal_observation_mask.final_mask)
        #print star_data.proj_x,star_data.proj_y
        pop_of_interest = pops[2] # can be 'MP', 'IM', 'MR'
        colormap='b' if pop_of_interest == 'MP' else 'g' if pop_of_interest == 'IM' else 'r'
        ylims = (-60,60) #ylimits from Walker plot
        xlims=(60,-60) #xlimits from Walker plot
        fig=plt.figure(1,figsize=(9,7))
        proj_x_one_pop = star_data[2].proj_x
        proj_y_one_pop = star_data[2].proj_y
        print ('NStars is ' + str(len(corr_dec_one_pop)) )
        plt.scatter(corr_ra_one_pop * 60.0,corr_dec_one_pop * 60.0 ,s=4.,color=colormap)
        zmesh,Rmesh=np.meshgrid(best_fit_brightness_profile.zlikeaxis,best_fit_brightness_profile.Rlikeaxis)
        zmesh = zmesh*180.0/math.pi*rs_best/dist_best*60.0
        Rmesh = Rmesh*180.0/math.pi*rs_best/dist_best*60.0
        best_fit_mesh = best_fit_brightness_profile.surfaceProbDensity
        log_levels=logList(np.max(np.abs(best_fit_mesh)) * 0.001, np.max(np.abs(best_fit_mesh)),20)
        CS=plt.contour(Rmesh,zmesh,best_fit_mesh,norm=LogNorm(),levels=log_levels)
        CB_contour=plt.colorbar(shrink=0.8,extend='both',format='%.2e')

        matplotlib.rcParams['xtick.direction'] = 'out'
        matplotlib.rcParams['ytick.direction'] = 'out'
        #plt.ylim(ylims)
        #plt.xlim(xlims)
        plt.xlabel('ra sep')
        plt.ylabel('dec sep')
        plt.title('Fornax Stars Position on Sky')
        pltFileDir=computational_params_archive.getPlotDir()
        #Uncomment the following line if you want to save the file
        #plt.savefig(pltFileDir + gal_of_interest + '_' + pop_of_interest + '_metal_cuts_' + str((metallicity_cuts[gal_of_interest])[0]) + '_' + str((metallicity_cuts[gal_of_interest])[1]) + '.png')
        #plt.savefig(pltFileDir + gal_of_interest + '_' + pop_of_interest + '_full_pop_division' + '.png')
        plt.show()


    if saverun:
        save_dir = computational_params_archive.getMCMCOutDir() + extra_dir
        fixed_param_str = '_'.join(fixed_params.keys())
        fixed_disk_str = 'None' if fixed_disk_angle is None else str(np.around(fixed_disk_angle / math.pi, 5))
        fixed_halo_str = 'None' if fixed_halo_angle is None else str(np.around(fixed_halo_angle / math.pi, 5))
        file_name = 'MCMC_out_' + extra_save_str + h_type + '_' + d_type + '_WD_' + str(withDisk) + '_' + prolate_or_oblate + '_start_' + str(start_param_index) + '_' + galaxy + '_simul_' + 'PSM_' + str(pop_selection_method) + '_HEdge_' + str(halo_edge_on) + '_DEdge_' + str(disk_edge_on) + '_fixDA_' + fixed_disk_str + 'pi_fixHA_' + fixed_halo_str + 'pi_' + fixed_param_str + '_crowdMsk_' + str(apply_crowding_mask) + '_N_' + str(nIterations) + '.csv'
        if use_artificial_data and generation_params_storer_index != 'None':
            file_name = 'artificial_params_' + str(generation_params_storer_index) + '_' + file_name
        #header = 'dist, M, omega_phi, sigsqr ,el, rs , phi , theta , h_xhat, h_yhat, h_zhat, h_x_center, h_y_center, h_z_center, c, lam, Rd, eps, a, b, d_xhat, d_yhat, d_zhat, d_x_center, d_y_center, d_z_center, nVisits, logLikelihood'
        header = 'gal, halo, disk, pops, pop_select, mask?, n_mask_bins_R, n_mask_bins_z, sky_angle, target_n_sky_pix, n_los_bins, dist, c, M, rs, el, hsym_x, hsym_y, hsym_z, phi, theta, hc_x, hc_y, hc_z, ' + ','.join(cant.flattenListOfLists([[disp_str + '_' + pop for disp_str in ['sigsqr_RR', 'omega_phi', 'dispersion_b']] for pop in pops])) + ', lam, Rd, eps, a, b, d_xhat, d_yhat, d_zhat, d_x_center, d_y_center, d_z_center, nVisits, logLikelihood'
        print ('Saving run to file: ' + save_dir + file_name)
        print ('header = ' + str(header))
        print ('parameter_result_array = ' + str(parameter_result_array))
        np.savetxt(save_dir + file_name, parameter_result_array, delimiter=", ",header = header, fmt = '%s')

    else:
        return parameter_result_array



if __name__ == '__main__':
    population = ['fornax','dummy']
    n_iters = 10000 #19999
    halo_type = 'nfw'
    disk_type = 'sech_disk'
    pop_selection_method = 'metal_rigid'
    withDisk = 0
    halo_funct = readInPotentialInterpolator(halo_type)
    disk_funct = readInPotentialInterpolator(disk_type)
    start_param_index = 0
    mask_R_bins = 100 # 200
    mask_z_bins = 100 # 200
    max_sky_angle_deg = None
    n_iters_to_save = 100
    target_n_sky_pixels = 250 #1000
    n_los_bins = 11 #41
    include_morph_prob = 1
    include_kinem_prob = 1
    halo_edge_on = 1
    disk_edge_on = 1
    #n_los_vel_bins = 11 #16
    #max_v_los = 40 # ?
    #n_sig_rr_to_calc_v_los = 3.0 # ?
    apply_observation_mask = 1
    #fixed_params = {'omega_phi':0.0}
    #fixed_params = {'sigsqr_rr_inf_to_0_rat':1.0, 'r_sigsqr_rr0':1000.0, 'alpha_sigsqr_rr':1.0, 'r_beta0':1000.0, 'alpha_beta':1.0, 'gamma_for_beta_inf':0.0, 'omega_phi':0.0}
    fixed_params = {'omega_phi':0.0, 'el':1.0, }
    specified_param_ranges = {'el':[1.0, 10.0]}
    pop_independent_params = ['M', 'rs', 'el', 'halo_center', 'halo_sym_axis', 'eps', 'Rd', 'lam', 'disk_center', 'disk_sym_axis']
    extra_save_str = '_ellipsoidal' + '_includeMorph' + str(include_morph_prob) + '_includeKinemProb' + str(include_kinem_prob)
    #extra_save_str = ''
    runMCMCForDwarfGalaxyProfilesAllPops(population[0], withDisk = withDisk, disk_funct = disk_funct, halo_funct = halo_funct, h_type = halo_type, d_type = disk_type,
                                                  nIterations = n_iters, saverun = 1, show_best_fit = 0, start_param_index = start_param_index, pop_selection_method = pop_selection_method,
                                                  mask_R_bins = mask_R_bins, mask_z_bins = mask_z_bins, max_sky_angle_deg = max_sky_angle_deg, target_n_sky_pixels = target_n_sky_pixels,
                                                  n_los_bins = n_los_bins, #n_los_vel_bins = n_los_vel_bins, max_v_los = max_v_los,
                                                  halo_edge_on = halo_edge_on, disk_edge_on = disk_edge_on,
                                                  apply_observation_mask = apply_observation_mask, include_morph_prob = include_morph_prob, include_kinem_prob = include_kinem_prob,
                                                  specified_param_ranges = specified_param_ranges,
                                                  fixed_params = fixed_params, pop_independent_params = pop_independent_params, n_iters_to_save = n_iters_to_save, extra_save_str = extra_save_str  )
