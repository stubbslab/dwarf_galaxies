#Write MCMC type procedure for comparing stellar distributions 
#Due to computational constraints, assumes a fixed value of
#  halo ellipticity and disk height-to-radius scale.
#Once I get this working, perhaps we can generalize

from ObservedGalaxyStarData import ObservedGalaxyStarData
from PotentialFunctionArray import PotentialFunctionArray
#from generateGalaxyDensityProfile import GalaxyDensityProfile
from PotentialArchive import PotentialArchive 
from SurfaceBrightnessProfile import SurfaceBrightnessProfile 
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

def runMCMCForDwarfGalaxyProfilesAllPops(galaxy,
                                         disk_funct = 0, halo_funct = 0, withDisk=0, nIterations = 2, saverun = 0, show_best_fit = 0, pop_selection_method = 'none', start_param_index = 0,
                                         outer_step = 1.0, use_artificial_data = 0, generation_params_storer_index = 'None',generation_disk_funct = 0, generation_halo_funct = 0,
                                         max_R_bins = 200, max_z_bins = 200, apply_observation_mask = 1, h_type = '', d_type = '', halo_edge_on = 0, disk_edge_on = 0, apply_rotation = 1, 
                                         fixed_disk_angle = None, fixed_halo_angle = None, fixed_params = {}, apply_crowding_mask = 0, extra_save_str = '', star_data = 'real',
                                         start_param_draw_type = 'random', arcmin_limits_los = [-50.0, 50.0], extra_dir = '', specified_param_ranges = {}, prolate_or_oblate = 'either'):

    astro_archive = AstronomicalParameterArchive ()
    computational_params_archive = ComputationalArchive() 
    start_param_storer = StartMCMCParameterStorer(draw = start_param_draw_type )
    
    if apply_rotation:
        vphi = 'none'
    else:
        vphi = 0.0 
    start_params = start_param_storer.getStartParameters(start_param_index, withDisk)
    #allow user to look for only prolate or oblate halo, fixing start params accordingly
   # (it is up to the user to specify 'el' param range to keep it within the appropriate range)
    
    if prolate_or_oblate.lower() in ['p','prolate', 'prol']:
        if start_params['el'] < 1.0: start_params['el'] = 1.0 / start_params['el'] 
    elif prolate_or_oblate.lower() in ['o','oblate', 'obl']:
        if start_params['el'] > 1.0: start_params['el'] = 1.0 / start_params['el'] 
    print 'start_params = ' + str(start_params)
    
    for param in specified_param_ranges.keys(): 
        if start_params[param] < specified_param_ranges[param][0]: start_params[param] = specified_param_ranges[param][0]
        elif start_params[param] > specified_param_ranges[param][1]: start_params[param] = specified_param_ranges[param][0]
    
    dSph_archive = DwarfGalDataArchive(pop_selection_method = pop_selection_method)
    deg_to_rad = astro_archive.getDegToRad()
    arcmin_limits_R,arcmin_limits_z = dSph_archive.getObservationBounds([galaxy,'dummy_pop_var'],return_unit = 'arcmin')
    pops = dSph_archive.getPopulations(galaxy)

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
            DwarfGalaxyParametersStorer([galaxy,pop], withDisk=withDisk, el=start_params['el'], lam = start_params['lam'],
                                        halo_sym_axis = start_params['halo_sym_axis'], halo_center = start_params['halo_center'], halo_edge_on = halo_edge_on,
                                        M = start_params['M'], rs = start_params['rs'], Rd = start_params['Rd'], eps = start_params['eps'], vphi = vphi, 
                                        disk_sym_axis = start_params['disk_sym_axis'], disk_center = start_params['disk_center'], disk_edge_on = disk_edge_on, fixed_disk_angle = fixed_disk_angle, fixed_halo_angle = fixed_halo_angle, 
                                        sigsqr = star_data[i].sigSqr, specified_param_ranges = specified_param_ranges, stellar_data = star_data[i])
             ]
    
    if use_artificial_data:
        if generation_params_storer_index == 'None':
            print 'generation_params_storer was not entered.  Using real data.'
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
                print 'For population ' + pops[i] + ', we have ' + str(len(star_data[i].proj_x)) + ' artificial stars. '
        
    potential_archive = PotentialArchive()
    if not disk_funct: 
        print 'Creating the disk_funct from file: ' + disk_file
        disk_funct = PotentialFunctionArray(disk_file)
    if not halo_funct:
        print 'Creating the halo_funct from file: ' + halo_file
        halo_funct = PotentialFunctionArray(halo_file)
    print 'Halo and disk interpolating functions created. '

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

    #R,z,los_bins = compRZLos(Rd, Rd * outer_step, outer_step, arcmin_limits_R, arcmin_limits_z, arcmin_limits_los, rs, dist)
    R,z,los_bins = compRZLos(1.0, 1.0 * outer_step, outer_step, arcmin_limits_R, arcmin_limits_z, arcmin_limits_los, rs, dist)
    
    current_parameter_storers = start_dsph_parameter_storers

    current_brightness_profiles = []
    #current_observation_mask = dSph_archive.getObservationMask(popGalPairs[0],R * rs / dist * 1.0 / deg_to_rad, z * rs / dist * 1.0 / deg_to_rad )
    Rmesh_degrees, zmesh_degrees = np.meshgrid(R * (rs / dist) * (180.0 / math.pi), z * (rs / dist) * (180.0 / math.pi))
    universal_observation_mask = GalaxyMask(Rmesh_degrees, zmesh_degrees, galaxy, mask_types = ['n_vel_meas'])
    for storer in current_parameter_storers:
        new_profile = SurfaceBrightnessProfile(R, z, los_bins, gamma, storer, disk_interpolating_function = disk_funct,
                                               halo_interpolating_function = halo_funct, observation_mask = universal_observation_mask.final_mask)
        if apply_crowding_mask: 
            specific_crowding_mask = GalaxyMask(Rmesh_degrees, zmesh_degrees, 'fornax', mask_types = ['crowding'], probability_profile = new_profile)
            print 'crowding_mask min = ' + str(np.min(specific_crowding_mask.final_mask) )
            print 'crowding_mask max = ' + str(np.max(specific_crowding_mask.final_mask) )
            new_profile.normalizeSurfaceProfile(specific_crowding_mask.final_mask, renormalizing = 1)
        current_brightness_profiles = current_brightness_profiles + [new_profile]

    current_log_likelihood = 0.0
    for i in range(len(current_brightness_profiles)):
        #print 'current_likelihood for pop ' + popGalPairs[i][1] + ' ' + str(current_brightness_profiles[i].sumSurfaceBrightness(star_data[i].proj_x,star_data[i].proj_y))
 
        current_log_likelihood = current_log_likelihood + current_brightness_profiles[i].sumLogSurfaceBrightness(star_data[i].proj_x, star_data[i].proj_y)
    
    #Note that the parameters we vary as part of our MCMC algorithm of interest in varying MCMC 
    visited_points = [VisitedDsphMCMCPoint(current_parameter_storers[0],current_log_likelihood,1)]
    
    #Begin MCMC loop
    for i in range(nIterations):
        start = time.time()
        jump_multiplier = 1
        if i%1 == 100-1:
            print 'On iteration '  + str(i)
            jump_multiplier = 10
        current_parameter_sets = []
        test_parameter_sets = []
        for storer in current_parameter_storers:
            current_parameter_sets = current_parameter_sets + [storer.getParameterSet()]
            test_parameter_sets = test_parameter_sets + [storer.getParameterSet()[:]] 
        #rot_index = storer.whereEqualsOmegaphistorer.'omegaphi'
        rot_index = 2
        
        #Note that the indeces we vary are underlying properties of potential, and so must be the same for each population 
        for j in range(len(param_indeces_to_vary)):
            index = param_indeces_to_vary[j]
            varying_function = param_functions[j]
            #If parameter is in list of parameters to keep fixed, don't move it
            if index in param_indeces_to_fix.keys():
                new_parameter_value = param_indeces_to_fix[index]
            #Otherwise, vary it 
            else:
                new_parameter_value = varying_function(test_parameter_sets[0][index], jump_multiplier)
            
            for test_parameter_set in test_parameter_sets:
                test_parameter_set[index] = new_parameter_value
                if apply_rotation: 
                    test_parameter_set[rot_index] = 'none'
                else: 
                    test_parameter_set[rot_index] = 0.0 
                

        #test_parameter_set = [ random.gauss(current_parameter_set[j],gauss_sampling_width[j]) if j in param_indeces_to_vary else current_parameter_set[j] for j in range(len(current_parameter_set)) ] 

        test_parameter_storers = []
        for pair in range(len(popGalPairs)):
            test_parameter_storers = test_parameter_storers + [DwarfGalaxyParametersStorer(popGalPairs[pair], withDisk, *(test_parameter_sets[pair]), halo_edge_on = halo_edge_on, disk_edge_on = disk_edge_on, fixed_disk_angle = fixed_disk_angle, fixed_halo_angle = fixed_halo_angle, specified_param_ranges = specified_param_ranges, stellar_data = star_data[pair])]

        c_test = test_parameter_storers[0].c
        dist_test = test_parameter_storers[0].dist
        rs_test = test_parameter_storers[0].rs
        Rd_test = test_parameter_storers[0].Rd 

        #R_test,z_test,los_bins_test = compRZLos(Rd_test, Rd_test * outer_step, outer_step, arcmin_limits_R, arcmin_limits_z, arcmin_limits_los, rs_test, dist_test)
        R_test,z_test,los_bins_test = compRZLos(1.0, 1.0 * outer_step, outer_step, arcmin_limits_R, arcmin_limits_z, arcmin_limits_los, rs_test, dist_test)
        #if i % 50 == 50 - 1:
        #    print 'Rd_test = ' + str(Rd_test) + ' so, len(los_bins_test) = ' + str(len(los_bins_test)) 

        test_brightness_profiles = []

        for storer in test_parameter_storers:
            new_profile = SurfaceBrightnessProfile(R_test, z_test, los_bins_test, gamma, storer, disk_interpolating_function = disk_funct,
                                                   halo_interpolating_function = halo_funct, observation_mask = universal_observation_mask.final_mask)
            if apply_crowding_mask: 
                specific_crowding_mask = GalaxyMask(Rmesh_degrees, zmesh_degrees, 'fornax', mask_types = ['crowding'], probability_profile =  new_profile)
                specific_crowding_mask = GalaxyMask(Rmesh_degrees, zmesh_degrees, 'fornax', mask_types = ['crowding'], probability_profile = new_profile)
                print 'crowding_mask min = ' + str(np.min(specific_crowding_mask.final_mask) )
                print 'crowding_mask max = ' + str(np.max(specific_crowding_mask.final_mask) )
                new_profile.normalizeSurfaceProfile(specific_crowding_mask.final_mask, renormalizing = 1)
            test_brightness_profiles = test_brightness_profiles + [new_profile]
                                                                  
        test_log_likelihood = 0.0
        for k in range(len(test_brightness_profiles)):
            test_log_likelihood = test_log_likelihood + test_brightness_profiles[k].sumLogSurfaceBrightness(star_data[k].proj_x,star_data[k].proj_y)

        log_rel_likelihood = test_log_likelihood - current_log_likelihood
        #print 'The log_likelihood of the new configuration is: ' + str(test_log_likelihood)
        #print 'The log_relative likelihood for the new configuration is: ' + str(log_rel_likelihood) 
        
        if log_rel_likelihood > math.log(random.random()):
            #print 'Moving to new point. '
            visited_points = visited_points + [VisitedDsphMCMCPoint(test_parameter_storers[0],test_log_likelihood,1)]
            current_brightness_profiles = test_brightness_profiles 
            current_parameter_storers = test_parameter_storers
            current_log_likelihood = test_log_likelihood
        else:
            #print 'Not moving to another point.  '  
            visited_points[-1].n_visits = visited_points[-1].n_visits + 1

        end = time.time()
        if i % 100 == 100-1:
            print 'Took ' + str(end - start) + ' seconds for iteration ' + str(i) + '. '
    
    parameter_result_array = []
    max_visits_index = 0
    max_visits = 0
    for j in range(len(visited_points)):
        point = visited_points[j]
        storer = point.parameter_storer
    
        #For those parameters that are actually arrays (eg, the halo and disk symmetry axes),
        # we need to flatten them for writing to a file.
        # save_params =    'halo, disk, gal, pops, pop_select, step, WD, obs_mask, dist, omega_phi, sigsqr, c, M, rs, hsym_x, hsym_y, hsym_z, phi, theta, hc_x, hc_y, hc_z, Rd, el, eps, hsym_x, hsym_y, hsym_z, a, b, hc_x, hc_y, hc_z, lam, nVisits, logLikelihood'
        flattened_params_to_write = [h_type, d_type, galaxy] + ['/'.join(pops), pop_selection_method, outer_step, storer.withDisk, apply_observation_mask, storer.dist, storer.omegaphi, storer.sigsqr, storer.c, storer.M, storer.rs] + storer.halo_sym_axis + [storer.phi, storer.theta] + storer.halo_center + [storer.Rd, storer.el, storer.eps] + storer.disk_sym_axis + [storer.a, storer.b] + storer.disk_center + [storer.lam]
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
    max_parameter_storer = max_point.parameter_storer
    max_visit_parameters = max_parameter_storer.getParameterSet()
    
    dist_best = max_parameter_storer.dist
    rs_best = max_parameter_storer.rs
    Rd_best = max_parameter_storer.Rd

    #R_best,z_best,los_bins_best = compRZLos(Rd_best, Rd_best * outer_step, outer_step, arcmin_limits_R, arcmin_limits_z, arcmin_limits_los, rs_best, dist_best)
    R_best,z_best,los_bins_best = compRZLos(1.0, 1.0 * outer_step, outer_step, arcmin_limits_R, arcmin_limits_z, arcmin_limits_los, rs_best, dist_best)
    print 'max_visit_parameters = '
    print max_visit_parameters 
    best_fit_brightness_profile = SurfaceBrightnessProfile(R_best, z_best, los_bins, gamma, max_parameter_storer,
                                                                   disk_interpolating_function = disk_funct,
                                                                   halo_interpolating_function = halo_funct,
                                                                   observation_mask = universal_observation_mask.final_mask)
    if apply_crowding_mask: 
        best_fit_crowding_mask = GalaxyMask(Rmesh_degrees, zmesh_degrees, 'fornax', mask_types = ['crowding'], probability_profile =  best_fit_brightness_profile)
        print 'best_fit_mask min = ' + str(np.min(best_fit_crowding_mask.final_mask) )
        print 'best_fit_mask max = ' + str(np.max(best_fit_crowding_mask.final_mask) )
        best_fit_brightness_profile.normalizeSurfaceProfile(best_fit_crowding_mask.final_mask, renormalizing = 1)
    
    if show_best_fit:
        
        #print star_data.proj_x,star_data.proj_y
        pop_of_interest = pops[2] # can be 'MP', 'IM', 'MR'
        colormap='b' if pop_of_interest == 'MP' else 'g' if pop_of_interest == 'IM' else 'r'
        ylims = (-60,60) #ylimits from Walker plot
        xlims=(60,-60) #xlimits from Walker plot
        fig=plt.figure(1,figsize=(9,7))
        proj_x_one_pop = star_data[2].proj_x
        proj_y_one_pop = star_data[2].proj_y
        print 'NStars is ' + str(len(corr_dec_one_pop))
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
        #header = 'dist, M, omegaphi, sigsqr ,el, rs , phi , theta , h_xhat, h_yhat, h_zhat, h_x_center, h_y_center, h_z_center, c, lam, Rd, eps, a, b, d_xhat, d_yhat, d_zhat, d_x_center, d_y_center, d_z_center, nVisits, logLikelihood'
        header = 'halo, disk, gal, pops, pop_select, step, WD, obs_mask, dist, omega_phi, sigsqr, c, M, rs, hsym_x, hsym_y, hsym_z, phi, theta, hc_x, hc_y, hc_z, Rd, el, eps, dsym_x, dsym_y, dsym_z, a, b, dc_x, dc_y, dc_z, lam, nVisits, logLikelihood'
        print ('Saving run to file: ' + save_dir + file_name) 
        print ('header = ' + str(header)) 
        print ('parameter_result_array = ' + str(parameter_result_array))  
        np.savetxt(save_dir + file_name, parameter_result_array, delimiter=", ",header = header, fmt = '%s')

    else:
        return parameter_result_array 
        
        
    
if __name__ == '__main__':
    el = 0.5
    lam = 0.2
    population = ['fornax','MP']
    runMCMCForDwarfGalaxyProfiles(el, lam, population, disk_funct = 0, halo_funct = 0, withDisk=0, nIterations = 2, saverun = 0, show_best_fit = 0, pop_selection_method = 'none')
