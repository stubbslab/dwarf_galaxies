#Write MCMC type procedure for comparing stellar distributions 
#Due to computational constraints, assumes a fixed value of
#  halo ellipticity and disk height-to-radius scale.
#Once I get this working, perhaps we can generalize

from BoylanArtificialGalaxyStarData import BoylanArtificialGalaxyStarData
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

def runMCMCForDwarfGalaxyProfilesBoylanArtificial(halo_number,
                                                  viewer_phi = 0.0, viewer_theta = 0.0, dist_from_sun = None, disk_funct = 0, halo_funct = 0, withDisk = 0,
                                                  nIterations = 2, saverun = 0, show_best_fit = 0, start_param_draw_type = 'random', 
                                                  start_param_index = 0,
                                                  outer_step = 1.0, max_R_bins = 200, max_z_bins = 200,
                                                  arcmin_limits_R = None, arcmin_limits_z = None, arcmin_limits_los = [-50.0, 50.0],
                                                  generation_params_storer_index = 'None', generation_disk_funct = 0, generation_halo_funct = 0,  h_type = '', d_type = '',
                                                  halo_edge_on = 0, disk_edge_on = 0, fixed_disk_angle = None, fixed_halo_angle = None, fixed_params = {},
                                                  specified_param_ranges = {}, prolate_or_oblate = 'either', 
                                                  extra_save_str = '', extra_dir = ''):
    astro_archive = AstronomicalParameterArchive ()
    computational_params_archive = ComputationalArchive() 
    start_param_storer = StartMCMCParameterStorer(draw = start_param_draw_type )
    start_params = start_param_storer.getStartParameters(start_param_index, withDisk)
    
    if prolate_or_oblate.lower() in ['p','prolate', 'prol']:
        if start_params['el'] < 1.0: start_params['el'] = 1.0 / start_params['el'] 
    elif prolate_or_oblate.lower() in ['o','oblate', 'obl']:
        if start_params['el'] > 1.0: start_params['el'] = 1.0 / start_params['el'] 
    print 'start_params = ' + str(start_params)

    print 'specified_param_ranges = ' + str(specified_param_ranges) 
    for param in specified_param_ranges.keys(): 
        if start_params[param] < specified_param_ranges[param][0]: start_params[param] = specified_param_ranges[param][0]
        elif start_params[param] > specified_param_ranges[param][1]: start_params[param] = specified_param_ranges[param][0] 
    
    dSph_archive = DwarfGalDataArchive()
    deg_to_rad = astro_archive.getDegToRad()
    
    pops = ['all'] # Perhaps someday, I will want to include more than one population in these calculations 

    star_data = []
    for i in range(len(pops)):
        pop = pops[i]
        star_data = star_data + [BoylanArtificialGalaxyStarData(halo_number,
                                                                viewer_phi = viewer_phi, viewer_theta = viewer_theta, dist_from_sun = dist_from_sun,
                                                                arcmin_limits_R = arcmin_limits_R, arcmin_limits_z = arcmin_limits_z)] 


    arcmin_limits_R,arcmin_limits_z = [star_data[0].arcmin_limits_R, star_data[0].arcmin_limits_z]
    dist = star_data[0].dist
    print 'arcmin_limits_R = ' + str(arcmin_limits_R)
    print 'dist = ' + str(dist)
    popGalPairs = []
    start_dsph_parameter_storers  = []
    
    for i in range(len(pops)):
        pop = pops[i]
        haloNumPopPairs = popGalPairs + [['boylan' + str(halo_number),pop]]

        
        #We presently don't want to have the mass parameter varied
        #start_dsph_parameter_storers = start_dsph_parameter_storers + [ 
        #    DwarfGalaxyParametersStorer([galaxy,pop], withDisk=withDisk, el=start_params['el'], lam = start_params['lam'],
        #                                M = start_params['M'], rs = start_params['rs'], phi = start_params['phi'], theta = start_params['theta'],
        #                                zeta = start_params['zeta'], eps = start_params['eps'], a = start_params['a'], b= start_params['b'],
        #                                sigsqr = star_data[i].sigSqr)
        #     ]
        #So we do it without the mass
        start_dsph_parameter_storers = start_dsph_parameter_storers + [ 
            DwarfGalaxyParametersStorer(['boylan' + str(halo_number),pop], withDisk=withDisk, dist = star_data[i].dist, el=start_params['el'], lam = start_params['lam'], 
                                        halo_sym_axis = start_params['halo_sym_axis'], halo_center = start_params['halo_center'], halo_edge_on = halo_edge_on,
                                        M = start_params['M'], rs = start_params['rs'], zeta = start_params['zeta'], eps = start_params['eps'],
                                        disk_sym_axis = start_params['disk_sym_axis'], disk_center = start_params['disk_center'], disk_edge_on = disk_edge_on, fixed_disk_angle = fixed_disk_angle, fixed_halo_angle = fixed_halo_angle, 
                                        sigsqr = star_data[i].sigSqr, specified_param_ranges = specified_param_ranges)
             ]
        
    print 'len(start_dsph_parameter_storers = ' + str(len(start_dsph_parameter_storers)) 
        
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
    print 'dist = ' + str(dist) 
    rs = start_dsph_parameter_storers[0].rs
    zeta = start_dsph_parameter_storers[0].zeta 
    gamma = astro_archive.getGamma()
    (param_indeces_to_vary, param_functions) =  start_dsph_parameter_storers[0].getMCMCInfo()
    param_indeces_to_fix = {}
    for param in fixed_params.keys():
        index = start_dsph_parameter_storers[0].param_order_array.index(param)
        param_indeces_to_fix[index] = fixed_params[param] 

    #R,z,los_bins = compRZLos(zeta, zeta * outer_step, outer_step, arcmin_limits_R, arcmin_limits_z, arcmin_limits_los, rs, dist)
    R,z,los_bins = compRZLos(1.0, 1.0 * outer_step, outer_step, arcmin_limits_R, arcmin_limits_z, arcmin_limits_los, rs, dist)
    
    current_parameter_storers = start_dsph_parameter_storers

    current_brightness_profiles = []
    #current_observation_mask = dSph_archive.getObservationMask(popGalPairs[0],R * rs / dist * 1.0 / deg_to_rad, z * rs / dist * 1.0 / deg_to_rad )
    Rmesh_degrees, zmesh_degrees = np.meshgrid(R * (rs / dist) * (180.0 / math.pi), z * (rs / dist) * (180.0 / math.pi)) 
    for storer in current_parameter_storers:
        new_profile = SurfaceBrightnessProfile(R, z, los_bins, gamma, storer, disk_interpolating_function = disk_funct,
                                               halo_interpolating_function = halo_funct)

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
        if i%100 == 100-1:
            #print 'On iteration '  + str(i)
            jump_multiplier = 10
        current_parameter_sets = []
        test_parameter_sets = []
        for storer in current_parameter_storers:
            current_parameter_sets = current_parameter_sets + [storer.getParameterSet()]
            test_parameter_sets = test_parameter_sets + [storer.getParameterSet()[:]]
        
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
                

        #test_parameter_set = [ random.gauss(current_parameter_set[j],gauss_sampling_width[j]) if j in param_indeces_to_vary else current_parameter_set[j] for j in range(len(current_parameter_set)) ] 

        test_parameter_storers = []
        for pair in range(len(haloNumPopPairs)):
            test_parameter_storers = test_parameter_storers + [DwarfGalaxyParametersStorer(haloNumPopPairs [pair], withDisk, *(test_parameter_sets[pair]), 
                                                                                           halo_edge_on = halo_edge_on, disk_edge_on = disk_edge_on, fixed_disk_angle = fixed_disk_angle,fixed_halo_angle = fixed_halo_angle, specified_param_ranges = specified_param_ranges)]

        c_test = test_parameter_storers[0].c
        dist_test = test_parameter_storers[0].dist
        rs_test = test_parameter_storers[0].rs
        zeta_test = test_parameter_storers[0].zeta 

        #R_test,z_test,los_bins_test = compRZLos(zeta_test, zeta_test * outer_step, outer_step, arcmin_limits_R, arcmin_limits_z, arcmin_limits_los, rs_test, dist_test)
        R_test,z_test,los_bins_test = compRZLos(1.0, 1.0 * outer_step, outer_step, arcmin_limits_R, arcmin_limits_z, arcmin_limits_los, rs_test, dist_test)
        #if i % 50 == 50 - 1:
        #    print 'zeta_test = ' + str(zeta_test) + ' so, len(los_bins_test) = ' + str(len(los_bins_test)) 

        test_brightness_profiles = []

        for storer in test_parameter_storers:
            new_profile = SurfaceBrightnessProfile(R_test, z_test, los_bins_test, gamma, storer, disk_interpolating_function = disk_funct,
                                                   halo_interpolating_function = halo_funct)

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
        if i % 200 == 200-1:
            print 'Took ' + str(end - start) + ' seconds for iteration ' + str(i) + '. '
    
    parameter_result_array = []
    max_visits_index = 0
    max_visits = 0
    for j in range(len(visited_points)):
        point = visited_points[j]
    
        #For those parameters that are actually arrays (eg, the halo and disk symmetry axes),
        # we need to flatten them for writing to a file.
        flattened_params_to_write = []
        for param in point.parameter_storer.getParameterSet():
            #print 'type(param) = ' + str(type(param))
            if np.size(np.array(param)) == 1: flattened_params_to_write = flattened_params_to_write + [param]
            else:  flattened_params_to_write = flattened_params_to_write + param
            
        parameter_result_array = parameter_result_array + [flattened_params_to_write + [point.n_visits] + [point.log_likelihood] ]
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
    zeta_best = max_parameter_storer.zeta

    #R_best,z_best,los_bins_best = compRZLos(zeta_best, zeta_best * outer_step, outer_step, arcmin_limits_R, arcmin_limits_z, arcmin_limits_los, rs_best, dist_best)
    R_best,z_best,los_bins_best = compRZLos(1.0, 1.0 * outer_step, outer_step, arcmin_limits_R, arcmin_limits_z, arcmin_limits_los, rs_best, dist_best)
    print 'max_visit_parameters = '
    print max_visit_parameters 
    best_fit_brightness_profile = SurfaceBrightnessProfile(R_best, z_best, los_bins, gamma, max_parameter_storer,
                                                                   disk_interpolating_function = disk_funct,
                                                                   halo_interpolating_function = halo_funct)

    
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
        file_name = 'MCMC_output_' + extra_save_str + h_type + '_' + d_type + '_WD_' + str(withDisk) + '_' + prolate_or_oblate + '_start_' + str(start_param_index) + '_BoylanHalo' + str(halo_number)  + '_HEdge_' + str(halo_edge_on) + '_DEdge_' + str(disk_edge_on) + '_fixDAng_' + str(fixed_disk_angle) + '_fixHAng_'+ str(fixed_halo_angle) + fixed_param_str + '_N_' + str(nIterations) + '.csv'
        
        header = 'dist, M, vphi, sigsqr ,el, rs , phi , theta , h_xhat, h_yhat, h_zhat, h_x_center, h_y_center, h_z_center, c, lam, zeta, eps, a, b, d_xhat, d_yhat, d_zhat, d_x_center, d_y_center, d_z_center, nVisits, logLikelihood'
        print 'Saving run to file: ' + save_dir + file_name  
        np.savetxt(save_dir + file_name, parameter_result_array, delimiter=",",header = header)

    else:
        return parameter_result_array 
        
        
    
if __name__ == '__main__':
    el = 0.5
    lam = 0.2
    population = ['fornax','MP']
    runMCMCForDwarfGalaxyProfiles(el, lam, population, disk_funct = 0, halo_funct = 0, withDisk=0, nIterations = 2, saverun = 0, show_best_fit = 0, pop_selection_method = 'none')
