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

def getArtificialStarPositions(N_stars, generation_params_storer, pop_sigsqr, max_R_bins, max_z_bins, outer_step, arcmin_limits_R,arcmin_limits_z, arcmin_limits_los, generation_disk_funct, generation_halo_funct, withDisk, plot_artificial_star_positions = 0):
    print 'Creating artifical stars... '
    start_art = time.time()

    astro_archive = AstronomicalParameterArchive()
    pot_archive = PotentialArchive()

    deg_to_rad = astro_archive.getDegToRad()
    gamma = astro_archive.getGamma()
    
    el = generation_params_storer.el
    lam = generation_params_storer.lam
    rs = generation_params_storer.rs
    dist = generation_params_storer.dist
    zeta = generation_params_storer.zeta 
    
    #generation_disk_file = pot_archive.getDiskFile(lam) 
    #generation_halo_file = pot_archive.getHaloFile(el)

    if not generation_disk_funct:
        generation_disk_funct = PotentialFunctionArray(generation_disk_file)
    if not generation_halo_funct:
        generation_halo_funct = PotentialFunctionArray(generation_halo_file) 

    generation_params_storer.sigsqr = pop_sigsqr

    generation_R,generation_z,generation_los = compRZLos(zeta, zeta * outer_step, outer_step, arcmin_limits_R, arcmin_limits_z, arcmin_limts_los, rs, dist ) 

    print 'Generating surface_profile... '
    generation_surface_profile = SurfaceBrightnessProfile(generation_R,generation_z,generation_los, 
                                                          gamma, generation_params_storer,
                                                          generation_disk_file, generation_halo_file,
                                                          disk_interpolating_function = generation_disk_funct,
                                                          halo_interpolating_function = generation_halo_funct) 
    print 'Distributing random stars... '
    #print 'max_R_bins = '
    #print max_R_bins
    #print 'max_z_bins = '
    #print max_z_bins 
    artificial_star_RAs_in_scale_radii, artificial_star_Decs_in_scale_radii = distributeRandomStars(N_stars, generation_surface_profile, max_R_bins, max_z_bins, withDisk)
    print 'Rescaling RA and Decs of random stars...'
    artificial_star_RAs_in_deg = [RA * rs / dist * 1 / deg_to_rad for RA in artificial_star_RAs_in_scale_radii]
    artificial_star_Decs_in_deg = [Dec * rs / dist * 1 / deg_to_rad for Dec in artificial_star_Decs_in_scale_radii]

    end_art = time.time()
    print 'Took ' + str(end_art - start_art) + 's to generate artificial stars. '

    if plot_artificial_star_positions:
        plt.scatter(artificial_star_RAs_in_deg,artificial_star_Decs_in_deg)
        plt.xlim(-1.0,1.0)
        plt.ylim(-1.0,1.0) 
        plt.show()
    return artificial_star_RAs_in_deg,artificial_star_Decs_in_deg


def runSimulMCMCForDwarfGalaxyProfilesAllPops(galaxy, disk_funct = 0, halo_funct = 0, withDisk=0, pop_selection_method = 'none', start_param_indeces = [0], outer_step = 1.0, apply_observation_mask = 1,
                                              nIterations = 2, saverun = 0, show_best_fit = 0, h_type = '', d_type = '',
                                              use_artificial_data = 0, generation_params_storer_index = 'None', generation_disk_funct = 0, generation_halo_funct = 0, max_R_bins = 200, max_z_bins = 200 ):
    astro_archive = AstronomicalParameterArchive ()
    computational_params_archive = ComputationalArchive() 
    start_param_storer = StartMCMCParameterStorer()
    n_chains = len(start_param_indeces) 
    start_params_array = [start_param_storer.getStartParameters(start_param_index, withDisk) for start_param_index in start_param_indeces]
    dSph_archive = DwarfGalDataArchive()
    deg_to_rad = astro_archive.getDegToRad()
    arcmin_limits_R,arcmin_limits_z = dSph_archive.getObservationBounds([galaxy,'dummy_pop_var'],return_unit = 'arcmin')
    arcmin_limits_los = [-50.0,50.0] 
    pops = dSph_archive.getPopulations(galaxy)
    #print 'pops = '
    #print pops   
    star_data = []
    popGalPairs = []
    for pop in pops:
        popGalPairs = popGalPairs + [[galaxy,pop]]
        star_data = star_data + [ObservedGalaxyStarData([galaxy,pop], pop_selection_method = pop_selection_method)]
    
    start_dsph_parameter_storers_array = [] 
    for start_params in start_params_array:
        start_dsph_parameter_storers  = []
        for i in range(len(pops)):
            pop = pops[i]
            #We presently don't want to have the mass parameter varied
            #start_dsph_parameter_storers = start_dsph_parameter_storers + [ 
            #    DwarfGalaxyParametersStorer([galaxy,pop], withDisk=withDisk, el=start_params['el'], lam = start_params['lam'],
            #                                M = start_params['M'], rs = start_params['rs'], phi = start_params['phi'], theta = start_params['theta'],
            #                                zeta = start_params['zeta'], eps = start_params['eps'], a = start_params['a'], b= start_params['b'],
            #                                sigsqr = star_data[i].sigSqr)
            #     ]
            #So we do it without the mass
            start_dsph_parameter_storers = start_dsph_parameter_storers + [ 
                DwarfGalaxyParametersStorer([galaxy,pop], withDisk=withDisk, el=start_params['el'], lam = start_params['lam'],
                                            rs = start_params['rs'], halo_sym_axis = start_params['halo_sym_axis'], halo_center = start_params['halo_center'],
                                            zeta = start_params['zeta'], eps = start_params['eps'], disk_sym_axis = start_params['disk_sym_axis'], disk_center = start_params['disk_center'],
                                            sigsqr = star_data[i].sigSqr)
                 ]
        start_dsph_parameter_storers_array = start_dsph_parameter_storers_array + [start_dsph_parameter_storers]
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
                star_data[i].corrRa, star_data[i].corrDec = getArtificialStarPositions(len(star_data[i].corrRa), generation_params_storer,
                                                                                         start_dsph_parameter_storers[i].sigsqr,
                                                                                         max_R_bins, max_z_bins, outer_step,
                                                                                         arcmin_limits_R, arcmin_limits_z,arcmin_limits_los,
                                                                                         generation_disk_funct, generation_halo_funct,withDisk )
                star_data[i].corrRa = np.array(star_data[i].corrRa)
                star_data[i].corrDec = np.array(star_data[i].corrDec) 
                print 'For population ' + pops[i] + ', we have ' + str(len(star_data[i].corrRa)) + ' artificial stars. '
        
    potential_archive = PotentialArchive()
    if not disk_funct: 
        print 'Creating the disk_funct from file: ' + disk_file
        disk_funct = PotentialFunctionArray(disk_file)
    if not halo_funct:
        print 'Creating the halo_funct from file: ' + halo_file
        halo_funct = PotentialFunctionArray(halo_file)
    print 'Halo and disk interpolating functions created. '

    cs = [start_dsph_parameter_storers[0].c for start_dsph_parameter_storers in start_dsph_parameter_storers_array]
    dists = [start_dsph_parameter_storers[0].dist for start_dsph_parameter_storers in start_dsph_parameter_storers_array]
    rss = [start_dsph_parameter_storers[0].rs for start_dsph_parameter_storers in start_dsph_parameter_storers_array]
    zetas = [start_dsph_parameter_storers[0].zeta for start_dsph_parameter_storers in start_dsph_parameter_storers_array]
    gamma = astro_archive.getGamma()
    (param_indeces_to_vary, param_functions) =  start_dsph_parameter_storers_array[0][0].getMCMCInfo()

    coordinate_arrays = [compRZLos(zetas[i], zetas[i] * outer_step, outer_step, arcmin_limits_R, arcmin_limits_z, arcmin_limits_los, rss[i], dists[i]) for i in range(n_chains)]

    #R,z,los_bins = compRZLos(zeta, zeta * outer_step, outer_step, arcmin_limits_R, arcmin_limits_z, arcmin_limits_los, rs, dist)
    
    current_parameter_storers_array = start_dsph_parameter_storers_array 

    current_brightness_profiles_array = []
    #current_observation_mask = dSph_archive.getObservationMask(popGalPairs[0],R * rs / dist * 1.0 / deg_to_rad, z * rs / dist * 1.0 / deg_to_rad )
    current_observation_mask_array = [dSph_archive.getObservationMask(popGalPairs[0],coordinate_arrays[i][0] * rss[i] / dists[i] * 1.0 / deg_to_rad, coordinate_arrays[i][1] * rss[i] / dists[i] * 1.0 / deg_to_rad ) for i in range(n_chains)]
    for i in range(n_chains):
        current_parameter_storers = current_parameter_storers_array[i]
        R,z,los_bins = coordinate_arrays[i]
        current_observation_mask = current_observation_mask_array[i]
        current_brightness_profiles = []
        for storer in current_parameter_storers:
            current_brightness_profiles = current_brightness_profiles + [SurfaceBrightnessProfile(R, z, los_bins, gamma, storer,
                                                                                                  disk_interpolating_function = disk_funct,
                                                                                                  halo_interpolating_function = halo_funct,
                                                                                                  observation_mask = current_observation_mask ) ]
        current_brightness_profiles_array = current_brightness_profiles_array + [current_brightness_profiles]

    current_log_likelihood_array = []
    for current_brightness_profiles in current_brightness_profiles_array:
        current_log_likelihood = 0.0
        for i in range(len(current_brightness_profiles)):
            #print 'current_likelihood for pop ' + popGalPairs[i][1] + ' ' + str(current_brightness_profiles[i].sumSurfaceBrightness(star_data[i].corrRa,star_data[i].corrDec))
            current_log_likelihood = current_log_likelihood + current_brightness_profiles[i].sumLogSurfaceBrightness(star_data[i].corrRa,star_data[i].corrDec)
        current_log_likelihood_array = current_log_likelihood_array + [current_log_likelihood]
    
    #Note that the parameters we vary as part of our MCMC algorithm of interest in varying MCMC
    visited_points_array = []
    for i in range(n_chains):
        current_log_likelihood = current_log_likelihood_array[i]
        current_parameter_storers = current_parameter_storers_array[i]
        visited_points = [VisitedDsphMCMCPoint(current_parameter_storers[0],current_log_likelihood,1)]
        visited_points_array = visited_points_array + [visited_points]
    
    #Begin MCMC loop
    for i in range(nIterations):
        start = time.time()
        jump_multiplier = 1
        if i%100 == 100-1:
            print 'On iteration '  + str(i)
            jump_multiplier = 10
        #current_parameter_sets_array = []
        #test_parameter_sets_array = []
        for chain in range(n_chains):
            print 'i = ' + str(i) + ' and chain = ' + str(chain)
            current_parameter_storers = current_parameter_storers_array[chain]
            visited_points = visited_points_array[chain]
            current_log_likelihood = current_log_likelihood_array[chain]
            current_parameter_sets = []
            test_parameter_sets = []
            for storer in current_parameter_storers:
                #storer.printContents() 
                current_parameter_sets = current_parameter_sets + [storer.getParameterSet()]
                test_parameter_sets = test_parameter_sets + [storer.getParameterSet()[:]]
            #print 'current_log_likelihood = ' + str(current_log_likelihood) 
            #Note that the indeces we vary are underlying properties of potential, and so must be the same for each population 
            for j in range(len(param_indeces_to_vary)):
                index = param_indeces_to_vary[j]
                varying_function = param_functions[j]
                new_parameter_value = varying_function(test_parameter_sets[0][index], jump_multiplier)
            
                for test_parameter_set in test_parameter_sets:
                    test_parameter_set[index] = new_parameter_value 
                

            #test_parameter_set = [ random.gauss(current_parameter_set[j],gauss_sampling_width[j]) if j in param_indeces_to_vary else current_parameter_set[j] for j in range(len(current_parameter_set)) ] 

            test_parameter_storers = []
            #print 'len(test_parameter_sets) = ' + str(len(test_parameter_sets))
            #print 'len(popGalPairs) = ' + str(len(popGalPairs))
            for pair in range(len(popGalPairs)):
                #print 'pair = ' + str(pair) 
                test_parameter_storers = test_parameter_storers + [DwarfGalaxyParametersStorer(popGalPairs[pair], withDisk, *(test_parameter_sets[pair]))]

            c_test = test_parameter_storers[0].c
            dist_test = test_parameter_storers[0].dist
            rs_test = test_parameter_storers[0].rs
            zeta_test = test_parameter_storers[0].zeta 

            R_test,z_test,los_bins_test = compRZLos(zeta_test, zeta_test * outer_step, outer_step, arcmin_limits_R, arcmin_limits_z, arcmin_limits_los, rs_test, dist_test)
            #if i % 50 == 50 - 1:
            #    print 'zeta_test = ' + str(zeta_test) + ' so, len(los_bins_test) = ' + str(len(los_bins_test)) 

            test_brightness_profiles = []

            test_observation_mask = dSph_archive.getObservationMask(popGalPairs[0],R_test * rs_test / dist_test * 1.0 / deg_to_rad, z_test * rs_test / dist_test * 1.0 / deg_to_rad )
            for storer in test_parameter_storers:
                test_brightness_profiles = test_brightness_profiles + [SurfaceBrightnessProfile(R_test, z_test, los_bins_test, gamma, storer, 
                                                                                            disk_interpolating_function = disk_funct,
                                                                                            halo_interpolating_function = halo_funct,
                                                                                            observation_mask = test_observation_mask) ]
                                                                  
            test_log_likelihood = 0.0
            for k in range(len(test_brightness_profiles)):
                test_log_likelihood = test_log_likelihood + test_brightness_profiles[k].sumLogSurfaceBrightness(star_data[k].corrRa,star_data[k].corrDec)

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
            current_parameter_storers_array[chain] = current_parameter_storers 
            current_brightness_profiles_array[chain] = current_brightness_profiles
            current_log_likelihood_array[chain] = current_log_likelihood 
            visited_points_array[chain] = visited_points 
        end = time.time()
        if i % 100 == 100-1:
            print 'Took ' + str(end - start) + ' seconds for iteration ' + str(i) + '. '

    #Now we want to write each chain to a separate file
    for chain in range(n_chains):
    
        parameter_result_array = []
        max_visits_index = 0
        max_visits = 0
        visited_points = visited_points_array[chain] 
        for j in range(len(visited_points)):
            point = visited_points[j]
    
            #For those parameters that are actually arrays (eg, the halo and disk symmetry axes),
            # we need to flatten them for writing to a file.
            flattened_params_to_write = []
            for param in point.parameter_storer.getParameterSet():
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

        R_best,z_best,los_bins_best = compRZLos(zeta_best, zeta_best * outer_step, outer_step, arcmin_limits_R, arcmin_limits_z, arcmin_limits_los, rs_best, dist_best)
        print 'max_visit_parameters = '
        print max_visit_parameters 
        best_fit_observation_mask = dSph_archive.getObservationMask(popGalPairs[0],R_best * rs_best / dist_best * 1.0 / deg_to_rad, z_best * rs_best / dist_best * 1.0 / deg_to_rad )
        best_fit_brightness_profile = SurfaceBrightnessProfile(R_best, z_best, los_bins, gamma, max_parameter_storer,
                                                                   disk_interpolating_function = disk_funct,
                                                                   halo_interpolating_function = halo_funct,
                                                                   observation_mask = best_fit_observation_mask) 
    
    
        if show_best_fit:
        
            #print star_data.corrRa,star_data.corrDec
            pop_of_interest = pops[2] # can be 'MP', 'IM', 'MR'
            colormap='b' if pop_of_interest == 'MP' else 'g' if pop_of_interest == 'IM' else 'r'
            ylims = (-60,60) #ylimits from Walker plot
            xlims=(60,-60) #xlimits from Walker plot
            fig=plt.figure(1,figsize=(9,7))
            corr_ra_one_pop = star_data[2].corrRa
            corr_dec_one_pop = star_data[2].corrDec
            print 'NStars is ' + str(len(corr_dec_one_pop))
            plt.scatter(corr_ra_one_pop * 60.0,corr_dec_one_pop * 60.0 ,s=4.,color=colormap)
            zmesh,Rmesh=np.meshgrid(best_fit_brightness_profile.zlikeaxis,best_fit_brightness_profile.Rlikeaxis)
            zmesh = zmesh*180.0/math.pi*rs_best/dist_best*60.0
            Rmesh = Rmesh*180.0/math.pi*rs_best/dist_best*60.0
            best_fit_mesh = best_fit_brightness_profile.surfaceProbDensity 
            log_levels=logList(np.min(np.abs(best_fit_mesh)),np.max(np.abs(best_fit_mesh)),20)
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
            save_dir = computational_params_archive.getMCMCOutDir()
            start_param_index = start_param_indeces[chain] 
            file_name = 'MCMC_output_data_' + h_type + '_' + d_type + '_WD_' + str(withDisk) + '_' + 'start_' + str(start_param_index) + '_' + galaxy + '_simultaneous_' + 'PSM_' + str(pop_selection_method) + '_N_' + str(nIterations) + '.csv'
            if use_artificial_data and generation_params_storer_index != 'None':
                file_name = 'artificial_params_' + str(generation_params_storer_index) + '_' + file_name
            header = 'dist, M, vphi, sigsqr ,el, rs , phi , theta , h_xhat, h_yhat, h_zhat, h_x_center, h_y_center, h_z_center, c, lam, zeta, eps, a, b, d_xhat, d_yhat, d_zhat, d_x_center, d_y_center, d_z_center, nVisits, logLikelihood'

            np.savetxt(save_dir + file_name, parameter_result_array, delimiter=",",header = header)
        
        
    
if __name__ == '__main__':
    el = 0.5
    lam = 0.2
    population = ['fornax','MP']
    runMCMCForDwarfGalaxyProfiles(el, lam, population, disk_funct = 0, halo_funct = 0, withDisk=0, nIterations = 2, saverun = 0, show_best_fit = 0, pop_selection_method = 'none')
