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


def runMCMCForDwarfGalaxyProfilesAllPops(galaxy, disk_funct = 0, halo_funct = 0, withDisk=0, nIterations = 2, saverun = 0, show_best_fit = 0, pop_selection_method = 'none', start_param_index = 0, outer_step = 1.0, use_artificial_data = 0, generation_params_storer_index = 'None',generation_disk_funct = 0, generation_halo_funct = 0, max_R_bins = 200, max_z_bins = 200,apply_observation_mask = 1):
    astro_archive = AstronomicalParameterArchive ()
    computational_params_archive = ComputationalArchive() 
    start_param_storer = StartMCMCParameterStorer()
    start_params = start_param_storer.getStartParameters(start_param_index, withDisk) 
    dSph_archive = DwarfGalDataArchive()
    deg_to_rad = astro_archive.getDegToRad()
    arcmin_limits_R,arcmin_limits_z = dSph_archive.getObservationBounds([galaxy,'dummy_pop_var'],return_unit = 'arcmin')
    arcmin_limits_los = [-50.0,50.0] 
    pops = dSph_archive.getPopulations(galaxy)
    #print 'pops = '
    #print pops   
    star_data = []
    popGalPairs = []
    start_dsph_parameter_storers  = []
    for i in range(len(pops)):
        pop = pops[i]
        popGalPairs = popGalPairs + [[galaxy,pop]]
        star_data = star_data + [ObservedGalaxyStarData([galaxy,pop], pop_selection_method = pop_selection_method)]

        
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
                                        rs = start_params['rs'], phi = start_params['phi'], theta = start_params['theta'],
                                        zeta = start_params['zeta'], eps = start_params['eps'], a = start_params['a'], b= start_params['b'],
                                        sigsqr = star_data[i].sigSqr)
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
                star_data[i].corrRa, star_data[i].corrDec = getArtificialStarPositions(len(star_data[i].corrRa), generation_params_storer,
                                                                                         start_dsph_parameter_storers[i].sigsqr,
                                                                                         max_R_bins, max_z_bins, outer_step,
                                                                                         arcmin_limits_R, arcmin_limits_z,arcmin_limits_los,
                                                                                         generation_disk_funct, generation_halo_funct,withDisk )
                star_data[i].corrRa = np.array(star_data[i].corrRa)
                star_data[i].corrDec = np.array(star_data[i].corrDec) 
                print 'For population ' + pops[i] + ', we have ' + str(len(star_data[i].corrRa)) + ' artificial stars. '
                #print 'For population ' + pops[i] + ', the star data is: '
                #print 'star_data[i].corrRa = '
                #print star_data[i].corrRa
                #print 'star_data[i].corrDec = '
                #print star_data[i].corrDec 
        
    potential_archive = PotentialArchive()
    #disk_file = potential_archive.getDiskFile(lam)
    #halo_file = potential_archive.getHaloFile(el)
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
    zeta = start_dsph_parameter_storers[0].zeta 
    gamma = astro_archive.getGamma()
    (param_indeces_to_vary, param_ranges, gauss_sampling_width) =  start_dsph_parameter_storers[0].getMCMCInfo()
    #print 'start_dsph_parameter_storer.dist = ' + str(start_dsph_parameter_storers[0].dist)


    #arcmin_limits_R = 60.0
    #arcmin_limits_z = 60.0
    #R_inner=np.arange(0,arcmin_limits_R*zeta + outer_step_R * zeta,outer_step_R * zeta)
    #z_inner=np.arange(0,arcmin_limits_z*zeta + outer_step_z * zeta,outer_step_z * zeta)
    #R_outer=np.arange(0,arcmin_limits_R + outer_step_R, outer_step_R)
    #z_outer=np.arange(0,arcmin_limits_z + outer_step_z, outer_step_z)
    #R_arcmin = np.unique(np.around(np.concatenate( (-R_outer,-R_inner,R_inner,R_outer) ),5))
    #z_arcmin = np.unique(np.around(np.concatenate( (-z_outer,-z_inner,z_inner,z_outer) ),5))
    #R_arcmin = np.sort(R_arcmin)
    #z_arcmin = np.sort(z_arcmin)
    
    #print 'max, min R in arcminutes are: ' + str(max(R_outer)) + ', ' + str(min(R_outer))
    #print 'max, min z in arcminutes are: ' + str(max(z_outer)) + ', ' + str(min(z_outer))
    #los_bins_outer = np.arange(max(R_arcmin)*-1.0,max(R_arcmin)*1.0 + 1.0,1.0,dtype=np.float)
    #los_bins_inner = np.arange(max(R_arcmin)*-0.1,max(R_arcmin)*0.1+0.1,0.1,dtype=np.float)
    #los_bins_arcmin = np.unique(np.around(np.concatenate( (los_bins_outer,los_bins_inner) ),5))
    #los_bins_arcmin = np.sort(los_bins_arcmin)
    #los_bins_arcmin = R_arcmin 
    
    #R=R_arcmin*1/60.0*1/180.0*math.pi*dist/rs
    #z=z_arcmin*1/60.0*1/180.0*math.pi*dist/rs
    #los_bins = los_bins_arcmin * 1/60.0 * 1/180.0*math.pi*dist/rs
    R,z,los_bins = compRZLos(zeta, zeta * outer_step, outer_step, arcmin_limits_R, arcmin_limits_z, arcmin_limits_los, rs, dist)
    
    #print 'zeta = ' + str(zeta) + ' so, len(los_bins) = ' + str(len(los_bins)) 
    #print 'los_bins = ' 
    #print los_bins
    #print ' R = '
    #print R
    #print 'z = '
    #print z 
    
    current_parameter_storers = start_dsph_parameter_storers

    current_brightness_profiles = []
    current_observation_mask = dSph_archive.getObservationMask(popGalPairs[0],R * rs / dist * 1.0 / deg_to_rad, z * rs / dist * 1.0 / deg_to_rad )
    for storer in current_parameter_storers:
        #print 'storer.population = '
        #print 'len(R) = ' + str(len(R))
        #print 'len(z) = ' + str(len(z))
        #print storer.population
        current_brightness_profiles = current_brightness_profiles + [SurfaceBrightnessProfile(R, z, los_bins, gamma, storer,
                                                                                              disk_interpolating_function = disk_funct,
                                                                                              halo_interpolating_function = halo_funct,
                                                                                              observation_mask = current_observation_mask ) ]

    #for current_brightness_profile in current_brightness_profiles: 
        #print 'current_brightness_profile.Rlikeaxis = '
        #print current_brightness_profile.Rlikeaxis
        #print 'current_brightness_profile.zlikeaxis = '
        #print current_brightness_profile.zlikeaxis
        
    #The log likelihood for the WHOLE galaxy should be the sum of the log likelihoods of the individual populations (product of their likelihoods). 
    current_log_likelihood = 0.0
    for i in range(len(current_brightness_profiles)):
        #print 'current_likelihood for pop ' + popGalPairs[i][1] + ' ' + str(current_brightness_profiles[i].sumSurfaceBrightness(star_data[i].corrRa,star_data[i].corrDec))
        current_log_likelihood = current_log_likelihood + current_brightness_profiles[i].sumLogSurfaceBrightness(star_data[i].corrRa,star_data[i].corrDec)
    
    #Note that the parameters we vary as part of our MCMC algorithm of interest in varying MCMC 
    visited_points = [VisitedDsphMCMCPoint(current_parameter_storers[0],current_log_likelihood,1)]
    
    #Begin MCMC loop
    for i in range(nIterations):
        start = time.time()
        jump_multiplier = 1
        #print 'On iteration ' + str(i) 
        if i%100 == 100-1:
            print 'On iteration '  + str(i)
            jump_multiplier = 10
        #if i%50 == 50-1:
        #    print 'On iterations ' + str(i) 
        #print 'Have visited ' + str(len(visited_points)) + ' points with parameter sets: '
        #for point in visited_points:
        #        print point.parameter_storer.getParameterSet() 
        current_parameter_sets = []
        test_parameter_sets = []
        for storer in current_parameter_storers:
            current_parameter_sets = current_parameter_sets + [storer.getParameterSet()]
            test_parameter_sets = test_parameter_sets + [storer.getParameterSet()[:]]
        #print 'current parameter set 0 is: '
        #print current_parameter_sets[0]
        #test_parameter_set = current_parameter_set[:]
        #test_parameter_set = 
        #print 'And we are going to vary indeces ' 
        #print param_indeces_to_vary
        
        #Note that the indeces we vary are underlying properties of potential, and so must be the same for each population 
        for j in range(len(param_indeces_to_vary)):
            index = param_indeces_to_vary[j]
            new_parameter_value = random.gauss(test_parameter_sets[0][index],gauss_sampling_width[j]*jump_multiplier)
            if new_parameter_value < param_ranges[j][0]:
                new_parameter_value = param_ranges[j][0]
            if new_parameter_value > param_ranges[j][1]:
                new_parameter_value = param_ranges[j][1]
            for test_parameter_set in test_parameter_sets: 
                test_parameter_set[index] = new_parameter_value 
                

        #test_parameter_set = [ random.gauss(current_parameter_set[j],gauss_sampling_width[j]) if j in param_indeces_to_vary else current_parameter_set[j] for j in range(len(current_parameter_set)) ] 

        test_parameter_storers = []
        for pair in range(len(popGalPairs)):
            test_parameter_storers = test_parameter_storers + [DwarfGalaxyParametersStorer(popGalPairs[pair], withDisk, *(test_parameter_sets[pair]))]
        #for popGalPair in popGalPairs:
        #    test_parameter_storers = test_parameter_storers + [DwarfGalaxyParametersStorer(popGalPair, withDisk, *test_parameter_set)]
        #print 'Defining new test parameter 0 as: '
        #print test_parameter_storers[0].getParameterSet()

        c_test = test_parameter_storers[0].c
        dist_test = test_parameter_storers[0].dist
        rs_test = test_parameter_storers[0].rs
        zeta_test = test_parameter_storers[0].zeta 
        
        #R_inner_test=np.arange(0,arcmin_limits_R*zeta_test + outer_step_R * zeta_test, outer_step_R * zeta_test)
        #z_inner_test=np.arange(0,arcmin_limits_z*zeta_test + outer_step_z * zeta_test, outer_step_z * zeta_test)
        #R_outer_test=np.arange(0,arcmin_limits_R + outer_step_R, outer_step_R)
        #z_outer_test=np.arange(0,arcmin_limits_z + outer_step_z, outer_step_z)
        #R_arcmin_test = np.unique(np.around(np.concatenate( (-R_outer_test,-R_inner_test,R_inner_test,R_outer_test) ),5))
        #z_arcmin_test = np.unique(np.around(np.concatenate( (-z_outer_test,-z_inner_test,z_inner_test,z_outer_test) ),5))
        #R_arcmin_test = np.sort(R_arcmin_test)
        #z_arcmin_test = np.sort(z_arcmin_test)
        
        #R_test=R_arcmin_test*1/60.0*1/180.0*math.pi*dist_test/rs_test
        #z_test=z_arcmin_test*1/60.0*1/180.0*math.pi*dist_test/rs_test


        #los_bins_outer_test = np.arange(max(R_test)*-1.0,max(R_test)*1.0 + 1.0,1.0,dtype=np.float)
        #los_bins_inner_test = np.arange(max(R_test)*-0.1,max(R_test)*0.1+0.1,0.1,dtype=np.float)
        #los_bins_test = np.unique(np.around(np.concatenate( (los_bins_outer_test,los_bins_inner_test) ),5))
        #los_bins_arcmin_test = R_arcmin_test
        #los_bins_test = los_bins_arcmin_test*1/60.0*1/180.0*math.pi*dist_test/rs_test

        R_test,z_test,los_bins_test = compRZLos(zeta_test, zeta_test * outer_step, outer_step, arcmin_limits_R, arcmin_limits_z, arcmin_limits_los, rs_test, dist_test)
        if i % 50 == 50 - 1:
            print 'zeta_test = ' + str(zeta_test) + ' so, len(los_bins_test) = ' + str(len(los_bins_test)) 

        test_brightness_profiles = []

        #print 'rs_test = ' + str(rs_test)
        #print 'dist_test = ' + str(dist_test)
        test_observation_mask = dSph_archive.getObservationMask(popGalPairs[0],R_test * rs_test / dist_test * 1.0 / deg_to_rad, z_test * rs_test / dist * 1.0 / deg_to_rad )
        for storer in test_parameter_storers:
            #print 'storer.dist = ' + str(storer.dist)
            #print 'storer.rs = ' + str(storer.rs)
            #print 'storer.zeta = ' + str(storer.zeta)
            #print 'storer.a = ' + str(storer.a)
            #print 'storer.b = ' + str(storer.b)
            #print 'storer.el = ' + str(storer.el)
            #print 'storer.population = '
            #print storer.population 
            #print 'len(R_test) = ' + str(len(R_test))
            #print 'len(z_test) = ' + str(len(z_test))
            #storer.printContents()
            test_brightness_profiles = test_brightness_profiles + [SurfaceBrightnessProfile(R_test, z_test, los_bins_test, gamma, storer, 
                                                                                            disk_interpolating_function = disk_funct,
                                                                                            halo_interpolating_function = halo_funct,
                                                                                            observation_mask = test_observation_mask) ]
                                                                  
        test_log_likelihood = 0.0
        for k in range(len(test_brightness_profiles)):
            #print 'i = ' + str(i)
            #print 'test_brightness_profiles[i].dist = ' + str(test_brightness_profiles[i].dist)
            #print 'test_brightness_profiles[i].rs = ' + str(test_brightness_profiles[i].rs)
            #print 'test_brightness_profiles[i].Rlikeaxis = '
            #print test_brightness_profiles[i].Rlikeaxis
            #print 'test_brightness_profiles[i].Rlikeaxis = '
            #print test_brightness_profiles[i].Rlikeaxis
            #print 'star_data[i].corrRa = '
            #print star_data[i].corrRa
            #print 'star_data[i].corrDec = '
            #print star_data[i].corrDec            
            #print 'test_likelihood for pop ' + popGalPairs[i][1] + ' ' + str(test_brightness_profiles[i].sumLogSurfaceBrightness(star_data[i].corrRa,star_data[i].corrDec))
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

        end = time.time()
        if i % 100 == 100-1:
            print 'Took ' + str(end - start) + ' seconds for iteration ' + str(i) + '. '

            
    #print 'Visited a total of: ' + str(len(visited_points)) + ' points with following log likelihoods: '
    #print [point.log_likelihood for point in visited_points]
    #print 'And parameter sets: '
    
    parameter_result_array = np.zeros((len(visited_points),len(visited_points[0].parameter_storer.getParameterSet()) + 2 ))
    max_visits_index = 0
    max_visits = 0
    for j in range(len(visited_points)):
        point = visited_points[j]
        #print point.parameter_storer.getParameterSet()
        parameter_result_array[j,:] = point.parameter_storer.getParameterSet() + [point.n_visits] + [point.log_likelihood]
        if point.n_visits > max_visits:
            max_visits = point.n_visits
            max_visits_index = j
            
    #print 'max_visits = ' + str(max_visits)
    #print 'max_visits index = ' + str(max_visits_index)
    max_point = visited_points[max_visits_index]
    max_parameter_storer = max_point.parameter_storer
    max_visit_parameters = max_parameter_storer.getParameterSet()
    
    dist_best = max_parameter_storer.dist
    rs_best = max_parameter_storer.rs
    zeta_best = max_parameter_storer.zeta

    #R_inner_best=np.arange(0,arcmin_limits_R*zeta_best + outer_step_R * zeta_best, outer_step_R * zeta_best)
    #z_inner_best=np.arange(0,arcmin_limits_z*zeta_best + outer_step_z * zeta_best, outer_step_z * zeta_best)
    #R_outer_best=np.arange(0,arcmin_limits_R + outer_step_R, outer_step_R)
    #z_outer_best=np.arange(0,arcmin_limits_z + outer_step_z, outer_step_z)
    #R_arcmin_best = np.unique(np.around(np.concatenate( (-R_outer_best,-R_inner_best,R_inner_best,R_outer_best) ),5))
    #z_arcmin_best = np.unique(np.around(np.concatenate( (-z_outer_best,-z_inner_best,z_inner_best,z_outer_best) ),5))
    #R_arcmin_best = np.sort(R_arcmin_best)
    #z_arcmin_best = np.sort(z_arcmin_best)

    #los_bins_arcmin_best = R_arcmin_best
        
    #R_best=R_arcmin_best*1/60.0*1/180.0*math.pi*dist_best/rs_best
    #z_best=z_arcmin_best*1/60.0*1/180.0*math.pi*dist_best/rs_best
    #los_bins_best = los_bins_arcmin_best*1/60.0*1/180.0*math.pi*dist_best/rs_best
    R_best,z_best,los_bins_best = compRZLos(zeta_best, zeta_best * outer_step, outer_step, arcmin_limits_R, arcmin_limits_z, arcmin_limits_los, rs_best, dist_best)
    #R_best=R_arcmin*1/60.0*1/180.0*math.pi*dist_best/rs_best
    #z_best=z_arcmin*1/60.0*1/180.0*math.pi*dist_best/rs_best
    #los_bins_best = np.unique(np.around(np.concatenate((np.arange(max(R_best)*-1.0,max(R_best)*1.0 + 1.0,1.0,dtype=np.float),np.arange(max(R_best)*-0.5*0.1,max(R_best)*0.5*0.1+0.05,0.05,dtype=np.float))),5))
    print 'max_visit_parameters = '
    print max_visit_parameters 
    #print 'start_parameters = ' 
    #print start_dsph_parameter_storers[0].getParameterSet()
    best_fit_observation_mask = dSph_archive.getObservationMask(popGalPairs[0],R_best * rs_best / dist_best * 1.0 / deg_to_rad, z_best * rs_best / dist_best * 1.0 / deg_to_rad )
    best_fit_brightness_profile = SurfaceBrightnessProfile(R_best, z_best, los_bins, gamma, max_parameter_storer,
                                                                   disk_interpolating_function = disk_funct,
                                                                   halo_interpolating_function = halo_funct,
                                                                   observation_mask = best_fit_observation_mask) 
    
    #print parameter_result_array   
    
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
        file_name = 'MCMC_output_data_WD_' + str(withDisk) + '_' + 'start_' + str(start_param_index) + '_' + galaxy + '_simultaneous_' + 'PSM_' + str(pop_selection_method) + '_N_' + str(nIterations) + '.csv'
        if use_artificial_data and generation_params_storer_index != 'None':
            file_name = 'artificial_params_' + str(generation_params_storer_index) + '_' + file_name
        header = 'dist, M, vphi, sigsqr ,el, rs , phi , theta , c, lam, zeta, eps, a, b, nVisits, logLikelihood'

        np.savetxt(save_dir + file_name, parameter_result_array, delimiter=",",header = header)
        
        
    
if __name__ == '__main__':
    el = 0.5
    lam = 0.2
    population = ['fornax','MP']
    runMCMCForDwarfGalaxyProfiles(el, lam, population, disk_funct = 0, halo_funct = 0, withDisk=0, nIterations = 2, saverun = 0, show_best_fit = 0, pop_selection_method = 'none')
