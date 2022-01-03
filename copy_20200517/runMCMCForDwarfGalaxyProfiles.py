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
#from DwarfGalaxyParametersStorer import getMCMCInfo
from AstronomicalParameterArchive import AstronomicalParameterArchive
from ComputationalArchive import ComputationalArchive 
from VisitedDsphMCMCPoint import VisitedDsphMCMCPoint
from logList import logList 
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from StartMCMCParameterStorer import StartMCMCParameterStorer
import copy 
import random
import math 
import numpy as np 

def runMCMCForDwarfGalaxyProfiles(el, lam, population, disk_funct = 0, halo_funct = 0, withDisk=0, nIterations = 2, saverun = 0, show_best_fit = 0, pop_selection_method = 'none', start_param_index = 0):
    astro_archive = AstronomicalParameterArchive ()
    computational_params_archive = ComputationalArchive()
    start_param_storer = StartMCMCParameterStorer()
    start_params = getStartParameters(start_param_index, withDisk = withDisk) 
    star_data = ObservedGalaxyStarData(population, pop_selection_method = pop_selection_method)
    start_dsph_parameter_storer = DwarfGalaxyParametersStorer(population,withDisk=withDisk,el=el,lam=lam,
                                                              rs = start_params[0], phi = start_params[1],
                                                              eps = start_params[2], a = start_params[3] )
    potential_archive = PotentialArchive()
    disk_file = potential_archive.getDiskFile(lam)
    halo_file = potential_archive.getHaloFile(el)
    if not disk_funct: 
        print 'Creating the disk_funct from file: ' + disk_file
        disk_funct = PotentialFunctionArray(disk_file)
    if not halo_funct:
        print 'Creating the halo_funct from file: ' + halo_file
        halo_funct = PotentialFunctionArray(halo_file)
    print 'Halo and disk interpolating functions created. '

    c = start_dsph_parameter_storer.c
    dist = start_dsph_parameter_storer.dist
    rs = start_dsph_parameter_storer.rs 
    gamma = astro_archive.getGamma()
    (param_indeces_to_vary, param_ranges, gauss_sampling_width) =  start_dsph_parameter_storer.getMCMCInfo()
    print 'start_dsph_parameter_storer.dist = ' + str(start_dsph_parameter_storer.dist)
    
    arcmin_limits_R = 60.0
    arcmin_limits_z = 60.0
    R_inner=np.arange(0,arcmin_limits_R*0.1 + 0.5,0.5)
    z_inner=np.arange(0,arcmin_limits_z*0.1 + 0.5,0.5)
    R_outer=np.arange(0,arcmin_limits_R + 2.0,2.0)
    z_outer=np.arange(0,arcmin_limits_z + 2.0,2.0)
    R_arcmin = np.unique(np.concatenate((-R_inner,-R_outer,R_inner,R_outer)))
    z_arcmin = np.unique(np.concatenate((-z_inner,-z_outer,z_inner,z_outer)))
    print 'max, min R in arcminutes are: ' + str(max(R_outer)) + ', ' + str(min(R_outer))
    print 'max, min z in arcminutes are: ' + str(max(z_outer)) + ', ' + str(min(z_outer))
    R=R_arcmin*1/60.0*1/180.0*math.pi*dist/rs
    z=z_arcmin*1/60.0*1/180.0*math.pi*dist/rs
    
    los_bins = np.unique(np.around(np.concatenate((np.arange(max(R)*-1.0,max(R)*1.0 + 0.25,0.25,dtype=np.float),np.arange(max(R)*-0.5*0.1,max(R)*0.5*0.1+0.05,0.05,dtype=np.float))),5))
    #print 'los_bins = ' 
    #print los_bins
    #print ' R = '
    #print R
    #print 'z = '
    #print z 
    
    current_parameter_storer = start_dsph_parameter_storer
    
    current_brightness_profile = SurfaceBrightnessProfile(R, z, los_bins, gamma, current_parameter_storer,
                                                                  disk_file, halo_file, 
                                                                  disk_interpolating_function = disk_funct,
                                                                  halo_interpolating_function = halo_funct       ) 
    print 'current_brightness_profile.Rlikeaxis = '
    print current_brightness_profile.Rlikeaxis
    print 'current_brightness_profile.zlikeaxis = '
    print current_brightness_profile.zlikeaxis
    
    current_log_likelihood = current_brightness_profile.sumLogSurfaceBrightness(star_data.corrRa,star_data.corrDec)
    current_star_probabilities = current_brightness_profile.starProbabilityValues(star_data.corrRa,star_data.corrDec)
    current_star_probabilities_by_old = current_brightness_profile.starProbabilityValuesOld(star_data.corrRa,star_data.corrDec)
    print 'length of star probabilities new method is: ' + str(len(current_star_probabilities))
    print 'length of star probabilities old method is: ' + str(len(current_star_probabilities_by_old))  
    print 'and there should be a total of ' + str(len(star_data.corrRa)) + ' stars'
    print 'To compare old and new method, '
    print 'The shape old star probabilities are : '
    print np.shape(current_star_probabilities_by_old)
    print 'The shape new star probabilities are : '
    print np.shape(current_star_probabilities)

    visited_points = [VisitedDsphMCMCPoint(current_parameter_storer,current_log_likelihood,1)]
    
    #Begin MCMC loop
    for i in range(nIterations): 
        if i%100 == 0: print 'On iteration '  + str(i)
        #print 'Have visited ' + str(len(visited_points)) + ' points with parameter sets: '
        #for point in visited_points:
        #        print point.parameter_storer.getParameterSet() 
        current_parameter_set = current_parameter_storer.getParameterSet()
        print 'current parameter set is: '
        print current_parameter_set
        test_parameter_set = current_parameter_set[:]
        #print 'And we are going to vary indeces ' 
        #print param_indeces_to_vary 
        for j in range(len(param_indeces_to_vary)):
            index = param_indeces_to_vary[j]
            test_parameter_set[index] = random.gauss(test_parameter_set[index],gauss_sampling_width[j])
            if test_parameter_set[index] < param_ranges[j][0]:
                test_parameter_set[index] = param_ranges[j][0]
            if test_parameter_set[index] > param_ranges[j][1]:
                test_parameter_set[index] = param_ranges[j][1]

        #test_parameter_set = [ random.gauss(current_parameter_set[j],gauss_sampling_width[j]) if j in param_indeces_to_vary else current_parameter_set[j] for j in range(len(current_parameter_set)) ] 
        print 'Defining new test parameters as: '
        print test_parameter_set 
        test_parameter_storer = DwarfGalaxyParametersStorer(population, withDisk, *test_parameter_set)
        
        c_test = test_parameter_storer.c
        dist_test = test_parameter_storer.dist
        rs_test = test_parameter_storer.rs 
        R_test=R_arcmin*1/60.0*1/180.0*math.pi*dist_test/rs_test
        z_test=z_arcmin*1/60.0*1/180.0*math.pi*dist_test/rs_test
        los_bins_test = np.unique(np.around(np.concatenate((np.arange(max(R_test)*-1.0,max(R_test)*1.0 + 0.5,0.5,dtype=np.float),np.arange(max(R_test)*-0.5*0.1,max(R_test)*0.5*0.1+0.05,0.05,dtype=np.float))),5))
        
        test_brightness_profile = SurfaceBrightnessProfile(R_test, z_test, los_bins_test, gamma, test_parameter_storer, 
                                                                  disk_file,halo_file,
                                                                  disk_interpolating_function = disk_funct,
                                                                  halo_interpolating_function = halo_funct        )
        print 'test_brightness_profile.Rlikeaxis = '
        print test_brightness_profile.Rlikeaxis
        print 'test_brightness_profile.zlikeaxis = '
        print test_brightness_profile.zlikeaxis
        test_log_likelihood = test_brightness_profile.sumLogSurfaceBrightness(star_data.corrRa,star_data.corrDec)
        log_rel_likelihood = test_log_likelihood - current_log_likelihood
        print 'The log_likelihood of the new configuration is: ' + str(test_log_likelihood)
        print 'The log_relative likelihood for the new configuration is: ' + str(log_rel_likelihood) 
        
        if log_rel_likelihood > math.log(random.random()):
            print 'Moving to new point. '
            visited_points = visited_points + [VisitedDsphMCMCPoint(test_parameter_storer,test_log_likelihood,1)]
            current_brightness_profile = test_brightness_profile 
            current_parameter_storer = test_parameter_storer
            current_log_likelihood = test_log_likelihood
        else:
            print 'Not moving to another point.  '  
            visited_points[-1].n_visits = visited_points[-1].n_visits + 1

            
    #print 'Visited a total of: ' + str(len(visited_points)) + ' points with following log likelihoods: '
    #print [point.log_likelihood for point in visited_points]
    #print 'And parameter sets: '
    
    parameter_result_array = np.zeros((len(visited_points) ,len(visited_points[0].parameter_storer.getParameterSet()) + 2 ))
    max_visits_index = 0
    max_visits = 0
    for j in range(len(visited_points)):
        point = visited_points[j]
        #print point.parameter_storer.getParameterSet()
        parameter_result_array[j,:] = point.parameter_storer.getParameterSet() + [point.n_visits] + [point.log_likelihood]
        if point.n_visits > max_visits:
            max_visits = point.n_visits
            max_visits_index = j
            
    print 'max_visits = ' + str(max_visits)
    print 'max_visits index = ' + str(max_visits_index)
    max_point = visited_points[max_visits_index]
    max_parameter_storer = max_point.parameter_storer
    max_visit_parameters = max_parameter_storer.getParameterSet()
    
    dist_best = max_parameter_storer.dist
    rs_best = max_parameter_storer.rs 
    R_best=R_arcmin*1/60.0*1/180.0*math.pi*dist_best/rs_best
    z_best=z_arcmin*1/60.0*1/180.0*math.pi*dist_best/rs_best
    los_bins_best = np.unique(np.around(np.concatenate((np.arange(max(R_best)*-1.0,max(R_best)*1.0 + 0.5,0.5,dtype=np.float),np.arange(max(R_best)*-0.5*0.1,max(R_best)*0.5*0.1+0.05,0.05,dtype=np.float))),5))
    print 'max_visit_parameters = '
    print max_visit_parameters 
    print 'start_parameters = ' 
    print start_dsph_parameter_storer.getParameterSet()
    best_fit_brightness_profile = SurfaceBrightnessProfile(R_best, z_best, los_bins, gamma, max_parameter_storer,
                                                                  disk_file, halo_file, 
                                                                  disk_interpolating_function = disk_funct,
                                                                  halo_interpolating_function = halo_funct       ) 
    
    #print parameter_result_array   
    
    if show_best_fit:
        
        #print star_data.corrRa,star_data.corrDec
        pop_of_interest = population[1] # can be 'MP', 'IM', 'MR'
        colormap='b' if pop_of_interest == 'MP' else 'g' if pop_of_interest == 'IM' else 'r'
        ylims = (-80,80) #ylimits from Walker plot
        xlims=(80,-80) #xlimits from Walker plot
        fig=plt.figure(1,figsize=(9,7))
        corr_ra_one_pop = star_data.corrRa
        corr_dec_one_pop = star_data.corrDec
        print 'NStars is ' + str(len(corr_dec_one_pop))
        plt.scatter(corr_ra_one_pop * 60.0,corr_dec_one_pop * 60.0 ,s=4.,color=colormap)
        zmesh,Rmesh=np.meshgrid(best_fit_brightness_profile.zlikeaxis,best_fit_brightness_profile.Rlikeaxis)
        zmesh = zmesh*180.0/math.pi*rs_best/dist_best*60.0
        Rmesh = Rmesh*180.0/math.pi*rs_best/dist_best*60.0
        best_fit_mesh = best_fit_brightness_profile.surfaceBrightness
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

        pltFileDir = computational_params_archive.getPlotDir()
        #Uncomment the following line if you want to save the file
        #plt.savefig(pltFileDir + gal_of_interest + '_' + pop_of_interest + '_metal_cuts_' + str((metallicity_cuts[gal_of_interest])[0]) + '_' + str((metallicity_cuts[gal_of_interest])[1]) + '.png')
        #plt.savefig(pltFileDir + gal_of_interest + '_' + pop_of_interest + '_full_pop_division' + '.png')
        plt.show()
    
        
    if saverun:
        save_dir = computational_params_archive.getMCMCOutDir 
        file_name = 'MCMC_output_data_WD_' + str(withDisk) + '_' + str(population[0]) + str(population[1]) + '_PSM_' + str(pop_selection_method) + '_el_' + str(el) + '_lam_' + str(lam) + '_N_' + str(nIterations) + '.csv'
        header = 'dist, M, vphi, sigsqr ,el, rs , phi , theta , c, lam, zeta, eps, a, b'

        np.savetxt(save_dir + file_name, parameter_result_array, delimiter=",",header = header)
        
        
    
if __name__ == '__main__':
    el = 0.5
    lam = 0.2
    population = ['fornax','MP']
    runMCMCForDwarfGalaxyProfiles(el, lam, population, disk_funct = 0, halo_funct = 0, withDisk=0, nIterations = 2, saverun = 0, show_best_fit = 0, pop_selection_method = 'none')
