import math
import numpy as np
from AstronomicalParameterArchive import AstronomicalParameterArchive
from ObservedGalaxyStarData import ObservedGalaxyStarData
from logList import logList 
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

def testFunction(number):
    return number + 1.0 

def measureDSphStarVelocityCorrelationAlongAxis (population, rot_plane_angle, pop_selection_method = 'none', show_plot = 0, dist_n_bins = 20, fit_order = 1, correctVhel = 1):
    print testFunction(10) 
    astro_archive = AstronomicalParameterArchive ()
    star_data = ObservedGalaxyStarData(population, pop_selection_method = pop_selection_method)
    corr_ra = star_data.corrRa
    corr_dec = star_data.corrDec
    if correctVhel:
        #print 'Measuring corrected Vhel, using reported Fornax PM.' 
        Vhel = star_data.corr_Vhel
        VhelE = star_data.corr_VhelE
    else:
        #print 'Measuring raw Vhel.' 
        Vhel = star_data.Vhel
        VhelE = star_data.VhelE
    #VhelE = star_data.VhelE
    dist_along_rot_plane = [math.cos(rot_plane_angle) * corr_ra[i] + math.sin(rot_plane_angle) * corr_dec[i] for i in range(len(corr_ra)) ]
    dist_perp_rot_plane = [-math.sin(rot_plane_angle) * corr_ra[i] + math.cos(rot_plane_angle) * corr_dec[i] for i in range(len(corr_dec)) ]
    n_stars = len(dist_along_rot_plane)
    mean_dist_along_rot_plane = sum(dist_along_rot_plane) / n_stars
    mean_Vhel = sum(Vhel) / n_stars 
    pearson_correlation_statistic_numerator = sum([dist_along_rot_plane[i] * Vhel[i] for i in range(n_stars)]) - n_stars * mean_Vhel * mean_dist_along_rot_plane
    pearson_correlation_statistic_denominator = math.sqrt(sum([dist_along_rot_plane[i] ** 2 for i in range(n_stars)]) - n_stars * mean_dist_along_rot_plane ** 2) * math.sqrt(sum([Vhel[i] ** 2 for i in range(n_stars)]) - n_stars * mean_Vhel ** 2)
    pearson_correlation_statistic = (pearson_correlation_statistic_numerator) / (pearson_correlation_statistic_denominator)                          
    #print 'pearson_correlation_statistic = ' + str(pearson_correlation_statistic)
    bin_step = (max(dist_along_rot_plane) - min(dist_along_rot_plane)) / dist_n_bins
    dist_bins_start = np.arange(min(dist_along_rot_plane), max(dist_along_rot_plane), bin_step)
    dist_bins_centers = [start + bin_step / 2.0 for start in dist_bins_start]
    dist_bins_start[0] = dist_bins_start[0] - bin_step / 10.0
    dist_bins_start[-1] = dist_bins_start[0-1] + bin_step / 10.0    
    #print dist_bins_start 
    binned_vhels = [[0.0,0.0,0.0] for bin in dist_bins_start]
    for i in range(n_stars):
        for j in range(len(dist_bins_start)):
            #print 'j = ' + str(j)
            if dist_bins_start[j] < dist_along_rot_plane[i] and dist_bins_start[j] + bin_step >= dist_along_rot_plane[i]:
                binned_vhels[j][0] = binned_vhels[j][0] + Vhel[i]
                binned_vhels[j][1] = binned_vhels[j][1] + VhelE[i] ** 2.0
                binned_vhels[j][2] = binned_vhels[j][2] + 1
                break 
    mean_vhels = [binned_vhel[0] / binned_vhel[2] if binned_vhel[0] != 0.0 else 0.0 for binned_vhel in binned_vhels]
    mean_vhels_errs = [math.sqrt(binned_vhel[1]) / binned_vhel[2] if binned_vhel[0] != 0.0 else 0.0 for binned_vhel in binned_vhels]
    good_mean_vhels = [mean_vhels[i] for i in range(len(mean_vhels)) if mean_vhels[i] > 0.0]
    good_mean_vhels_errs = [mean_vhels_errs[i] for i in range(len(mean_vhels)) if mean_vhels[i] > 0.0]
    good_dist_bins_centers = [dist_bins_centers[i] for i in range(len(mean_vhels)) if mean_vhels[i] > 0.0]
    

    fit_params = np.polyfit(np.array(good_dist_bins_centers), np.array(good_mean_vhels), fit_order, w = 1/np.array(good_mean_vhels_errs))
    #print 'fit_params = '
    #print fit_params
    #np.flip(fit_params,0) 
    #print dist_bins_centers
    #print mean_vhels
    best_fit_curve = [ sum([(dist ** n) * np.flip(fit_params,0)[n] for n in range(fit_order + 1) ]) for dist in dist_bins_centers]
    #print best_fit_curve 

    if show_plot:
        f, axarr = plt.subplots(2, sharex = True)
        axarr[0].errorbar(dist_along_rot_plane, Vhel, yerr = VhelE, linestyle = 'None')
        axarr[1].errorbar(dist_bins_centers, mean_vhels, yerr = mean_vhels_errs, linestyle = 'None')
        axarr[1].plot(dist_bins_centers, best_fit_curve,color = 'red') 
        axarr[0].set_xlabel('Separation from galactic center (deg)')
        axarr[1].set_xlabel('Separation from galactic center (deg)')
        axarr[0].set_ylabel('Heliocentric velocity (km/s)')
        axarr[1].set_ylabel('Averaged heliocentric velocity (km/s)')
        best_fit_str = 'Best fit params, from highest to lowest order: '
        for i in range(fit_order + 1): best_fit_str = best_fit_str + ' ' + str(fit_params[i])[0:5]
        print dist_bins_centers[1]
        print mean_vhels[1]
        axarr[1].text(dist_bins_centers[0], mean_vhels[0],best_fit_str,color = 'red')
        #plt.errorbar(dist_along_rot_plane, Vhel, yerr = VhelE, linestyle = 'None')
        #plt.show() 
        #plt.scatter(dist_bins_centers, mean_vhels)
        plt.show()

    return pearson_correlation_statistic, fit_params
