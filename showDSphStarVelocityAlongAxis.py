import math
import numpy as np
from AstronomicalParameterArchive import AstronomicalParameterArchive
from ObservedGalaxyStarData import ObservedGalaxyStarData
from logList import logList 
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np

def showDSphStarVelocityAlongAxis (population, rot_plane_angle, pop_selection_method = 'none', show_plot = 0):
    astro_archive = AstronomicalParameterArchive ()
    star_data = ObservedGalaxyStarData(population, pop_selection_method = pop_selection_method)
    corr_ra = star_data.corrRa
    corr_dec = star_data.corrDec
    Vhel = star_data.Vhel
    VhelE = star_data.VhelE
    dist_along_rot_plane = [math.cos(rot_plane_angle) * corr_ra[i] + math.sin(rot_plane_angle) * corr_dec[i] for i in range(len(corr_ra)) ]
    dist_perp_rot_plane = [-math.sin(rot_plane_angle) * corr_ra[i] + math.cos(rot_plane_angle) * corr_dec[i] for i in range(len(corr_dec)) ]
    n_stars = len(dist_along_rot_plane)
    mean_dist_along_rot_plane = sum(dist_along_rot_plane) / n_stars
    mean_Vhel = sum(Vhel) / n_stars 
    pearson_correlation_statistic_numerator = sum([dist_along_rot_plane[i] * Vhel[i] for i in range(n_stars)]) - n_stars * mean_Vhel * mean_dist_along_rot_plane
    pearson_correlation_statistic_denominator = math.sqrt(sum([dist_along_rot_plane[i] ** 2 for i in range(n_stars)]) - n_stars * mean_dist_along_rot_plane ** 2) * math.sqrt(sum([Vhel[i] ** 2 for i in range(n_stars)]) - n_stars * mean_Vhel ** 2)
    pearson_correlation_statistic = (pearson_correlation_statistic_numerator) / (pearson_correlation_statistic_denominator)                          
    if show_plot:
        plt.errorbar(dist_along_rot_plane, Vhel, yerr = VhelE, linestyle = 'None')
        plt.show() 
    print 'pearson_correlation_statistic = ' + str(pearson_correlation_statistic)
   
    return pearson_correlation_statistic
