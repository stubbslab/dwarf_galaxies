import numpy as np
import math
from SurfaceBrightnessProfile import SurfaceBrightnessProfile
from DwarfGalaxyParametersStorer import DwarfGalaxyParametersStorer
from AstronomicalParameterArchive import AstronomicalParameterArchive
from logList import logList 
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from compRZLos import compRZLos
import time

def measureSurfaceProfDifference(storer1, storer2, inner_lim1, inner_lim2, inner_step1, inner_step2, outer_step1, outer_step2, arcmin_limits_R, arcmin_limits_z, halo_file1, disk_file1, halo_file2, disk_file2, halo_funct1=None, disk_funct1=None, halo_funct2=None, disk_funct2=None, plot_diff = 0, plot_individual = 0):
    astro_archive = AstronomicalParameterArchive() 
    gamma = astro_archive.getGamma() 
    R1,z1,los1 = compRZLos(inner_lim1, inner_step1, outer_step1, arcmin_limits_R, arcmin_limits_z, storer1.rs, storer1.dist)
    R2,z2,los2 = compRZLos(inner_lim2, inner_step2, outer_step2, arcmin_limits_R, arcmin_limits_z, storer2.rs, storer2.dist)
    print 'R1 = '
    print R1
    print 'R2 = '
    print R2
    print 'los1 = '
    print los1
    print 'los2 = '
    print los2
    start1 = time.time()
    surfaceProfile1 = SurfaceBrightnessProfile(R1, z1, los1, gamma, storer1, disk_file1, halo_file1, disk_interpolating_function = disk_funct1, halo_interpolating_function = halo_funct1)
    surface_prob_interpolator1 = surfaceProfile1.onSkyInterpolator
    end1 = time.time()
    print 'Took ' + str(end1-start1) + ' s to create and get surface interpolator from storer1. '
    start2=time.time()
    surfaceProfile2 = SurfaceBrightnessProfile(R2, z2, los2, gamma, storer2, disk_file2, halo_file2, disk_interpolating_function = disk_funct2, halo_interpolating_function = halo_funct2)
    surface_prob_interpolator2 = surfaceProfile2.onSkyInterpolator
    end2 = time.time() 
    print 'Took ' + str(end2-start2) + ' s to create and get surface interpolator from storer1. '
    combined_R = np.sort(np.array(list(set(R1) | set(R2))))[1:-1]
    combined_z = np.sort(np.array(list(set(z1) | set(z2))))[1:-1]
    samplingPoints = np.dstack(( combined_R.tolist(), combined_z.tolist() ))[0]

    #print 'combined_R = '
    #print combined_R
    #print 'combined_z = '
    #print combined_z
    #print 'samplingPoints = '
    #print samplingPoints 
    surface_prob_values1 = surface_prob_interpolator1(samplingPoints)
    surface_prob_values2 = surface_prob_interpolator2(samplingPoints)
    #print 'surface_prob_values1 = '
    #print surface_prob_values1
    #print 'surface_prob_values2 = '
    #print surface_prob_values2
    #print 'surface_prob_values2 - surface_prob_values1 = '
    #print [surface_prob_values2[i] - surface_prob_values1[i] for i in range(len(surface_prob_values1))]
    interpolated_points_1 = surface_prob_interpolator1(np.dstack(np.meshgrid(combined_R,combined_z))) 
    interpolated_points_2 = surface_prob_interpolator2(np.dstack(np.meshgrid(combined_R,combined_z)))
    dif_mesh = interpolated_points_2 - interpolated_points_1
    scaled_dif_mesh = dif_mesh / interpolated_points_2
    scaled_dif_mesh = np.abs(scaled_dif_mesh)
    print 'interpolated_points_1 = '
    print interpolated_points_1
    print 'interpolated_points_2 = '
    print interpolated_points_2

    if plot_individual:
        Rmesh,zmesh = np.meshgrid(combined_R,combined_z)
        fig = plt.figure()
        log_levels_indiv = logList(np.min(np.abs(interpolated_points_2)),np.max(np.abs(interpolated_points_2)),20)
        plt.subplot(1,2,1)
        CS=plt.contour(Rmesh,zmesh,interpolated_points_1,norm=LogNorm(),levels=log_levels_indiv)
        CB_contour=plt.colorbar(shrink=0.8,extend='both',format='%.2e')
        plt.subplot(1,2,2)
        CS=plt.contour(Rmesh,zmesh,interpolated_points_2,norm=LogNorm(),levels=log_levels_indiv)
        CB_contour=plt.colorbar(shrink=0.8,extend='both',format='%.2e')
        plt.show()
        #plt.close('all') 
    elif plot_diff:
        fig = plt.figure() 
        Rmesh,zmesh = np.meshgrid(combined_R,combined_z)
        log_levels=logList(np.min(np.abs(scaled_dif_mesh)),np.max(np.abs(scaled_dif_mesh)),20)
        CS=plt.contour(Rmesh,zmesh,scaled_dif_mesh,norm=LogNorm(),levels=log_levels)
        plt.title('Scaled difference in surface probability') 
        CB_contour=plt.colorbar(shrink=0.8,extend='both',format='%.2e')
        plt.show()

    #raw_input() 
    return 0
