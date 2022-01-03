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

#measure the difference in surface profiles for various outer_step_widths 
def measureErrorOfSurfaceProfs(fine_surface_prof, storer, inner_lim, inner_step, outer_steps, arcmin_limits_R, arcmin_limits_z, halo_file, disk_file, halo_funct=None, disk_funct=None, plot = 0):
    astro_archive = AstronomicalParameterArchive() 
    gamma = astro_archive.getGamma()
    surfaceProfilesToCompare = []
    scale_dif_mesh_array = []
    interpolated_points_course_array = []
    combined_R_array = []
    combined_z_array = []
    for outer_step in outer_steps: 
        R,z,los = compRZLos(inner_lim, inner_step, outer_step, arcmin_limits_R, arcmin_limits_z, storer.rs, storer.dist)
        #print 'R = '
        #print R
        #print 'los = '
        #print los
        start = time.time()
        surfaceProfileToCompare = SurfaceBrightnessProfile(R, z, los, gamma, storer, disk_file, halo_file, disk_interpolating_function = disk_funct, halo_interpolating_function = halo_funct)
        surface_prob_interpolator_course = surfaceProfileToCompare.onSkyInterpolator
        end = time.time()
        print 'Took ' + str(end-start) + ' s to create and get surface interpolator from storer1. '
        combined_R = np.sort(np.array(list(set(R) | set(fine_surface_prof.Rlikeaxis))))[1:-1]
        combined_z = np.sort(np.array(list(set(z) | set(fine_surface_prof.zlikeaxis))))[1:-1]
        samplingPoints = np.dstack(( combined_R.tolist(), combined_z.tolist() ))[0]

        #print 'combined_R = '
        #print combined_R
        #print 'combined_z = '
        #print combined_z
        #print 'samplingPoints = '
        #print samplingPoints 
        surface_prob_values_course = surface_prob_interpolator_course(samplingPoints)
        surface_prob_values_fine = surface_prob_interpolator_fine(samplingPoints)
        #print 'surface_prob_values1 = '
        #print surface_prob_values1
        #print 'surface_prob_values2 = '
        #print surface_prob_values2
        #print 'surface_prob_values2 - surface_prob_values1 = '
        #print [surface_prob_values2[i] - surface_prob_values1[i] for i in range(len(surface_prob_values1))]
        interpolated_points_course = surface_prob_interpolator_course(np.dstack(np.meshgrid(combined_R,combined_z))) 
        interpolated_points_fine = surface_prob_interpolator_fine(np.dstack(np.meshgrid(combined_R,combined_z)))
        dif_mesh = interpolated_points_fine - interpolated_points_course
        scaled_dif_mesh = dif_mesh / interpolated_points_fine
        scaled_dif_mesh = np.abs(scaled_dif_mesh)
        print 'interpolated_points_course = '
        print interpolated_points_course
        combined_R_array = Rmesh_array + [combined_R]
        combined_z_array = zmesh_array + [combined_z]
        interpolated_points_course_array = interpolated_points_course_array + [interpolated_points_course]
        scaled_dif_mesh_array = scaled_dif_mesh_array + [scaled_dif_mesh] 
        
    print 'interpolated_points_fine = '
    print interpolated_points_fine

    if plot:
        f, axarr = plt.subplots(2,len(scale_dif_mesh_array) + 1,figsize = (10, 5 * len(combined_R) )) 
        for i in range(len(scale_dif_mesh_array)):
            combined_R = combined_R_array[i]
            combined_z = combined_z_array[i]
            interpolated_points_course = interpolated_points_course_array[i]
            scaled_dif_mesh = scaled_dif_mesh_array[i] 
            Rmesh,zmesh = np.meshgrid(combined_R,combined_z)
            log_levels_full = logList(np.min(np.abs(interpolated_points_fine)),np.max(np.abs(interpolated_points_fine)),20)
            CS1 = axarr[1,i+1].contour(Rmesh,zmesh,interpolated_points_course,norm=LogNorm(),levels=log_levels_full)
            log_levels_diff = logList(np.min(np.abs(scaled_dif_mesh,norm=LogNorm(),levels=log_levels_full)))
            CS2 = axarr[2,i+1].contour(Rmesh,zmesh,scaled_dif_mesh,norm=LogNorm(),levels=log_levels_diff)
            CB_contour=plt.colorbar(shrink=0.8,extend='both',format='%.2e')
        plt.show()
    #elif plot_diff:
    #    fig = plt.figure() 
    #    Rmesh,zmesh = np.meshgrid(combined_R,combined_z)
    #    log_levels=logList(np.min(np.abs(scaled_dif_mesh)),np.max(np.abs(scaled_dif_mesh)),20)
    #    CS=plt.contour(Rmesh,zmesh,scaled_dif_mesh,norm=LogNorm(),levels=log_levels)
    #    plt.title('Scaled difference in surface probability') 
    #    CB_contour=plt.colorbar(shrink=0.8,extend='both',format='%.2e')
    #    plt.show()
    #    if return_plot: return fig 

    #raw_input() 
    return 0
