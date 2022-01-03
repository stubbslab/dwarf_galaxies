import numpy as np
import math
from SurfaceBrightnessProfile import SurfaceBrightnessProfile
from DwarfGalaxyParametersStorer import DwarfGalaxyParametersStorer
from AstronomicalParameterArchive import AstronomicalParameterArchive
from ComputationalArchive import ComputationalArchive 
from logList import logList 
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from compRZLos import compRZLos
import time

#measure the difference in surface profiles for various outer_step_widths 
def measureErrorOfSurfaceProfs(fine_surface_prof, fine_step, storer, inner_lim, inner_step, outer_steps, arcmin_limits_R, arcmin_limits_z, halo_file, disk_file, halo_funct=None, disk_funct=None, plot = 1, save_fig = 0):
    surface_prob_interpolator_fine = fine_surface_prof.onSkyInterpolator  
    astro_archive = AstronomicalParameterArchive()
    computational_params_archive = ComputationalArchive()
    gamma = astro_archive.getGamma()
    #surfaceProfilesToCompare = []
    #scaled_dif_mesh_array = []
    #interpolated_points_course_array = []
    #combined_R_array = []
    #combined_z_array = []
    if plot:
        large_font = 9
        med_font = 7
        small_font = 5
        matplotlib.rcParams.update({'font.size':small_font})
        #fig, axarr = plt.subplots(ncols = 2,nrows = len(scaled_dif_mesh_array) + 1,figsize = (6, 3 * len(outer_steps) ))
        fig = plt.figure(figsize = (7.2,2.0 * (len(outer_steps) + 1)))
        Rmesh,zmesh = np.meshgrid(fine_surface_prof.Rlikeaxis * storer.rs / storer.dist * 180 / math.pi * 60.0,fine_surface_prof.zlikeaxis * storer.rs / storer.dist * 180 / math.pi * 60.0)
        samplingPoints = np.dstack(( fine_surface_prof.Rlikeaxis.tolist(), fine_surface_prof.Rlikeaxis.tolist() ))[0]
        interpolated_points_fine = surface_prob_interpolator_fine(np.dstack(np.meshgrid(fine_surface_prof.Rlikeaxis,fine_surface_prof.zlikeaxis)))
        log_levels_full = logList(np.min(np.abs(interpolated_points_fine)),np.max(np.abs(interpolated_points_fine)),20)
        plt.subplot(len(outer_steps)+1, 3, 1)
        plt.contour(Rmesh,zmesh,interpolated_points_fine,norm=LogNorm(),levels=log_levels_full)
        plt.colorbar(shrink=0.8,extend='both',format='%.1e')
        plt.ylabel('Dist from center(arcmin)', fontsize = med_font) 
        plt.title('Finely Sampled Normalized Surface Probability',fontsize = med_font)
    
    for i in range(len(outer_steps)):
        outer_step = outer_steps[i]
        print 'outer_step[' + str(i) + '] = ' + str(outer_step) 
        R,z,los = compRZLos(inner_lim, inner_step, outer_step, arcmin_limits_R, arcmin_limits_z, storer.rs, storer.dist)
        #print 'R = '
        #print R
        #print 'los = '
        #print los
        start = time.time()
        surfaceProfileToCompare = SurfaceBrightnessProfile(R, z, los, gamma, storer, disk_file, halo_file, disk_interpolating_function = disk_funct, halo_interpolating_function = halo_funct)
        surface_prob_interpolator_course = surfaceProfileToCompare.onSkyInterpolator
        end = time.time()
        print 'Took ' + str(end-start) + ' s to create and get surface interpolator for outer_step = ' + str(outer_step)
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
        #combined_R_array = combined_R_array + [combined_R]
        #combined_z_array = combined_z_array + [combined_z]
        #interpolated_points_course_array = interpolated_points_course_array + [interpolated_points_course]
        #scaled_dif_mesh_array = scaled_dif_mesh_array + [scaled_dif_mesh]

        if plot:
            Rmesh,zmesh = np.meshgrid(combined_R * storer.rs / storer.dist * 180 / math.pi * 60.0, combined_z * storer.rs / storer.dist * 180 / math.pi * 60.0)
            log_levels_full = logList(np.min(np.abs(interpolated_points_fine)),np.max(np.abs(interpolated_points_fine)),20)
            #axarr[i+1,0].contour(Rmesh,zmesh,interpolated_points_course,norm=LogNorm(),levels=log_levels_full)
            #axarr[i+1,0].colorbar(shrink = 0.8, extend = 'both', format = '%.1e')
            plt.subplot(len(outer_steps)+1 ,3, 3 * i + 2)
            CS1 = plt.contour(Rmesh,zmesh,interpolated_points_course,norm=LogNorm(),levels=log_levels_full)
            plt.colorbar(shrink=0.8,extend='both',format='%.1e')
            plt.ylabel('Outer step ' + str(outer_step) + '\n Dist from center (arcmin)', fontsize = med_font)
            if i == 0:
                plt.title('Normalized Surface Probability Density',fontsize = med_font) 
            if i == len(outer_steps)-1:
                plt.xlabel('Dist from center (arcmin)',fontsize = small_font)
                
            
            log_levels_diff = logList(np.min(np.abs(scaled_dif_mesh)),np.max(np.abs(scaled_dif_mesh)),20)
            #axarr[i+1,1].contour(Rmesh,zmesh,scaled_dif_mesh,norm=LogNorm(),levels=log_levels_diff)
            #axarr[i+1,1].colorbar(shrink = 0.8, extend = 'both', format = '%.1e')
            plt.subplot(len(outer_steps)+1 ,3, 3 * i + 3)
            CS2 = plt.contour(Rmesh,zmesh,scaled_dif_mesh,norm=LogNorm(),levels=log_levels_diff)
            plt.colorbar(shrink=0.8,extend='both',format='%.1e')
            #CB_contour=plt.colorbar(shrink=0.8,extend='both',format='%.1e')
            if i == 0:
                plt.title('Normalized Probability Density Error',fontsize = med_font) 
            if i == len(outer_steps)-1:
                plt.xlabel('Dist from center (arcmin)', fontsize = med_font)

    if plot:
        #plt.suptitle('Prob error for halo with rs=' + str(storer.rs) + ', el = ' + str(storer.el),fontsize = large_font) 
        plt.show()
        
    #print 'interpolated_points_fine = '
    #print interpolated_points_fine

    #if plot:
    #    #fig, axarr = plt.subplots(ncols = 2,nrows = len(scaled_dif_mesh_array) + 1,figsize = (6, 3 * len(outer_steps) ))
    #    fig = plt.figure(figsize = (6,3 * len(outer_steps))) 
    #    for i in range(len(scaled_dif_mesh_array)):
    #        print ' i = ' + str(i) 
    #        combined_R = combined_R_array[i]
    #        combined_z = combined_z_array[i]
    #        interpolated_points_course = interpolated_points_course_array[i]
    #        scaled_dif_mesh = scaled_dif_mesh_array[i] 
    #        Rmesh,zmesh = np.meshgrid(combined_R,combined_z)
    #        log_levels_full = logList(np.min(np.abs(interpolated_points_fine)),np.max(np.abs(interpolated_points_fine)),20)
    #        #axarr[i+1,0].contour(Rmesh,zmesh,interpolated_points_course,norm=LogNorm(),levels=log_levels_full)
    #        #axarr[i+1,0].colorbar(shrink = 0.8, extend = 'both', format = '%.2e')
    #        plt.subplot(len(outer_steps) ,2, 2 * i -1 )
    #        CS1 = plt.contour(Rmesh,zmesh,interpolated_points_course,norm=LogNorm(),levels=log_levels_full)
    #        plt.colorbar(shrink=0.8,extend='both',format='%.2e')
    #        
    #        log_levels_diff = logList(np.min(np.abs(scaled_dif_mesh)),np.max(np.abs(scaled_dif_mesh)),20)
    #        #axarr[i+1,1].contour(Rmesh,zmesh,scaled_dif_mesh,norm=LogNorm(),levels=log_levels_diff)
    #        #axarr[i+1,1].colorbar(shrink = 0.8, extend = 'both', format = '%.2e')
    #        plt.subplot(len(outer_steps) ,2, 2 * i)
    #        CS2 = plt.contour(Rmesh,zmesh,scaled_dif_mesh,norm=LogNorm(),levels=log_levels_diff)
    #        plt.colorbar(shrink=0.8,extend='both',format='%.2e')
    #        #CB_contour=plt.colorbar(shrink=0.8,extend='both',format='%.2e')
    #    plt.show()
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
    if save_fig:
        plot_save_dir = computational_params_archive.getPlotDir() 
        fig_name = 'surfProbErr_rs_' + str(storer.rs) + '_el_' + str(storer.el) + 'fstep_' + str(fine_step)
        fig.savefig(plot_save_dir + fig_name + '.pdf',bbox_inches = 'tight') 
    return 0
