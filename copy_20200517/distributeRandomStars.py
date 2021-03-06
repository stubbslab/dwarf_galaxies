import numpy as np
import math
import random 
from DwarfGalaxyParametersStorer import DwarfGalaxyParametersStorer
from SurfaceBrightnessProfile import SurfaceBrightnessProfile
from AstronomicalParameterArchive import AstronomicalParameterArchive
import itertools
import time

#Rbin_walls and zbin_walls should be given in scale radii (not angles) 
def distributeRandomStars(NStars, surface_profile,max_R_bins, max_z_bins,withDisk):
    #WD = surface_profile.WithDisk
    zeta = surface_profile.zeta
    Raxis = surface_profile.Rlikeaxis
    zaxis = surface_profile.zlikeaxis
    Rlims = [min(Raxis),max(Raxis)]
    zlims = [min(zaxis),max(zaxis)]
    if withDisk:
        max_z_bins = max_z_bins / 2
        max_R_bins = max_R_bins / 2
        Rbin_walls_inner = np.arange(Rlims[0] * zeta, Rlims[1] * zeta + (Rlims[1] - Rlims[0]) * zeta / (2.0 * max_R_bins), (Rlims[1] - Rlims[0]) * zeta / max_R_bins)
        zbin_walls_inner = np.arange(zlims[0] * zeta, zlims[1] * zeta + (zlims[1] - zlims[0]) * zeta / (2.0 * max_z_bins), (zlims[1] - zlims[0]) * zeta / max_z_bins)
        Rbin_walls_outer = np.arange(Rlims[0], Rlims[1] + (Rlims[1] - Rlims[0]) / (2.0 * max_R_bins), (Rlims[1] - Rlims[0]) / max_R_bins)
        zbin_walls_outer = np.arange(zlims[0], zlims[1] + (zlims[1] - zlims[0]) / (2.0 * max_z_bins), (zlims[1] - zlims[0]) / max_z_bins)
        Rbin_walls = np.unique(np.around(np.concatenate( (Rbin_walls_inner,Rbin_walls_outer) ),10))
        zbin_walls = np.unique(np.around(np.concatenate( (zbin_walls_inner,zbin_walls_outer) ),10))
    else:
        Rbin_walls = np.arange(Rlims[0], Rlims[1] + (Rlims[1] - Rlims[0]) / (2.0 * max_R_bins), (Rlims[1] - Rlims[0]) / max_R_bins)
        zbin_walls= np.arange(zlims[0], zlims[1] + (zlims[1] - zlims[0]) / (2.0 * max_z_bins), (zlims[1] - zlims[0]) / max_z_bins)        

    start1 = time.time()
    star_RAs = []
    star_Decs = []
    zleft_bin_wall = zbin_walls[0]
    Rleft_bin_wall = Rbin_walls[0]
    prob_array = []
    prob_interpolator = surface_profile.scaleRadiusInterpolator
    Rbin_centers = [(Rbin_walls[i+1] + Rbin_walls[i])/2.0 for i in range(len(Rbin_walls)-1) ]
    zbin_centers = [(zbin_walls[i+1] + zbin_walls[i])/2.0 for i in range(len(zbin_walls)-1) ]
    Rbin_centers_mesh,zbin_centers_mesh = np.meshgrid(Rbin_centers,zbin_centers)
    surface_probabilities = prob_interpolator (np.dstack((Rbin_centers_mesh,zbin_centers_mesh)))
    R_bin_widths = [(Rbin_walls[i+1] - Rbin_walls[i]) for i in range(len(Rbin_walls)-1) ]
    z_bin_widths = [(zbin_walls[i+1] - zbin_walls[i]) for i in range(len(zbin_walls)-1) ]
    R_bin_widths_mesh,z_bin_widths_mesh = np.meshgrid(R_bin_widths,z_bin_widths)
    area_mesh = R_bin_widths_mesh * z_bin_widths_mesh
    prob_array = surface_probabilities * area_mesh
    prob_array = prob_array.flatten()
    total_prob = np.sum(prob_array)

    #print 'prob_array = '
    #print prob_array
    #for zright_bin_wall in zbin_walls[1:]:
    #    for Rright_bin_wall in Rbin_walls[1:]:
    #        #print (Rleft_bin_wall,zleft_bin_wall)
    #        #print  'prob_interpolator(( (Rright_bin_wall + Rleft_bin_wall) / 2.0, (zright_bin_wall + zleft_bin_wall) / 2.0)) = '
    #        #print prob_interpolator(( (Rright_bin_wall + Rleft_bin_wall) / 2.0, (zright_bin_wall + zleft_bin_wall) / 2.0 ))
    #        new_probability = prob_interpolator(( (Rright_bin_wall + Rleft_bin_wall) / 2.0, (zright_bin_wall + zleft_bin_wall) / 2.0)) * (Rright_bin_wall - Rleft_bin_wall) * (zright_bin_wall - zleft_bin_wall)
    #        prob_array = prob_array + [new_probability]
    #        Rleft_bin_wall = Rright_bin_wall
    #    zleft_bin_wall = zright_bin_wall
    #    Rleft_bin_wall = Rbin_walls[0] 
            
    #print 'len(prob_array) = ' + str(len(prob_array))
    #print 'prob_array = '
    #print prob_array
    #total_prob = sum(prob_array)
    #print  'total_prob = ' + str(total_prob)
    prob_array = prob_array / total_prob
    prob_array = prob_array.tolist() 
    #integrated_prob_array_old = [ sum(prob_array[0:i+1]) for i in range(len(prob_array)) ]
    
    integrated_prob_array = np.zeros(len(prob_array))
    for i in range(len(prob_array)):
        if i > 0:
            integrated_prob_array[i] = prob_array[i] + integrated_prob_array[i-1]
        else:
            integrated_prob_array[0] = prob_array[i]
    #print 'Finished creating integrated_prob_array'

    Rz_points = (np.vstack(([Rbin_centers_mesh.T],[zbin_centers_mesh.T])).T).reshape(np.size(Rbin_centers_mesh),2)
    
    for star in range(NStars):
        #start = time.time()
        random_val = random.random()
        #print 'random_val = ' + str(random_val) 
        for i in range(len(integrated_prob_array)):
            if random_val <= integrated_prob_array[i]:
                star_RAs = star_RAs + [Rz_points[i][0]]
                star_Decs = star_Decs + [Rz_points[i][1]]
                #print 'Rz_points[i] = '
                #print Rz_points[i] 
                break
        #end = time.time()
        #print 'took ' + str(end-start) + ' seconds for star ' + str(star)

    return star_RAs, star_Decs 
    

    
