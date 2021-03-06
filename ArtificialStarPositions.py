import numpy as np
import math
import time 
from AstronomicalParameterArchive import AstronomicalParameterArchive
from PotentialArchive import PotentialArchive
from distributeRandomStars import distributeRandomStars
import matplotlib.pyplot as plt
from compRZLos import compRZLos
from SurfaceBrightnessProfile import SurfaceBrightnessProfile
from DwarfGalDataArchive import DwarfGalDataArchive 

class ArtificialStarPositions:

    def getStarPositions(self):
        return self.artificial_star_RAs_in_deg, self.artificial_star_Decs_in_deg

    def plotStarPositions(self): 
        plt.scatter(self.artificial_star_RAs_in_deg, self.artificial_star_Decs_in_deg)
        x_lims = [x_lim / 60.0 for x_lim in self.arcmin_limits_R]
        y_lims = [y_lim / 60.0 for y_lim in self.arcmin_limits_z]
        plt.ylim = y_lims 
        plt.xlim = x_lims 

        plt.xlabel('ra sep')
        plt.ylabel('dec sep')
        plt.title('Artificial Stars Position on Sky')
        #pltFileDir = computational_params_archive.getPlotDir() 
        #Uncomment the following line if you want to save the file
        #plt.savefig(pltFileDir + gal_of_interest + '_' + pop_of_interest + '_metal_cuts_' + str((metallicity_cuts[gal_of_interest])[0]) + '_' + str((metallicity_cuts[gal_of_interest])[1]) + '.png')
        #plt.savefig(pltFileDir + gal_of_interest + '_' + pop_of_interest + '_full_pop_division' + '.png')
        plt.axis('scaled')
        plt.axis(x_lims + y_lims) #combined limits from Walker plot 
        plt.show()
        return 0
        
    def __init__(self, N_stars, generation_params_storer, max_R_bins, max_z_bins, outer_step,
                 arcmin_limits_R, arcmin_limits_z, arcmin_limits_los, generation_disk_funct, generation_halo_funct, withDisk):
        self.N_stars = N_stars
        self.max_R_bins = max_R_bins
        self.max_z_bins = max_z_bins
        self.outer_step = outer_step
        self.arcmin_limits_R =  arcmin_limits_R
        self.arcmin_limits_z =  arcmin_limits_z
        self.arcmin_limits_los =  arcmin_limits_los
        
        print ('Creating artificial stars... ')
        start_art = time.time() 


        astro_archive = AstronomicalParameterArchive()
        pot_archive = PotentialArchive()
        dSph_archive = DwarfGalDataArchive()

        deg_to_rad = astro_archive.getDegToRad()
        gamma = astro_archive.getGamma()
    
        el = generation_params_storer.el
        lam = generation_params_storer.lam
        rs = generation_params_storer.rs
        dist = generation_params_storer.dist
        zeta = generation_params_storer.zeta 

        print ('sigsqr = ' + str(generation_params_storer.sigsqr))

        #generation_R,generation_z,generation_los = compRZLos(zeta, zeta * outer_step, outer_step, arcmin_limits_R, arcmin_limits_z, arcmin_limts_los, rs, dist )
        generation_R,generation_z,generation_los = compRZLos(1.0, 1.0 * outer_step, outer_step, arcmin_limits_R, arcmin_limits_z, arcmin_limits_los, rs, dist ) 

        current_observation_mask = dSph_archive.getObservationMask(generation_params_storer.population, generation_R * rs / dist * 1.0 / deg_to_rad, generation_z * rs / dist * 1.0 / deg_to_rad ) 
        generation_surface_profile = SurfaceBrightnessProfile(generation_R,generation_z,generation_los, 
                                                              gamma, generation_params_storer,
                                                              disk_interpolating_function = generation_disk_funct,
                                                              halo_interpolating_function = generation_halo_funct,
                                                              observation_mask = current_observation_mask ) 
        #print 'Distributing random stars... '
        #print 'max_R_bins = '
        #print max_R_bins
        #print 'max_z_bins = '
        #print max_z_bins 
        artificial_star_RAs_in_scale_radii, artificial_star_Decs_in_scale_radii = distributeRandomStars(N_stars, generation_surface_profile, max_R_bins, max_z_bins, withDisk)
        #print 'Rescaling RA and Decs of random stars...'
        artificial_star_RAs_in_deg = [RA * rs / dist * 1 / deg_to_rad for RA in artificial_star_RAs_in_scale_radii]
        artificial_star_Decs_in_deg = [Dec * rs / dist * 1 / deg_to_rad for Dec in artificial_star_Decs_in_scale_radii]

        end_art = time.time()
        print ('Took ' + str(end_art - start_art) + 's to generate artificial stars. ') 
        self.artificial_star_RAs_in_deg = artificial_star_RAs_in_deg
        self.artificial_star_Decs_in_deg = artificial_star_Decs_in_deg
    
