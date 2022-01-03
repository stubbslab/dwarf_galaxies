#Here we define a python object that basically includes all of the information on a dSph observation.
#It is a repository class; the data is actually read in using another function. 

from readInBoylanArtificialGalaxyData import readInBoylanArtificialStarData
import math

class BoylanArtificialGalaxyStarData:
    def __init__(self, halo_number, viewer_phi = 0.0, viewer_theta = 0.0, dist_from_sun = None, arcmin_limits_R = None, arcmin_limits_z = None, inclusion_range = None):
        #data_archive = DwarfGalDataArchive()
        #self.dataFile = dataArchive.getFile(population)
        #self.dist = dataArchive.getDistance(population)
        #self.dist = data_archive.getDistanceFromSun(population) #distance from sun  in pc {"carina":105000,"fornax":147000,"sculptor":86000,"sextans":86000}
        #self.M = data_archive.getTotalMass(population) # total mass in M_sun {"carina":1.28*10**8,"fornax":1.28*10**8,"sculptor":1.28*10**8,"sextans":1.28*10**8} 
        #self.gal_center_tuple = data_archive.getRaDec(population) 
        #(self.dist, self.M) = ExtraGalDataArchive.readInExtraGalData(population, pop_selection_method)
        (self.xs,              self.ys,          self.zs,
         self.viewed_xs,       self.viewed_ys,   self.viewed_zs,
         self.viewed_RAs,      self.viewed_Decs,
         self.corr_ra,         self.corr_dec,
         self.all_proj_x,      self.all_proj_y, 
         self.proj_x,          self.proj_y,
         self.los_vals,        self.sigSqr,
         self.arcmin_limits_R, self.arcmin_limits_z, self.dist) = readInBoylanArtificialStarData(halo_number,
                                                                                                 viewer_phi = viewer_phi, viewer_theta = viewer_theta,
                                                                                                 dist_from_sun = dist_from_sun,
                                                                                                 arcmin_limits_R = arcmin_limits_R,  arcmin_limits_z = arcmin_limits_z,
                                                                                                 inclusion_range = inclusion_range) 
