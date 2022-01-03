#Here we define a python object that basically includes all of the information on a dSph observation.
#It is a repository class; the data is actually read in using another function.

from readInGalaxyStarData import readInGalaxyStarData
from AstronomicalParameterArchive import AstronomicalParameterArchive
import math
import numpy as np
from DwarfGalDataArchive import DwarfGalDataArchive 


class ObservedGalaxyStarData:

    #The geometric position of a star determines how it's los velocity (the velocity I have) relates to v_R, v_z, v_\phi
    #So far, we just assume \sig_RR ^ 2 = \sig_los ^ 2 .  I don't think that's right, but I don't know a better strategy right now...
    def measureDispersionSignal(self, population, dist, halo_center, halo_sym_axis, dispersion_b, gamma0, beta_exp, beta_R1, beta_R2 ):
        #proj_pos = [np.array((self.proj_x[i] * deg_to_rad * dist - halo_center[0], 0.0, self.proj_y[i] * deg_to_rad * dist - halo_center[2])) for i in range(len(self.proj_x))]
        data_archive = DwarfGalDataArchive()
        sigsqr_los = data_archive.getVelocityDispersion(population)
        sigsqr_RR = sigsqr_los #This probably isn't right, but I have no better idea right now.
        return sigsqr_RR

    def measureRotationSignal(self, dist, halo_center, halo_sym_axis):
        astro_arch = AstronomicalParameterArchive()
        deg_to_rad = astro_arch.getDegToRad()
        pc_to_m = astro_arch.getParsecToM()
        rot_axis = np.array(halo_sym_axis)
        #print 'dist = ' + str(dist)
        #print 'halo_center = ' + str(halo_center)
        #print 'halo_sym_axis = ' + str(halo_sym_axis)
        #print 'self.proj_x = ' + str(self.proj_x)
        #print 'self.proj_y = ' + str(self.proj_y)
        #print      [np.array((self.proj_x[i] * deg_to_rad * dist - halo_center[0], 0.0, self.proj_y[i] * deg_to_rad * dist - halo_center[2])) for i in range(len(self.proj_x))]
        proj_pos = [np.array((self.proj_x[i] * deg_to_rad * dist - halo_center[0], 0.0, self.proj_y[i] * deg_to_rad * dist - halo_center[2])) for i in range(len(self.proj_x))]
        #print 'proj_pos = '  str(proj_pos)
        parallel_vecs = [rot_axis * sum(rot_axis * vec) for vec in proj_pos]
        #print 'parallel_vecs = ' + str(parallel_vecs)
        orthog_vecs = [proj_pos[i] - parallel_vecs[i] for i in range(len(parallel_vecs))]
        orthog_seps = [np.cross(vec, rot_axis)[1] for vec in orthog_vecs]
        orthog_seps_km = [sep * pc_to_m * 10.0 ** -3.0 for sep in orthog_seps]
        #print 'orthog_seps = ' + str(orthog_seps)
        star_vlos = self.corr_Vhel
        star_vlos_errs = self.corr_VhelE
        #print 'star_vlos = ' + str(star_vlos)
        line_fit = np.polyfit(orthog_seps_km, star_vlos, 1, w = [1.0 / err for err in star_vlos_errs])
        #print 'line_fit = ' + str(line_fit)
        #print 'sum(np.array(star_vlos * orthog_seps)) = ' + str(sum(np.array(star_vlos) * orthog_seps))
        return line_fit[0]

    def __init__(self, population, pop_selection_method='none'):
        #data_archive = DwarfGalDataArchive()
        #self.dataFile = dataArchive.getFile(population)
        #self.dist = dataArchive.getDistance(population)
        #self.dist = data_archive.getDistanceFromSun(population) #distance from sun  in pc {"carina":105000,"fornax":147000,"sculptor":86000,"sextans":86000}
        #self.M = data_archive.getTotalMass(population) # total mass in M_sun {"carina":1.28*10**8,"fornax":1.28*10**8,"sculptor":1.28*10**8,"sextans":1.28*10**8}
        #self.gal_center_tuple = data_archive.getRaDec(population)
        #(self.dist, self.M) = ExtraGalDataArchive.readInExtraGalData(population, pop_selection_method)
        ( self.Target,  self.Field,  self.Date,   self.RA,             self.DEC,        self.Vmag,
          self.VImag,   self.Vhel,   self.VhelE,  self.corr_Vhel,      self.corr_VhelE, self.SigFe,
          self.SigFeE,  self.SigMg,  self.SigMgE, self.V_Intensities,  self.SigMg_corr, self.corrRa,
          self.corrDec, self.proj_x, self.proj_y, self.half_light_rad, self.sigSqr,     self.sigSqrE
          ) = readInGalaxyStarData(population,pop_selection_method)
        self.membership_fraction = 1.0
