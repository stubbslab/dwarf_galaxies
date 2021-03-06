#Define the SurfaceBrightnessProfile class.
#This takes in a MassDensity class, integrates it along a los, and then creates an interpolator for the observed surface brightness.

from MassDensity import MassDensity
import numpy as np
import math
import time
from scipy.interpolate import RegularGridInterpolator
from DwarfGalaxyParametersStorer import DwarfGalaxyParametersStorer
from DwarfGalDataArchive import DwarfGalDataArchive
from AstronomicalParameterArchive import AstronomicalParameterArchive
from cantrips import getAreaMesh
#from generatePotentialInterplator import generatePotentialInterpolator

class SurfaceBrightnessProfile:

    #Sum the log of the value of the dSph surface brightness profile at every observed star position.
    #Must have corrRa and corrDec given in degrees from center of galaxy.
    def sumLogSurfaceBrightness(self, proj_x, proj_y):
        #startSum = time.time()
        astro_archive = AstronomicalParameterArchive()
        deg_to_rad = astro_archive.getDegToRad()
        proj_x_in_scale_radii = proj_x*deg_to_rad*self.dist/self.rs
        proj_y_in_scale_radii = proj_y*deg_to_rad*self.dist/self.rs

        log_likelihoods = np.log(self.onSkyInterpolator(np.dstack((proj_x_in_scale_radii,proj_y_in_scale_radii))[0]))
        #endSum = time.time()
        #print 'Took ' + str(endSum - startSum) + 's to do log sum.'
        return np.sum(log_likelihoods)

    #Sum the value of the dSph surface brightness profile at every observed star position
    #Must have corrRa and corrDec given in degrees from center of galaxy
    #Note that this quantity CANNOT serve as a good proxy for the likelihood of an array, since we're adding together values where a true probability would demand multiplication.
    def sumSurfaceBrightness(self,corrRa,corrDec):
        #print np.dstack(np.meshgrid(corrRa,corrDec))
        #print np.shape(np.dstack(np.meshgrid(corrRa,corrDec)))
        #print np.meshgrid(corrRa,corrDec)
        astro_archive = AstronomicalParameterArchive()
        deg_to_rad = astro_archive.getDegToRad()
        corr_Ra_in_scale_radii = corrRa*deg_to_rad*self.dist/self.rs
        corr_Dec_in_scale_radii = corrDec*deg_to_rad*self.dist/self.rs
        return np.sum(self.onSkyInterpolator(np.dstack((corr_Ra_in_scale_radii,corr_Dec_in_scale_radii))[0]))

    #Multiply the value of the dSph surface brightness profile at every observed star position
    #Must have corrRa and corrDec given in degrees from center of galaxy
    #Note that this quantity CAN serve as a proxy for probability, but this value can quickly grow too small for python to effectively handle well.
    def multiplySurfaceBrightness(self,corrRA,corrDec):
        astro_archive = AstronomicalParameterArchive()
        deg_to_rad = astro_archive.getDegToRad()
        corr_Ra_in_scale_radii = corrRa*deg_to_rad*self.dist/self.rs
        corr_Dec_in_scale_radii = corrDec*deg_to_rad*self.dist/self.rs
        return np.prod(self.onSkyInterpolator(np.dstack((corr_Ra_in_scale_radii,corr_Dec_in_scale_radii))[0]))

    #Return the value of the dSph surface brightness profile at every observed star position
    #Must have corrRa and corrDec given in degrees from center of galaxy
    def starProbabilityValues(self, corrRa, corrDec):
        astro_archive = AstronomicalParameterArchive()
        deg_to_rad = astro_archive.getDegToRad()
        corr_Ra_in_scale_radii = corrRa*deg_to_rad*self.dist/self.rs
        corr_Dec_in_scale_radii = corrDec*deg_to_rad*self.dist/self.rs
        return self.onSkyInterpolator(np.dstack((corr_Ra_in_scale_radii,corr_Dec_in_scale_radii))[0])

    #Normalized the surface probability, allowing for the possibility of appling a new mask and renormalizing
    def normalizeSurfaceProfile(self, mask, renormalizing = 0):

        if renormalizing:
            surface_prob = self.surfaceProbDensity_perSqrScaleRadius * mask.transpose()
        else:
            surface_prob = self.surfaceBrightness * mask.transpose()

        self.areaMesh_scaleRadii = getAreaMesh([self.projXMesh_scaleRadii, self.projYMesh_scaleRadii]).transpose()

        #RA_bin_lengths_scaleRadii = [RA_scaleRadii[1]-RA_scaleRadii[0]] + [(RA_scaleRadii[i+1]-RA_scaleRadii[i-1])/2.0 for i in range(len(RA_scaleRadii)-1)[1:]] + [RA_scaleRadii[len(RA_scaleRadii)-1] - RA_scaleRadii[(len(RA_scaleRadii)-2)] ]
        #Dec_bin_lengths_scaleRadii = [Dec_scaleRadii[1]-Dec_scaleRadii[0]] + [(Dec_scaleRadii[i+1]-Dec_scaleRadii[i-1])/2.0 for i in range(len(Dec_scaleRadii)-1)[1:]] + [Dec_scaleRadii[len(Dec_scaleRadii)-1] - Dec_scaleRadii[(len(Dec_scaleRadii)-2)] ]

        #self.areaMeshScaleRadii = (np.meshgrid(RA_bin_lengths_scaleRadii,Dec_bin_lengths_scaleRadii)[0] * np.meshgrid(RA_bin_lengths_scaleRadii,Dec_bin_lengths_scaleRadii)[1]).transpose()
        #This method turned out to be much slower.
        #areaMeshScaleRadii = np.array([[RA_length * Dec_length for Dec_length in Dec_bin_lengths_scaleRadii] for RA_length in RA_bin_lengths_scaleRadii])
        self.normalization_constant = np.sum(surface_prob * self.areaMesh_scaleRadii)
        #print 'normalization_constant = ' + str(normalization_constant)
        #for i in range(np.shape(surface_prob)[0]):
        #    for j in range(np.shape(surface_prob)[1]):
        #        prob_at_point = surface_prob[i][j]
        #        #print 'surface_prob[i][j] = ' + str(surface_prob[i][j])
        #        if i == 0:
        #            prob_at_point = prob_at_point * (Decmesh_scaleRadii[i+1][j]-Decmesh_scaleRadii[i][j])
        #        elif i == np.shape(surface_prob)[0]-1:
        #            prob_at_point = prob_at_point * (Decmesh_scaleRadii[i][j]-Decmesh_scaleRadii[i-1][j])
        #        else:
        #            prob_at_point = prob_at_point * (Decmesh_scaleRadii[i+1][j]-Decmesh_scaleRadii[i-1][j]) / 2.0
        #        if j == 0:
        #            prob_at_point = prob_at_point * (RAmesh_scaleRadii[i][j+1]-RAmesh_scaleRadii[i][j])
        #        elif j == np.shape(surface_prob)[1]-1:
        #            prob_at_point = prob_at_point * (RAmesh_scaleRadii[i][j]-RAmesh_scaleRadii[i][j-1])
        #        else:
        #            prob_at_point = prob_at_point * (RAmesh_scaleRadii[i][j+1]-RAmesh_scaleRadii[i][j-1]) / 2.0
        #        normalization_constant = normalization_constant + prob_at_point
        self.unnormalized_surface_prob = surface_prob
        self.surface_prob = surface_prob / self.normalization_constant
        self.surfaceProbDensity_perSqrScaleRadius = self.surface_prob
        self.surfaceProbDensity_perSqrRadians = self.surfaceProbDensity_perSqrScaleRadius * (self.dist / self.rs) ** 2.0 # surface probability per square radian

        #This gives an interpolator for surface probability density per square scale radius
        self.scaleRadiusInterpolator = RegularGridInterpolator((self.Rlikeaxis,self.zlikeaxis),self.surfaceProbDensity_perSqrScaleRadius,method = 'linear')
        #This gives an interpolator for surface probability density per square radian (not per square degree or scale radii, but square radian)

        self.onSkyInterpolator = RegularGridInterpolator((self.Rlikeaxis,self.zlikeaxis), self.surfaceProbDensity_perSqrRadians, method = 'linear')


    #The old method I was using the match up the stars on the sky.
    # Obviosly wrong, since it was taking the inner product of the RA and Dec lists, rather than just matching them up elementwise.
    #def starProbabilityValuesOld(self,corrRa,corrDec ):
    #    astro_archive = AstronomicalParameterArchive()
    #    deg_to_rad = astro_archive.getDegToRad()
    #    corr_Ra_in_scale_radii = corrRa*deg_to_rad*self.dist/self.rs
    #    corr_Dec_in_scale_radii = corrDec*deg_to_rad*self.dist/self.rs
    #    #print 'np.dstack(np.meshgrid(corr_Ra_in_scale_radii,corr_Dec_in_scale_radii))'
    #    #print np.dstack(np.meshgrid(corr_Ra_in_scale_radii,corr_Dec_in_scale_radii))
    #    print 'np.shape(np.dstack(np.meshgrid(corr_Ra_in_scale_radii,corr_Dec_in_scale_radii)))'
    #    print np.shape(np.dstack(np.meshgrid(corr_Ra_in_scale_radii,corr_Dec_in_scale_radii)))
    #    return self.onSkyInterpolator(np.dstack(np.meshgrid(corr_Ra_in_scale_radii,corr_Dec_in_scale_radii)))


    #A surface brightness profile is defined by an underlying massDensity, along with some angles that allow one to integrate the density along the los.
    #We then create an interpolator over the grid generated after integration along los
    # (note this is different than the graviational potential interpolator used in computing the mass density profile).
    #Here, R, z, and los_bins are all given in number of scale radii
    def __init__(self, R, z, los_bins, gamma, parameter_storer, disk_file=None, halo_file=None, C=1.0, disk_interpolating_function=None, halo_interpolating_function=None, population = None, observation_mask = None, masking_fields = [], halo_type = 'nfw', disk_type = 'sech_disk'):
        #parameter_storer.printContents()
        dwarf_archive = DwarfGalDataArchive()
        self.Rlikeaxis = R
        self.zlikeaxis = z
        self.los_bins = los_bins
        self.gamma = gamma
        self.C = C
        self.dist = parameter_storer.dist
        self.rs = parameter_storer.rs
        self.el = parameter_storer.el
        self.lam = parameter_storer.lam
        self.eps = parameter_storer.eps
        self.zeta = parameter_storer.zeta
        self.rs = parameter_storer.rs
        self.phi = parameter_storer.phi
        self.theta = parameter_storer.theta
        self.halo_sym_axis = parameter_storer.halo_sym_axis
        self.a = parameter_storer.a
        self.b = parameter_storer.b
        self.disk_sym_axis = parameter_storer.disk_sym_axis
        self.M = parameter_storer.M
        self.c = parameter_storer.c
        self.omegaphi = parameter_storer.omegaphi
        self.sigsqr_RR = parameter_storer.sigsqr_RR
        self.dispersion_b = parameter_storer.dispersion_b
        self.beta_int_dispersion_funct = parameter_storer.beta_int_dispersion_funct
        self.population = population
        self.observation_mask = observation_mask
        self.disk_center = parameter_storer.disk_center
        self.halo_center = parameter_storer.halo_center

        self.massDensity = MassDensity(self.Rlikeaxis, self.zlikeaxis, self.gamma, self.M, self.rs, self.eps, self.c, self.lam, self.zeta, self.el, self.omegaphi, self.C,
                                       self.sigsqr_RR, 
                                       dispersion_b = self.dispersion_b, beta_funct_over_param_int = self.beta_int_dispersion_funct, dist = self.dist, disk_file=disk_file, halo_file=halo_file,
                                       disk_interpolating_function=disk_interpolating_function, halo_interpolating_function=halo_interpolating_function,
                                       halo_type = halo_type, disk_type = disk_type, halo_center = self.halo_center, disk_center = self.disk_center)
        #endDefMass = time.time()
        #print 'Took ' + str(endDefMass - startDefMass) + 's to define massDensity.'
        #print 'In surface prob density:'
        #print 'self.massDensity.integrate_los(self.phi,self.theta,self.a,self.b,self.los_bins) = '
        #print self.massDensity.integrate_los(self.phi,self.theta,self.a,self.b,self.los_bins)

        #surfaceBrightness is the physical surface brightness expected, based on the PHYSICAL mass profile generated.
        #startIntegrate = time.time()
        self.surfaceBrightness = self.massDensity.integrate_los(self.halo_sym_axis, self.disk_sym_axis, self.los_bins)

        #endIntegrate = time.time()
        #print 'Took ' + str(endIntegrate - startIntegrate) + 's to integrate los.'
        #print 'np.shape(self.surfaceBrightness) = '
        #print np.shape(self.surfaceBrightness)
        #Normalize surface brightness to unity?  Now is this truly a normalization, since I'm only summing over points I've measured?  I don't think so...
        #I should normalize somehow.
        #Since this is a SURFACE PROBABILITY DENSITY, I need to multiply by the surface area of the portion of the sky that each point describes.  Yes?
        # Here, we normalize the surface density to get probability density per square scale radius. That normalization is convenient for ...
        # But for comparison of two profiles on the sky, I need to normalize by angular size on sky
        # (that way, I'm guarenteed to have each differential area element the same).
        #Rmesh, zmesh = np.meshgrid(self.Rlikeaxis,self.zlikeaxis)
        #RAmesh,Decmesh are defined here in radians on sky
        #startDefInterp = time.time()
        self.projXMesh_radians, self.projYMesh_radians = np.meshgrid(self.Rlikeaxis * self.rs / self.dist, self.zlikeaxis * self.rs / self.dist)
        projX_scaleRadii = self.Rlikeaxis
        projY_scaleRadii = self.zlikeaxis
        self.projXMesh_scaleRadii, self.projYMesh_scaleRadii = np.meshgrid(projX_scaleRadii, projY_scaleRadii)
        normalization_constant = 0.0

        #surface_probability is the OBSERVED probability of finding a star at a given position, which must account for the incomplete sky sampling region
        if self.observation_mask is None:
            if population is None:
                self.observation_mask = np.zeros(np.shape(self.surfaceBrightness)) + 1.0
            else:
                self.observation_mask = GalaxyMask(self.projXMesh_radians, self.projYMesh_radians, population[0], mask_types = ['n_vel_meas'], masking_fields = masking_fields, units_of_coord_meshes = 'radians')

        self.normalizeSurfaceProfile(self.observation_mask, renormalizing = 0)
        #print self.normalization_constant
