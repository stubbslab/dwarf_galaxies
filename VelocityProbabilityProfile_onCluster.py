#Define the VelocityProbabilityProfile class.
#This takes in a MassDensity class for a dwarf spheroidal galaxy and observations of star positions and line of sight velocites.
#From this, it calculates the likelihood that the stars would have the line of sight velocities that we measure,
# given the likelihood that they would be at a particular location on the sky and the likelihood that they would
# have the measured los velocity.

from MassDensitySpherical import SMassDensity
import numpy as np
import math
import time
from scipy.interpolate import RegularGridInterpolator
from SphericalGalaxyParametersStorer import DwarfGalaxyParametersStorer
from AstronomicalParameterArchive import AstronomicalParameterArchive
from DwarfGalDataArchive import DwarfGalDataArchive
import AstronomicalParameterArchive as apa
from cantrips import getAreaMesh
import sphericalVelocityDispersionFuncts as sphVDF
import scipy.stats as stats
import scipy.integrate as integrate
import scipy.special as special

class VelocityProbabilityCalculator:

    #Computes the likelihood of measuring a stellar line of sight velocity
    # given the star's positions in the galaxy.
    #We'll need to integrate along the line of sight since we don't
    # actually know where the star is.
    #Assumes that the 3d vectors (x, y, z) are given in the same
    # units as the scale radii and distance (typically pc) and that
    # they are with reference to the center of the dark matter halo,
    # offset from the light center of the galaxy by my halo center variable.
    #i.e. The x, y, and z are already given the galaxy's natural coordinate system
    #
    #For the purposes of calculating the intersections of the velocity plane with the equi-probability ellipsoidals, we should place teh vectors in r, x, y coordinates.
    #Returns in units of prob / (km/s), assuming the kinematic parameters are all specified in km/s
    def computeProbabilityOfLOSVel(self, los_velocity, x, y, z, observer_vec, vel_int_lower_bound = 0.0, vel_int_upper_bound = np.inf ):
        #los_velocity = abs(los_velocity)
        r_vec = np.array([x, y, z]) #- self.halo_center
        center_theta, center_offset, center_phi = self.halo_center
        #observer_vec = np.array([0.0, - center_offset - self.dist, 0.0])
        observer_vec = np.array(observer_vec)
        los_vec = r_vec - observer_vec
        #los_vec = np.array([0.0, self.dist, 0.0]) + np.array(self.halo_center) + np.array([x, y, z]) * self.rs #the line of sight vector, in cartesian units in the light-center of the dSph reference frame
        #los_vec = np.array([0.0, self.dist, 0.0]) + r_vec
        los_mag = np.sqrt(np.sum([elem ** 2.0 for elem in los_vec]))
        if los_mag > 0.0:
            los_unit = los_vec / los_mag
        else:
            los_unit = [0.0, 0.0, 1.0]
        r_mag = np.sqrt(np.sum([elem ** 2.0 for elem in r_vec]))

        #If the r magnitude is 0.0, there is a degenereate definition of the
        # r direction.  So just assume it's along the line of sight.

        if r_mag > 0.0:
            los_hat_dot_r_hat = np.sum(r_vec * los_vec) / (los_mag * r_mag)
        else:
            los_hat_dot_r_hat = 1.0
        los_unit_r = los_hat_dot_r_hat
        los_unit_t = np.sqrt(1.0 - los_unit_r ** 2.0 )
        los_r = los_unit_r * los_velocity
        los_t = los_unit_t * los_velocity

        #print ('[r_vec, observer_vec, los_vec, los_mag, los_unit, r_mag, los_hat_dot_r_hat, los_unit_r, los_unit_t, los_r, los_t] = ' + str([r_vec, observer_vec, los_vec, los_mag, los_unit, r_mag, los_hat_dot_r_hat, los_unit_r, los_unit_t, los_r, los_t] ))
        if np.isnan(los_t):
            print ('[los_t, los_r, observer_vec, los_vec, r_vec] = ' + str([los_t, los_r, observer_vec, los_vec, r_vec]))
        #los_r = np.sum(np.array(los_vec) * r_vec) / (r_mag) # the dot_product of r_hat and los
        #los_t = los_mag ** 2.0 - los_r ** 2.0
        sigsqr = self.sigsqr
        omega_phi = self.omega_phi
        ellipsoidal_axis_ratio = 1.0
        ellipsoidal_eccentricity = ellipsoidal_axis_ratio ** 2.0
        #R = np.sqrt(x ** 2.0 + y ** 2.0)

        pc_to_m = self.astro_arch.getParsecToM()
        expected_velocity_vector = np.array([-y, x, 0.0]) * pc_to_m * 10.0 ** (-3.0) * omega_phi   #the rotation direction of the star : rot_unit_vector * v_rot = [-y,x,0] / R * omega_phi * R = [-y,x,0] * omega_phi
        expected_v_dot_los = np.sum(expected_velocity_vector * los_unit)
        kappa = abs(los_velocity - expected_v_dot_los)

        #print ('[self.sigsqr, self.omega_phi, ellipsoidal_axis_ratio')
        scaled_d = 1.0 / (ellipsoidal_axis_ratio ** 2.0 * los_unit_t ** 2.0 + los_unit_r ** 2.0)
        beta_part1 = los_unit_t ** 2.0 * (1.0 + 1.0 / (ellipsoidal_axis_ratio ** 2.0)) + los_unit_r ** 2.0 * 2.0 / (ellipsoidal_axis_ratio ** 2.0)
        beta_part2 = los_unit_t ** 2.0 / (ellipsoidal_axis_ratio ** 2.0) + los_unit_r ** 2.0 / (ellipsoidal_axis_ratio ** 4.0)
        scaled_beta_plus = (beta_part1 + np.sqrt(max(0.0, beta_part1 ** 2.0  - 4.0 * beta_part2))) / 2.0
        scaled_beta_minus = (beta_part1 - np.sqrt(max(0.0, beta_part1 ** 2.0  - 4.0 * beta_part2)))/ 2.0
        #print ('[self.sigsqr, self.omega_phi, scaled_d, beta_part1, beta_part2, scaled_beta_plus, scaled_beta_minus] = ' + str([self.sigsqr, self.omega_phi, scaled_d, beta_part1, beta_part2, scaled_beta_plus, scaled_beta_minus] ))

        #The minimum velocity vector that flags an ellipse that just touches the plane of possible velocity vectors
        #vel_intersection_min = los_velocity ** 2.0 / np.sqrt(sigsqr_rr * los_unit_r ** 2.0 + sigsqr_tangent * los_unit_t ** 2.0)
        vel_intersection_min = scaled_d

        prob_scaling = 1.0 / np.sqrt(scaled_beta_plus * scaled_beta_minus) * np.sqrt( sigsqr  / (2.0 * np.pi * sigsqr ** 2.0))
        #if np.isnan(prob_scaling): print ('[scaled_beta_plus, scaled_beta_minus, sigsqr_tangent] = ' + str([scaled_beta_plus, scaled_beta_minus, sigsqr]))
        prob_term1 = np.exp(-scaled_d * (kappa) ** 2.0 / (2.0 * sigsqr))
        #prob_term2 = -np.sqrt(np.pi) * scaled_d / 2.0 * np.sqrt(sigsqr_rr) / np.abs(los_velocity) * special.erfc(scaled_d * los_velocity / np.sqrt(sigsqr_rr))
        prob_term2 = 0.0
        total_prob_of_star_having_losv_at_position  = prob_scaling * (prob_term1 + prob_term2)
        #print('[prob_scaling, prob_term1, prob_term2 ] = ' + str([prob_scaling, prob_term1, prob_term2 ]))
        #print ('[prob_scaling, prob_term1, prob_term2] = ' + str([prob_scaling, prob_term1, prob_term2]))
        return total_prob_of_star_having_losv_at_position

    #def __init__(self, parameter_storer, halo_type, disk_file=None, halo_file=None, C=1.0, observation_mask = None, masking_fields = [], population = None):
    def __init__(self, parameter_storer, C = 1.0, observation_mask = None, masking_fields = [], population = None):
        dwarf_archive = DwarfGalDataArchive()
        astro_archive = apa.AstronomicalParameterArchive()
        self.astro_arch = astro_archive
        #self.x_sky_axis = x
        #self.z_sky_axis = z
        #self.los_bins = los_bins
        sigsqr = parameter_storer.sigsqr_RR
        omega_phi = parameter_storer.omega_phi
        self.gamma = astro_archive.getGamma()
        self.C = C
        self.dist = parameter_storer.dist
        self.rs = parameter_storer.rs
        self.sigsqr = sigsqr
        self.omega_phi = omega_phi
        #self.population = population
        #self.halo_type = halo_type
        self.halo_center = parameter_storer.halo_center

        #self.massDensity = SMassDensity(self.x_sky_axis, self.z_sky_axis, self.M, self.rs, self.dist, self.c, self.C, self.halo_type ,
        #                                dispersion_rr_params = self.dispersion_rr_params, dispersion_beta_params = self.dispersion_beta_params,
        #                                halo_center = self.halo_center)
