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
        #print ('[los_velocity, x, y, z] = ' + str([los_velocity, x, y, z] ))
        r_vec = np.array([x, y, z]) #- self.halo_center
        center_theta, center_offset, center_phi = self.halo_center
        #observer_vec = np.array([0.0, - center_offset - self.dist, 0.0])
        observer_vec = np.array(observer_vec)
        los_vec = r_vec - observer_vec
        #los_vec = np.array([0.0, self.dist, 0.0]) + np.array(self.halo_center) + np.array([x, y, z]) * self.rs #the line of sight vector, in cartesian units in the light-center of the dSph reference frame
        #los_vec = np.array([0.0, self.dist, 0.0]) + r_vec
        los_mag = np.sqrt(np.sum([elem ** 2.0 for elem in los_vec]))
        r_mag = np.sqrt(np.sum([elem ** 2.0 for elem in r_vec]))
        #If the r magnitude is 0.0, there is a degenereate definition of the
        # r direction.  So just assume it's along the line of sight.
        if r_mag > 0.0:
            los_hat_dot_r_hat = np.sum(r_vec * los_vec) / (los_mag * r_mag)
        else:
            los_hat_dot_r_hat = 1.0
        los_unit_r = los_hat_dot_r_hat
        los_unit_t = np.sqrt(1.0 - los_unit_r ** 2.0 )
        if np.isnan(los_unit_t):
            print ('[r_vec, observer_vec] = ' + str([r_vec, observer_vec])) 
            print ('[los_vec, r_vec, los_mag, r_mag] = ' + str([los_vec, r_vec, los_mag, r_mag]))
            print ('[los_unit_t, los_unit_r, (1.0 - los_unit_r ** 2.0)] = ' + str([los_unit_t, los_unit_r, (1.0 - los_unit_r ** 2.0)]))
        los_r = los_unit_r * los_velocity
        los_t = los_unit_t * los_velocity
        #los_r = np.sum(np.array(los_vec) * r_vec) / (r_mag) # the dot_product of r_hat and los
        #los_t = los_mag ** 2.0 - los_r ** 2.0
        sigsqr_rr = self.radialDispersionFunct (r_mag)
        beta_dispersion = self.betaDispersionFunct (r_mag)
        sigsqr_tangent = sigsqr_rr * (1.0 - beta_dispersion)
        ellipsoidal_axis_ratio = np.sqrt(sigsqr_tangent / sigsqr_rr)
        ellipsoidal_eccentricity = ellipsoidal_axis_ratio ** 2.0 #my \xi
        #print ('[los_velocity, sigsqr_rr, sigsqr_tangent] = ' + str([los_velocity, sigsqr_rr, sigsqr_tangent] ))

        scaled_d = 1.0 / (ellipsoidal_axis_ratio ** 2.0 * los_unit_t ** 2.0 + los_unit_r ** 2.0)
        beta_part1 = los_unit_t ** 2.0 * (1.0 + 1.0 / (ellipsoidal_axis_ratio ** 2.0)) + los_unit_r ** 2.0 * 2.0 / (ellipsoidal_axis_ratio ** 2.0)
        beta_part2 = los_unit_t ** 2.0 / (ellipsoidal_axis_ratio ** 2.0) + los_unit_r ** 2.0 / (ellipsoidal_axis_ratio ** 4.0)
        scaled_beta_plus = (beta_part1 + np.sqrt(max(0.0, beta_part1 ** 2.0  - 4.0 * beta_part2))) / 2.0
        scaled_beta_minus = (beta_part1 - np.sqrt(max(0.0, beta_part1 ** 2.0  - 4.0 * beta_part2)))/ 2.0
        if np.isnan(scaled_beta_plus * scaled_beta_minus): print ('[scaled_beta_plus, scaled_beta_minus, beta_part1 ** 2.0, 4.0 * beta_part2] = ' + str([scaled_beta_plus, scaled_beta_minus, beta_part1 ** 2.0, 4.0 * beta_part2]))
        #print ('[scaled_d, beta_part1, beta_part2, scaled_beta_plus, scaled_beta_minus] = ' + str([scaled_d, beta_part1, beta_part2, scaled_beta_plus, scaled_beta_minus]))

        #The minimum velocity vector that flags an ellipse that just touches the plane of possible velocity vectors
        #vel_intersection_min = los_velocity ** 2.0 / np.sqrt(sigsqr_rr * los_unit_r ** 2.0 + sigsqr_tangent * los_unit_t ** 2.0)
        vel_intersection_min = scaled_d
        #vel_intersection_min = 0.0
        #print ('vel_intersection_min = ' + str(vel_intersection_min))

        #print ('[r_vec, los_vec, los_mag, r_mag, los_hat_dot_r_hat, los_r, los_t, sig_rr, beta_dispersion, sig_tangent, ellipsoidal_eccentricity] = ' + str([r_vec, los_vec, los_mag, r_mag, los_hat_dot_r_hat, los_r, los_t, sig_rr, beta_dispersion, sig_tangent, ellipsoidal_eccentricity]))

        #The following derivation of the ellipse from the intersection of a plane and an ellipsoid is based on:
        # Section 6 of https://file.scirp.org/pdf/AM20121100009_89014420.pdf

        #scaled_ellipsoidal_a = np.sqrt((1.0 - los_hat_dot_r_hat ** 2.0) + ellipsoidal_axis_ratio ** 2.0 * los_hat_dot_r_hat ** 2.0 )
        #scaled_ellipsoidal_c = scaled_ellipsoidal_a / ellipsoidal_axis_ratio
        #print ('[scaled_ellipsoidal_a, scaled_ellipsoidal_c] = ' + str([scaled_ellipsoidal_a, scaled_ellipsoidal_c]))
        #ellipsoidal_a_funct = lambda dist_perp_vel: scaled_ellipsoidal_a * dist_perp_vel
        #ellipsoidal_c_funct = lambda dist_perp_vel: scaled_ellipsoidal_c * dist_perp_vel
        #ellipsoidal_a_funct = lambda dist_perp_vel: dist_perp_vel * np.sqrt(1.0 - eccentricity ** 2.0 * los_hat_dor_r_hat ** 2.0)  #???
        #ellipsoidal_c_funct = lambda dist_perp_vel: ellipsoidal_a_funct(dist_perp_vel) * eccentricity #???

        #scaled_d = los_velocity ** 4.0 / ( (los_t * scaled_ellipsoidal_a) ** 2.0 + (los_r * scaled_ellipsoidal_c) ** 2.0) if abs(los_velocity) > 0.0 else 0.0
        #d_funct = lambda dist_perp_vel: 1.0 / dist_perp_vel ** 2.0 * scaled_d if abs(dist_perp_vel) > 0.0 else 0.0

        #scaled_beta_term1 = 1.0
        #scaled_beta_term2 = -(los_unit_t ** 2.0 * (1.0 / scaled_ellipsoidal_a ** 2.0 + 1.0 / scaled_ellipsoidal_c ** 2.0) + los_unit_r ** 2.0 * (2.0 / scaled_ellipsoidal_a** 2.0))
        #scaled_beta_term3 = (los_unit_r ** 2.0 / (scaled_ellipsoidal_a ** 4.0) + los_unit_t ** 2.0 / (scaled_ellipsoidal_a ** 2.0 * scaled_ellipsoidal_c ** 2.0))
        #print ('[scaled_d, scaled_beta_term1, scaled_beta_term2,scaled_beta_term3 ] = ' + str([scaled_d, scaled_beta_term1, scaled_beta_term2,scaled_beta_term3 ]))

        #scaled_beta_plus = (- scaled_beta_term2 + np.sqrt(scaled_beta_term2 ** 2.0 - 4.0 * scaled_beta_term1 * scaled_beta_term3) ) / (2.0 * scaled_beta_term1)
        #scaled_beta_minus = (- scaled_beta_term2 - np.sqrt(scaled_beta_term2 ** 2.0 - 4.0 * scaled_beta_term1 * scaled_beta_term3) ) / (2.0 * scaled_beta_term1)
        #print ('[scaled_d, scaled_beta_term1, scaled_beta_term2,scaled_beta_term3, scaled_beta_plus, scaled_beta_minus] = ' + str([scaled_d, scaled_beta_term1, scaled_beta_term2,scaled_beta_term3, scaled_beta_plus, scaled_beta_minus ]))
        #beta_plus_funct = lambda dist_perp_vel: 1.0 / dist_perp_vel ** 2.0 * scaled_beta_plus if abs(dist_perp_vel) > 0.0 else 0.0
        #beta_minus_funct = lambda dist_perp_vel: 1.0 / dist_perp_vel ** 2.0 * scaled_beta_minus if abs(dist_perp_vel) > 0.0 else 0.0

        #The eccentricity of the ellipse formed by intersecting the the plane with the ellipsoid
        #ellipse_eccentricity = np.sqrt(1 - scaled_beta_minus / scaled_beta_plus)

        #ellipse_length_integral = integrate.quad(lambda t_int: np.sqrt(1.0 - ellipse_eccentricity ** 2.0 * np.sin(t_int) ** 2.0), 0.0, np.pi / 2.0)[0]


        #print ('[d_funct(2 * los_velocity), beta_plus_funct(2 * los_velocity)] = ' + str([d_funct(2 * los_velocity), beta_plus_funct(2 * los_velocity)]))
        #ellipse_length_funct = lambda dist_perp_vel: (4.0 * np.sqrt((1.0 - d_funct(dist_perp_vel) ) / beta_plus_funct(dist_perp_vel)) * ellipse_length_integral
        #                                              if (abs(d_funct(dist_perp_vel)) < 1.0 and abs(dist_perp_vel) > 0.0) else 0.0 )
        #print('[ellipse_length_funct(2 * los_velocity), d_funct(2 * los_velocity), beta_plus_funct(2 * los_velocity), ellipse_length_integral] = ' + str([ellipse_length_funct(2 * los_velocity), d_funct(2 * los_velocity), beta_plus_funct(2 * los_velocity), ellipse_length_integral]))

        #radial_prob_funct = lambda v_r: stats.norm.pdf(v_r, 0.0, sigsqr_rr)  #Probability the radial velocity vector is v_r
        #tangential_prob_funct = lambda v_phi:  1 / (np.sqrt(2.0 * np.pi) * sigsqr_tangent) * stats.norm.pdf(v_phi, 0.0, sigsqr_tangent) * v_phi

        #prob_density_funct = lambda dist_perp_vel: radial_prob_funct(ellipsoidal_c_funct(dist_perp_vel)) * tangential_prob_funct(ellipsoidal_a_funct(dist_perp_vel))

        #prob_on_ellipse_funct = lambda los_vel_len: prob_density_funct((1 if los_velocity >= 0.0  else -1) * los_vel_len) * ellipse_length_funct((1 if los_velocity >= 0.0 else -1) * los_vel_len)
        #total_prob_of_star_having_losv_at_position =  (1 if los_velocity >= 0.0  else -1) * integrate.quad(prob_on_ellipse_funct, max(vel_int_lower_bound, vel_intersection_min), vel_int_upper_bound) [0]
        #total_prob_of_star_having_losv_at_position = np.sum(int_funct (los_x_vels_mesh, los_y_vels_mesh)) * (int_mesh_area)

        prob_scaling = 1.0 / np.sqrt(scaled_beta_plus * scaled_beta_minus) * np.sqrt( sigsqr_rr  / (2.0 * np.pi * sigsqr_tangent ** 2.0))
        if np.isnan(prob_scaling): print ('[scaled_beta_plus, scaled_beta_minus, sigsqr_tangent, prob_scaling] = ' + str([scaled_beta_plus, scaled_beta_minus, sigsqr_tangent, prob_scaling]))
        prob_term1 = np.exp(-scaled_d * (los_velocity) ** 2.0 / (2.0 * sigsqr_rr))
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
        #self.x_sky_axis = x
        #self.z_sky_axis = z
        #self.los_bins = los_bins
        self.gamma = astro_archive.getGamma()
        self.C = C
        self.dist = parameter_storer.dist
        self.rs = parameter_storer.rs
        self.sigsqr_rr_0 = parameter_storer.sigsqr_rr_0
        self.sigsqr_rr_inf = parameter_storer.sigsqr_rr_inf
        self.r_sigsqr_rr0 = parameter_storer.r_sigsqr_rr0
        self.alpha_sigsqr_rr = parameter_storer.alpha_sigsqr_rr
        self.gamma_for_beta_inf = parameter_storer.gamma_for_beta_inf
        self.r_beta0 = parameter_storer.r_beta0
        self.alpha_beta = parameter_storer.alpha_beta
        #self.population = population
        #self.halo_type = halo_type
        self.halo_center = parameter_storer.halo_center
        self.dispersion_rr_params = [self.sigsqr_rr_0, self.sigsqr_rr_inf, self.r_sigsqr_rr0, self.alpha_sigsqr_rr]
        self.dispersion_beta_params = [self.gamma_for_beta_inf, self.r_beta0, self.alpha_beta]

        #self.massDensity = SMassDensity(self.x_sky_axis, self.z_sky_axis, self.M, self.rs, self.dist, self.c, self.C, self.halo_type ,
        #                                dispersion_rr_params = self.dispersion_rr_params, dispersion_beta_params = self.dispersion_beta_params,
        #                                halo_center = self.halo_center)

        self.radialDispersionFunct = sphVDF.returnSigsqrRRFunct(self.dispersion_rr_params)
        self.betaDispersionFunct = sphVDF.returnBetaFunct(self.dispersion_beta_params)
