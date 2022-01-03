#Define the SurfaceBrightnessProfile class.
#This takes in a MassDensity class, integrates it along a los, and then creates an interpolator for the observed surface brightness.

from MassDensitySpherical import SMassDensity
import numpy as np
import math
import time
import scipy.interpolate as interpolate
from SphericalGalaxyParametersStorer import DwarfGalaxyParametersStorer
from AstronomicalParameterArchive import AstronomicalParameterArchive
from DwarfGalDataArchive import DwarfGalDataArchive
import AstronomicalParameterArchive as apa
import cantrips as cant
import SphericalVelocityProbabilityProfile as svpp
import SashasAstronomyTools as atools
import scipy.integrate as integrate
import GalaxyMask as gm
import matplotlib.pyplot as plt
#from generatePotentialInterplator import generatePotentialInterpolator




#These are x_sky, z_sky as angular coordinates ALREADY corrected offset from the UNMOVED dSph center
def getStarCartesianPositions(sky_thetas, sky_phis, dist_los, halo_center, gal_dist):
    center_ra_ang, center_offset, center_dec_ang = halo_center
    deg_to_rad = np.pi / 180.0
    center_phi = cant.goodArctan(center_ra_ang, center_dec_ang)
    center_theta = np.sqrt(center_ra_ang ** 2.0 + center_dec_ang ** 2.0) * deg_to_rad
    #sky_x, sky_y, sky_z = [np.sin(sky_phi) * np.sin(sky_theta) * dist_los, np.cos(sky_theta) * dist_los, np.cos(sky_phi) * np.sin(sky_theta) * dist_los]
    sky_xs, sky_ys, sky_zs = np.array([(np.sin(sky_phis) * np.sin(sky_thetas) * dist_los).tolist(), (np.cos(sky_thetas) * dist_los).tolist(), (np.cos(sky_phis) * np.sin(sky_thetas) * dist_los).tolist()])
    center_x, center_y, center_z = [np.sin(center_phi) * np.sin(center_theta) * (gal_dist + center_offset), np.cos(center_theta) * (gal_dist + center_offset), np.cos(center_phi) * np.sin(center_theta) * (gal_dist + center_offset)]
    #total_ang_offset = c.measureAngularSeparationOnSky([cor_x_sky, cor_z_sky], [0.0, 0.0], return_radian = 1)
    #y_dist = dist_los  * np.cos(total_ang_offset)
    #y = y_dist - halo_center[1] - gal_dist
    #x_sqr_p_z_sqr = dist_los ** 2.0 - y_dist ** 2.0
    #offset_in_sky_coords = np.array([sky_xs, sky_ys, sky_zs]) - np.array([center_x, center_y, center_z])
    offsets_in_sky_coords = [sky_xs - center_x, sky_ys - center_y, sky_zs - center_z]
    offsets_in_gal_coords = offsets_in_sky_coords[:] #For a spherical galaxy, no rotations I need worry about
    #if np.sqrt(cor_x_sky ** 2.0 + cor_z_sky ** 2.0) > 0.0:
    #    cos_angle_on_sky = cor_x_sky / np.sqrt(cor_x_sky ** 2.0 + cor_z_sky ** 2.0) # north off east
    #    sin_angle_on_sky = cor_z_sky / np.sqrt(cor_x_sky ** 2.0 + cor_z_sky ** 2.0)
    #    x = np.sqrt(x_sqr_p_z_sqr) * cos_angle_on_sky
    #    z = np.sqrt(x_sqr_p_z_sqr) * sin_angle_on_sky
    #else:
    #    x, z = [0.0, 0.0]
    return offsets_in_gal_coords

class SurfaceBrightnessProfile:

#    #Sum the log of the value of the dSph surface brightness profile at every observed star position.
#    #Must have corrRa and corrDec given in degrees from center of galaxy.
#    #proj_x and proj_z should be given in degrees.
#    def sumLogSurfaceBrightness(self, proj_x, proj_z):
#        #startSum = time.time()
#
#        deg_to_rad = self.astro_arch.getDegToRad()
#        #Because we assume the system is spherically symmetric, I only calculate the interpolator for the upper left quadrant of the sky.
#        #Then I can compute the full probability by just taking the absolute value of the star positions - should cut down on computation
#        #  time because I only need to compute 1/4th the sky projection.
#        proj_x_in_scale_radii = np.abs(proj_x)*deg_to_rad*self.dist/self.rs
#        proj_y_in_scale_radii = np.abs(proj_z)*deg_to_rad*self.dist/self.rs
#
#        print ('[np.min(proj_x_in_scale_radii), np.max(proj_x_in_scale_radii)] = ' + str([np.min(proj_x_in_scale_radii), np.max(proj_x_in_scale_radii)]))
#        print ('[np.min(proj_y_in_scale_radii), np.max(proj_y_in_scale_radii)] = ' + str([np.min(proj_x_in_scale_radii), np.max(proj_y_in_scale_radii)]))
#
#        log_likelihoods = np.log(self.onSkyInterpolator(np.dstack((proj_x_in_scale_radii,proj_y_in_scale_radii))[0]))
#        #endSum = time.time()
#        #print 'Took ' + str(endSum - startSum) + 's to do log sum.'
#        return np.sum(log_likelihoods)

    #Sum the value of the dSph surface brightness profile at every observed star position
    #Must have corrRa and corrDec given in degrees from center of galaxy
    #Note that this quantity CANNOT serve as a good proxy for the likelihood of an array, since we're adding together values where a true probability would demand multiplication.
#    def sumSurfaceBrightness(self,corrRa,corrDec):
#        #print np.dstack(np.meshgrid(corrRa,corrDec))
#        #print np.shape(np.dstack(np.meshgrid(corrRa,corrDec)))
#        #print np.meshgrid(corrRa,corrDec)
#        deg_to_rad = self.astro_arch.getDegToRad()
#        corr_Ra_in_scale_radii = corrRa*deg_to_rad*self.dist/self.rs
#        corr_Dec_in_scale_radii = corrDec*deg_to_rad*self.dist/self.rs
#        return np.sum(self.onSkyInterpolator(np.dstack((corr_Ra_in_scale_radii,corr_Dec_in_scale_radii))[0]))

#    #Multiply the value of the dSph surface brightness profile at every observed star position
#    #Must have corrRa and corrDec given in degrees from center of galaxy
#    #Note that this quantity CAN serve as a proxy for probability, but this value can quickly grow too small for python to effectively handle well.
#    def multiplySurfaceBrightness(self,corrRA,corrDec):
#        deg_to_rad = self.astro_arch.getDegToRad()
#        corr_Ra_in_scale_radii = corrRa*deg_to_rad*self.dist/self.rs
#        corr_Dec_in_scale_radii = corrDec*deg_to_rad*self.dist/self.rs
#        return np.prod(self.onSkyInterpolator(np.dstack((corr_Ra_in_scale_radii,corr_Dec_in_scale_radii))[0]))

    #Return the value of the dSph surface brightness profile at every observed star position
    #Must have corrRa and corrDec given in degrees from center of galaxy
    def getStarProbabilityValues(self, corrRa, corrDec, v_losses):
        deg_to_rad = self.astro_arch.getDegToRad()
        sky_coords = [[np.sqrt(corrRa[i] ** 2.0 + corrDec[i] ** 2.0) * deg_to_rad, cant.goodArctan(corrRa[i], corrDec[i])]  for i in range(len(corrRa)) ]
        #sky_coords = [[np.sqrt(corrRa[i] ** 2.0 + corrDec[i] ** 2.0) * deg_to_rad, (np.arccos(corrRa[i] / np.sqrt(corrRa[i] ** 2.0 + corrDec[i] ** 2.0)) + np.pi * (corrDec[i] < 0.0)) * deg_to_rad]
        #               for i in range(len(corrRa)) ] # Convert from RA and Dec separation to traditional theta-phi offsets, and puts them in radians

        probabilities = self.calculateProbOnPhaseSpaceGrid(sky_coords, v_losses, solid_angles = [0.1], v_los_bin_widths = [1.0])
        log_probs = [np.log10(probability) for probability in probabilities]
        return np.sum(log_probs)

    def getProbAtSeveralPositions(self, sky_coords, dist_loses, v_loses, compute_vel_prob = 1, compute_pos_prob = 1):
            start = time.time()
            theta_skys, phi_skys = [[sky_coord[0] for sky_coord in sky_coords], [sky_coord[1] for sky_coord in sky_coords]]
            individual_prob_funct = lambda target_index: 1.0 / self.normalization_constant * solid_angles[target_index] * v_los_bin_widths[target_index] * np.sum([los_dist_integrand_elems[i] * self.getProbAtSinglePosition(sky_coords[target_index], los_dist_centers[i], v_los_to_calc[target_index], compute_pos_prob = compute_pos_prob, compute_vel_prob = compute_vel_prob) for i in range(n_los_bins)])

            #x_halo_offsets, y_halo_offsets, z_halo_offsets = [getStarCartesianPositions(theta_skys, phi_skys, dist_los, self.halo_center, self.dist, self.halo_sym_axis) for dist_los in dist_loses] # star position in galaxy halo coordinates
            halo_offsets = np.array([getStarCartesianPositions(theta_skys, phi_skys, dist_los, self.halo_center, self.dist) for dist_los in dist_loses])
            #halo_offsets = np.transpose(halo_offsets)
            x_offsets = np.array(cant.flattenListOfLists([halo_offset_at_los[0] for halo_offset_at_los in halo_offsets]))
            y_offsets = np.array(cant.flattenListOfLists([halo_offset_at_los[1] for halo_offset_at_los in halo_offsets]))
            z_offsets = np.array(cant.flattenListOfLists([halo_offset_at_los[2] for halo_offset_at_los in halo_offsets]))
            r_s = np.sqrt(x_offsets ** 2.0 + y_offsets ** 2.0 + z_offsets ** 2.0)
            #print ('r_s = ' + str(r_s))
            if compute_vel_prob and self.include_kinem_prob:
                #print ('Here 1 B')
                #start = time.time()
                vel_probs = cant.flattenListOfLists([[self.velProbCalc.computeProbabilityOfLOSVel(v_loses[v_los_index], x_offsets[i * len(v_loses) + v_los_index], y_offsets[i * len(v_loses) + v_los_index], z_offsets[i * len(v_loses) + v_los_index], self.observer_offset_in_gal_coords ) for v_los_index in range(len(v_loses)) ] for i in range(len(dist_loses))])
                #end = time.time()
                #print ('took ' + str(end - start) + 's for velocity probability at single star position. ')
            else:
                #print ('Here 2 B')
                vel_probs = cant.flattenListOfLists([[1.0 for i in range(len(dist_loses))] for v_los_index in range(len(v_loses))])
            #print ('np.sqrt(x_offset ** 2.0 + y_offset ** 2.0 + z_offset ** 2.0) / self.rs = ' + str(np.sqrt(x_offset ** 2.0 + y_offset ** 2.0 + z_offset ** 2.0)  / self.rs))
            if compute_pos_prob and self.include_morph_prob:
                #print ('Here 3 B')
                #start = time.time()
                pos_probs = self.massDensity.get_density_profile(r_s / self.rs)
                #end = time.time( )
                #print ('took ' + str(end - start) + 's for morphological probability at single star position. ')
            else:
                #print ('Here 4 B')
                pos_probs = 0.0 * r_s + 1.0
            #print ('[vel_probs, pos_probs] = ' + str([vel_probs, pos_probs]))

            #print ('pos_probs = ' + str(pos_probs))
            overall_prob_array = vel_probs * pos_probs
            #print ('[overall_prob_array, dist_loses, sky_coords] = ' + str([overall_prob_array, dist_loses, sky_coords] ))
            reshaped_prob_array = np.reshape(overall_prob_array, (len(dist_loses), len(sky_coords)))
            #print ('reshaped_prob_array = ' + str(reshaped_prob_array))
            #print ('[sky_coord, dist_los, x_halo_offset, y_halo_offset, z_halo_offset, R_halo / self.rs, z_halo / self.rs, v_los, vel_prob, pos_prob] = ' + str([sky_coord, dist_los, x_halo_offset, y_halo_offset, z_halo_offset, R_halo / self.rs, z_halo / self.rs, v_los, vel_prob, pos_prob] ))
            #if np.isnan(overall_prob): print ('[sky_coord, dist_los, x_halo_offset, y_halo_offset, z_halo_offset, R_halo / self.rs, z_halo / self.rs, v_los, vel_prob, pos_prob] = ' + str([sky_coord, dist_los, x_halo_offset, y_halo_offset, z_halo_offset, R_halo / self.rs, z_halo / self.rs, v_los, vel_prob, pos_prob]))
            #if np.isinf(overall_prob): print ('[sky_coord, dist_los, x_halo_offset, y_halo_offset, z_halo_offset, R_halo / self.rs, z_halo / self.rs, v_los, vel_prob, pos_prob] = ' + str([sky_coord, dist_los, x_halo_offset, y_halo_offset, z_halo_offset, R_halo / self.rs, z_halo / self.rs, v_los, vel_prob, pos_prob]))
            #print ('reshaped_prob_array = ' + str(reshaped_prob_array))
            return reshaped_prob_array

    def getProbAtSinglePosition(self, sky_coord, dist_los, v_los, compute_vel_prob = 1, compute_pos_prob = 1):
        #start = time.time()
        theta_sky, phi_sky = sky_coord
        x_offset, y_offset, z_offset = getStarCartesianPositions([theta_sky], [phi_sky], dist_los, self.halo_center, self.dist)
        x_offset, y_offset, z_offset = [float(x_offset[0]), float(y_offset[0]), float(z_offset[0])]
        #print ('[x_offset, y_offset, z_offset] = ' + str([x_offset, y_offset, z_offset]))
        r = np.sqrt(x_offset ** 2.0 + y_offset ** 2.0 + z_offset ** 2.0)
        #print ('r = ' + str(r))
        if compute_vel_prob and self.include_kinem_prob:
            #print ('Here 1 A')
            #start = time.time()
            vel_prob = self.velProbCalc.computeProbabilityOfLOSVel(v_los, x_offset, y_offset, z_offset, self.observer_offset_in_gal_coords )
            #end = time.time()
            #print ('took ' + str(end - start) + 's for velocity probability at single star position. ')
        #    print ('[Here1, v_los, vel_prob] = ' + str([v_los, vel_prob] ))
        else:
            #print ('Here 2 A')
            vel_prob = 1.0
        #    print ('[Here2, v_los, vel_prob] = ' + str([v_los, vel_prob] ))
        #print ('np.sqrt(x_offset ** 2.0 + y_offset ** 2.0 + z_offset ** 2.0) / self.rs = ' + str(np.sqrt(x_offset ** 2.0 + y_offset ** 2.0 + z_offset ** 2.0)  / self.rs))
        if compute_pos_prob and self.include_morph_prob:
            #print ('Here 3 A')
            #start = time.time()
            pos_prob = self.massDensity.get_density_profile(r / self.rs)
            #end = time.time()
            #print ('took ' + str(end - start) + 's for morphological probability at single star position. ')
        else:
            #print ('Here 4 A')
            pos_prob = 1.0
        #print ('pos_prob = ' + str(pos_prob))
        overall_prob = vel_prob * pos_prob
        #print ('[sky_coord, dist_los, x_offset, y_offset, z_offset, r / self.rs, v_los, vel_prob, pos_prob] = ' + str([sky_coord, dist_los, x_offset, y_offset, z_offset, r / self.rs, v_los, vel_prob, pos_prob]))
        #if np.isnan(overall_prob): print ('[sky_coord, dist_los, v_los, x_offset, y_offset, z_offset, r / self.rs, vel_prob, pos_prob, overall_prob] = ' + str([sky_coord, dist_los, v_los, x_offset, y_offset, z_offset, r / self.rs, vel_prob, pos_prob, overall_prob]))
        #print ('overall_prob = ' + str(overall_prob))
        return overall_prob

    def getNormalizationConstant(self):
        sky_coords = self.sky_pixel_coords
        #v_losses = self.v_los_to_calc
        v_losses = [0.0]

        sky_coords_indeces_mesh, star_v_los_indeces_mesh = np.meshgrid(range(len(sky_coords)), range(len(v_losses)))
        flat_sky_coords_indeces = cant.flattenListOfLists(sky_coords_indeces_mesh)
        flat_v_los_indeces = cant.flattenListOfLists(star_v_los_indeces_mesh)
        sky_coords_vec = [sky_coords[sky_coord_index] for sky_coord_index in flat_sky_coords_indeces]
        v_losses_vec = [v_losses[v_los_index] for v_los_index in flat_v_los_indeces]
        solid_angles_vec = [self.sky_pixel_solid_angles[sky_coord_index] for sky_coord_index in flat_sky_coords_indeces]
        #v_los_bin_widths_vec = [self.v_los_bin_widths[v_los_index] for v_los_index in flat_v_los_indeces]
        v_los_bin_widths_vec = [1.0]

        unnormalized_surface_prob_array = self.calculateProbOnPhaseSpaceGrid(sky_coords_vec, v_losses_vec, solid_angles = solid_angles_vec, v_los_bin_widths = v_los_bin_widths_vec, compute_pos_prob = 1, compute_vel_prob = 0)
        #unnormalized_surface_prob_array = [prob * 2.0 for prob in unnormalized_surface_prob_array] #We multiply by 2 becuase we use only positive los velocities to normalize, but a star could actually have a positive or negative los velocity
        unnormalized_surface_prob_array = np.reshape(unnormalized_surface_prob_array, [len(v_losses), len(sky_coords)])
        masked_unnormalized_surface_prob_arrays = self.mask_vals * unnormalized_surface_prob_array
        #This cube tells me the unnormalized probability of observing a star at a position with a velocity.
        #Now I need to normalize the entire cube to unity, and then the cube tells me the normalization parameter.
        #I don't even need to interpolate over this - I can just repeat this analysis on the individual stars,
        # with the calculated normalization constant
        prob_normalization_constant = np.sum(masked_unnormalized_surface_prob_arrays)
        #print ('New normalization constant is: ' + str(prob_normalization_constant ))
        #f, axarr = plt.subplots(2,2)
        #axarr[0,0].scatter(self.sky_pixel_xs, self.sky_pixel_zs, c = self.mask_vals, marker = '.')
        #axarr[1,0].scatter(self.sky_pixel_xs, self.sky_pixel_zs, c = unnormalized_surface_prob_array[0], marker = '.')
        #axarr[0,1].scatter(self.sky_pixel_xs, self.sky_pixel_zs, c = masked_unnormalized_surface_prob_arrays[0], marker = '.')
        #plt.draw()
        #plt.pause(2.0)
        #plt.close('all')
        return prob_normalization_constant

    #We need to measure the total probability of observing a star at a given position
    # with a given observed line of sight velocity.  This is achieved by multiplying
    # the probability of each star.  That is, the probability that I would observe a
    # star at a location on the sky with a given line of sight velocity.
    #We make this measurement with sky coordinates, sky_coords.
    #The sky_coords should be the on-sky positions of the stars, expressed
    #  as spherical theta, phi with respect to the reported dSph center.
    #This is different from the convention in which the degrees are
    # reported as xs and ys on the sky.  I will need to redo my
    # calculation of star positions.
    #sky_coords should be given in units of radians v_los should be given in units km/s
    def calculateProbOnPhaseSpaceGrid(self, sky_coords, v_los_to_calc, solid_angles = [0.1], v_los_bin_widths = [1.0], compute_pos_prob = 1, compute_vel_prob = 1):
        #sky_coords = self.sky_pixel_coords
        #solid_angles = self.sky_pixel_solid_angles
        if len(solid_angles) == 1:
            solid_angles = [solid_angles[0] for coord in sky_coords]
        if len(v_los_bin_widths) == 1:
            v_los_bin_widths = [v_los_bin_widths[0] for v_los in v_los_to_calc]

        los_dist_edges = self.los_bin_edges + self.dist
        n_los_bins = len(los_dist_edges) - 1
        los_dist_centers = [(los_dist_edges[i] + los_dist_edges[i-1]) / 2.0 for i in range(1, n_los_bins+1 )]
        los_dist_integrand_elems = np.array([(los_dist_edges[i+1] - los_dist_edges[i]) * los_dist_centers[i] ** 2.0 for i in range(n_los_bins)])
        individual_prob_funct = lambda sky_coord, coord_solid_angle, v_los, v_los_width: 1.0 / self.normalization_constant * coord_solid_angle * v_los_width * np.sum([(los_dist_edges[i+1] - los_dist_edges[i]) * los_dist_centers[i] ** 2.0 * self.getProbAtSinglePosition(sky_coord, los_dist_centers[i], v_los, compute_pos_prob = compute_pos_prob, compute_vel_prob = compute_vel_prob) for i in range(n_los_bins)])

        individual_prob_funct = lambda target_index: 1.0 / self.normalization_constant * solid_angles[target_index] * v_los_bin_widths[target_index] * np.sum([los_dist_integrand_elems[i] * self.getProbAtSinglePosition(sky_coords[target_index], los_dist_centers[i], v_los_to_calc[target_index], compute_pos_prob = compute_pos_prob, compute_vel_prob = compute_vel_prob) for i in range(n_los_bins)])
        #v_prob_funct = np.vectorize(individual_prob_funct)
        #prob_array_vecA = v_prob_funct(range(len(sky_coords)))
        start = time.time()
        #prob_array_vecA = [individual_prob_funct(index) for index in range(len(sky_coords))]
        raw_prob_array = self.getProbAtSeveralPositions(sky_coords, los_dist_centers, v_los_to_calc, compute_pos_prob = compute_pos_prob, compute_vel_prob = compute_vel_prob)
        prob_array_vec = 1.0 / self.normalization_constant * np.array(solid_angles) * np.array(v_los_bin_widths) * los_dist_integrand_elems.dot(raw_prob_array)
        #print ('[sky_coords, los_dist_centers, v_los_to_calc] = ' + str([sky_coords, los_dist_centers, v_los_to_calc]))
        #print('[[sky_coords[target_index], los_dist_centers[i], v_los_to_calc[target_index] for i in range(n_los_bins)] for target_index in range(len(sky_coords))] = ' + str([[[sky_coords[target_index], los_dist_centers[i], v_los_to_calc[target_index]] for i in range(n_los_bins)] for target_index in range(len(sky_coords))]))
        #print ('prob_array_vecA - prob_array_vecB = ' + str(prob_array_vecA - prob_array_vecB))
        #print ('[np.sum(prob_array_vecA ), np.sum(prob_array_vecB), np.sum(prob_array_vecA - prob_array_vecB )] = ' + str([np.sum(prob_array_vecA ), np.sum(prob_array_vecB), np.sum(prob_array_vecA - prob_array_vecB )]))
        #print ('[sky_coords[0:10], v_los_to_calc[0:10], prob_array_vec[0:10], solid_angles[0:10], self.normalization_constant] = ' + str([sky_coords[0:10], v_los_to_calc[0:10], prob_array_vec[0:10], solid_angles[0:10], self.normalization_constant]))
        end = time.time()
        #print ('Took ' + str(end - start) + 's to compute prob on these phase space parameters ')

        #print (breaknow)
        return prob_array_vec

    #A surface brightness profile is defined by an underlying massDensity, along with some angles that allow one to integrate the density along the los.
    #We then create an interpolator over the grid generated after integration along los
    # (note this is different than the graviational potential interpolator used in computing the mass density profile).
    #Here, x, z, and los_bins are all given in number of scale radii
    #Here, x, z are given in units of scale radii
    #los_vel_bins are in km/s
    #Can I integrate the product along the line of sight using integrate.quad and still have it take a reasonable amount of time?
    def __init__(self, parameter_storer, sharedQuantityHolder, halo_type, disk_file=None, halo_file=None, C=1.0, masking_fields = [], include_morph_prob = 1, include_kinem_prob = 1):

        dwarf_archive = DwarfGalDataArchive()
        self.astro_arch = apa.AstronomicalParameterArchive()
        deg_to_rad = self.astro_arch.getDegToRad()

        self.include_morph_prob = include_morph_prob
        self.include_kinem_prob = include_kinem_prob

        #self.sky_pixel_coords, self.sky_pixel_solid_angles, self.sky_pixel_xs, self.sky_pixel_zs, self.v_los_to_calc, self.los_bin_edges, self.observation_mask = [sharedQuantityHolder.sky_pixel_coords, sharedQuantityHolder.sky_pixel_solid_angles, sharedQuantityHolder.sky_pixel_xs, sharedQuantityHolder.sky_pixel_zs, sharedQuantityHolder.v_los_to_calc, sharedQuantityHolder.los_bin_edges, sharedQuantityHolder.mask]
        self.sky_pixel_coords, self.sky_pixel_solid_angles, self.sky_pixel_xs, self.sky_pixel_zs, self.los_bin_edges, self.observation_mask = [sharedQuantityHolder.sky_pixel_coords, sharedQuantityHolder.sky_pixel_solid_angles, sharedQuantityHolder.sky_pixel_xs, sharedQuantityHolder.sky_pixel_zs, sharedQuantityHolder.los_bin_edges, sharedQuantityHolder.mask]
        #v_los_to_calc = self.v_los_to_calc
        #v_los_bin_edges = [(v_los_to_calc[i] + v_los_to_calc[i-1]) / 2.0 for i in range(1, len(v_los_to_calc)) ]
        #v_los_bin_edges = [v_los_to_calc[0] - (v_los_bin_edges[0] - v_los_to_calc[0])] + v_los_bin_edges + [v_los_to_calc[-1] + (v_los_to_calc[-1] - v_los_bin_edges[-1])]
        #self.v_los_bin_widths = [v_los_bin_edges[i] - v_los_bin_edges[i-1] for i in range(1, len(v_los_bin_edges))]

        self.n_sky_pixels = len(self.sky_pixel_coords)
        self.mask_vals = self.observation_mask.final_mask_interp(np.array([[self.sky_pixel_xs[i], self.sky_pixel_zs[i]] for i in range(self.n_sky_pixels) ]) / deg_to_rad)

        #parameter_storer.printContents()

        #self.los_bin_edges = np.array(self.los_bin_edges)
        self.gamma = self.astro_arch.getGamma()
        self.C = C
        self.dist = parameter_storer.dist
        self.rs = parameter_storer.rs
        self.sigsqr_rr_0 = parameter_storer.sigsqr_rr_0
        self.sigsqr_rr_inf_to_0_rat = parameter_storer.sigsqr_rr_inf_to_0_rat
        self.sigsqr_rr_inf = self.sigsqr_rr_0 * self.sigsqr_rr_inf_to_0_rat
        self.r_sigsqr_rr0 = parameter_storer.r_sigsqr_rr0
        self.alpha_sigsqr_rr = parameter_storer.alpha_sigsqr_rr
        self.gamma_for_beta_inf = parameter_storer.gamma_for_beta_inf
        self.r_beta0 = parameter_storer.r_beta0
        self.alpha_beta = parameter_storer.alpha_beta
        self.rs = parameter_storer.rs
        self.M = parameter_storer.M
        self.c = parameter_storer.c
        self.halo_type = halo_type
        self.halo_center = parameter_storer.halo_center
        self.dispersion_rr_params = [self.sigsqr_rr_0, self.sigsqr_rr_inf, self.r_sigsqr_rr0, self.alpha_sigsqr_rr]
        self.dispersion_beta_params = [self.gamma_for_beta_inf, self.r_beta0, self.alpha_beta]

        self.observer_offset_in_gal_coords = getStarCartesianPositions([0.0], [0.0], 0.0, self.halo_center, self.dist)
        self.observer_offset_in_gal_coords = [float(self.observer_offset_in_gal_coords[0]), float(self.observer_offset_in_gal_coords[1]), float(self.observer_offset_in_gal_coords[2])]

        self.massDensity = SMassDensity(self.M, self.rs, self.dist, self.c, self.C, self.halo_type ,
                                        dispersion_rr_params = self.dispersion_rr_params, dispersion_beta_params = self.dispersion_beta_params,
                                        halo_center = self.halo_center)
        self.velProbCalc = svpp.VelocityProbabilityCalculator(parameter_storer)

        self.normalization_constant = 1.0 #We need a stand in so we can calculate the normalization constant
        self.normalization_constant = self.getNormalizationConstant()
