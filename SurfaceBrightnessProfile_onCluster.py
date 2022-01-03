#Define the SurfaceBrightnessProfile class.
#This takes in a MassDensity class, integrates it along a los, and then creates an interpolator for the observed surface brightness.

from MassDensity import MassDensity
import numpy as np
import math
import time
from scipy.interpolate import RegularGridInterpolator
import DwarfGalaxyParametersStorer as dgps
from DwarfGalDataArchive import DwarfGalDataArchive
import AstronomicalParameterArchive as apa
from cantrips import getAreaMesh
import SharedObjectHolder as soh
import cantrips as cant
import VelocityProbabilityProfile as vpp
import cantrips as cant
import matplotlib.pyplot as plt
#from generatePotentialInterplator import generatePotentialInterpolator

#These are x_sky, z_sky as angular coordinates ALREADY corrected offset from the UNMOVED dSph center
def getStarCartesianPositions(sky_thetas, sky_phis, dist_los, center, gal_dist, sym_axis):
    center_ra_ang, center_offset, center_dec_ang = center
    deg_to_rad = np.pi / 180.0
    center_phi = cant.goodArctan(center_ra_ang, center_dec_ang)
    center_theta = np.sqrt(center_ra_ang ** 2.0 + center_dec_ang ** 2.0) * deg_to_rad
    #sky_x, sky_y, sky_z = np.array([np.sin(sky_phis) * np.sin(sky_thetas) * dist_los, np.cos(sky_thetas) * dist_los, np.cos(sky_phis) * np.sin(sky_thetas) * dist_los])
    sky_xs, sky_ys, sky_zs = np.array([(np.sin(sky_phis) * np.sin(sky_thetas) * dist_los).tolist(), (np.cos(sky_thetas) * dist_los).tolist(), (np.cos(sky_phis) * np.sin(sky_thetas) * dist_los).tolist()])
    center_x, center_y, center_z = np.array([np.sin(center_phi) * np.sin(center_theta) * (gal_dist + center_offset), np.cos(center_theta) * (gal_dist + center_offset), np.cos(center_phi) * np.sin(center_theta) * (gal_dist + center_offset)])
    #total_ang_offset = c.measureAngularSeparationOnSky([cor_x_sky, cor_z_sky], [0.0, 0.0], return_radian = 1)
    #y_dist = dist_los  * np.cos(total_ang_offset)
    #y = y_dist - halo_center[1] - gal_dist
    #x_sqr_p_z_sqr = dist_los ** 2.0 - y_dist ** 2.0

    offsets_in_sky_coords = [sky_xs - center_x, sky_ys - center_y, sky_zs - center_z]
    gal_theta, gal_phi = [np.arccos(sym_axis[2]), cant.goodArctan(sym_axis[0], sym_axis[1])]
    offsets_in_gal_coords = [(offsets_in_sky_coords[0] * np.cos(gal_phi) + offsets_in_sky_coords[1] * np.sin(gal_phi)), (-offsets_in_sky_coords[0] * np.sin(gal_phi) + offsets_in_sky_coords[1] * np.cos(gal_phi)), offsets_in_sky_coords[2]]
    offsets_in_gal_coords = [offsets_in_gal_coords[0] * np.cos(gal_theta) - np.sin(gal_theta) * offsets_in_gal_coords[2], offsets_in_gal_coords[1], np.sin(gal_theta) * offsets_in_gal_coords[0] + np.cos(gal_theta) * offsets_in_gal_coords[2]]
    #For a spherical galaxy, no rotations I need worry about
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
#    def sumLogSurfaceBrightness(self, proj_x, proj_y):
#        #startSum = time.time()
#        astro_archive = AstronomicalParameterArchive()
#        deg_to_rad = astro_archive.getDegToRad()
#        proj_x_in_scale_radii = proj_x*deg_to_rad*self.dist/self.rs
#        proj_y_in_scale_radii = proj_y*deg_to_rad*self.dist/self.rs
#
#        log_likelihoods = np.log(self.onSkyInterpolator(np.dstack((proj_x_in_scale_radii,proj_y_in_scale_radii))[0]))
#        #endSum = time.time()
#        #print 'Took ' + str(endSum - startSum) + 's to do log sum.'
#        return np.sum(log_likelihoods)

#    #Sum the value of the dSph surface brightness profile at every observed star position
#    #Must have corrRa and corrDec given in degrees from center of galaxy
#    #Note that this quantity CANNOT serve as a good proxy for the likelihood of an array, since we're adding together values where a true probability would demand multiplication.
#    def sumSurfaceBrightness(self,corrRa,corrDec):
#        #print np.dstack(np.meshgrid(corrRa,corrDec))
#        #print np.shape(np.dstack(np.meshgrid(corrRa,corrDec)))
#        #print np.meshgrid(corrRa,corrDec)
#        astro_archive = AstronomicalParameterArchive()
#        deg_to_rad = astro_archive.getDegToRad()
#        corr_Ra_in_scale_radii = corrRa*deg_to_rad*self.dist/self.rs
#        corr_Dec_in_scale_radii = corrDec*deg_to_rad*self.dist/self.rs
#        return np.sum(self.onSkyInterpolator(np.dstack((corr_Ra_in_scale_radii,corr_Dec_in_scale_radii))[0]))

#    #Multiply the value of the dSph surface brightness profile at every observed star position
#    #Must have corrRa and corrDec given in degrees from center of galaxy
#    #Note that this quantity CAN serve as a proxy for probability, but this value can quickly grow too small for python to effectively handle well.
#    def multiplySurfaceBrightness(self,corrRA,corrDec):
#        astro_archive = AstronomicalParameterArchive()
#        deg_to_rad = astro_archive.getDegToRad()
#        corr_Ra_in_scale_radii = corrRa*deg_to_rad*self.dist/self.rs
#        corr_Dec_in_scale_radii = corrDec*deg_to_rad*self.dist/self.rs
#        return np.prod(self.onSkyInterpolator(np.dstack((corr_Ra_in_scale_radii,corr_Dec_in_scale_radii))[0]))

    #Return the value of the dSph surface brightness profile at every observed star position
    #Must have corrRa and corrDec given in degrees from center of galaxy
    def getStarProbabilityValues(self, corrRa, corrDec, v_losses):
        deg_to_rad = self.astro_archive.getDegToRad()
        sky_coords = [[np.sqrt(corrRa[i] ** 2.0 + corrDec[i] ** 2.0) * deg_to_rad, cant.goodArctan(corrRa[i], corrDec[i])]  for i in range(len(corrRa)) ]

        log_probs = np.log10(self.calculateProbOnPhaseSpaceGrid(sky_coords, v_losses, solid_angles = [0.1], v_los_bin_widths = [1.0]))
        #log_probs = [np.log10(probability) for probability in probabilities]
        return np.sum(log_probs)


    def getProbAtSeveralPositions(self, sky_coords, dist_loses, v_loses, compute_vel_prob = 1, compute_pos_prob = 1):
        start = time.time()
        theta_skys, phi_skys = [[sky_coord[0] for sky_coord in sky_coords], [sky_coord[1] for sky_coord in sky_coords]]

        #x_halo_offsets, y_halo_offsets, z_halo_offsets = [getStarCartesianPositions(theta_skys, phi_skys, dist_los, self.halo_center, self.dist, self.halo_sym_axis) for dist_los in dist_loses] # star position in galaxy halo coordinates
        halo_offsets = np.array([getStarCartesianPositions(theta_skys, phi_skys, dist_los, self.halo_center, self.dist, self.halo_sym_axis) for dist_los in dist_loses])
        x_halo_offsets = np.array(cant.flattenListOfLists([halo_offset_at_los[0] for halo_offset_at_los in halo_offsets]))
        y_halo_offsets = np.array(cant.flattenListOfLists([halo_offset_at_los[1] for halo_offset_at_los in halo_offsets]))
        z_halo_offsets = np.array(cant.flattenListOfLists([halo_offset_at_los[2] for halo_offset_at_los in halo_offsets]))
        #x_disk_offsets, y_disk_offsets, z_disk_offsets = [getStarCartesianPositions(theta_skys, phi_skys, dist_los, self.disk_center, self.dist, self.disk_sym_axis) for dist_los in dist_loses]
        disk_offsets = np.array([getStarCartesianPositions(theta_skys, phi_skys, dist_los, self.halo_center, self.dist, self.halo_sym_axis) for dist_los in dist_loses])
        x_disk_offsets = np.array(cant.flattenListOfLists([disk_offset_at_los[0] for disk_offset_at_los in disk_offsets]))
        y_disk_offsets = np.array(cant.flattenListOfLists([disk_offset_at_los[1] for disk_offset_at_los in disk_offsets]))
        z_disk_offsets = np.array(cant.flattenListOfLists([disk_offset_at_los[2] for disk_offset_at_los in disk_offsets]))
        R_halos = np.sqrt(x_halo_offsets ** 2.0 + y_halo_offsets ** 2.0)
        z_halos = z_halo_offsets
        R_disks = np.sqrt(x_halo_offsets ** 2.0 + x_halo_offsets ** 2.0)
        z_disks = z_disk_offsets
        if compute_vel_prob and self.include_kinem_prob:
            #start = time.time()
            vel_probs = cant.flattenListOfLists([[self.velProbCalc.computeProbabilityOfLOSVel(v_loses[v_los_index], x_halo_offsets[i * len(v_loses) + v_los_index], y_halo_offsets[i * len(v_loses) + v_los_index], z_halo_offsets[i * len(v_loses) + v_los_index], self.observer_offset_in_gal_coords) for v_los_index in range(len(v_loses))] for i in range(len(dist_loses))])
            #end = time.time()
            #print ('took ' + str(end - start) + 's for velocity probability at single star position. ')
        else:
            vel_probs = cant.flattenListOfLists([[1.0 for i in range(len(dist_loses))] for v_los_index in range(len(v_loses))])
        #print ('np.sqrt(x_offset ** 2.0 + y_offset ** 2.0 + z_offset ** 2.0) / self.rs = ' + str(np.sqrt(x_offset ** 2.0 + y_offset ** 2.0 + z_offset ** 2.0)  / self.rs))
        mass_density_exponent_scaling = self.gamma * self.M / (self.rs * self.el)
        mass_density_exponent_scaling = 1.0
        if compute_pos_prob and self.include_morph_prob:
            start = time.time()
            pos_probs = self.massDensity.get_density_cylindrical(R_halos / self.rs, z_halos / self.rs, R_disks / self.rs, z_disks / self.rs, mass_density_exponent_scaling)
            end = time.time()
            #print ('took ' + str(end - start) + 's for morphological probability at single star position. ')
            #print ('took ' + str(end - start) + 's for morphological probability at single star position. ')
        else:
            pos_probs = 0.0 * R_halos + 1.0
        #print ('[vel_probs, pos_probs] = ' + str([vel_probs, pos_probs]))

        overall_prob_array = vel_probs * pos_probs

        reshaped_prob_array = np.reshape(overall_prob_array, (len(dist_loses), len(sky_coords))) #TO BE FILLED IN
        #print ('[sky_coord, dist_los, x_halo_offset, y_halo_offset, z_halo_offset, R_halo / self.rs, z_halo / self.rs, v_los, vel_prob, pos_prob] = ' + str([sky_coord, dist_los, x_halo_offset, y_halo_offset, z_halo_offset, R_halo / self.rs, z_halo / self.rs, v_los, vel_prob, pos_prob] ))
        #if np.isnan(overall_prob): print ('[sky_coord, dist_los, x_halo_offset, y_halo_offset, z_halo_offset, R_halo / self.rs, z_halo / self.rs, v_los, vel_prob, pos_prob] = ' + str([sky_coord, dist_los, x_halo_offset, y_halo_offset, z_halo_offset, R_halo / self.rs, z_halo / self.rs, v_los, vel_prob, pos_prob]))
        #if np.isinf(overall_prob): print ('[sky_coord, dist_los, x_halo_offset, y_halo_offset, z_halo_offset, R_halo / self.rs, z_halo / self.rs, v_los, vel_prob, pos_prob] = ' + str([sky_coord, dist_los, x_halo_offset, y_halo_offset, z_halo_offset, R_halo / self.rs, z_halo / self.rs, v_los, vel_prob, pos_prob]))
        #print ('reshaped_prob_array = ' + str(reshaped_prob_array))
        return reshaped_prob_array


    def getProbAtSinglePosition(self, sky_coord, dist_los, v_los, compute_vel_prob = 1, compute_pos_prob = 1):
        start = time.time()
        theta_sky, phi_sky = sky_coord

        x_halo_offset, y_halo_offset, z_halo_offset = getStarCartesianPositions([theta_sky], [phi_sky], dist_los, self.halo_center, self.dist, self.halo_sym_axis)  # star position in galaxy halo coordinates
        x_halo_offset, y_halo_offset, z_halo_offset = [float(x_halo_offset), float(y_halo_offset), float(z_halo_offset)]
        x_disk_offset, y_disk_offset, z_disk_offset = getStarCartesianPositions([theta_sky], [phi_sky], dist_los, self.disk_center, self.dist, self.disk_sym_axis)
        x_disk_offset, y_disk_offset, z_disk_offset = [float(x_disk_offset), float(y_disk_offset), float(z_disk_offset)]
        R_halo = np.sqrt(x_halo_offset ** 2.0 + y_halo_offset ** 2.0)
        z_halo = z_halo_offset
        R_disk = np.sqrt(x_disk_offset ** 2.0 + y_disk_offset ** 2.0)
        z_disk = z_disk_offset
        #print ('[x_offset, y_offset, z_offset] = ' + str([x_offset, y_offset, z_offset] ))
        if compute_vel_prob and self.include_kinem_prob:
            #start = time.time()
            vel_prob = self.velProbCalc.computeProbabilityOfLOSVel(v_los, x_halo_offset, y_halo_offset, z_halo_offset, self.observer_offset_in_gal_coords)
            #end = time.time()
            #print ('took ' + str(end - start) + 's for velocity probability at single star position. ')
        else:
            vel_prob = 1.0
        #print ('np.sqrt(x_offset ** 2.0 + y_offset ** 2.0 + z_offset ** 2.0) / self.rs = ' + str(np.sqrt(x_offset ** 2.0 + y_offset ** 2.0 + z_offset ** 2.0)  / self.rs))
        mass_density_exponent_scaling = self.gamma * self.M / (self.rs * self.el)
        mass_density_exponent_scaling = 1.0
        if compute_pos_prob and self.include_morph_prob:
            start = time.time()
            pos_prob = self.massDensity.get_density_cylindrical(R_halo / self.rs, z_halo / self.rs, R_disk / self.rs, z_disk / self.rs, mass_density_exponent_scaling)
            end = time.time()
            print ('took ' + str(end - start) + 's for morphological probability at single star position. ')
        else:
            pos_prob = 1.0
        overall_prob = vel_prob * pos_prob
        #print ('[sky_coord, dist_los, x_halo_offset, y_halo_offset, z_halo_offset, R_halo / self.rs, z_halo / self.rs, v_los, vel_prob, pos_prob] = ' + str([sky_coord, dist_los, x_halo_offset, y_halo_offset, z_halo_offset, R_halo / self.rs, z_halo / self.rs, v_los, vel_prob, pos_prob] ))
        #if np.isnan(overall_prob): print ('[sky_coord, dist_los, x_halo_offset, y_halo_offset, z_halo_offset, R_halo / self.rs, z_halo / self.rs, v_los, vel_prob, pos_prob] = ' + str([sky_coord, dist_los, x_halo_offset, y_halo_offset, z_halo_offset, R_halo / self.rs, z_halo / self.rs, v_los, vel_prob, pos_prob]))
        #if np.isinf(overall_prob): print ('[sky_coord, dist_los, x_halo_offset, y_halo_offset, z_halo_offset, R_halo / self.rs, z_halo / self.rs, v_los, vel_prob, pos_prob] = ' + str([sky_coord, dist_los, x_halo_offset, y_halo_offset, z_halo_offset, R_halo / self.rs, z_halo / self.rs, v_los, vel_prob, pos_prob]))
        #print ('overall_prob = ' + str(overall_prob))
        return overall_prob

    #Normalized the surface probability, allowing for the possibility of appling a new mask and renormalizing
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

        #flat_v_los_indeces = c.flattenListOfLists(star_v_los_indeces_mesh)
        #sky_coords_vec = [sky_coords[sky_coord_index] for sky_coord_index in flat_sky_coords_indeces]
        #v_losses_vec = [v_losses[v_los_index] for v_los_index in flat_v_los_indeces]
        #solid_angles_vec = [self.sky_pixel_solid_angles[sky_coord_index] for sky_coord_index in flat_sky_coords_indeces]
        #v_los_bin_widths_vec = [self.v_los_bin_widths[v_los_index] for v_los_index in flat_v_los_indeces]

        unnormalized_surface_prob_array = self.calculateProbOnPhaseSpaceGrid(sky_coords_vec, v_losses_vec, solid_angles = solid_angles_vec, v_los_bin_widths = v_los_bin_widths_vec, compute_pos_prob = 1, compute_vel_prob = 0)
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
        individual_prob_funct = lambda sky_coord, coord_solid_angle, v_los, v_los_width: 1.0 / self.normalization_constant * coord_solid_angle * v_los_width * np.sum([(los_dist_edges[i+1] - los_dist_edges[i]) * los_dist_centers[i] ** 2.0 * self.getProbAtSinglePosition(sky_coord, los_dist_centers[i], v_los) for i in range(n_los_bins)])
        individual_prob_funct = lambda target_index: 1.0 / self.normalization_constant * solid_angles[target_index] * v_los_bin_widths[target_index] * np.sum([los_dist_integrand_elems[i] * self.getProbAtSinglePosition(sky_coords[target_index], los_dist_centers[i], v_los_to_calc[target_index]) for i in range(n_los_bins)])
        #v_prob_funct = np.vectorize(individual_prob_funct)
        #prob_array_vec = v_prob_funct(range(len(sky_coords)))
        start = time.time()

        #prob_array_vec = [individual_prob_funct(index) for index in range(len(sky_coords))]
        #print ('[np.shape(np.array(los_dist_integrand_elems)), np.shape(self.getProbAtSeveralPositions(sky_coords, los_dist_centers, v_los_to_calc))]  = ' + str([np.shape(np.array(los_dist_integrand_elems)), np.shape(self.getProbAtSeveralPositions(sky_coords, los_dist_centers, v_los_to_calc))] ))

        raw_prob_array = self.getProbAtSeveralPositions(sky_coords, los_dist_centers, v_los_to_calc, compute_pos_prob = compute_pos_prob, compute_vel_prob = compute_vel_prob)
        #print ('[len(sky_coords), len(v_los_to_calc), len(los_dist_integrand_elems), np.shape(raw_prob_array)] = ' + str([len(sky_coords), len(v_los_to_calc), len(los_dist_integrand_elems), np.shape(raw_prob_array)]))
        prob_array_vec = 1.0 / self.normalization_constant * np.array(solid_angles) * np.array(v_los_bin_widths) * los_dist_integrand_elems.dot(raw_prob_array)
        #print ('[sky_coords[0:10], v_los_to_calc[0:10], prob_array_vec[0:10], solid_angles[0:10], self.normalization_constant] = ' + str([sky_coords[0:10], v_los_to_calc[0:10], prob_array_vec[0:10], solid_angles[0:10], self.normalization_constant]))
        end = time.time()
        #print ('Took ' + str(end - start) + 's')

        return prob_array_vec


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
    #We then create an interpolator over the grid generated after integration along los.
    #We also need to normalize this surface brightness profile.  This is new.
    #Then we can calculate the likelihood at every position where a star is defined.  We actually do not
    #  need to define an interpolator.
    # (note this is different than the graviational potential interpolator used in computing the mass density profile).
    #Here, R, z, and los_bins are all given in number of scale radii
    def __init__(self, parameter_storer, sharedQuantityHolder, halo_type, disk_type, disk_file=None, halo_file=None, C=1.0, disk_interpolating_function=None, halo_interpolating_function=None, masking_fields = [], include_morph_prob = 1, include_kinem_prob = 1):
        #parameter_storer.printContents()
        self.include_morph_prob = include_morph_prob
        self.include_kinem_prob = include_kinem_prob
        dwarf_archive = DwarfGalDataArchive()
        self.astro_archive = apa.AstronomicalParameterArchive()
        deg_to_rad = self.astro_archive.getDegToRad()
        self.gamma = self.astro_archive.getGamma()
        #self.sky_pixel_coords, self.sky_pixel_solid_angles, self.sky_pixel_xs, self.sky_pixel_zs, self.v_los_to_calc, self.los_bin_edges, self.observation_mask = [sharedQuantityHolder.sky_pixel_coords, sharedQuantityHolder.sky_pixel_solid_angles, sharedQuantityHolder.sky_pixel_xs, sharedQuantityHolder.sky_pixel_zs, sharedQuantityHolder.v_los_to_calc, sharedQuantityHolder.los_bin_edges, sharedQuantityHolder.mask]
        self.sky_pixel_coords, self.sky_pixel_solid_angles, self.sky_pixel_xs, self.sky_pixel_zs, self.los_bin_edges, self.observation_mask = [sharedQuantityHolder.sky_pixel_coords, sharedQuantityHolder.sky_pixel_solid_angles, sharedQuantityHolder.sky_pixel_xs, sharedQuantityHolder.sky_pixel_zs, sharedQuantityHolder.los_bin_edges, sharedQuantityHolder.mask]

        self.n_sky_pixels = len(self.sky_pixel_coords)
        self.mask_vals = self.observation_mask.final_mask_interp(np.array([[self.sky_pixel_xs[i], self.sky_pixel_zs[i]] for i in range(self.n_sky_pixels) ]) / deg_to_rad)

        #parameter_storer.printContents()

        self.los_bin_edges = np.array(self.los_bin_edges)
        #self.Rlikeaxis = R
        #self.zlikeaxis = z
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
        self.omega_phi = parameter_storer.omega_phi
        self.sigsqr_RR = parameter_storer.sigsqr_RR
        self.dispersion_b = parameter_storer.dispersion_b
        #self.beta_int_dispersion_funct = parameter_storer.beta_int_dispersion_funct
        self.disk_center = parameter_storer.disk_center
        self.halo_center = parameter_storer.halo_center
        self.halo_type = halo_type
        self.disk_type = disk_type


        self.observer_offset_in_gal_coords = getStarCartesianPositions([0.0], [0.0], 0.0, self.halo_center, self.dist, self.halo_sym_axis)
        self.observer_offset_in_gal_coords = [float(self.observer_offset_in_gal_coords[0]), float(self.observer_offset_in_gal_coords[1]), float(self.observer_offset_in_gal_coords[2])]

        self.massDensity = MassDensity(self.gamma, self.M, self.rs, self.eps, self.c, self.lam, self.zeta, self.el, self.omega_phi, self.C, self.sigsqr_RR,
                                       dispersion_b = self.dispersion_b, dist = self.dist, disk_file=disk_file, halo_file=halo_file,
                                       disk_interpolating_function=disk_interpolating_function, halo_interpolating_function=halo_interpolating_function,
                                       halo_type = halo_type, disk_type = disk_type, halo_center = self.halo_center, disk_center = self.disk_center)

        self.velProbCalc = vpp.VelocityProbabilityCalculator(parameter_storer)
        self.normalization_constant = 1.0
        self.normalization_constant = self.getNormalizationConstant()
        #endDefMass = time.time()
        #print 'Took ' + str(endDefMass - startDefMass) + 's to define massDensity.'
        #print 'In surface prob density:'
        #print 'self.massDensity.integrate_los(self.phi,self.theta,self.a,self.b,self.los_bins) = '
        #print self.massDensity.integrate_los(self.phi,self.theta,self.a,self.b,self.los_bins)

        #surfaceBrightness is the physical surface brightness expected, based on the PHYSICAL mass profile generated.
        #startIntegrate = time.time()
        #self.surfaceBrightness = self.massDensity.integrate_los(self.halo_sym_axis, self.disk_sym_axis, self.los_bins)

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
        #self.projXMesh_radians, self.projYMesh_radians = np.meshgrid(self.Rlikeaxis * self.rs / self.dist, self.zlikeaxis * self.rs / self.dist)
        #projX_scaleRadii = self.Rlikeaxis
        #projY_scaleRadii = self.zlikeaxis
        #self.projXMesh_scaleRadii, self.projYMesh_scaleRadii = np.meshgrid(projX_scaleRadii, projY_scaleRadii)
        #normalization_constant = 0.0

        #surface_probability is the OBSERVED probability of finding a star at a given position, which must account for the incomplete sky sampling region
        #if self.observation_mask is None:
        #    if population is None:
        #        self.observation_mask = np.zeros(np.shape(self.surfaceBrightness)) + 1.0
        #    else:
        #        self.observation_mask = GalaxyMask(self.projXMesh_radians, self.projYMesh_radians, population[0], mask_types = ['n_vel_meas'], masking_fields = masking_fields, units_of_coord_meshes = 'radians')

        #self.normalizeSurfaceProfile(self.observation_mask, renormalizing = 0)
        #print self.normalization_constant
