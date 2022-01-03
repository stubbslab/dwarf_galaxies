import SashasAstronomyTools as atools
import GalaxyMask as gm
import numpy as np
#Holds the objects that loops through the MCMC share
class SharedObjects:

    def computeSkyPixels(self):
        #print('[self.approx_n_sky_pixels, self.maximum_sky_angle_rad, self.sky_pix_n_start_angles] = ' + str([self.approx_n_sky_pixels, self.maximum_sky_angle_rad, self.sky_pix_n_start_angles]))
        self.sky_pixel_coords, self.sky_pixel_solid_angles = atools.computeSkyPixels(self.approx_n_sky_pixels, self.maximum_sky_angle_rad, n_start_angles = self.sky_pix_n_start_angles)
        self.sky_pixel_xs, self.sky_pixel_zs = [np.array([coord[0] * np.cos(coord[1]) for coord in self.sky_pixel_coords]), np.array([coord[0] * np.sin(coord[1]) for coord in self.sky_pixel_coords])]
        return 1

    def generateGalaxyMask(self, compute_mask):
        Rmesh_degrees, zmesh_degrees = np.meshgrid(self.R_mask_deg, self.z_mask_deg)
        if compute_mask:
            mask = gm.GalaxyMask(Rmesh_degrees, zmesh_degrees, self.gal, mask_types = self.mask_types)
        else:
            mask = np.zeros(np.shape(Rmesh_degrees)) + 1.0
        return mask

    def __init__(self, dist, R_mask_deg, z_mask_deg, los_bin_edges_as_fracs, approx_n_sky_pixels, maximum_sky_angle_rad, #v_los_to_calc,
                       sky_pix_n_start_angles = 16, mask_types = ['n_vel_meas'], gal = 'fornax',compute_mask = 1 ):
        self.dist = dist
        #self.storer = storer
        self.R_mask_deg = R_mask_deg
        self.z_mask_deg = z_mask_deg
        self.mask_types = mask_types
        self.gal = gal
        #self.v_los_to_calc  = v_los_to_calc
        self.mask = self.generateGalaxyMask(compute_mask)
        dist_of_solid_angle = np.tan(maximum_sky_angle_rad) * self.dist
        self.los_bin_edges = dist_of_solid_angle * np.array(los_bin_edges_as_fracs)
        self.approx_n_sky_pixels = approx_n_sky_pixels
        self.maximum_sky_angle_rad = maximum_sky_angle_rad
        self.sky_pix_n_start_angles = sky_pix_n_start_angles
        sky_pixels_computed = self.computeSkyPixels()
