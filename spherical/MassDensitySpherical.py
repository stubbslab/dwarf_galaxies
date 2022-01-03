import SphericalDMPotentials as sdmp
import AstronomicalParameterArchive as apa
import sphericalVelocityDispersionFuncts as svdf
import numpy as np
import math

class SMassDensity:

    def get_density_profile(self, r_halo):
        #print ('r_halo = ' + str(r_halo))
        km_to_pc = 1.0 / self.astro_archive.getParsecToM() * 10.0 ** 3.0
        term1 = self.sigsqr_rr_funct(0.0) / self.sigsqr_rr_funct(r_halo)
        term2_exp = - (2.0 * (self.beta_over_r_int_funct(r_halo) - self.beta_over_r_int_funct(0.0)))
        term2 = np.exp(term2_exp)
        term3_exp = - (  (self.gamma*self.M)/(self.rs) * self.dphi_over_sigsqr_rr_Nint_funct(r_halo) )
        term3 = np.exp(term3_exp)

        #print ('(self.gamma*self.M)/(self.rs) = ' + str((self.gamma*self.M)/(self.rs)))
        #print ('[term1, term2, term3] = ' + str([term1, term2, term3]))

        return term1 * term2 * term3

        #Turn a 3-dimensional dwarf galaxy distribution into a 2-dimensional distribution by summing a series of cross sections perpendicular to the viewing line of site.
        #The integration_rect_centers tell you the positions along los where you wish to measure a cross section, measured in scale radii from the center of the galaxy.
        #I think I should change this to be the bin borders
        #Each cross section is weighted by the portion of the viewing axis that lies between the bounds of the two bins.

        #|   |   |   |   | | | |   |   |   |
        #  x   x   x   x  x x x  x   x   x
    def integrate_los(self, integration_bin_centers,
                      apply_tidal_constraint = 0, los_bin_step = 1.0, exponent_scaling = 1.0, overflow_exponent_rescaling = -100.0) :
        if apply_tidal_constraint:
            integration_bin_limits = self.getLosTidalLimits(tide_source = 'MW',sampling_density = los_bin_step)
            integration_bin_borders = np.arange(integration_bin_limits[0], integration_bin_limits[1] + los_bin_step / 2.0, los_bin_step)
        else:
            n_centers = len(integration_bin_centers)
            integration_bin_borders = ( [integration_bin_centers[0] - (integration_bin_centers[1] - integration_bin_centers[0]) / 2.0]
                                        + [(integration_bin_centers[i+1] + integration_bin_centers[i]) / 2.0 for i in range(n_centers - 1)]
                                        + [integration_bin_centers[n_centers-1] + (integration_bin_centers[n_centers-1] - integration_bin_centers[n_centers - 2]) /2.0 ] )

        integration_bin_borders.sort()
        integrated_array=np.zeros((len(self.x_values),len(self.z_values)))
        #next_x=integration_rect_centers[1]
        #To have a well defined calculation for first step, have to effectively add an extra point
        #prev_x=integration_rect_centers[0] - (next_x - integration_rect_centers[0])
        #Note that we have one more bin border than the number of bins, so the for loop goes up to the length of the array minus 1
        integrated_array = sum([ self.get_cross_section(integration_bin_centers[i], exponent_scaling = exponent_scaling)
                                  * (integration_bin_borders[i + 1] - integration_bin_borders[i]) for i in range(len(integration_bin_centers))])
        if math.isinf(np.sum(integrated_array)):
            print ('np.sum(integrated_array) is infinite.  Recycling with new overall scaling...' )
            return self.integrate_los(integration_bin_centers,
                                      apply_tidal_constraint = apply_tidal_constraint, los_bin_step = los_bin_step, exponent_scaling = overflow_exponent_rescaling + exponent_scaling)
        else:
            return integrated_array

    #We want to mass density on some plane that intersects our distribution at some arbitrary angle and depth
    #The basic way we do the rotations is to define 6, 2d meshgrids.
    #Each of these corresponds to the value of the natural x, y, and z coordinates of the DM disk and halos on the 2d cross-section of interest.
    def get_cross_section(self, los_depth, exponent_scaling = 1.0 ):

        #x goes right on sky, z goes up on sky, y goes along los AWAY from observer
        y_sky = los_depth # the distance PAST gal center, alogn los
        xmesh_sky ,zmesh_sky = np.meshgrid(self.x_values, self.z_values)

        #The galaxy profile is spherical, and so we need only one mesh.
        rmesh_halo = np.sqrt((xmesh_sky - self.halo_center[0] / self.rs) ** 2.0 + (y_sky - self.halo_center[1] / self.rs) ** 2.0 + (zmesh_sky - self.halo_center[2] / self.rs) ** 2.0)

        #Any issues with the coordinate mesh calculations?
        if (np.any(np.isnan(xmesh_sky))):
            print ('!!!!! xmesh_sky contains np.nan !!!!!')
        if (np.any(np.isnan(y_sky))):
            print ('!!!!! y_sky contains np.nan !!!!!')
        if (np.any(np.isnan(zmesh_sky))):
            print ('!!!!! zmesh_sky contains np.nan !!!!!')
        if (np.any(np.isnan(rmesh_halo))):
            print ('!!!!! rmesh_halo contains np.nan !!!!!')
        if (np.any(0.0 > (np.around(rmesh_halo, 7)))):
            print ('!!!!! rmesh_halo value possibly invalid (negative) !!!!!')
            print ('(np.min( (np.around(rmesh_halo, 7)) )) = ' + str(np.min( (np.around(rmesh_halo, 7)) )))

        proj_array = self.get_density_profile(rmesh_halo)

        #NOTE WE RETURN THE TRANSPOSE SINCE THAT IS THE THING THAT WE CAN FEED INTO AN INTERPOLATOR
        #print 'Applying correction.'
        #return np.transpose(proj_array * proj_array_where_valid + proj_array_correction )
        #return np.transpose(proj_array * proj_array_correction) # if you impose exponential cut off past some limit
        return np.transpose(proj_array) # if you don't want to impose any exponential cut off

    #We define the spherical velocity dispersion profile as: sigsqr_rr(r) = sigsqr_rr_inf - (sigsqr_rr_inf - sigsqr_rr_0)  * 1 / (1 + r/r_sig0 ** alpha_rr), alpha_rr > 0.0
    # dispersion_rr_params = [sigsqr_rr_0, sigsqr_rr_inf, r0, alpha_rr]
    # NOTE: sigsqr_rr_0 and sigsqr_rr_inf are unitful quantities and should be given in m/s.  A typical value for Fornax is ~100 km/s
    #We define the velocity beta dispersion profile as: beta(r) = beta_inf - beta_inf * 1 / (1 + r/r_beta0 ** alpha_beta) = beta_inf * (r/r_beta0 ** alpha_beta) / (1 + r/r_beta0 ** alpha_beta)
    # dispersion_beta_params = [gamma_inf, r_beta0, alpha_beta]


    def __init__(self, M, rs, dist, c, C, halo_type,
                 dispersion_rr_params = [100.0, 100.0, 1.0, 2.0], dispersion_beta_params = [0.0, 1.0, 1.0], halo_center = [0.0, 0.0, 0.0],
                 n_interp_rs_for_velocity_dispersion_int = 101):
         self.astro_archive = apa.AstronomicalParameterArchive()
         gamma = self.astro_archive.getGamma()
         scaled_spherical_potential_funct = sdmp.getPotentialFunction(halo_type, c)
         self.gamma = gamma
         self.M = M
         self.rs = rs
         self.dist = dist
         self.spherical_potential_funct = lambda rs: (self.gamma*self.M)/(self.rs) * scaled_spherical_potential_funct(rs)

         #self.x_values = xs
         #self.z_values = zs
         #self.x_range=[np.min(xs),np.max(xs)]
         #self.z_range=[np.min(zs),np.max(zs)]

         self.dispersion_rr_params = dispersion_rr_params
         self.dispersion_beta_params = dispersion_beta_params
         self.sigsqr_rr_funct = svdf.returnSigsqrRRFunct(dispersion_rr_params)  # int (beta_inf - beta_inf * 1 / (1 + r/r_beta0 ** alpha_beta) / r)
         self.beta_over_r_int_funct = svdf.returnBetaOverRIntFunct(dispersion_beta_params)
         self.c=c

         self.dphi_over_sigsqr_rr_Nint_funct = svdf.returnDPotOverSig_rr_Nint_Funct(dispersion_rr_params, halo_type, n_interp_rs_for_velocity_dispersion_int, self.c)
         self.M_MW = self.astro_archive.getMWMass()
         self.C=C
         self.halo_center = halo_center[:]
