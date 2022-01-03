

#Do a sequence of the MCMC Series algorithms, with varying run conditions.
#You can change various parameters (present set listed in commented line right above parameter_serires assignment), and have program execute in sequence
#This script can be run from command line without having to be in python environment (bash$ python doMCMCSeriesAllPops.py)
#This method tries to be as efficient as possible, reading in the various potentialFunctionArrays only once, as that can be a time consuming part of the process.

from runMCMCForDwarfGalaxyProfilesAllPops import runMCMCForDwarfGalaxyProfilesAllPops
from PotentialFunctionArray import PotentialFunctionArray
from PotentialArchive import PotentialArchive
from ComputationalArchive import ComputationalArchive
from SurfaceBrightnessProfile import SurfaceBrightnessProfile
import SphericalGalaxyParametersStorer as sgps
from AstronomicalParameterArchive import AstronomicalParameterArchive
from DwarfGalDataArchive import DwarfGalDataArchive
from ObservedGalaxyStarData import ObservedGalaxyStarData
from readInPotentialInterpolator import readInPotentialInterpolator
from compRZLos import compRZLos
from GalaxyMask import GalaxyMask
import SphericalSharedObjectHolder as ssoh
import SphericalSurfaceBrightnessProfile as ssbp
import math
import numpy as np
import cantrips as c

if __name__ == '__main__':
    #parameters at present are: [el, lam, galaxy, number of iterations, population selection method, start parameter index]
    print ('Starting computation of likelihood series. ')
    astro_archive = AstronomicalParameterArchive ()
    deg_to_rad = astro_archive.getDegToRad()
    compute_params_archive = ComputationalArchive()
    dSph_archive = DwarfGalDataArchive()
    r_half_light = 791.0
    M_star = 10.0 ** 7.39

    saverun = 1
    save_array = []
    results_number = 'spherical_best_fit_'
    gamma = astro_archive.getGamma()

    #                     halo       disk          gal       pops               pop_select      step,  WD, apply_observation_mask,  M                rs                      halo_sym_axis                  halo_center                        Rd                     el     eps         disk_sym_axis                disk_center                        lam

    halo_index = 0
    gal_index = halo_index + 1
    pops_index = gal_index + 1
    pop_select_index = pops_index + 1
    n_sky_pixels_index = pop_select_index + 1
    n_los_bins_index = n_sky_pixels_index + 1
    mask_R_index = n_los_bins_index + 1
    mask_z_index = mask_R_index + 1
    include_morph_index = mask_z_index + 1
    include_kinem_index = include_morph_index + 1
    apply_observation_mask_index = include_kinem_index + 1
    M_index = apply_observation_mask_index + 1
    rs_index = M_index + 1
    h_center_index = rs_index + 1
    sigsqr_rr_index = h_center_index + 1
    kinem_beta_index = sigsqr_rr_index + 1

    #                     halo       gal      pops              pop_select     n_sky_pix, n_los_bins, mask_R_bins, mark_z_bins, incl_morph?, incl_kin?, mask?  M         rs         halo_center              sigsqr_rr                                                                      kin beta,
    #Gaussian-on-one-variable peaks
    parameter_series = [
                        ['nfw',     'fornax', ['MP','IM','MR'], 'metal_rigid', 500,      21,          100,          100,         1,           1,         1,    14.4e+9, 7.5e+03, [-0.0186, 0.0, -0.0653], [[159.0, 1.0, 1000.0, 1.0], [127.5, 1.0, 1000.0, 1.0], [107.4, 1.0, 1000.0, 1.0]], [[0.0, 1000.0, 1.0], [0.0, 1000.0, 1.0], [0.0, 1000.0, 1.0]] ],#Sph, Iso, NFW best
                        ['nfw',     'fornax', ['MP','IM','MR'], 'metal_rigid', 500,      21,          100,          100,         1,           1,         1,    14.5e+9, 7.3e+03, [-0.0186, 0.0, -0.0652], [[168.9, 1.0, 1000.0, 1.0], [132.4, 1.0, 1000.0, 1.0], [112.0, 1.0, 1000.0, 1.0]], [[0.86, 2900.0, 1.0], [0.05, 3000.0, 1.0], [0.24, 2400.0, 1.0]] ],#Sph, Aniso, NFW best
                        #['nfw',     'fornax', ['MP','IM','MR'], 'metal_rigid', 500,      21,          100,          100,         1,           1,         1,    14.25e+9, 7.484e+03, [-0.0281, 0.0, -0.0576], [[158.6, 1.0, 1000.0, 1.0], [127.3, 1.0, 1000.0, 1.0], [107.3, 1.0, 1000.0, 1.0]], [[0.0, 1000.0, 1.0], [0.0, 1000.0, 1.0], [0.0, 1000.0, 1.0]] ],#Sph, Iso, NFW best
                        #['nfw',     'fornax', ['MP','IM','MR'], 'metal_rigid', 500,      21,          100,          100,         1,           1,         1,    14.31e+9, 7.205e+03, [-0.0275, 0.0, -0.0573], [[169.3, 1.0, 1000.0, 1.0], [133.0, 1.0, 1000.0, 1.0], [112.5, 1.0, 1000.0, 1.0]], [[0.57, 2553, 1.0], [0.11, 3114, 1.0], [0.18, 1628, 1.0]] ],#Sph, Aniso, NFW best
                        #['nfw',     'fornax', ['MP','IM','MR'], 'metal_rigid', 500,      21,          100,          100,         1,           1,         1,    1.4315e+10, 7.205e+03, [-0.0187, 0.0, -0.0652], [[169.2, 1.0, 1000.0, 1.0], [133.0, 1.0, 1000.0, 1.0], [112.6, 1.0, 1000.0, 1.0]], [[0.56, 2356, 1.0], [0.12, 2121, 1.0], [0.19, 891.0, 1.0]] ],#Sph, Aniso, NFW Other
                        #['nfw',     'fornax', ['MP','IM','MR'], 'metal_rigid', 500,      21,          100,          100,         1,           1,         1,    1.4315e+10, 7.205e+03, [-0.0189, 0.0, -0.0653], [[169.2, 1.0, 1000.0, 1.0], [133.0, 1.0, 1000.0, 1.0], [112.6, 1.0, 1000.0, 1.0]], [[0.0, 2441, 1.0], [0.0, 2836, 1.0], [0.0, 1446, 1.0]] ], #Sph Iso NFW with Aniso params
                        #['nfw',     'fornax', ['MP','IM','MR'], 'metal_rigid', 500,      21,          100,          100,         1,           1,         1,    1.4315e+10, 7.205e+03, [-0.019, 0.0, -0.065],   [[169.2, 1.0, 1000.0, 1.0], [133.0, 1.0, 1000.0, 1.0], [112.6, 1.0, 1000.0, 1.0]], [[0.56, 2356, 1.0], [0.0, 2121, 1.0], [0.0, 79.1, 1.0]] ], #Sph Iso NFW with MP Aniso params
                        ['burkert', 'fornax', ['MP','IM','MR'], 'metal_rigid', 500,      21,          100,          100,         1,           1,         1,    1.15e+9, 7.3e+02, [-0.0188, 0.0, -0.0657], [[160.7, 1.0, 1000.0, 1.0], [128.2, 1.0, 1000.0, 1.0], [106.6, 1.0, 1000.0, 1.0]], [[0.0, 1000.0, 1.0], [0.0, 1000.0, 1.0], [0.0, 1000.0, 1.0]] ],#Sph, Iso, Burkert best
                        ['burkert', 'fornax', ['MP','IM','MR'], 'metal_rigid', 500,      21,          100,          100,         1,           1,         1,    1.35e+9, 7.7e+02 , [-0.0184, 0.0, -0.0656], [[175.9, 1.0, 1000.0, 1.0], [135.6, 1.0, 1000.0, 1.0], [113.1, 1.0, 1000.0, 1.0]], [[0.74, 1600.0, 1.0], [0.09, 1700.0, 1.0], [0.28, 2400.0, 1.0]] ],#Sph, Anso, Burkert best
                        #['burkert', 'fornax', ['MP','IM','MR'], 'metal_rigid', 500,      21,          100,          100,         1,           1,         1,    1.24e+9, 7.66e+02, [-0.0289, 0.0, -0.0578], [[161.6, 1.0, 1000.0, 1.0], [128.3, 1.0, 1000.0, 1.0], [106.3, 1.0, 1000.0, 1.0]], [[0.0, 1000.0, 1.0], [0.0, 1000.0, 1.0], [0.0, 1000.0, 1.0]] ],#Sph, Iso, Burkert best
                        #['burkert', 'fornax', ['MP','IM','MR'], 'metal_rigid', 500,      21,          100,          100,         1,           1,         1,    1.26e+9, 7.42e+02 , [-0.0284, 0.0, -0.0571], [[174.2, 1.0, 1000.0, 1.0], [135.1, 1.0, 1000.0, 1.0], [112.3, 1.0, 1000.0, 1.0]], [[0.57, 1770.0, 1.0], [0.12, 2184.0, 1.0], [0.21, 1832.0, 1.0]] ],#Sph, Anso, Burkert best
                        #['burkert', 'fornax', ['MP','IM','MR'], 'metal_rigid', 500,      21,          100,          100,         1,           1,         1,    1.311e+9, 7.60e+02 , [-0.0182, 0.0, -0.0659], [[174.2, 1.0, 1000.0, 1.0], [136.0, 1.0, 1000.0, 1.0], [112.9, 1.0, 1000.0, 1.0]], [[0.0, 1883.4, 1.0], [0.0, 1030.8, 1.0], [0.0, 1796.3, 1.0]] ],#Sph, Iso Burkert with Aniso params
                        #['burkert', 'fornax', ['MP','IM','MR'], 'metal_rigid', 500,      21,          100,          100,         1,           1,         1,    1.311e+9, 7.60e+02 , [-0.0182, 0.0, -0.0659], [[174.2, 1.0, 1000.0, 1.0], [136.0, 1.0, 1000.0, 1.0], [112.9, 1.0, 1000.0, 1.0]], [[0.61, 1883.4, 1.0], [0.0, 1030.8, 1.0], [0.0, 1796.3, 1.0]] ],#Sph, Iso Burkert with MP Aniso params
                        #['nfw',     'fornax', ['MP','IM','MR'], 'metal_rigid', 500,      21,          100,          100,         1,           1,         1,    1.4315e+10, 7.205e+03, [-0.019, 0.0, -0.065], [[169.2, 1.0, 1000.0, 1.0], [133.0, 1.0, 1000.0, 1.0], [112.6, 1.0, 1000.0, 1.0]], [[0.0, 2356, 1.0], [0.0, 2121, 1.0], [0.0, 79.1, 1.0]] ], #Sph Iso NFW with MP Aniso params
                        #['burkert', 'fornax', ['MP','IM','MR'], 'metal_rigid', 500,      21,          100,          100,         1,           1,         1,    1.164e+9, 7.36e+02, [-0.0180, 0.0, -0.0663], [[161.6, 1.0, 1000.0, 1.0], [128.3, 1.0, 1000.0, 1.0], [106.3, 1.0, 1000.0, 1.0]], [[0.0, 1000.0, 1.0], [0.0, 1000.0, 1.0], [0.0, 1000.0, 1.0]] ],#Sph, Iso, Burkert best
                         ]

    #arcmin_limits_z =50.0
    likelihood_array = []
    save_array = []
    #Now actually compute the likelihood for each of the specified point.
    for parameter_set in parameter_series:
        print ('Computing likelihood for parameter_set: ' + str(parameter_set))
        parameter_storers = []
        surface_profiles = []


        halo_type = parameter_set[halo_index]
        galaxy = parameter_set[gal_index]
        M = parameter_set[M_index]
        rs = parameter_set[rs_index]
        dist = dSph_archive.getDistanceFromSun([galaxy,'dummy_pop_variable'])
        #phi = parameter_set[phi_index]
        #theta = parameter_set[theta_index]
        halo_center = parameter_set[h_center_index]
        pops = parameter_set[pops_index]
        if len(pops) == 0:
            pops = dSph_archive.getPopulations(galaxy)
        pop_selection_method = parameter_set[pop_select_index]
        target_n_sky_pixels = parameter_set[n_sky_pixels_index]
        n_los_steps = parameter_set[n_los_bins_index]
        mask_R_bins = parameter_set[mask_R_index]
        mask_z_bins = parameter_set[mask_z_index]
        include_morph_prob = parameter_set[include_morph_index]
        include_kinem_prob = parameter_set[include_kinem_index]
        apply_observation_mask = parameter_set[apply_observation_mask_index]
        dispersion_rr_params = parameter_set[sigsqr_rr_index]
        kinem_beta_params = parameter_set[kinem_beta_index]
        print ('[dispersion_rr_params, kinem_beta_params] = ' + str([dispersion_rr_params, kinem_beta_params]))

        deg_limits_R,deg_limits_z = dSph_archive.getObservationBounds([galaxy,'dummy_pop_var'],return_unit = 'degrees')

        max_sky_angle_deg = np.max(np.abs([deg_limits_R ,deg_limits_z]))
        mask_R = np.linspace(-max_sky_angle_deg, max_sky_angle_deg, mask_R_bins)
        mask_z = np.linspace(-max_sky_angle_deg, max_sky_angle_deg, mask_z_bins)

        star_data = []
        for i in range(len(pops)):
            pop = pops[i]
            star_data = star_data + [ObservedGalaxyStarData([galaxy,pop], pop_selection_method = pop_selection_method)]
        dSphParamStorers = [ sgps.DwarfGalaxyParametersStorer( [galaxy, pops[i]], dist = dist, M = M, rs = rs, halo_center = halo_center,
                                                              sigsqr_rr_0 = dispersion_rr_params[i][0], sigsqr_rr_inf_to_0_rat = dispersion_rr_params[i][1], r_sigsqr_rr0 = dispersion_rr_params[i][2], alpha_sigsqr_rr = dispersion_rr_params[i][3],
                                                              gamma_for_beta_inf = kinem_beta_params[i][0], r_beta0 = kinem_beta_params[i][1], alpha_beta = kinem_beta_params[i][2], stellar_data = star_data[i] )
                            for i in range(len(pops)) ]


        print ('Generating shared objects (including sky mask)...')
        shared_object_holder = ssoh.SharedObjects(dSphParamStorers[0].dist, mask_R, mask_z, np.linspace(-1.0, 1.0, n_los_steps), target_n_sky_pixels, max_sky_angle_deg * deg_to_rad, compute_mask = apply_observation_mask)
        print ('Computing surface probability of chosen params...')
        #[storer.printContents() for storer in dSphParamStorers]
        s_probs = [ssbp.SurfaceBrightnessProfile(storer, shared_object_holder, halo_type, include_morph_prob = include_morph_prob, include_kinem_prob = include_kinem_prob ) for storer in dSphParamStorers]
        log_likelihood = 0.0


        for i in range(len(pops)):
            population = [galaxy,pops[i]]
            log_probability_of_stars = s_probs[i].getStarProbabilityValues(star_data[i].proj_x, star_data[i].proj_y, star_data[i].corr_Vhel)
            log_likelihood = log_likelihood + log_probability_of_stars

            #plt.scatter(np.array(star_data.corrRa) * 60.0 ,np.array(star_data.corrDec) * 60.0)
        pop_string = ' '.join(pops)
        halo_sym_string = ''
        halo_center_string = ''
        if halo_center is 'none':
            halo_center_string = 'none'
        else:
            halo_center_string = ' '.join([str(elem) for elem in halo_center])

        #                  halo     gal    pops    pop_select   n_sky_pix, n_los_bins,  mask_R_bins, mark_z_bins, incl_morph?,      incl_kin?,     mask?         M   rs       halo_center                    sigsqr_rr                 kin beta,
        header = ','.join(['halo', 'gal', 'pops', 'popSelect', 'nSkyPix', 'nLosBins',  'maskRBins', 'maskZBins', 'InclMorphProb?', 'InclKinProb', 'ApplyMask?', 'M','rs', 'c', 'haloc', ' ',' ', 'sigsqr_rrParams_MP', ' ',' ', ' ', 'sigsqr_rrParams_IM', ' ',' ', ' ', 'sigsqr_rrParams_MR', ' ',' ', ' ', 'kinemBetaParams_MP', ' ',' ', 'kinemBetaParams_IM', ' ',' ', 'kinemBetaParams_MR', ' ', ' ', 'logLikelihood'])
        save_array = save_array + [[halo_type, galaxy, '/'.join(pops), pop_selection_method, target_n_sky_pixels, n_los_steps, mask_R_bins, mask_z_bins, include_morph_prob, include_kinem_prob, apply_observation_mask]
                                    + [dSphParamStorers[0].M, dSphParamStorers[0].rs, dSphParamStorers[0].c]
                                    + halo_center
                                    + [c.flattenListOfLists(dispersion_rr_params)]
                                    + [c.flattenListOfLists(kinem_beta_params)]
                                    + [log_likelihood]
                                   ]

                                    #+ parameter_set[0:halo_sym_axis_index] + [halo_sym_string]
                                    #+ parameter_set[halo_sym_axis_index+1:halo_center_index] + [halo_center_string]
                                    #+ parameter_set[halo_center_index+1:disk_sym_axis_index] + [disk_sym_string]
                                    #+ parameter_set[disk_sym_axis_index+1:disk_center_index] + [disk_center_string]
                                    #+ parameter_set[disk_center_index+1:pops_index] + [pop_string] + parameter_set[pops_index+1:] + [log_likelihood]]
        likelihood_array = likelihood_array + [log_likelihood]


    #print likelihood_array
    #save_array = [parameter_series[i] + [likelihood_array[i]] for i in range(len(likelihood_array)) ]
    #print save_array
    #print np.array(save_array)

    if saverun:
        save_dir = compute_params_archive.getSpotProbabilityDir()
        print ('save_dir = ' + save_dir )
        file_name = 'best_fit_probabilities_number_' + str(results_number) + '.csv'
        print ('Saving series to ' + save_dir + file_name)
        #header = 'halo, disk, gal, pops, pop_select, step, WD, obs_mask, dist, omega_phi, sigsqr, c, M, rs, hsym_x, hsym_y, hsym_z, phi, theta, hc_x, hc_y, hc_z, Rd, el, eps, hsym_x, hsym_y, hsym_z, a, b, hc_x, hc_y, hc_z, lam, logLikelihood'
        np.savetxt(save_dir + file_name, np.array(save_array), delimiter = ',',fmt = '%10s',header = header)
