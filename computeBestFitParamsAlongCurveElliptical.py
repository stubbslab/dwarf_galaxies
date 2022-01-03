

#Do a sequence of the MCMC Series algorithms, with varying run conditions.
#You can change various parameters (present set listed in commented line right above parameter_serires assignment), and have program execute in sequence
#This script can be run from command line without having to be in python environment (bash$ python doMCMCSeriesAllPops.py)
#This method tries to be as efficient as possible, reading in the various potentialFunctionArrays only once, as that can be a time consuming part of the process.

from runMCMCForDwarfGalaxyProfilesAllPops import runMCMCForDwarfGalaxyProfilesAllPops
from PotentialFunctionArray import PotentialFunctionArray
from PotentialArchive import PotentialArchive
from ComputationalArchive import ComputationalArchive
from SurfaceBrightnessProfile import SurfaceBrightnessProfile
import DwarfGalaxyParametersStorer as dgps
from AstronomicalParameterArchive import AstronomicalParameterArchive
from DwarfGalDataArchive import DwarfGalDataArchive
from ObservedGalaxyStarData import ObservedGalaxyStarData
from readInPotentialInterpolator import readInPotentialInterpolator
from compRZLos import compRZLos
from GalaxyMask import GalaxyMask
import SharedObjectHolder as soh
import SurfaceBrightnessProfile as sbp
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
    Gyr_to_s = 365.2425 * 24 * 3600 * 10.0 ** 9.0

    saverun = 1
    save_array = []
    results_number = 'elliptical_best_fit_'
    gamma = astro_archive.getGamma()

    #                     halo       disk          gal       pops               pop_select      step,  WD, apply_observation_mask,  M                rs                      halo_sym_axis                  halo_center                        Rd                     el     eps         disk_sym_axis                disk_center                        lam

    halo_index = 0
    disk_index = halo_index + 1
    gal_index = disk_index + 1
    pops_index = gal_index + 1
    pop_select_index = pops_index + 1
    withDisk_index = pop_select_index + 1
    n_sky_pixels_index = withDisk_index + 1
    n_los_bins_index = n_sky_pixels_index + 1
    mask_R_index = n_los_bins_index + 1
    mask_z_index = mask_R_index + 1
    include_morph_index = mask_z_index + 1
    include_kinem_index = include_morph_index + 1
    apply_observation_mask_index = include_kinem_index + 1
    M_index = apply_observation_mask_index + 1
    rs_index = M_index + 1
    h_center_index = rs_index + 1
    el_index = h_center_index + 1
    phi_index = el_index + 1
    theta_index = phi_index + 1
    kinem_index = theta_index + 1

    #                     halo       disk         gal      pops              pop_select     withDisk  n_sky_pix, n_los_bins, mask_R_bins, mark_z_bins, incl_morph?, incl_kin?, mask? M           rs         halo_center              el,    phi,   theta,  sigsqr_RR                                                                      kin beta,
    #Gaussian-on-one-variable peaks
    parameter_series = [
                        ['nfw',     'sech_disk',  'fornax', ['MP','IM','MR'], 'metal_rigid', 0,        800,      25,          100,         100,         1,           1,         1,    1.4702e+10, 5.079e+03, [-0.0272, 0.0, -0.0570], 2.140,  0.816, 0.0,    [[158.6, 0.42 / Gyr_to_s], [127.7, 0.084 / Gyr_to_s], [113.2, 0.23 / Gyr_to_s]] ], #Elliptical, Prolate, NFW best
                        ['nfw',     'sech_disk',  'fornax', ['MP','IM','MR'], 'metal_rigid', 0,        800,      25,          100,         100,         1,           1,         1,    1.1402e+10, 8.016e+03, [-0.0260, 0.0, -0.0576], 0.412,  2.390, 0.0,    [[158.3, 0.20 / Gyr_to_s], [127.2, -1.34 / Gyr_to_s], [114.0, -1.56 / Gyr_to_s]] ], #Elliptical, Oblate, NFW best #Sph Iso NFW with Aniso params
                        ['burkert',  'sech_disk',  'fornax', ['MP','IM','MR'], 'metal_rigid', 0,        800,      25,          100,         100,         1,           1,         1,    2.249e+9,   0.7624e+03, [-0.0279, 0.0, -0.0571], 1.930,  0.820, 0.0,    [[163.4, 0.30 / Gyr_to_s], [128.5, -0.051 / Gyr_to_s], [110.5, 0.093 / Gyr_to_s]] ],#Elliptical, Prolate, Burkert best
                        ['burkert',  'sech_disk',  'fornax', ['MP','IM','MR'], 'metal_rigid', 0,        800,      25,          100,         100,         1,           1,         1,    2.489e+9,   1.388e+03, [-0.0284, 0.0, -0.0571], 0.466,  2.390, 0.0,    [[161.4, 0.14 / Gyr_to_s], [127.7, -1.33 / Gyr_to_s], [111.1, -1.50 / Gyr_to_s]]  ],#Elliptical, Oblate, Burkert best
                        #['burkert', 'fornax', ['MP','IM','MR'], 'metal_rigid', 500,      21,          100,         100,         1,           1,         1,    1.311e+9, 7.60e+02 , [-0.0182, 0.0, -0.0659], [[174.2, 1.0, 1000.0, 1.0], [136.0, 1.0, 1000.0, 1.0], [112.9, 1.0, 1000.0, 1.0]], [[0.0, 1883.4, 1.0], [0.0, 1030.8, 1.0], [0.0, 1796.3, 1.0]] ],#Sph, Iso Burkert with Aniso params
                        #['burkert', 'fornax', ['MP','IM','MR'], 'metal_rigid', 500,      21,          100,         100,         1,           1,         1,    1.311e+9, 7.60e+02 , [-0.0182, 0.0, -0.0659], [[174.2, 1.0, 1000.0, 1.0], [136.0, 1.0, 1000.0, 1.0], [112.9, 1.0, 1000.0, 1.0]], [[0.61, 1883.4, 1.0], [0.0, 1030.8, 1.0], [0.0, 1796.3, 1.0]] ],#Sph, Iso Burkert with MP Aniso params
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
        halo_funct = readInPotentialInterpolator(halo_type)
        disk_type = parameter_set[disk_index]
        disk_funct = readInPotentialInterpolator(disk_type)
        galaxy = parameter_set[gal_index]
        M = parameter_set[M_index]
        rs = parameter_set[rs_index]
        dist = dSph_archive.getDistanceFromSun([galaxy,'dummy_pop_variable'])
        #phi = parameter_set[phi_index]
        #theta = parameter_set[theta_index]
        halo_center = parameter_set[h_center_index]
        el = parameter_set[el_index ]
        phi = parameter_set[phi_index]
        theta = parameter_set[theta_index]
        pops = parameter_set[pops_index]
        if len(pops) == 0:
            pops = dSph_archive.getPopulations(galaxy)
        pop_selection_method = parameter_set[pop_select_index]
        withDisk = parameter_set[withDisk_index]
        target_n_sky_pixels = parameter_set[n_sky_pixels_index]
        n_los_steps = parameter_set[n_los_bins_index]
        mask_R_bins = parameter_set[mask_R_index]
        mask_z_bins = parameter_set[mask_z_index]
        include_morph_prob = parameter_set[include_morph_index]
        include_kinem_prob = parameter_set[include_kinem_index]
        apply_observation_mask = parameter_set[apply_observation_mask_index]
        kinem_params = parameter_set[kinem_index]
        print ('kinem_params = ' + str(kinem_params))

        deg_limits_R,deg_limits_z = dSph_archive.getObservationBounds([galaxy,'dummy_pop_var'],return_unit = 'degrees')

        max_sky_angle_deg = np.max(np.abs([deg_limits_R ,deg_limits_z]))
        mask_R = np.linspace(-max_sky_angle_deg, max_sky_angle_deg, mask_R_bins)
        mask_z = np.linspace(-max_sky_angle_deg, max_sky_angle_deg, mask_z_bins)

        star_data = []
        for i in range(len(pops)):
            pop = pops[i]
            star_data = star_data + [ObservedGalaxyStarData([galaxy,pop], pop_selection_method = pop_selection_method)]
        dSphParamStorers = [ dgps.DwarfGalaxyParametersStorer( [galaxy, pops[i]], withDisk = withDisk, dist = dist, M = M, rs = rs, halo_center = halo_center, el = el, phi = phi, theta = theta,
                                                              sigsqr_RR = kinem_params[i][0], omega_phi = kinem_params[i][1], stellar_data = star_data[i] )
                            for i in range(len(pops)) ]


        print ('Generating shared objects (including sky mask)...')
        shared_object_holder = soh.SharedObjects(dSphParamStorers[0].dist, mask_R, mask_z, np.linspace(-1.0, 1.0, n_los_steps), target_n_sky_pixels, max_sky_angle_deg * deg_to_rad, compute_mask = apply_observation_mask)
        print ('Computing surface probability of chosen params...')
        #[storer.printContents() for storer in dSphParamStorers]
        s_probs = [sbp.SurfaceBrightnessProfile(storer, shared_object_holder, halo_type, disk_type, disk_interpolating_function=disk_funct, halo_interpolating_function=halo_funct, include_morph_prob = include_morph_prob, include_kinem_prob = include_kinem_prob ) for storer in dSphParamStorers]

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

        #                  halo       disk   gal      pops  pop_select   withDisk  n_sky_pix, n_los_bins,  mask_R_bins, mark_z_bins, incl_morph?,        incl_kin?,     mask?       M    rs    'c'    halo_center      el,    phi,   theta,    sigsqr_RR       omega_phi
        header = ','.join(['halo', 'disk',   'gal', 'pops', 'popSelect', 'WD'  ,    'nSkyPix', 'nLosBins',  'maskRBins', 'maskZBins', 'InclMorphProb?', 'InclKinProb', 'ApplyMask?', 'M', 'rs', 'c',  'haloc', ' ',' ', 'el', 'phi',  'theta', 'sigsqr_RR_MP', 'omega_phi_MP','sigsqr_RR_IM', 'omega_phi_IM','sigsqr_RR_IM', 'omega_phi_IM', 'logLikelihood'])
        save_array = save_array + [[halo_type, disk_type, galaxy, '/'.join(pops), pop_selection_method, withDisk, target_n_sky_pixels, n_los_steps, mask_R_bins, mask_z_bins, include_morph_prob, include_kinem_prob, apply_observation_mask]
                                    + [dSphParamStorers[0].M, dSphParamStorers[0].rs, dSphParamStorers[0].c]
                                    + halo_center
                                    + [dSphParamStorers[0].el, dSphParamStorers[0].phi, dSphParamStorers[0].theta]
                                    + [c.flattenListOfLists(kinem_params)]
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
