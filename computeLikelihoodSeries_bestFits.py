

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
    Rd_index = theta_index + 1
    lam_index = Rd_index + 1
    eps_index = lam_index + 1
    a_index = eps_index + 1
    b_index = a_index + 1
    d_center_index = b_index + 1
    kinem_index = d_center_index + 1

    #                     halo       disk          gal      pops              pop_select     withDisk  n_sky_pix, n_los_bins, mask_R_bins, mark_z_bins, incl_morph?, incl_kin?, mask? M          rs         halo_center              el,    phi,   theta,  Rd,     lam,  eps, a,                  b,   disk_center,             sigsqr_RR                                                                      kin beta,
    #Gaussian-on-one-variable peaks
    parameter_series = [
                        ['nfw',      'sech_disk',  'fornax', ['MP','IM','MR'], 'metal_rigid', 0,        500,      21,          100,         100,         1,           1,         1,    1592548141.18399, 5400.31435994469, [-0.0968933428949802,	0, -0.220285937663979], 8.02010126564584,  1.82403771653697, 0.0,   1000.0, 0.1,  0.0, 178.0 * deg_to_rad, 0.0, [0.0,   0.0,   0.0   ], [[317.594020814284,	3.16887385068114E-16 ], [308.181329447375,	3.16887385068114E-16], [304.633452428443,	3.16887385068114E-16]] ], #Verify I get same results as chains: should be 1401.33508154923
                        ['nfw',      'sech_disk',  'fornax', ['MP','IM','MR'], 'metal_rigid', 0,        500,      21,          100,         100,         1,           1,         1,    10175852629.8668, 3654.1501386456, [-0.278486566088987,	0,	-0.422721469532029], 0.129432961171195,  1.62679870205106, 0.0,   1000.0, 0.1,  0.0, 178.0 * deg_to_rad, 0.0, [0.0,   0.0,   0.0   ], [[306.204379324607,	3.16887385068114E-16], [308.531456519294,	3.16887385068114E-16], [311.242007141175,	3.16887385068114E-16]] ], #Verify I get same results as chains: should be -527.7101908
                        ['burkert',      'sech_disk',  'fornax', ['MP','IM','MR'], 'metal_rigid', 0,        500,      21,          100,         100,         1,           1,         1,    10214468954.6258, 3677.10992978586, [-0.27686400289225,	0,	-0.433039947463569], 0.129907733337948,  1.5808116555335, 0.0,   1000.0, 0.1,  0.0, 178.0 * deg_to_rad, 0.0, [0.0,   0.0,   0.0   ], [[305.761366879201,	3.16887385068114E-16], [312.791391965791,	3.16887385068114E-16], [305.970708075373,	3.16887385068114E-16]] ], #Verify I get same results as chains: should be 869.2788919
                        ['burkert',      'sech_disk',  'fornax', ['MP','IM','MR'], 'metal_rigid', 0,        500,      21,          100,         100,         1,           1,         1,    10402629093.017216, 3893.763856755854, [-0.2854811847350313, 0.0, -0.42129771963630824], 0.12257976197336011,  1.583549927482053, 0.0,   1000.0, 0.1,  0.0, 178.0 * deg_to_rad, 0.0, [0.0,   0.0,   0.0   ], [[305.52977305196777, 3.143705904145926e-16], [301.3333048905686, 3.1617281233632034e-16], [ 306.63212675457714, 3.1519314358109243e-16]] ], #Verify I get same results as chains: should be 932.0868788141978
                        ['nfw',      'sech_disk',  'fornax', ['MP','IM','MR'], 'metal_rigid', 1,        500,      21,          100,         100,         1,           1,         1,    11820887925.1032, 5421.92333547376, [-0.196547489952411, 0.0, 0.268462045561764], 0.426,  2.385, 0.0,   861.942751775546, 0.1,  0.00615526190771908, 178.0 * deg_to_rad, 0.0, [-0.0047969679219269,   0.0,   0.0089316709258229  ], [[295.997985418317, -1.8808556838281E-16], [290.582477976986, -1.85755582668944E-16], [ 289.554345668677, -2.05789512995068E-16]] ], #Verify I get same results as chains: should be 911.794754440306
                        #['burkert',  'sech_disk',  'fornax', ['MP','IM','MR'], 'metal_rigid', 0,        500,      21,          100,         100,         1,           1,         1,     2.11e+9, 0.729e+03, [-0.0289, 0.0, -0.0578], 1.980,  0.817, 0.0,   1000.0, 0.1,  0.0, 178.0 * deg_to_rad, 0.0, [0.0,   0.0,   0.0   ], [[162.1, 0.38 / Gyr_to_s], [128.5, 0.02 / Gyr_to_s], [111.3, 0.33 / Gyr_to_s]] ],#Verify I get same results as chains: should be XXX
                        ['nfw',      'sech_disk',  'fornax', ['MP','IM','MR'], 'metal_rigid', 0,        500,      21,          100,         100,         1,           1,         1,    14.6e+9, 5.1e+03, [-0.0289, 0.0, -0.0570], 2.130,  0.815, 0.0,   1000.0, 0.1,  0.0, 178.0 * deg_to_rad, 0.0, [0.0,   0.0,   0.0   ], [[156.9, 0.56 / Gyr_to_s], [127.1, -0.03 / Gyr_to_s], [109.4, -0.06 / Gyr_to_s]] ], #Verify I get same results as chains: should be XXX
                        ['nfw',      'sech_disk',  'fornax', ['MP','IM','MR'], 'metal_rigid', 0,        500,      21,          100,         100,         1,           1,         1,    10.8e+9, 7.8e+03, [-0.0282, 0.0, -0.0571], 0.426,  2.385, 0.0,   1000.0, 0.1,  0.0, 178.0 * deg_to_rad, 0.0, [0.0,   0.0,   0.0   ], [[157.9, -0.095 / Gyr_to_s], [127.6, -1.36 / Gyr_to_s], [109.7, -0.97 / Gyr_to_s]] ], #Elliptical, Oblate, NFW best
                        ['burkert',  'sech_disk',  'fornax', ['MP','IM','MR'], 'metal_rigid', 0,        500,      21,          100,         100,         1,           1,         1,     2.3e+9, 0.76e+03, [-0.0302, 0.0, -0.0572], 1.950,  0.817, 0.0,   1000.0, 0.1,  0.0, 178.0 * deg_to_rad, 0.0, [0.0,   0.0,   0.0   ], [[162.8, 0.14 / Gyr_to_s], [127.9, -0.069 / Gyr_to_s], [106.3, -0.35 / Gyr_to_s]] ],#Elliptical, Prolate, Burkert best
                        ['burkert',  'sech_disk',  'fornax', ['MP','IM','MR'], 'metal_rigid', 0,        500,      21,          100,         100,         1,           1,         1,     2.6e+9, 1.42e+03, [-0.0289, 0.0, -0.0570], 0.463,  2.384, 0.0,   1000.0, 0.1,  0.0, 178.0 * deg_to_rad, 0.0, [0.0,   0.0,   0.0   ], [[161.6, -0.04 / Gyr_to_s], [127.5, -1.37 / Gyr_to_s], [106.5, -0.94 / Gyr_to_s]]  ],#Elliptical, Oblate, Burkert best
                        ['nfw',      'sech_disk',  'fornax', ['MP','IM','MR'], 'metal_rigid', 1,        500,      21,          100,         100,         1,           1,         1,    14.6e+9, 5.11e+03, [-0.0287, 0.0, -0.0572], 2.130,  0.815, 0.0,   1558.0, 0.1,  0.012, 178.0 * deg_to_rad, 0.0, [0.0034,   0.0,   0.036   ], [[154.6, 0.61 / Gyr_to_s], [124.9, 0.092 / Gyr_to_s], [107.1, -0.115 / Gyr_to_s]] ], #Verify I get same results as chains: should be XXX
                        ['nfw',      'sech_disk',  'fornax', ['MP','IM','MR'], 'metal_rigid', 1,        500,      21,          100,         100,         1,           1,         1,    10.8e+9, 7.8e+03, [-0.0283, 0.0, -0.0571], 0.426,  2.385, 0.0,     79.0, 0.1,  0.015, 178.0 * deg_to_rad, 0.0, [0.158,   0.0,   -0.14   ], [[154.3, -0.14 / Gyr_to_s], [124.2, -1.41 / Gyr_to_s], [106.3, -0.99 / Gyr_to_s]] ], #Elliptical, Oblate, NFW best
                        ['burkert',  'sech_disk',  'fornax', ['MP','IM','MR'], 'metal_rigid', 1,        500,      21,          100,         100,         1,           1,         1,     2.3e+9, 0.76e+03, [-0.0293, 0.0, -0.0575], 1.950,  0.817, 0.0,   2398.0, 0.1,  0.037, 178.0 * deg_to_rad, 0.0, [-0.034,   0.0,   -0.261   ], [[158.3, 0.56 / Gyr_to_s], [126.3, 0.05 / Gyr_to_s], [106.7, -0.17 / Gyr_to_s]] ],#Elliptical, Prolate, Burkert best
                        ['burkert',  'sech_disk',  'fornax', ['MP','IM','MR'], 'metal_rigid', 1,        500,      21,          100,         100,         1,           1,         1,     2.6e+9, 1.42e+03, [-0.0294, 0.0, -0.0573], 0.463,  2.384, 0.0,   2208.0, 0.1,  0.022, 178.0 * deg_to_rad, 0.0, [0.050,   0.0,   0.029   ], [[157.9, -0.10 / Gyr_to_s], [125.3, -1.39 / Gyr_to_s], [105.4, -0.99 / Gyr_to_s]]  ],#Elliptical, Oblate, Burkert best
                        ['nfw',      'sech_disk',  'fornax', ['MP','IM','MR'], 'metal_rigid', 1,        500,      21,          100,         100,         1,           1,         1,    14.6e+9, 5.11e+03, [-0.0287, 0.0, -0.0572], 2.130,  0.815, 0.0,   1558.0, 0.1,  0.0, 178.0 * deg_to_rad, 0.0, [0.0034,   0.0,   0.036   ], [[154.6, 0.61 / Gyr_to_s], [124.9, 0.092 / Gyr_to_s], [107.1, -0.115 / Gyr_to_s]] ], #Verify I get same results as chains: should be XXX
                        ['nfw',      'sech_disk',  'fornax', ['MP','IM','MR'], 'metal_rigid', 1,        500,      21,          100,         100,         1,           1,         1,    10.8e+9, 7.8e+03, [-0.0283, 0.0, -0.0571], 0.426,  2.385, 0.0,     79.0, 0.1,  0.0, 178.0 * deg_to_rad, 0.0, [0.158,   0.0,   -0.14   ], [[154.3, -0.14 / Gyr_to_s], [124.2, -1.41 / Gyr_to_s], [106.3, -0.99 / Gyr_to_s]] ], #Elliptical, Oblate, NFW best
                        ['burkert',  'sech_disk',  'fornax', ['MP','IM','MR'], 'metal_rigid', 1,        500,      21,          100,         100,         1,           1,         1,     2.3e+9, 0.76e+03, [-0.0293, 0.0, -0.0575], 1.950,  0.817, 0.0,   2398.0, 0.1,  0.0, 178.0 * deg_to_rad, 0.0, [-0.034,   0.0,   -0.261   ], [[158.3, 0.56 / Gyr_to_s], [126.3, 0.05 / Gyr_to_s], [106.7, -0.17 / Gyr_to_s]] ],#Elliptical, Prolate, Burkert best
                        ['burkert',  'sech_disk',  'fornax', ['MP','IM','MR'], 'metal_rigid', 1,        500,      21,          100,         100,         1,           1,         1,     2.6e+9, 1.42e+03, [-0.0294, 0.0, -0.0573], 0.463,  2.384, 0.0,   2208.0, 0.1,  0.0, 178.0 * deg_to_rad, 0.0, [0.050,   0.0,   0.029   ], [[157.9, -0.10 / Gyr_to_s], [125.3, -1.39 / Gyr_to_s], [105.4, -0.99 / Gyr_to_s]]  ],#Elliptical, Oblate, Burkert best
                        #['nfw',      'sech_disk',  'fornax', ['MP','IM','MR'], 'metal_rigid', 1,        500,      21,          100,         100,         1,           1,         1,    10.92e+9, 7.735e+03, [-0.0275, 0.0, -0.0573], 0.427,  2.390, 0.0,   1000.0, 0.1,  0.0, 178.0 * deg_to_rad, 0.0, [0.0,   0.0,   0.0   ], [[158.5, 0.12 / Gyr_to_s], [128.1, -1.35 / Gyr_to_s], [113.3, -1.50 / Gyr_to_s]] ], #Do we get the same result with eps = 0, WD = 1 as WD = 0?
                        #['burkert',  'sech_disk',  'fornax', ['MP','IM','MR'], 'metal_rigid', 1,        500,      21,          100,         100,         1,           1,         1,     2.46e+9, 1.382e+03, [-0.0289, 0.0, -0.0578], 0.463,  2.390, 0.0,   1000.0, 0.1,  0.0, 178.0 * deg_to_rad, 0.0, [0.0,   0.0,   0.0   ], [[161.8, 0.13 / Gyr_to_s], [127.7, -1.34 / Gyr_to_s], [110.8, -1.54 / Gyr_to_s]]  ],#Do we get the same result with eps = 0, WD = 1 as WD = 0?
                        #['nfw',      'sech_disk',  'fornax', ['MP','IM','MR'], 'metal_rigid', 1,        500,      21,          100,         100,         1,           1,         1,     433350823.97961, 714.468497662, [-0.0104605328534637,0 ,	-0.00583343655105964], 0.427,  3.1066860685499, 0.0,   776.0803630, 0.1,  0.001724408, 2.39, 0.0, [0.00625200033559256,	0,	0.0070418581106660], [[107.306207695863,	-2.57589996011495E-18], [96.0947776901623,	-5.82198000577824E-18], [110.583916119766,	-2.9544295730564E-18] ]], #Verify I get same results as chains: should be 1955.92474286345
                        #['nfw',      'sech_disk',  'fornax', ['MP','IM','MR'], 'metal_rigid', 1,        500,      21,          100,         100,         1,           1,         1,     494970089.95193, 788.876661781, [-0.0030719177319025,0 ,	0.000739272190392563], 0.427,  3.1066860685499, 0.0,   729.3645518, 0.1,  0.000000000, 2.39, 0.0, [0.00631711132748324,	0,	0.0101654473807818], [[101.759838365943,	-9.46930888617616E-18], [98.3746195977642,	-4.60090573800163E-18], [105.76651678909,	-7.23728325538744E-19]]  ],#erify I get same results as chains: should be 1916.353584
                        #['nfw',      'sech_disk',  'fornax', ['MP','IM','MR'], 'metal_rigid', 1,        500,      21,          100,         100,         1,           1,         1,     494970089.95193, 788.876661781, [-0.0030719177319025,0 ,	0.000739272190392563], 0.427,  3.1066860685499, 0.0,   729.3645518, 0.1,  0.000000000, 2.39, 0.0, [0.00631711132748324,	0,	0.0101654473807818], [[101.759838365943,	0.0], [98.3746195977642,	0.0], [105.76651678909,	0.0]]  ],#erify I get same results as chains: should be 1916.353584
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
        Rd = parameter_set[Rd_index]
        lam = parameter_set[lam_index]
        eps = parameter_set[eps_index]
        a = parameter_set[a_index]
        b = parameter_set[b_index]
        disk_center = parameter_set[d_center_index]
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

        print('[M, rs, dist, halo_center, el, phi, theta, Rd ,  eps, lam , disk_center, a,  b, kinem_params, pops] = ' + str([M, rs, dist, halo_center, el, phi, theta, Rd ,  eps, lam , disk_center, a,  b, kinem_params, pops] ))

        deg_limits_R,deg_limits_z = dSph_archive.getObservationBounds([galaxy,'dummy_pop_var'],return_unit = 'degrees')

        max_sky_angle_deg = np.max(np.abs([deg_limits_R ,deg_limits_z]))
        mask_R = np.linspace(-max_sky_angle_deg, max_sky_angle_deg, mask_R_bins)
        mask_z = np.linspace(-max_sky_angle_deg, max_sky_angle_deg, mask_z_bins)

        star_data = []
        for i in range(len(pops)):
            pop = pops[i]
            star_data = star_data + [ObservedGalaxyStarData([galaxy,pop], pop_selection_method = pop_selection_method)]
        dSphParamStorers = [ dgps.DwarfGalaxyParametersStorer( [galaxy, pops[i]], withDisk = withDisk, dist = dist, M = M, rs = rs, halo_center = halo_center, el = el,
                                                              halo_sym_axis = 'none', disk_sym_axis = 'none',
                                                              phi = phi, theta = theta,
                                                              Rd = Rd, eps = eps, lam = lam, disk_center = disk_center, a = a, b = b,
                                                              sigsqr_RR = kinem_params[i][0], omega_phi = kinem_params[i][1], stellar_data = star_data[i] )
                            for i in range(len(pops)) ]
        print ('[storer.phi for storer in dSphParamStorers ] = ' + str([storer.phi for storer in dSphParamStorers ]))


        print ('Generating shared objects (including sky mask)...')
        shared_object_holder = soh.SharedObjects(dSphParamStorers[0].dist, mask_R, mask_z, np.linspace(-1.0, 1.0, n_los_steps), target_n_sky_pixels, max_sky_angle_deg * deg_to_rad, compute_mask = apply_observation_mask)
        print ('Computing surface probability of chosen params...')
        #[storer.printContents() for storer in dSphParamStorers]
        s_probs = [sbp.SurfaceBrightnessProfile(storer, shared_object_holder, halo_type, disk_type, disk_interpolating_function=disk_funct, halo_interpolating_function=halo_funct, include_morph_prob = include_morph_prob, include_kinem_prob = include_kinem_prob ) for storer in dSphParamStorers]

        log_likelihood = 0.0


        for i in range(len(pops)):
            population = [galaxy,pops[i]]
            dSphParamStorers[i].printContents()
            print('s_probs[i].normalization_constant = ' + str(s_probs[i].normalization_constant ))
            log_probability_of_stars = s_probs[i].getStarProbabilityValues(star_data[i].proj_x, star_data[i].proj_y, star_data[i].corr_Vhel)
            print ('log_probability_of_stars = ' + str(log_probability_of_stars))
            log_likelihood = log_likelihood + log_probability_of_stars

            #plt.scatter(np.array(star_data.corrRa) * 60.0 ,np.array(star_data.corrDec) * 60.0)
        pop_string = ' '.join(pops)
        halo_sym_string = ''
        halo_center_string = ''
        if halo_center is 'none':
            halo_center_string = 'none'
        else:
            halo_center_string = ' '.join([str(elem) for elem in halo_center])
        #                  # gal	 halo	 disk	 pops	 pop_select	 mask?	 WD, n_mask_bins_R	 n_mask_bins_z	 sky_angle	 target_n_sky_pix	 n_los_bins	 dist	 c	 M	 rs	 el	 hsym_x	 hsym_y	 hsym_z	 phi	 theta	 hc_x	 hc_y	 hc_z	 sigsqr_RR_MP	omega_phi_MP	dispersion_b_MP	sigsqr_RR_IM	omega_phi_IM	dispersion_b_IM	sigsqr_RR_MR	omega_phi_MR	dispersion_b_MR	 lam	 Rd	 eps	 a	 b	 d_xhat	 d_yhat	 d_zhat	 d_x_center	 d_y_center	 d_z_center	 nVisits	 logLikelihood
        #                  halo       disk   gal      pops  pop_select   mask?       withDisk  mask_R_bins, mark_z_bins, n_sky_pix, n_los_bins,  inclu_morph?   incl_kin?,         dist    c      M    rs    halo_center      el,    phi,   theta,      Rd    eps   lam    a   b    disk_center   sigsqr_RR       omega_phi
        header = ','.join(['gal', 'halo', 'disk', 'pops', 'popSelect', 'ApplyMask?', 'WD'  ,  'maskRBins', 'maskZBins',   'nSkyPix', 'nLosBins',  'InclMorphProb?', 'InclKinProb', 'dist', 'c',  'M', 'rs',  'haloc', ' ',' ', 'el', 'phi',  'theta', 'Rd',    'eps' ,  'lam' ,   'a' ,  'b' ,   'disk_center', ' ', ' ', 'sigsqr_RR_MP', 'omega_phi_MP','sigsqr_RR_IM', 'omega_phi_IM','sigsqr_RR_IM', 'omega_phi_IM', 'logLikelihood'])
        save_array = save_array + [[galaxy, halo_type, disk_type, '/'.join(pops), pop_selection_method, apply_observation_mask, withDisk, target_n_sky_pixels, n_los_steps, mask_R_bins, mask_z_bins, include_morph_prob, include_kinem_prob]
                                    + [dSphParamStorers[0].dist, dSphParamStorers[0].c, dSphParamStorers[0].M, dSphParamStorers[0].rs]
                                    + halo_center
                                    + [dSphParamStorers[0].el, dSphParamStorers[0].phi, dSphParamStorers[0].theta]
                                    + [dSphParamStorers[0].Rd, dSphParamStorers[0].eps, dSphParamStorers[0].lam, dSphParamStorers[0].a, dSphParamStorers[0].b]
                                    + disk_center
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
