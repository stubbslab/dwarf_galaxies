import numpy as np
import math
from DwarfGalDataArchive import DwarfGalDataArchive
from ObservedGalaxyStarData import ObservedGalaxyStarData
import random
from SurfaceBrightnessProfile import SurfaceBrightnessProfile
from DwarfGalaxyParametersStorer import DwarfGalaxyParametersStorer
from compRZLos import compRZLos
from readInPotentialInterpolator import readInPotentialInterpolator
from AstronomicalParameterArchive import AstronomicalParameterArchive
from runMCMCForSphericalGalaxyProfilesAllPops import runMCMCForDwarfGalaxyProfilesAllPopsSpherical
from cantrips import removeListElement
from GalaxyMask import GalaxyMask
from manageMCMCResults import measureMCMCResults
import StartMCMCSphericalParameterStorer as smcmc
import manageSphericalMCMCResults as msmr
import time
import computeBestFitParamsAlongCurveSpherical as cbfpacs
from ComputationalArchive import ComputationalArchive

def randomlySortStars(galaxy, w_replace = 1, pop_selection_method = 'metal_rigid', n_resamples = 'full', save = 0, file_name = ''):
    dSph_archive = DwarfGalDataArchive()
    pops = dSph_archive.getPopulations(galaxy)
    popGalPairs = []
    orig_star_data = []
    star_data = []
    resampled_star_data = []
    flat_star_data = []
    for i in range(len(pops)):
        pop = pops[i]
        popGalPairs = popGalPairs + [[galaxy,pop]]
        new_stars = ObservedGalaxyStarData([galaxy,pop], pop_selection_method = pop_selection_method)
        star_data = star_data + [new_stars]
        resampled_star_data = resampled_star_data + [new_stars]
        #flat_star_data = flat_star_data + new_stars
    n_stars_per_pop = [len(stars.RA) for stars in star_data]

    n_stars = sum(n_stars_per_pop)

    if n_resamples == 'full':
        n_resamples = n_stars

    original_sample_RAs = [stars.corrRa for stars in star_data]
    original_sample_Decs = [stars.corrDec for stars in star_data]
    original_sample_SigMg_corr = [stars.SigMg_corr for stars in star_data]
    original_sample_corrVHels = [stars.corr_Vhel for stars in star_data]
    original_sample_corrVHelEs = [stars.corr_VhelE for stars in star_data]
    original_sample_proj_x = [stars.proj_x for stars in star_data]
    original_sample_proj_y = [stars.proj_y for stars in star_data]

    resampled_RAs = [[] for stars in star_data]
    resampled_Decs = [[] for stars in star_data]
    resampled_SigMg_corr = [[] for stars in star_data]
    resampled_corrVHels = [[] for stars in star_data]
    resampled_corrVHelEs = [[] for stars in star_data]
    resampled_proj_x = [[] for stars in star_data]
    resampled_proj_y = [[] for stars in star_data]

    n_remaining_stars = n_stars
    sample_RAs = original_sample_RAs[:]
    sample_Decs = original_sample_Decs[:]
    sample_SigMg_corr = original_sample_SigMg_corr[:]
    sample_corrVHels = original_sample_corrVHels[:]
    sample_corrVHelEs = original_sample_corrVHelEs[:]
    sample_proj_x = original_sample_proj_x[:]
    sample_proj_y = original_sample_proj_y[:]


    for sample in range(n_resamples):
        star_index = random.randint(0, n_remaining_stars - 1)

        for j in range(len(sample_RAs)):
            n_stars_in_pop = len(sample_RAs[j])
            if star_index < n_stars_in_pop:
                ras = sample_RAs[j]
                decs = sample_Decs[j]
                SigMg_corr = sample_SigMg_corr[j]
                corrVHels = sample_corrVHels[j]
                corrVHelEs = sample_corrVHelEs[j]
                proj_x = sample_proj_x[j]
                proj_y = sample_proj_y[j]

                resampled_RAs[j] = resampled_RAs[j] + [ras[star_index]]
                resampled_Decs[j] = resampled_Decs[j] + [decs[star_index]]
                resampled_SigMg_corr[j] = resampled_SigMg_corr[j] + [SigMg_corr[star_index]]
                resampled_corrVHels[j] = resampled_corrVHels[j] + [corrVHels[star_index]]
                resampled_corrVHelEs[j] = resampled_corrVHelEs[j] + [corrVHelEs[star_index]]
                resampled_proj_x[j] = resampled_proj_x[j] + [proj_x[star_index]]
                resampled_proj_y[j] = resampled_proj_y[j] + [proj_y[star_index]]

                if w_replace is 0:
                    sample_RAs[j] = np.array(removeListElement((sample_RAs[j]).tolist()  , star_index))
                    sample_Decs[j] = np.array(removeListElement((sample_Decs[j]).tolist() , star_index))
                    sample_SigMg_corr[j] = np.array(removeListElement((sample_SigMg_corr[j]).tolist(), star_index))
                    sample_corrVHels[j] = np.array(removeListElement((sample_corrVHels[j]).tolist(), star_index))
                    sample_corrVHelEs[j] = np.array(removeListElement((sample_corrVHelEs[j]).tolist(), star_index))
                    sample_proj_x[j] = np.array(removeListElement((sample_proj_x[j]).tolist(), star_index))
                    sample_proj_y[j] = np.array(removeListElement((sample_proj_y[j]).tolist(), star_index))
                    n_remaining_stars = n_remaining_stars - 1
                break
            else:
                star_index = star_index - n_stars_in_pop
    for i in range(len(resampled_star_data)):
        resampled_star_data[i].corrRA = resampled_RAs[i]
        resampled_star_data[i].corrDec = resampled_Decs[i]
        resampled_star_data[i].SigMg_corr = resampled_SigMg_corr[i]
        resampled_star_data[i].corr_Vhel = resampled_corrVHels[i]
        resampled_star_data[i].corr_VhelE = resampled_corrVHelEs[i]
        resampled_star_data[i].proj_x = np.array(resampled_proj_x[i])
        resampled_star_data[i].proj_y = np.array(resampled_proj_y[i])
        mean_Vhel = sum(resampled_star_data[i].corr_Vhel) / float(len(resampled_star_data[i].corr_Vhel))
        sigsqr = sum([(vel - mean_Vhel) ** 2.0 for vel in resampled_star_data[i].corr_Vhel]) / float(len(resampled_star_data[i].corr_Vhel))
        print ('sigsqr = ' + str(sigsqr))
        resampled_star_data[i].sigSqr = sigsqr
    if save:
        computational_params_archive = ComputationalArchive()
        header = 'corrRA, corrDec, corr_Vhel, proj_x, proj_y'
        save_dir = computational_params_archive.getMCMCOutDir() + 'stellarRandomizations/'
        file_name = file_name + '.csv'
        np.savetxt(save_dir + file_name, [resampled_RAs, resampled_Decs, resampled_proj_x, resampled_proj_y], delimiter=",", header = header)

    return resampled_star_data

def doBootstrapOnMCMC(galaxy, model_types, halo_types, param_sets_to_measure,
                      n_resamples = 'full', w_replace = 1, pop_selection_method = 'metal_rigid', save_runs = 1, disk_MCMC_multiplier = 1,
                      fixed_params_set= [{}, {}, {}, {}] , save_res = 1, n_interim_iters_to_save = 500,
                      target_n_sky_pix = 500, mask_R_bins = 100, mask_z_bins = 100, max_sky_angle_deg = None, n_los_bins = 41, #41
                      n_randomizations = 1, n_MCMC_iterations = 3, n_MCMC_chains = 1, include_morph_prob = 1, include_kinem_prob = 1,
                      extra_iterator_save_number = 0, extra_save_str = '', prolate_or_oblate_set = ['either'], n_ignore_frac = 0.2,
                      curve_param_box_size = [3000.0, 1.0 * 10.0 ** 9.0], n_points_along_curve_to_find_max = 50,
                      smallest_max_val_to_be_on_curve = 50.0, show_results = 0, randomization_file_name = '', pop_independent_params = ['M', 'rs', 'halo_center'],
                      do_randomization = 1, n_single_param_bins = 30, pop_key_strs = ['_MP','_IM','_MR']
                      ):

    dSph_archive = DwarfGalDataArchive()
    pops = dSph_archive.getPopulations(galaxy)
    disk_MCMC_multiplier = int(disk_MCMC_multiplier)
    res_sets = []
    model_to_start_index = {'nfw_iso':0,'nfw_aniso':1,'burkert_iso':2,'burkert_aniso':3}

    params_to_curve_fit = ['rs','M']

    n_successful_computations = 0
    n_attempted_computations = 0

    MCMC_results_sets = []
    MCMC_results_file_names = []
    best_fit_vals_dict_set = []

    while n_successful_computations < n_randomizations:
        if do_randomization:
            resampled_stars = randomlySortStars('fornax', w_replace = w_replace, file_name = 'bootstrapSelectedStars_' + str(n_successful_computations))
            #resampled_stars = [star_set[0:50] for star_set in randomlySortStars('fornax', w_replace = w_replace, file_name = 'bootstrapSelectedStars_' + str(n_successful_computations))]
        else:
            resampled_stars = []
            for i in range(len(pops)):
                pop = pops[i]
                #new_stars = ObservedGalaxyStarData([galaxy,pop], pop_selection_method = pop_selection_method)
                new_stars = ObservedGalaxyStarData([galaxy,pop], pop_selection_method = pop_selection_method)
                resampled_stars = resampled_stars + [new_stars]

        n_successful_models = 0
        for j in range(len(model_types)):
            model_type = model_types[j]
            halo_type = halo_types[j]
            fixed_params = fixed_params_set[j]
            params_to_measure = param_sets_to_measure[j]
            central_param_storer = smcmc.StartMCMCParameterStorer( draw = 'best_fit' )
            ref_best_fit_params = central_param_storer.getStartParameters(model_to_start_index[model_type], 0)
            print ('Our reference best fit parameters are: ' + str(ref_best_fit_params))

            print ('Running MCMC for model ' + model_type + ' without a disk.')

            print ('n_los_bins = ' + str(n_los_bins))
            [runMCMCForDwarfGalaxyProfilesAllPopsSpherical(galaxy, nIterations = n_MCMC_iterations, pop_selection_method = pop_selection_method,
                                                           start_param_index = model_to_start_index[model_type], saverun = save_runs, star_data = resampled_stars,
                                                           start_param_draw_type = 'best_fit_random', fixed_params = fixed_params,
                                                           extra_save_str = extra_save_str + 'randomization' + str(n_successful_computations+extra_iterator_save_number) + '_chain' + str(chain) + '_' + model_type + '_',
                                                           extra_dir = 'MCMCResultsFromResampledStars/',
                                                           mask_R_bins = mask_R_bins, mask_z_bins = mask_z_bins, max_sky_angle_deg = max_sky_angle_deg, target_n_sky_pixels = target_n_sky_pix,
                                                           n_los_bins = n_los_bins, halo_type = halo_type,
                                                           include_morph_prob = include_morph_prob, include_kinem_prob = include_kinem_prob,
                                                           pop_independent_params = pop_independent_params, n_iters_to_save = n_interim_iters_to_save, )
                    for chain in range(n_MCMC_chains)]

            MCMC_files = ['MCMCResultsFromResampledStars/MCMC_out_' + extra_save_str + 'randomization' + str(n_successful_computations+extra_iterator_save_number) + '_chain' + str(chain) + '_' + model_type + '_' + halo_type + '_start_' + str(model_to_start_index[model_type] )
                              + '_fornax_simul_PSM_metal_rigidmask_0_N_'
                              + str(n_MCMC_iterations) + '.csv' for chain in range(n_MCMC_chains)]

            MCMC_results_file_names = MCMC_results_file_names + [MCMC_files]

            try:
                if show_results:
                    msmr.showMCMCResults(MCMC_files, ['M','rs','h_x_center','h_z_center', 'sigsqr_rr_0'], n_var_per_plot = 2, fig_size_unit = 1.0, params_to_fit = {'rs':'buffer','M':'buffer'}, smallest_max_val_to_be_on_curve = smallest_max_val_to_be_on_curve, n_fitted_points = n_points_along_curve_to_find_max)
                print ('Trying to measure best fit for files ' + str(MCMC_files))
                best_fit_vals_dict = {}
                for param_to_measure in params_to_measure:
                    full_fit_results_for_pop_dependent_param = msmr.measureMCMCResults(MCMC_files, [param_to_measure], ['single_hump_shift'], ['single_hump_shift'], n_bins_set = [n_single_param_bins], show_fit = show_results, n_ignore = int(n_MCMC_iterations * n_ignore_frac))
                    print ('param_to_measure = ' + str(param_to_measure))
                    print ('full_fit_results_for_pop_dependent_param = ' + str(full_fit_results_for_pop_dependent_param))
                    if len(full_fit_results_for_pop_dependent_param) == 1:
                        best_fit_vals_dict[param_to_measure] =full_fit_results_for_pop_dependent_param[0][1]
                    else:
                        for pop_index in range(len(full_fit_results_for_pop_dependent_param)):
                            best_fit_vals_dict[param_to_measure + pop_key_strs[pop_index]] = full_fit_results_for_pop_dependent_param[pop_index][1]
                #best_fit_vals_dict = {param_to_measure: msmr.measureMCMCResults(MCMC_files, [param_to_measure], ['single_hump_shift'], ['single_hump_shift'], n_bins_set = [15], show_fit = show_results)[0][0][1]
                #                      for param_to_measure in params_to_measure}
                best_fit_vals_dict_set = best_fit_vals_dict_set + [best_fit_vals_dict]
                print ('For model ' + str(model_type) + ', found the following best fit values: ' + str(best_fit_vals_dict))
                print ('Commencing curve fitting for parametes: ' + str(params_to_curve_fit))
                ref_best_fit_x_param = ref_best_fit_params[params_to_curve_fit[0]]
                ref_best_fit_y_param = ref_best_fit_params[params_to_curve_fit[1]]
                #print 'x_range_to_check = ' + str([ref_best_fit_x_param - curve_param_box_size[0]/2, ref_best_fit_x_param + curve_param_box_size[0]/2 ])
                #print 'y_range_to_check = ' + str([ref_best_fit_y_param - curve_param_box_size[1]/2, ref_best_fit_y_param + curve_param_box_size[1]/2 ])
                #{'h_x_center':-0.019, 'h_z_center':-0.065, 'sigsqr_rr_0_MP':169.6, 'sigsqr_rr_0_IM':132.8, 'sigsqr_rr_0_MR':112.5, 'gamma_for_beta_inf_MP':0.5, 'gamma_for_beta_inf_IM':0.1, 'gamma_for_beta_inf_IM':0.16, 'r_beta0_MP':2441, 'r_beta0_IM':2836, 'r_beta0_MR':1446}
                print ('best_fit_vals_dict = ' + str(best_fit_vals_dict))
                curve_fitting_results = cbfpacs.computeBestFitParamsAlongCurve(MCMC_files, 'fornax',
                                                                         best_fit_vals_dict, halo_type, n_los_bins, target_n_sky_pix, mask_R_bins,
                                                                         n_fitted_points_along_curve = n_points_along_curve_to_find_max, smallest_max_val_to_be_on_curve = smallest_max_val_to_be_on_curve, show = show_results,
                                                                         curve_params = params_to_curve_fit, include_morph_prob = include_morph_prob, include_kinem_prob = include_kinem_prob,)
                if len(curve_fitting_results) == 0:
                    print ('Failed to identify enough significant points along the curve for M-rs curve fitting.  We will discard this attempt. ')
                    print ('Have completed ' + str(n_successful_computations) + ' of ' + str(n_randomizations) + ' needed for completion. ')
                    break
                #print 'curve_fitting_results = ' + str(curve_fitting_results)
                #print 'curve_fitting_results[-1] = ' + str(curve_fitting_results[-1])
                #print " curve_fitting_results[-1]['x'] = " + str( curve_fitting_results[-1]['x'])

                best_fit_vals_dict[params_to_curve_fit[0]] = curve_fitting_results[-1]['x'][0]
                best_fit_vals_dict[params_to_curve_fit[1]] = curve_fitting_results[-1]['x'][1]
                MCMC_results_sets = MCMC_results_sets + [model_type, best_fit_vals_dict, -1.0 * curve_fitting_results[-1]['fun']]
                print ('MCMC_results_sets now: ' + str(MCMC_results_sets))
                n_successful_models = n_successful_models + 1
            except RuntimeError:
                print ('On attempt ' + str(n_attempted_computations))
                print ('Measuring MCMC results failed for at least one parameter for model ' + model_type + '.  This randomization will be reattempted. ')
                print ('Have completed ' + str(n_successful_computations) + ' of ' + str(n_randomizations) + ' needed for completion. ')
                break

        print ('MCMC_results_sets now: ' + str(MCMC_results_sets ))
        #if len(best_fit_vals_dict_set) == len(model_types):
        #    print ('Measuring MCMC results succeeded for all parameters for all models.  This randomization will be used. ')

        #    MCMC_WD_results_sets = []
        #    for j in range(len(model_types)):
        #        doWithDisk = doWithDisksSet[j]
        #        model_type = model_types[j]
        #        disk_funct = disk_functs[j]
        #        halo_type = halo_types[j]
        #        fixed_params = fixed_params_set[j]

        #            model_type = model_types[j]
        #            halo_type = halo_types[j]
        #            params_to_measure = param_sets_to_measure[j]
        #        print ('For model ' + model_types[j] + ' doWithDisk = ' + str(doWithDisk))

        #    n_successful_computations = n_successful_computations + 1
        #    print ('Have completed ' + str(n_successful_computations) + ' of ' + str(n_randomizations) + ' needed for completion. ')
        print ('[n_successful_models, len(model_types)] = ' + str([n_successful_models, len(model_types)]))
        if n_successful_models == len(model_types):
            print ('Successfully finished attempt ' + str(n_attempted_computations) + ' (total target is ' + str(n_randomizations) + ')')
            n_successful_computations = n_successful_computations + 1
        else:
            print ('Failed somewhere on attempt ' + str(n_attempted_computations) + ' (total target is ' + str(n_randomizations) + ')')
        print ('[n_successful_models, n_successful_computations, len(model_types)] = ' + str([n_successful_models, n_successful_computations, len(model_types)]))
        n_attempted_computations = n_attempted_computations + 1
        print ('We have attempted ' + str(n_attempted_computations) )

    if save_res:
        computational_params_archive = ComputationalArchive()
        if len(randomization_file_name) == 0:
            randomization_file_name = 'bootstrap_randomization_' + str(extra_iterator_save_number) + '_' + extra_save_str
        np.save(computational_params_archive.getRandomizationDir() + randomization_file_name, MCMC_results_sets )

    return MCMC_results_sets


def doBootstrapSampling(galaxy,
                        dist1, M1, vphi1, el1, rs1, phi1, theta1, halo_type1, halo_sym_axis1, halo_center1, c1, lam1, zeta1, eps1, a1, b1, disk_type1, disk_sym_axis1, disk_center1,
                        dist2, M2, vhpi2, el2, rs2, phi2, theta2, halo_type2, halo_sym_axis2, halo_center2, c2, lam2, zeta2, eps2, a2, b2, disk_type2, disk_sym_axis2, disk_center2,
                        n_resamples = 'full', w_replace = 1, pop_selection_method = 'metal_rigid',
                        outer_step = 2.0, arcmin_limits_los = [-60.0, 60.0], n_randomizations = 1):

    astro_archive = AstronomicalParameterArchive()
    gamma = astro_archive.getGamma()
    dSph_archive = DwarfGalDataArchive()
    pops = dSph_archive.getPopulations(galaxy)
    popGalPairs = []
    for i in range(len(pops)):
        pop = pops[i]
        popGalPairs = popGalPairs + [[galaxy,pop]]

    arcmin_limits_R,arcmin_limits_z = dSph_archive.getObservationBounds([galaxy,'dummy_var'],return_unit = 'arcmin')

    R1, z1, los_bins1 = compRZLos(1.0, 1.0 * outer_step, outer_step, arcmin_limits_R, arcmin_limits_z, arcmin_limits_los, rs1, dist1)
    R2, z2, los_bins2 = compRZLos(1.0, 1.0 * outer_step, outer_step, arcmin_limits_R, arcmin_limits_z, arcmin_limits_los, rs2, dist2)
    Rmesh_degrees1, zmesh_degrees1 = np.meshgrid(R1 * (rs1 / dist1) * (180.0 / math.pi), z1 * (rs1 / dist1) * (180.0 / math.pi))
    Rmesh_degrees2, zmesh_degrees2 = np.meshgrid(R2 * (rs2 / dist2) * (180.0 / math.pi), z2 * (rs2 / dist2) * (180.0 / math.pi))

    observation_mask1 = GalaxyMask(Rmesh_degrees1, zmesh_degrees1, galaxy, mask_types = ['n_vel_meas'])
    observation_mask2 = GalaxyMask(Rmesh_degrees2, zmesh_degrees2, galaxy, mask_types = ['n_vel_meas'])

    halo_funct1 = readInPotentialInterpolator(halo_type1)
    disk_funct1 = readInPotentialInterpolator(disk_type1)
    halo_funct2 = readInPotentialInterpolator(halo_type2)
    disk_funct2 = readInPotentialInterpolator(disk_type2)

    rs = []
    for rand_index in range(n_randomizations):
        start = time.time()

        resampled_star_data = randomlySortStars(galaxy, w_replace = w_replace, pop_selection_method = pop_selection_method, n_resamples = 'full')
        resampled_proj_x = [resampled_star_data[i].proj_x for i in range(len(resampled_star_data))]
        resampled_proj_y = [resampled_star_data[i].proj_y for i in range(len(resampled_star_data))]
        resampled_RAs = [resampled_star_data[i].corrRA for i in range(len(resampled_star_data)) ]
        resampled_Decs = [resampled_star_data[i].corrDec for i in range(len(resampled_star_data)) ]
        resampled_corrVHels = [resampled_star_data[i].corr_Vhel for i in range(len(resampled_star_data)) ]
        resampled_sigsqrs = [resampled_star_data[i].sigSqr for i in range(len(resampled_star_data)) ]
        #resampled_RAs, resampled_Decs, resampled_corrVHels = randomlySortStars(galaxy, w_replace = w_replace, pop_selection_method = pop_selection_method, n_resamples = 'full')
        log_likelihood1 = 0.0
        log_likelihood2 = 0.0

        for i in range(len(pops)):
            population = popGalPairs[i]
            random_select_Vhels = resampled_corrVHels[i]
            mean_Vhel = sum(random_select_Vhels) / float(len(random_select_Vhels))
            sigsqr = sum([(vel - mean_Vhel) ** 2.0 for vel in random_select_Vhels]) / float(len(random_select_Vhels))
            print ('sigsqr = ' + str(sigsqr))
            sigsqr = resampled_sigsqrs[i]
            print ('sigsqr = ' + str(sigsqr))
            storer1 = DwarfGalaxyParametersStorer(population, 0, sigsqr = sigsqr,
                                                 el = el1, M = M1, rs = rs1, phi = phi1, theta = theta1, halo_sym_axis = halo_sym_axis1, c = c1, halo_center = halo_center1,
                                                 lam = lam1, zeta = zeta1, eps = eps1, a = a1, b = b1, disk_sym_axis = disk_sym_axis1, disk_center = disk_center1)
            storer2 = DwarfGalaxyParametersStorer(population, 0, sigsqr = sigsqr,
                                                 el = el2, M = M2, rs = rs2, phi = phi2, theta = theta2, halo_sym_axis = halo_sym_axis2, c = c2, halo_center = halo_center2,
                                                 lam = lam2, zeta = zeta2, eps = eps2, a = a2, b = b2, disk_sym_axis = disk_sym_axis2, disk_center = disk_center2)


            surface_prof1 = SurfaceBrightnessProfile(R1, z1, los_bins1, gamma, storer1,
                                                     disk_interpolating_function = disk_funct1,
                                                     halo_interpolating_function = halo_funct1,
                                                     observation_mask = observation_mask1.final_mask)
            surface_prof2 = SurfaceBrightnessProfile(R2, z2, los_bins2, gamma, storer2,
                                                     disk_interpolating_function = disk_funct2,
                                                     halo_interpolating_function = halo_funct2,
                                                     observation_mask = observation_mask2.final_mask)

            log_likelihood1 = log_likelihood1 + surface_prof1.sumLogSurfaceBrightness( np.array(resampled_proj_x[i]), np.array(resampled_proj_y[i]) )

            log_likelihood2 = log_likelihood2 + surface_prof2.sumLogSurfaceBrightness( np.array(resampled_proj_x[i]), np.array(resampled_proj_y[i]) )

        print ('log_likelihood1 = ' + str(log_likelihood1))
        print ('log_likelihood2 = ' + str(log_likelihood2))
        r = math.e ** (log_likelihood1 - log_likelihood2)
        rs = rs + [r]
        end = time.time()
        print ('Iteration ' + str(rand_index+1) + ' of ' + str(n_randomizations) + ' took ' + str(end - start) + 's.')

    return rs
