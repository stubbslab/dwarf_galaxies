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
from runMCMCForDwarfGalaxyProfilesAllPops import runMCMCForDwarfGalaxyProfilesAllPops 
from cantrips import removeListElement
from GalaxyMask import GalaxyMask
from manageMCMCResults import measureMCMCResults
from StartMCMCParameterStorer import StartMCMCParameterStorer
from manageMCMCResults import showMCMCResults 
import time
from computeBestFitParamsAlongCurve import computeBestFitParamsAlongCurve
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
        print 'sigsqr = ' + str(sigsqr) 
        resampled_star_data[i].sigSqr = sigsqr 
    if save:
        computational_params_archive = ComputationalArchive() 
        header = 'corrRA, corrDec, corr_Vhel, proj_x, proj_y'
        save_dir = computational_params_archive.getMCMCOutDir() + 'stellarRandomizations/'
        file_name = file_name + '.csv'
        np.savetxt(save_dir + file_name, [resampled_RAs, resampled_Decs, resampled_proj_x, resampled_proj_y], delimiter=",", header = header)

    return resampled_star_data 

def doBootstrapOnMCMC(galaxy, model_types, halo_functs, disk_functs, doWithDisksSet,
                      n_resamples = 'full', w_replace = 1, pop_selection_method = 'metal_rigid', save_runs = 1, disk_MCMC_multiplier = 1,
                      specified_param_ranges_set = [{'a':0.0}], save_res = 1, 
                      outer_step = 2.0, arcmin_limits_los = [-50.0, 50.0], n_randomizations = 1, n_MCMC_iterations = 3, n_MCMC_chains = 1, 
                      extra_iterator_save_number = 0, extra_save_str = '', prolate_or_oblate_set = ['either'], 
                      curve_param_box_size = [3000.0, 1.0 * 10.0 ** 9.0], n_points_along_curve_to_find_max = 50,
                      smallest_max_val_to_be_on_curve = 50.0, show_results = 0, randomization_file_name = ''):
    
    disk_MCMC_multiplier = int(disk_MCMC_multiplier) 
    res_sets = []
    model_to_start_index = {'nfw_obl':0, 'nfw_pro':1, 'cored_obl':2, 'cored_pro':3, 'burkert_obl':4, 'burkert_pro':5}
    params_to_measure = ['el', 'phi', 'h_x_center', 'h_z_center']
    params_to_curve_fit = ['rs','M']
    
    n_successful_computations = 0
    n_attempted_computations = 0
    while n_successful_computations < n_randomizations:
        resampled_stars = randomlySortStars('fornax', w_replace = w_replace, file_name = 'bootstrapSelectedStars_' + str(n_successful_computations))
        MCMC_WOD_results_sets = []
        MCMC_results_file_names = []
        best_fit_vals_dict_set = []
        
        for j in range(len(model_types)):
            model_type = model_types[j]
            disk_funct = disk_functs[j]
            halo_funct = halo_functs[j]
            if len(prolate_or_oblate_set) == 1:
                prolate_or_oblate = prolate_or_oblate_set[0]
            else:
                prolate_or_oblate = prolate_or_oblate_set[j]
            if len(specified_param_ranges_set) == 1:
                specified_param_ranges = specified_param_ranges_set[0]
            else:
                specified_param_ranges = specified_param_ranges_set[j]
            central_param_storer = StartMCMCParameterStorer(draw = 'best_fit' )
            ref_best_fit_params = central_param_storer.getStartParameters(model_to_start_index[model_type], 0)
            print 'Our reference best fit parameters are: ' + str(ref_best_fit_params)
            
            print 'Running MCMC for model ' + model_type + ' without a disk.'
            
            [runMCMCForDwarfGalaxyProfilesAllPops(galaxy, withDisk = 0, disk_funct = disk_funct, halo_funct = halo_funct,
                                                  nIterations = n_MCMC_iterations, pop_selection_method = pop_selection_method,
                                                  start_param_index = model_to_start_index[model_type], halo_edge_on = 1,
                                                  disk_edge_on = 1, arcmin_limits_los = arcmin_limits_los,
                                                  saverun = save_runs, outer_step = outer_step, star_data = resampled_stars,
                                                  start_param_draw_type = 'best_fit_random',
                                                  specified_param_ranges = specified_param_ranges, 
                                                  extra_save_str = extra_save_str + 'randomization' + str(n_successful_computations+extra_iterator_save_number) + '_chain' + str(chain),
                                                  extra_dir = 'MCMCResultsFromResampledStars/', prolate_or_oblate = prolate_or_oblate)
             for chain in range(n_MCMC_chains)]
            WOD_MCMC_files = ['MCMCResultsFromResampledStars/MCMC_output_' + extra_save_str + 'randomization' + str(n_successful_computations+extra_iterator_save_number) + '_chain' + str(chain) +'__WD_0_' + prolate_or_oblate + '_start_'
                              + str(model_to_start_index[model_type]) + '_fornax_simul_PSM_metal_rigid_HEdge_1_DEdge_1_fixDAng_None_fixHAng_None_crowdMsk_0_N_'
                              + str(n_MCMC_iterations) + '.csv' for chain in range(n_MCMC_chains)]
            
            MCMC_results_file_names = MCMC_results_file_names + [WOD_MCMC_files]

            try:
                if show_results: 
                    showMCMCResults(WOD_MCMC_files, ['rs','M','el','phi','h_x_center','h_z_center'], n_var_per_plot = 2, fig_size_unit = 2.5, params_to_fit = {'rs':'buffer','M':'buffer'}, smallest_max_val_to_be_on_curve = smallest_max_val_to_be_on_curve, n_fitted_points = n_points_along_curve_to_find_max) 
                print 'Trying to measure best fit for files ' + str(WOD_MCMC_files) 
                best_fit_vals_dict = {param_to_measure: measureMCMCResults(WOD_MCMC_files, [param_to_measure], ['single_hump'], ['single_hump'], n_bins_set = [15], show_fit = show_results)[0][0][1]
                                      for param_to_measure in params_to_measure}
                best_fit_vals_dict_set = best_fit_vals_dict_set + [best_fit_vals_dict]
                print 'For model ' + str(model_type) + ', found the following best fit values: ' + str(best_fit_vals_dict)
                print 'Commencing curve fitting for parametes: ' + str(params_to_curve_fit)
                ref_best_fit_x_param = ref_best_fit_params[params_to_curve_fit[0]]
                ref_best_fit_y_param = ref_best_fit_params[params_to_curve_fit[1]]
                #print 'x_range_to_check = ' + str([ref_best_fit_x_param - curve_param_box_size[0]/2, ref_best_fit_x_param + curve_param_box_size[0]/2 ])
                #print 'y_range_to_check = ' + str([ref_best_fit_y_param - curve_param_box_size[1]/2, ref_best_fit_y_param + curve_param_box_size[1]/2 ]) 
                curve_fitting_results = computeBestFitParamsAlongCurve(WOD_MCMC_files, 'fornax', 0, best_fit_vals_dict, model_type, 'sech_disk',
                                                                       halo_funct = halo_funct, disk_funct = disk_funct, curve_params = params_to_curve_fit, 
                                                                       outer_step = outer_step, n_fitted_points_along_curve = n_points_along_curve_to_find_max,
                                                                       smallest_max_val_to_be_on_curve = smallest_max_val_to_be_on_curve, show = show_results
                                                                       #x_range_to_check = [], y_range_to_check = 'buffer'
                                                                       #x_range_to_check = [ref_best_fit_x_param - curve_param_box_size[0]/2, ref_best_fit_x_param + curve_param_box_size[0]/2 ],
                                                                       #y_range_to_check = [ref_best_fit_y_param - curve_param_box_size[1]/2, ref_best_fit_y_param + curve_param_box_size[1]/2 ]
                                                                       )
                #print 'curve_fitting_results = ' + str(curve_fitting_results)
                #print 'curve_fitting_results[-1] = ' + str(curve_fitting_results[-1])
                #print " curve_fitting_results[-1]['x'] = " + str( curve_fitting_results[-1]['x'])
                
                best_fit_vals_dict[params_to_curve_fit[0]] = curve_fitting_results[-1]['x'][0]
                best_fit_vals_dict[params_to_curve_fit[1]] = curve_fitting_results[-1]['x'][1]
                MCMC_WOD_results_sets = MCMC_WOD_results_sets + [model_type, best_fit_vals_dict, -1.0 * curve_fitting_results[-1]['fun']]
            except RuntimeError:
                print 'On attempt ' + str(n_attempted_computations)
                print 'Measuring MCMC results failed for at least one parameter for model ' + model_type + '.  This randomization will be reattempted. '
                print 'Have completed ' + str(n_successful_computations) + ' of ' + str(n_randomizations) + ' needed for completion. '
                break
        print 'On attempt ' + str(n_attempted_computations)
        if len(best_fit_vals_dict_set) == len(model_types):
            print 'Measuring MCMC results succeeded for all parameters for all models.  This randomization will be used. '
            
            MCMC_WD_results_sets = []
            for j in range(len(model_types)): 
                doWithDisk = doWithDisksSet[j]
                model_type = model_types[j]
                disk_funct = disk_functs[j]
                halo_funct = halo_functs[j] 
                print 'For model ' + model_types[j] + ' doWithDisk = ' + str(doWithDisk) 
                if doWithDisk:
                    print 'Working on MCMC with variable disk params.'
                    WOD_MCMC_file = MCMC_results_file_names[j]
                    fixed_params_from_WOD = best_fit_vals_dict_set[j]
                    if ('h_x_center' in params_to_measure) or ('h_z_center' in params_to_measure):
                        fixed_params_from_WOD['halo_center'] = [fixed_params_from_WOD['h_x_center'], 0.0, fixed_params_from_WOD['h_z_center']]
                        fixed_params_from_WOD.pop('h_x_center', None)
                        fixed_params_from_WOD.pop('h_z_center', None)
                    if ('phi' in params_to_measure) or ('theta' in params_to_measure):
                        phi = fixed_params_from_WOD['phi']
                        if 'theta' in fixed_params_from_WOD.keys():
                            theta = fixed_params_from_WOD['theta']
                        else:
                            theta = 0.0 
                        fixed_params_from_WOD['halo_sym_axis'] = [np.sin(phi) * np.cos(theta), np.sin(phi) * np.sin(theta), np.cos(phi)]
                    fixed_params_from_WOD['a'] = 0.0
                    
                    print 'Running MCMC for model ' + model_type + ' with a disk.' 
                    MCMC_WD_results_sets = MCMC_WD_results_sets + [runMCMCForDwarfGalaxyProfilesAllPops(galaxy, withDisk = 1, disk_funct = disk_funct, halo_funct = halo_funct,
                                                                       nIterations = n_MCMC_iterations*disk_MCMC_multiplier, pop_selection_method = pop_selection_method,
                                                                       start_param_index = model_to_start_index[model_type], halo_edge_on = 1, disk_edge_on = 1,
                                                                       fixed_params = fixed_params_from_WOD, arcmin_limits_los = arcmin_limits_los, saverun = save_runs, 
                                                                       outer_step = outer_step, star_data = resampled_stars, start_param_draw_type = 'best_fit',
                                                                       extra_save_str = 'randomization' + str(n_successful_computations+extra_iterator_save_number),
                                                                       extra_dir = 'MCMCResultsFromResampledStars/')]
            
            n_successful_computations = n_successful_computations + 1
            print 'Have completed ' + str(n_successful_computations) + ' of ' + str(n_randomizations) + ' needed for completion. '
                
                

            res_sets = res_sets + [MCMC_WOD_results_sets + MCMC_WD_results_sets]

    if save_res:
        computational_params_archive = ComputationalArchive()
        if len(randomization_file_name) == 0:
            randomization_file_name = 'bootstrap_randomization_' + str(extra_iterator_save_number) + '_' + extra_save_str 
        np.save(computational_params_archive.getRandomizationDir() + randomization_file_name, res_sets)

    return res_sets 
    

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
            print 'sigsqr = ' + str(sigsqr)
            sigsqr = resampled_sigsqrs[i]
            print 'sigsqr = ' + str(sigsqr) 
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
        
        print 'log_likelihood1 = ' + str(log_likelihood1)
        print 'log_likelihood2 = ' + str(log_likelihood2)
        r = math.e ** (log_likelihood1 - log_likelihood2)
        rs = rs + [r]
        end = time.time()
        print 'Iteration ' + str(rand_index+1) + ' of ' + str(n_randomizations) + ' took ' + str(end - start) + 's.'

    return rs
