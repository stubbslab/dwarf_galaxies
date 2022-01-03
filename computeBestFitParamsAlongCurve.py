

#Do a sequence of the MCMC Series algorithms, with varying run conditions.
#You can change various parameters (present set listed in commented line right above parameter_serires assignment), and have program execute in sequence
#This script can be run from command line without having to be in python environment (bash$ python doMCMCSeriesAllPops.py)
#This method tries to be as efficient as possible, reading in the various potentialFunctionArrays only once, as that can be a time consuming part of the process.

from runMCMCForDwarfGalaxyProfilesAllPops import runMCMCForDwarfGalaxyProfilesAllPops
from PotentialFunctionArray import PotentialFunctionArray
from PotentialArchive import PotentialArchive
from ComputationalArchive import ComputationalArchive
from SurfaceBrightnessProfile import SurfaceBrightnessProfile
from DwarfGalaxyParametersStorer import DwarfGalaxyParametersStorer
from AstronomicalParameterArchive import AstronomicalParameterArchive
from DwarfGalDataArchive import DwarfGalDataArchive
from ObservedGalaxyStarData import ObservedGalaxyStarData
from readInPotentialInterpolator import readInPotentialInterpolator
#from BoylanArtificialGalaxyStarData import BoylanArtificialGalaxyStarData
from compRZLos import compRZLos
from GalaxyMask import GalaxyMask
import scipy.optimize as optimize
import manageEllipticalMCMCResults as memr
from cantrips import safeSortOneListByAnother
import math
import SharedObjectHolder as soh
import scipy.interpolate as interpolate
import numpy as np
import matplotlib.pyplot as plt
from manageMCMCResults import getParamDisplayVals
from manageMCMCResults import getVarDisplayScalings

def computeLikelihoodGivenSeedAndFixedParams(galaxy, withDisk, seed_params, seed_x, seed_y,
                                             param_initializer, pops, pop_selection_method, apply_observation_mask, halo_funct, disk_funct, halo_type, disk_type,
                                             #dist, arcmin_limits_R, arcmin_limits_z, arcmin_limits_los, gamma, outer_step,
                                             dist, gamma, n_los_bins, target_n_sky_pixels, mask_R_bins,
                                             fixed_hx_rs_val, fixed_hz_rs_val, fixed_dx_rs_val, fixed_dz_rs_val, star_data_sets,
                                             pop_connector = '_', include_kinem_prob = 1, include_morph_prob = 1):

    astro_arch = AstronomicalParameterArchive()
    deg_to_rad = astro_arch.getDegToRad()
    dSph_archive = DwarfGalDataArchive(pop_selection_method = pop_selection_method)
    deg_limits_R,deg_limits_z = dSph_archive.getObservationBounds([galaxy,'dummy_pop_var'], return_unit = 'degrees')
    max_sky_angle_deg = np.max(np.abs([deg_limits_R ,deg_limits_z]))
    mask_R = np.linspace(-max_sky_angle_deg, max_sky_angle_deg, mask_R_bins)
    mask_z = mask_R[:]
    print ('n_los_bins = ' + str(n_los_bins))
    shared_object_holder = soh.SharedObjects(dist, mask_R, mask_z, np.linspace(-1.0, 1.0, n_los_bins), target_n_sky_pixels, max_sky_angle_deg * deg_to_rad, compute_mask = apply_observation_mask)
    Gyr_to_s = 365.2425 * 24 * 3600 * 10.0 ** 9.0

    param_initializer[seed_params[0]] = seed_x
    param_initializer[seed_params[1]] = seed_y
    param_initializer['halo_center'] = [fixed_hx_rs_val, 0.0, fixed_hz_rs_val]
    param_initializer['disk_center'] = [fixed_dx_rs_val, 0.0, fixed_dz_rs_val]

    #if seed_params[0] in ['rs']:
    #    rs_seed = seed_x
    #elif seed_params[1] in ['rs']:
    #    rs_seed = seed_y
    #else:
    #    rs_seed = None
    #if not rs_seed is None:
    #    current_h_center = param_initializer['halo_center']
    #    new_h_center = current_h_center[:]
    #    if 'halo_x_center' in param_initializer:
    #        new_h_center[0] = param_initializer['halo_x_center']
    #    elif 'h_x_center' in param_initializer:
    #        new_h_center[0] = param_initializer['h_x_center']
    #    elif not(fixed_hx_rs_val is None):
    #        new_h_center[0] = fixed_hx_rs_val
    #    if 'halo_z_center' in param_initializer:
    #        new_h_center[2] = param_initializer['halo_z_center']
    #    elif 'h_z_center' in param_initializer:
    #        new_h_center[2] = param_initializer['h_z_center']
    #    elif  not(fixed_hz_rs_val is None):
    #        new_h_center[2] = fixed_hz_rs_val
    #    param_initializer['halo_center'] = new_h_center

    #if not rs_seed is None:
    #    current_d_center = param_initializer['disk_center']
    #    new_d_center = current_d_center[:]
    #    if 'disk_x_center' in param_initializer:
    #        new_d_center[0] = param_initializer['disk_x_center']
    #    elif 'h_x_center' in param_initializer:
    #        new_d_center[0] = param_initializer['d_x_center']
    #    elif not(fixed_dx_rs_val is None):
    #        new_d_center[0] = fixed_dx_rs_val
    #    if 'disk_z_center' in param_initializer:
    #        new_d_center[2] = param_initializer['disk_z_center']
    #    elif 'd_z_center' in param_initializer:
    #        new_d_center[2] = param_initializer['d_z_center']
    #    elif  not(fixed_dz_rs_val is None):
    #        new_d_center[2] = fixed_dz_rs_val
    #    param_initializer['disk_center'] = new_d_center

    parameter_storers = []
    log_likelihood = 0.0

    for i in range(len(pops)):
        pop = pops[i]
        star_data = star_data_sets[i]
        population = [galaxy, pop]
        storer =  DwarfGalaxyParametersStorer(population, withDisk = withDisk, sigsqr_RR = param_initializer['sigsqr_RR' + pop_connector + pop],
                                              omega_phi = param_initializer['omega_phi' + pop_connector + pop]  ,
                                              el = param_initializer['el'], M = param_initializer['M'], rs = param_initializer['rs'],
                                              phi = param_initializer['phi'], theta = param_initializer['theta'],
                                              halo_sym_axis = param_initializer['halo_sym_axis'], halo_center = [param_initializer['h_x_center'], 0.0, param_initializer['h_z_center']],
                                              lam = param_initializer['lam'], Rd = param_initializer['Rd'], eps = param_initializer['eps'],
                                              a = param_initializer['a'], b = param_initializer['b'], #stellar_data = star_data,
                                              disk_sym_axis = param_initializer['disk_sym_axis'], disk_center = [param_initializer['d_x_center'], 0.0, param_initializer['d_z_center']], stellar_data = star_data,)

        parameter_storers = parameter_storers + [storer]
        surface_profile = SurfaceBrightnessProfile(storer, shared_object_holder, halo_type, disk_type,
                                                   disk_interpolating_function = disk_funct, halo_interpolating_function = halo_funct,
                                                   include_morph_prob = include_morph_prob, include_kinem_prob = include_kinem_prob)
        single_pop_log_likelihood = surface_profile.getStarProbabilityValues(star_data.proj_x, star_data.proj_y, star_data.corr_Vhel)
        #print ('single_pop_log_likelihood = ' + str(single_pop_log_likelihood) )
        log_likelihood = log_likelihood + single_pop_log_likelihood

    print ('Log likelihood = ' + str(log_likelihood) + ' with [seed_x, seed_y] = ' + str([seed_x, seed_y]))
    return log_likelihood


def computeBestFitParamsAlongCurve(results_files, galaxy, withDisk, fixed_params, halo_type, disk_type, n_los_bins, target_n_sky_pixels, mask_R_bins, #galaxy, withDisk, fixed_params, halo_type, disk_type,
                                   results_dir = '/Users/sashabrownsberger/Documents/Harvard/physics/randall/MCMCOutputs/', halo_funct = None, disk_funct = None,
                                   include_kinem_prob = 1, include_morph_prob = 1,
                                   pops = ['MP','IM','MR'], pop_connector = '_', pop_selection_method = 'metal_rigid', fit_funct_type = 'spline', data_type = 'real', apply_observation_mask = 1, dist = None,
                                   curve_params = ['rs','M'], n_xs_to_check = 10, x_range_to_check = [791.0 * 0.1, 791.0 * 10.0], y_range_to_check = [10.0 ** 7.39, 600 * 10.0 ** 7.39], max_refine_steps = 6, max_final_steps = 100, outer_step = 1.0, show = 1, apply_rotation = 1,
                                   x_min_width = 100.0, y_min_width = 5.0 * 10.0 ** 7.0, fixed_hx_rs_val = None, fixed_hz_rs_val = None, fixed_dx_rs_val = None, fixed_dz_rs_val = None, n_ignore = 0, fit_lines_slope = 0.0, go_beyond_MCMC = 0,
                                   n_fitted_points_along_curve = 200, n_xy_points_to_find_curve = 1000, n_points_past_max_to_break = 10, smallest_max_val_to_be_on_curve = 50.0, fixed_slope = None):

    dSph_archive = DwarfGalDataArchive()
    astro_archive = AstronomicalParameterArchive ()
    if len(pops) == 0:
        pops = dSph_archive.getPopulations(galaxy)
    if dist is None:
        dist = dSph_archive.getDistanceFromSun([galaxy,'dummy_pop_variable'])
    gamma = astro_archive.getGamma()

    param_display_vals = memr.getParamDisplayVals(pops = pops)
    astro_archive = AstronomicalParameterArchive ()
    var_display_scalings = memr.getVarDisplayScalings(pops = pops)
    x_scaling = var_display_scalings[curve_params[0]]
    y_scaling = var_display_scalings[curve_params[1]]
    #compute_params_archive = ComputationalArchive()

    if halo_funct is None:
        halo_funct = readInPotentialInterpolator(halo_type)
    if disk_funct is None:
        disk_funct = readInPotentialInterpolator(disk_type)

    if 'boy' in galaxy:
        arcmin_limits_R = [-60.0, 60.0]
        arcmin_limits_z = [-60.0, 60.0]
        apply_observation_mask = 0
    else:
        arcmin_limits_R,arcmin_limits_z = dSph_archive.getObservationBounds([galaxy,'dummy_var'],return_unit = 'arcmin')
    arcmin_limits_los = [-50.0,50.0]

    default_params = {'halo_sym_axis':None, 'h_y_center':0.0, 'theta':0.0,
                      'eps':0.0, 'lam':0.1, 'Rd':1000.0, 'a':0.0, 'b':0.0, 'disk_sym_axis':None,
                      'd_x_center':0.0, 'd_y_center':0.0, 'd_z_center':0.0, 'sigsqr_RR':100.0, 'omega_phi':0.0  }

    pop_dependent_vars = memr.getVarPopDependence()
    init_param_keys = list(default_params.keys())
    for var in init_param_keys:
        if var in pop_dependent_vars:
            for pop in pops:
                default_params[var + pop_connector + pop] = default_params[var]
            default_params.pop(var)

    #if fit_funct_type.lower() in ['poly','polynomial','p']:
    #    fit_funct = np.poly1d(fit_params)
    #elif fit_funct_type.lower() in ['spline', 'splines', 'sp', 's']:
    print ('param_ranges_to_fit = ' + str([x_range_to_check, y_range_to_check]))
    fit_xs, fit_ys = memr.fitMCMCResults(results_files, curve_params, measured_arrays = None, n_visits_array = None, param_ranges_to_fit = [x_range_to_check, y_range_to_check], withDisk = withDisk,
                                         results_dir = results_dir, n_ignore = n_ignore, theta_shift = None, n_fitted_points = n_fitted_points_along_curve, n_xy_points_to_fit = n_xy_points_to_find_curve, smallest_max_val_for_fit = smallest_max_val_to_be_on_curve, fixed_slope = fixed_slope)
    if len(fit_xs) != len(fit_ys) or len(fit_xs) < 1:
        print ('Not enough well located points for fitting the curve.  Returning an empty array. ')
        return []
    fit_funct = interpolate.interp1d( *safeSortOneListByAnother(fit_xs, [fit_xs, fit_ys]), kind = 'linear')
    #print 'fit_xs = ' + str(fit_xs)
    #print 'fit_ys = ' + str(fit_ys)
    if x_range_to_check in [None, 'buffer', 'all']:
        seed_xs = np.linspace(min(fit_xs), max(fit_xs), n_xs_to_check)
    else:
        seed_xs = np.linspace(max(x_range_to_check[0], min(fit_xs)), min(x_range_to_check[1], max(fit_xs)), n_xs_to_check)
    min_xs = []
    min_ys = []
    seed_ys = [fit_funct(x) for x in seed_xs]
    used_xs = []
    used_ys = []
    log_probs_along_line = []


    #read in default parameters and then overwrite with whichever ones the user gives
    param_initializer = default_params.copy()
    for param in fixed_params.keys():
        param_initializer[param] = fixed_params[param]
    if not withDisk:
        param_initializer['Rd'] = 1000.0
        param_initializer['eps'] = 0.0

    star_data_sets = []
    for i in range(len(pops)):
        pop = pops[i]
        population = [galaxy, pop]
        if 'boy' in population[0]:
            viewer_phi = 0.0
            viewer_theta = 0.0
            halo_number = population[0][-1]
            print ('halo_number = ' + str(halo_number) )
            star_data = BoylanArtificialGalaxyStarData(halo_number, viewer_phi = viewer_phi, viewer_theta = viewer_theta,
                                                       arcmin_limits_R = arcmin_limits_R, arcmin_limits_z = arcmin_limits_z)
        else:
            star_data = ObservedGalaxyStarData(population, pop_selection_method = pop_selection_method)
        star_data_sets = star_data_sets + [star_data]

    vals_along_curve = []
    n_points_past_max = 0
    max_log_likelihood = 0.0
    max_found = 0
    n_measured_points = 0

    while not(max_found) and ((n_measured_points < len(seed_xs) ) or go_beyond_MCMC):
        i = n_measured_points
        if i < len(seed_xs):
            seed_x = seed_xs[i]
            seed_y = seed_ys[i]
        else:
            prev_length_for_averaging = 5
            prev_x = min_xs[-1]
            avg_recent_x_step = np.median(min_xs[-1 * prev_length_for_averaging:-1]) / prev_length_for_averaging
            seed_x = prev_x + avg_recent_x_step
            prev_y = min_ys[-1]
            avg_recent_y_step = np.median(min_ys[-1 * prev_length_for_averaging:-1]) / prev_length_for_averaging
            seed_y = prev_y + avg_recent_y_step
            print ('We are extrapolating to the next point. ')
        param_initializer[curve_params[0]] = seed_x
        param_initializer[curve_params[1]] = seed_y
        print ('[seed_x, seed_y] = ' + str([seed_x, seed_y]))
        print ('param_initializer = ' + str( param_initializer))
        start_likelihood = computeLikelihoodGivenSeedAndFixedParams(galaxy, withDisk, curve_params, seed_x, seed_y,
                                                                    param_initializer, pops, pop_selection_method, apply_observation_mask, halo_funct, disk_funct, halo_type, disk_type,
                                                                    dist, gamma, n_los_bins, target_n_sky_pixels, mask_R_bins,
                                                                    fixed_hx_rs_val, fixed_hz_rs_val, fixed_dx_rs_val, fixed_dz_rs_val, star_data_sets, pop_connector = pop_connector,
                                                                    include_kinem_prob = include_kinem_prob, include_morph_prob = include_morph_prob )
        print ('Working on computation ' + str(i+1) + ' of ' + str(len(seed_xs)) + ': seed_x = ' + str(seed_x) + ' => start seed_y = ' + str(seed_y) + ' which has log_likelihood = ' + str(start_likelihood))

        #minimizing -1.0 * our function means maximizing our function (since it is always positive)
        # Also, minimize orthogonal to best fit polynomial
        #funct_to_minimize = lambda best_fit_x: -1.0 * computeLikelihoodGivenSeedAndFixedParams(galaxy, withDisk, curve_params,
        #                                                                                       - (np.sqrt(best_fit_x ** 2.0 + fit_funct(best_fit_x) ** 2.0)
        #                                                                                          * fit_funct(best_fit_x)) / (best_fit_x * np.sqrt(1 + (fit_funct(best_fit_x) / best_fit_x) ** 2.0 )),
        #                                                                                       (np.sqrt(best_fit_x ** 2.0 + fit_funct(best_fit_x) ** 2.0)
        #                                                                                          * 1.0) / (best_fit_x * np.sqrt(1 + (fit_funct(best_fit_x) / best_fit_x) ** 2.0 )),
        #                                                                                       param_initializer, pops, pop_selection_method, apply_observation_mask, halo_funct, disk_funct,
        #                                                                                       dist, arcmin_limits_R, arcmin_limits_z, arcmin_limits_los, gamma, outer_step,
        #                                                                                       fixed_hx_rs_val, fixed_hz_rs_val)

        if fixed_slope is None:
            fit_lines_slope = -1.0 * ((param_display_vals[curve_params[1]][1] - param_display_vals[curve_params[1]][0])
                                      / (param_display_vals[curve_params[0]][1] - param_display_vals[curve_params[0]][0]))
        else:
            fit_lines_slope = fixed_slope
        print ('fit_lines_slope = ' + str(fit_lines_slope) )
        funct_to_minimize = lambda best_fit_x: -1.0 * computeLikelihoodGivenSeedAndFixedParams(galaxy, withDisk, curve_params, best_fit_x, seed_y + fit_lines_slope * (best_fit_x - seed_x),
                                                                                               param_initializer, pops, pop_selection_method, apply_observation_mask, halo_funct, disk_funct, halo_type, disk_type,
                                                                                               dist, gamma, n_los_bins, target_n_sky_pixels, mask_R_bins,
                                                                                               fixed_hx_rs_val, fixed_hz_rs_val, fixed_dx_rs_val, fixed_dz_rs_val, star_data_sets, pop_connector = pop_connector,
                                                                                               include_kinem_prob = include_kinem_prob, include_morph_prob = include_morph_prob )
        min_results = optimize.minimize_scalar(funct_to_minimize,
                                               bounds = [max(seed_x - x_min_width, x_range_to_check[0]), min(seed_x + x_min_width, x_range_to_check[1])],
                                               method = 'bounded', options = {'maxiter':max_refine_steps})
        print ('min_results = ')
        print (min_results )
        min_x = min_results['x']
        min_y = seed_y + fit_lines_slope * (min_x - seed_x)
        log_prob = min_results['fun'] * -1.0
        min_xs = min_xs + [min_x]
        min_ys = min_ys + [min_y]
        #print 'log_prob = ' + str(log_prob)
        #print 'np.isnan(log_prob) = ' + str(np.isnan(log_prob))
        if np.isnan(log_prob):
            print ('log prob returned np.nan.  This probability will therefore be listed as 0. ')
            log_prob = 0.0
        used_xs = used_xs + [seed_x]
        used_ys = used_ys + [seed_y]
        log_probs_along_line = log_probs_along_line + [log_prob]
        n_measured_points = n_measured_points + 1
        if max_log_likelihood > log_prob:
            n_points_past_max = n_points_past_max + 1
        else:
            max_log_likelihood = log_prob
            n_points_past_max = 0
        if n_points_past_max >= n_points_past_max_to_break:
            print ('The maximum likelihood was encountered ' + str(n_points_past_max) + ' iterations ago.  Calling that sufficient. ')
            max_found = 1
    print ('min_xs = ' + str(min_xs))
    print ('min_ys = ' + str(min_ys) )
    print ('log_probs_along_line = ' + str(log_probs_along_line))
    if show:
        plt.scatter(min_xs, log_probs_along_line)
        plt.show()
    rough_best_fit_index = np.argmax(log_probs_along_line)
    print ('rough_best_fit_index = ' + str(rough_best_fit_index) )
    rough_best_fit_point = [min_xs[rough_best_fit_index], min_ys[rough_best_fit_index]]
    print ('rough_best_fit_point = ' + str(rough_best_fit_point) )

    y_scaling_for_minim = 10.0 ** (-6.0) #Brings M and rs to same order, which is necessary for effective minimization
    funct_to_minimize = lambda x_y_point: -1.0 * computeLikelihoodGivenSeedAndFixedParams(galaxy, withDisk, curve_params, x_y_point[0], x_y_point[1] * 1.0 / y_scaling_for_minim,
                                                                                          param_initializer, pops, pop_selection_method, apply_observation_mask, halo_funct, disk_funct, halo_type, disk_type,
                                                                                          dist, gamma, n_los_bins, target_n_sky_pixels, mask_R_bins,
                                                                                          fixed_hx_rs_val, fixed_hz_rs_val, fixed_dx_rs_val, fixed_dz_rs_val, star_data_sets, pop_connector = pop_connector,
                                                                                          include_kinem_prob = include_kinem_prob, include_morph_prob = include_morph_prob  )
    bounds = [(max(1.0, rough_best_fit_point[0] - x_min_width), rough_best_fit_point[0] + x_min_width),
              ((max(1.0, rough_best_fit_point[1] - y_min_width)) * y_scaling_for_minim, (rough_best_fit_point[1] + y_min_width) * y_scaling_for_minim)]
    #print 'bounds = ' + str(bounds)
    global_best_fit = optimize.minimize(funct_to_minimize, [np.around(rough_best_fit_point[0],0), np.around(rough_best_fit_point[1] * y_scaling_for_minim, 0)],
                                        bounds = bounds,
                                        method = 'Powell', tol = 10.0 ** (-2.0), options = {'maxiter':max_final_steps})
    print ('global_best_fit = ' + str(global_best_fit))
    print ("global_best_fit['x'] = " + str(global_best_fit['x']))


    return [used_xs, used_ys, log_probs_along_line, global_best_fit]
