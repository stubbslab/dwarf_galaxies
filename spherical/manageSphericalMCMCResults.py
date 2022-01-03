import csv
import numpy as np
import matplotlib
import matplotlib.cm as cm
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import scipy.integrate as integrate
import scipy.optimize as optimize
import scipy.special as special
import scipy.interpolate as interpolate
from logList import logList
from matplotlib.colors import LogNorm
from ModulateCyclicParam import ModulateCyclicParam
from cantrips import calculateFunctionAlongCurve
import math
import time
from ComputationalArchive import ComputationalArchive
from binData import binData
from scipy.interpolate import RegularGridInterpolator
from cantrips import safeSortOneListByAnother
from matplotlib.ticker import NullFormatter
from SphericalBackgroundMCMCInformationStorer import BackgroundStorer

#function object must looks like: {[(var_name1, var_name2, ..., var_namen), string_descriptor]: function using variables }


def updateVarDictForPopDependentVars(dict_to_update, pops = [], pop_connector = '_', by_pop_change_functs = None):
    #print ('dict_to_update = ' + str(dict_to_update))
    old_var_keys = list(dict_to_update.keys())
    pop_dependent_vars = getVarPopDependence()
    for var in old_var_keys:
        if var in pop_dependent_vars:
            for i in range(len(pops)):
                pop = pops[i]
                #print ('[var, pop] = ' + str([var, pop]))
                if by_pop_change_functs == None:
                    dict_to_update[var + pop_connector + pop] = dict_to_update[var]
                else:
                    update_funct = by_pop_change_functs[i]
                    dict_to_update[var + pop_connector + pop] = update_funct(dict_to_update[var])
            dict_to_update.pop(var)
    return dict_to_update

def getNBinsArray(pops = [], pop_connector = '_'):
    n_bins_array = {'M':200, 'rs':200, 'h_x_center':200, 'h_y_center':200, 'h_z_center':200,
                    'sigsqr_rr_0':100, 'sigsqr_rr_inf_to_0_rat':100, 'r_sigsqr_rr0':100, 'alpha_sigsqr_rr':100, 'gamma_for_beta_inf':100, 'r_beta0':50, 'alpha_beta':100, }
    if len(pops) > 0:
        n_bins_array = updateVarDictForPopDependentVars(n_bins_array, pops = pops, pop_connector = pop_connector)
    return n_bins_array

def getParamDisplayVals(pops = [], pop_connector = '_'):
    n_bins_array = getNBinsArray()
    bg_storer = BackgroundStorer()
    r_half_light = 791.0
    M_star = 10.0 ** 7.39
    param_display_vals = {'M':[*bg_storer.param_ranges['M'], n_bins_array['M']],
                          'rs':[*bg_storer.param_ranges['rs'], n_bins_array['rs']],
                          'h_x_center':[*bg_storer.param_ranges['h_x_center'], n_bins_array['h_x_center']],
                          'h_y_center':[*bg_storer.param_ranges['h_y_center'], n_bins_array['h_y_center']],
                          'h_z_center':[*bg_storer.param_ranges['h_z_center'], n_bins_array['h_z_center']],
                          'sigsqr_rr_0':[*bg_storer.param_ranges['sigsqr_rr_0'], n_bins_array['sigsqr_rr_0']],
                          'sigsqr_rr_inf_to_0_rat':[*bg_storer.param_ranges['sigsqr_rr_inf_to_0_rat'], n_bins_array['sigsqr_rr_inf_to_0_rat']],
                          'r_sigsqr_rr0':[*bg_storer.param_ranges['r_sigsqr_rr0'], n_bins_array['r_sigsqr_rr0']],
                          'alpha_sigsqr_rr':[*bg_storer.param_ranges['alpha_sigsqr_rr'], n_bins_array['alpha_sigsqr_rr']],
                          'gamma_for_beta_inf':[*bg_storer.param_ranges['gamma_for_beta_inf'], n_bins_array['gamma_for_beta_inf']],
                          'r_beta0':[*bg_storer.param_ranges['r_beta0'], n_bins_array['r_beta0']],
                          'alpha_beta':[*bg_storer.param_ranges['alpha_beta'], n_bins_array['alpha_beta']],
                          }
    param_display_vals = updateVarDictForPopDependentVars(param_display_vals, pops = pops, pop_connector = pop_connector)
    return param_display_vals

def getVarTicks(pops = [], pop_connector = '_'):
    #var_ticks = {'M' : [1000,2000,3000,4000, 5000, 6000, 7000, 8000, 9000], 'el':[10.0 ** -0.8, 10.0 ** -0.4, 10.0 ** -0.0, 10.0 ** 0.4, 10.0 ** 0.8], 'rs':[1000, 2000, 3000, 4000, 5000, 6000, 7000],
    #             'phi':[math.pi / 6.0, math.pi / 3.0, math.pi / 2.0, 2.0 * math.pi / 3.0, 5.0 * math.pi / 6.0], 'theta':[math.pi / 6.0, math.pi / 3.0, math.pi / 2.0, 2.0 * math.pi / 3.0, 5.0 * math.pi / 6.0],
    #             'h_x_center':[-650, -325, 0, 325, 650], 'h_y_center':[-650, -325, 0, 325, 650], 'h_z_center':[-650, -325, 0, 325, 650],
    #             'lam':[0.11, 0.13, 0.15, 0.17, 0.19], 'Rd':[1400, 2800, 4200, 5600, 7000], 'eps':[0.05, 0.15, 0.25, 0.35, 0.45],
    #             'a':np.around([math.pi / 6.0, math.pi / 3.0, math.pi / 2.0, 2.0 * math.pi / 3.0, 5.0 * math.pi / 6.0],3), 'b':np.around([math.pi / 6.0, math.pi / 3.0, math.pi / 2.0, 2.0 * math.pi / 3.0, 5.0 * math.pi / 6.0],3), 'd_x_center':[-650, -325, 0, 325, 650], 'd_y_center':[-650, -325, 0, 325, 650], 'd_z_center':[-650, -325, 0, 325, 650] }
    var_ticks = {'M' : [2000, 7000, 12000], 'rs':[2000, 4000, 6000],
                 'h_x_center':[-0.2, 0.2], 'h_y_center':[-0.2, 0.2], 'h_z_center':[-0.2, 0.2],
                 'sigsqr_rr_0':[100.0, 300.0],  'sigsqr_rr_inf_to_0_rat':[0.2, 1.0, 5.0], 'r_sigsqr_rr0':[2000, 4000, 6000], 'alpha_sigsqr_rr':[1.5, 2.5, 3.5, 4.5],
                 'gamma_for_beta_inf':[-0.5, 0.5], 'r_beta0':[2000, 4000, 6000], 'alpha_beta':[1.5, 3.0, 4.5]}
    var_ticks = updateVarDictForPopDependentVars(var_ticks, pops = pops, pop_connector = pop_connector)
    return var_ticks

def getVarTickLabels(pops = [], pop_connector = '_'):
    var_ticklabels = {'M' : ['2', '7', '12'], 'rs':['2', '4', '6'],
                 'h_x_center':[r'$-0.2$', '0.2'], 'h_y_center':[r'$-0.2$', '0.2'], 'h_z_center':[r'$-0.2$', '0.2'],
                 'sigsqr_rr_0':['100','300'],  'sigsqr_rr_inf_to_0_rat':['0.2', '1', '5'], 'r_sigsqr_rr0':['2', '4', '6'], 'alpha_sigsqr_rr':['1.5', '2.5', '3.5', '4.5'],
                 'gamma_for_beta_inf':[r'$-0.5$',  '0.5'], 'r_beta0':['2', '4', '6'], 'alpha_beta':['1.5', '3', '4.5']}
    var_ticklabels = updateVarDictForPopDependentVars(var_ticklabels, pops = pops, pop_connector = pop_connector)
    return var_ticklabels

def getBinTypes(pops = [], pop_connector = '_'):
    bin_types = {'M' : 'lin', 'rs':'lin', 'h_x_center':'lin', 'h_y_center':'lin', 'h_z_center':'lin',
                 'sigsqr_rr_0':'lin', 'sigsqr_rr_inf_to_0_rat':'log', 'r_sigsqr_rr0':'lon', 'alpha_sigsqr_rr':'lin',
                 'gamma_for_beta_inf':'lin', 'r_beta0':'lin', 'alpha_beta':'lin'}
    bin_types = updateVarDictForPopDependentVars(bin_types, pops = pops, pop_connector = pop_connector)
    return bin_types

def getDisplayLabels(pops = [], pop_connector = '_', pop_additional_labels = {r'$MP$', r'$IM$', r'$MR$'}):
    disp_labels = {'M':r'$M$ ', 'el':r'$\mathrm{log}_{10} (Q)$', 'rs':'$r_s$ ' ,
                      'h_x_center':r'$h_c$' + r'$_,$' + r'$_x$' , 'h_y_center':r'$h_c$' + r'$_,$' + r'$_y$' ,  'h_z_center':r'$h_c$' + r'$_,$' + r'$_z$' ,
                      'sigsqr_rr_0':r'$\sigma_{rr,0}$', 'sigsqr_rr_inf_to_0_rat':r'$\sigma_{rr,0}/\sigma_{rr,infty}$', 'r_sigsqr_rr0':r'$r_{\sigma, rr}$','alpha_sigsqr_rr':r'$\alpha_{\sigma, rr}$',
                      'gamma_for_beta_inf':r'$\gamma_{\infty}$', 'r_beta0':r'$r_{\beta,0}$', 'alpha_beta':r'$\alpha_{\beta}$',}
    disp_labels = updateVarDictForPopDependentVars(disp_labels, pops = pops, pop_connector = pop_connector, by_pop_change_functs = [lambda old_disp_label, i = i: old_disp_label + pop_additional_labels[i] for i in range(len(pop_additional_labels)) ] )

    return disp_labels

def getDispUnits(pops = [], pop_connector = '_'):
    var_units = {'M' : r'$10^9$' + ' ' + r'$M_\odot$', 'rs':'kpc', 'h_x_center':'deg', 'h_y_center':'deg', 'h_z_center':'deg',
                 'sigsqr_rr_0':'km/s', 'sigsqr_rr_inf_to_0_rat':'', 'r_sigsqr_rr0':'kpc', 'alpha_sigsqr_rr':'',
                 'gamma_for_beta_inf':'', 'r_beta0':'kpc', 'alpha_beta':''}
    var_units = updateVarDictForPopDependentVars(var_units, pops = pops, pop_connector = pop_connector)
    return var_units

def getVarCyclicity(pops = [], pop_connector = '_'):
    var_cyclicity = {'M' : 'none', 'rs':'none', 'h_x_center':'none', 'h_y_center':'none', 'h_z_center':'none',
                     'sigsqr_rr_0':'none', 'sigsqr_rr_inf_to_0_rat':'none', 'r_sigsqr_rr0':'none', 'alpha_sigsqr_rr':'none',
                     'gamma_for_beta_inf':'none', 'r_beta0':'none', 'alpha_beta':'none', }
    var_cyclicity = updateVarDictForPopDependentVars(var_cyclicity, pops = pops, pop_connector = pop_connector)
    return var_cyclicity

def getVarIndeces():
    var_indeces = {'M' : 14, 'rs':15, 'h_x_center':16, 'h_y_center':17, 'h_z_center':18,
                   'sigsqr_rr_0':{'MP':19,'IM':26,'MR':33}, 'sigsqr_rr_inf_to_0_rat':{'MP':20,'IM':27,'MR':34}, 'r_sigsqr_rr0':{'MP':21,'IM':28,'MR':35}, 'alpha_sigsqr_rr':{'MP':22,'IM':29,'MR':36},
                   'gamma_for_beta_inf':{'MP':23,'IM':30,'MR':37}, 'r_beta0':{'MP':24,'IM':31,'MR':38}, 'alpha_beta':{'MP':25,'IM':32,'MR':39}
                   }
    return var_indeces

def getVarPopDependence():
    pop_dependent_vars = ['sigsqr_rr_0', 'sigsqr_rr_inf_to_0_rat', 'r_sigsqr_rr0', 'alpha_sigsqr_rr', 'gamma_for_beta_inf', 'r_beta0', 'alpha_beta']
    return pop_dependent_vars

def getVarDisplayScalings(pops = [], pop_connector = '_'):
    var_display_scalings = {'M' : 10.0 ** (-6.0),'rs':1.0, 'phi':1.0, 'theta':1.0, 'h_x_center':1.0, 'h_y_center':1.0, 'h_z_center':1.0,
                            'sigsqr_rr_0':1.0, 'sigsqr_rr_inf_to_0_rat':1.0, 'r_sigsqr_rr0':1.0, 'alpha_sigsqr_rr':1.0,
                            'gamma_for_beta_inf':1.0, 'r_beta0':1.0, 'alpha_beta':1.0,}
    var_display_scalings = updateVarDictForPopDependentVars(var_display_scalings, pops = pops, pop_connector = pop_connector)
    return var_display_scalings

def readInMCMCResults(results_files,
                      functions_to_compute = [], vars_to_load = ['M', 'rs', 'h_x_center', 'h_z_center',], pops = ['MP', 'IM', 'MR'], pop_connector = '_',
                      results_dir = '/Users/sashabrownsberger/Documents/Harvard/physics/randall/MCMCOutputs/', n_ignore = 0):

    possible_operations = {'+','-','*','/'}
    var_indeces = getVarIndeces()
    pop_dependent_vars = getVarPopDependence()
    #y_ticks_by_var = {'M' : [0.2], 'el':[0.1, 0.5, 1.0, 2.0, 9.0], 'rs': [], 'phi': 6, 'theta':7, 'h_x_center':11, 'h_y_center':12, 'h_z_center':13,
    #                  'lam':[0., 10.15, 0.2}, 'Rd':{}, 'eps':17, 'a': 18, 'b':19, 'd_x_center':23, 'd_y_center':24, 'd_z_center':25}
    measured_arrays = {}
    likelihood_array = []
    n_visits_array = []
    for var in vars_to_load:
        if var in pop_dependent_vars:
            for pop in pops:
                measured_arrays[var + pop_connector + pop] = []
        else:
            measured_arrays[var] = []
    for function_list in functions_to_compute:
        measured_arrays[function_list[0]] = []

    n_visits_index = 40
    likelihood_index = 41

    for results_file in results_files:
        data_file = results_dir + results_file
        print ('Reading in ' + data_file )
        skip_first_line = True
        n_visits = 0
        with open(data_file) as csvfile:
            myreader = csv.reader(csvfile,delimiter=',')
            row_number = 0
            new_measured_arrays = {}
            new_likelihood_array = []
            new_n_visits_array = []
            for var in measured_arrays.keys():
                if var in pop_dependent_vars:
                    for pop in pops:
                        new_measured_arrays[var + pop_connector + pop] = []
                else:
                    new_measured_arrays[var] = []
            for row in myreader:
                if skip_first_line:
                    skip_first_line = False
                else:
                    row_number = row_number + 1
                    new_n_visits = int(row[n_visits_index])
                    if n_visits >= n_ignore:
                        for var in vars_to_load:
                            if var in pop_dependent_vars:
                                for pop in pops:
                                    new_measured_arrays[var + pop_connector + pop] = new_measured_arrays[var + pop_connector + pop]  + [float( row[ var_indeces[var][pop] ] )]
                            else:
                                new_measured_arrays[var] = new_measured_arrays[var] + [float( row[ var_indeces[var] ] )]
                        new_likelihood_array = new_likelihood_array + [float(row[likelihood_index])]
                        new_n_visits_array = new_n_visits_array + [new_n_visits]
                    n_visits = n_visits + new_n_visits

            likelihood_array = likelihood_array + new_likelihood_array
            n_visits_array = n_visits_array + new_n_visits_array
            n_visits = n_visits + float(row[n_visits_index])
            for var in measured_arrays.keys():
                measured_arrays[var] = measured_arrays[var] + new_measured_arrays[var]
    for function_list in functions_to_compute:
        measured_arrays[function_list[0]] = [ function_list[2]( *[measured_arrays[var][i] for var in function_list[1]]) for i in range(len(measured_arrays[vars_to_load[0]])) ]


    return [measured_arrays, n_visits_array, likelihood_array ]

def measureMCMCResults(results_files, vars_to_fit, fitting_functions, guess_params_set, function_to_fit = [],
                       pops = ['MP','IM','MR'], pop_connector = '_', pop_additional_labels = [r'$_{,MP}$', r'$_{,IM}$', r'$_{,MR}$'],
                       results_dir = '/Users/sashabrownsberger/Documents/Harvard/physics/randall/MCMCOutputs/', n_ignore = 0, n_bins_set = [200],
                       show_fit = 1, ranges_to_fit = ['all'], bin_buffer = 5, maxfev = 2000, bounds = [(-np.inf, np.inf)], max_n_sigma_right = 0.1):
    #print 'fitting_functions = ' + str(fitting_functions)
    #print 'guess_params_set = ' + str(guess_params_set)
    MCMC_background_storer = BackgroundStorer()
    #param_sampling_ranges = MCMC_background_storer.param_ranges
    param_display_vals = getParamDisplayVals(pops = pops)
    param_sampling_ranges = {key:[param_display_vals[key][0], param_display_vals[key][1]] for key in param_display_vals.keys()}
    print ('param_sampling_ranges = ' + str(param_sampling_ranges) )
    var_indeces = getVarIndeces()
    n_visits_index = 40
    likelihood_index = 41

    var_cyclicity = getVarCyclicity()
    n_bins_array = getNBinsArray(pops = pops)

    vars_to_read_in = vars_to_fit
    for funct in function_to_fit:
        vars_to_read_in = vars_to_read_in + funct[1]
    print ('vars_to_read_in = ' + str(vars_to_read_in))
    vars_to_read_in = list(set(vars_to_read_in))
    vars_to_read_in = [var for var in vars_to_read_in if var in var_indeces.keys()]
    #print 'vars_to_read_in = ' + str(vars_to_read_in)

    measured_arrays, n_visits_array, likelihood_array = readInMCMCResults(results_files, functions_to_compute = function_to_fit, vars_to_load = vars_to_read_in,
                                                                          results_dir = results_dir, n_ignore = n_ignore )


    print ('measured_arrays.keys() = ' + str(measured_arrays.keys() ))
    pop_dependent_vars = getVarPopDependence()
    replace_vars_to_disp = []
    for i  in range(len(vars_to_fit)):
        var = vars_to_fit[i]
        print ('[var, var in pop_dependent_vars] = ' + str([var, var in pop_dependent_vars] ))
        if var in pop_dependent_vars:
            new_vars_to_disp = [var + pop_connector + pop for pop in pops]
            replace_vars_to_disp = replace_vars_to_disp +  new_vars_to_disp
        else:
            replace_vars_to_disp = replace_vars_to_disp + [var]
    vars_to_fit = replace_vars_to_disp[:]

    if n_bins_set is None:
        n_bins_set = n_bins_array
    else:
        if  len(n_bins_set) == 1:
            n_bins_set = [n_bins_set[0] for var in vars_to_fit]
        n_bins_set = {vars_to_fit[i]:n_bins_set[i] for i in range(len(vars_to_fit))}
    if len(ranges_to_fit) == 1:
        ranges_to_fit = [ranges_to_fit[0] for var in vars_to_fit]
    if len(fitting_functions) == 1:
        fitting_functions = [fitting_functions[0] for var in vars_to_fit]
    if len(guess_params_set) == 1:
        guess_params_set = [guess_params_set[0] for var in vars_to_fit]
    if len(bounds) == 1:
        bounds = [bounds[0] for var in vars_to_fit]

    print ('vars_to_fit = ' + str(vars_to_fit))

    fit_results_set = []
    for i in range(len(vars_to_fit)):
        var_to_fit = vars_to_fit[i]
        print ('Working on fit for variable: ' + str(var_to_fit) )
        param_range = param_sampling_ranges[var_to_fit]
        fitting_function = fitting_functions[i]
        guess_params = guess_params_set[i]
        range_to_fit = ranges_to_fit[i]
        n_bins = n_bins_set[var_to_fit]
        print ('Binning data...' )
        x_vals, full_y_vals = binData(measured_arrays[var_to_fit], n_visits_array, n_bins = n_bins)
        print ('Done binning data...' )
        y_vals = [full_y_vals[0][j] * full_y_vals[2][j] for j in range(len(full_y_vals[0]))]
        max_index = np.argmax(y_vals)
        single_bounds = bounds[i]
        recompute_param_function = lambda params: params
        if range_to_fit is 'all':
            x_vals = x_vals[bin_buffer: len(x_vals) - bin_buffer]
            y_vals = y_vals[bin_buffer: len(y_vals) - bin_buffer]
        else:
            x_vals = x_vals[range_to_fit[0]:range_to_fit[1]]
            y_vals = y_vals[range_to_fit[0]:range_to_fit[1]]
        if fitting_function is 'single_hump':
            fitting_function = lambda xs, A, mu, sigma: A * np.exp(-(xs - mu)**2.0 /(2.0 * sigma ** 2.0))
            single_bounds = ([0.0, param_sampling_ranges[var_to_fit][0], 0.0], [np.inf, param_sampling_ranges[var_to_fit][1], (param_sampling_ranges[var_to_fit][1] - param_sampling_ranges[var_to_fit][0]) / 2.0])
        elif fitting_function is 'single_hump_shift':
            fitting_function = lambda xs, A, mu, sigma, power, shift: A * np.exp(-np.abs(xs - mu)**power /((sigma) ** power)) + shift
            #fitting_function = lambda xs, A, mu, sigma, shift: A * np.exp(-np.abs(xs - mu)**2.0 /((np.sqrt(2.0) * sigma) ** 2.0)) + shift
            single_bounds = ([0.0, param_sampling_ranges[var_to_fit][0], 0.0, 0.0, -np.inf], [np.inf, param_sampling_ranges[var_to_fit][1], (param_sampling_ranges[var_to_fit][1] - param_sampling_ranges[var_to_fit][0]) / 2.0, 10.0, np.inf])
            recompute_param_function = lambda params: [params[0], params[1], np.sqrt(params[2] ** 2.0 * special.gamma(3.0 / params[3]) / special.gamma(1.0 / params[3])), params[3], params[4]]
            #single_bounds = ([0.0, param_sampling_ranges[var_to_fit][0], 0.0,  -np.inf], [np.inf, param_sampling_ranges[var_to_fit][1], (param_sampling_ranges[var_to_fit][1] - param_sampling_ranges[var_to_fit][0]) / 2.0, np.inf])
        elif fitting_function is 'flat_line':
            fitting_function = lambda xs, A: [A for x in xs]
        elif fitting_function is 'gamma_shift':
            lower_bound = param_range[0]
            upper_bound = param_range[1]
            interval = upper_bound - lower_bound
            fitting_function = lambda xs, A, k, theta, shift: A * theta ** (-k) / math.gamma(k) * ((xs - lower_bound) / interval) ** (k - 1.0) * np.exp(-(xs - lower_bound) / (interval * theta) )
            recompute_param_function = lambda params: [params[0], (params[1]-1.0) * params[2] * interval + lower_bound, np.sqrt(params[1] * (params[2] * interval)**2.0), params[3]]
        elif fitting_function is 'beta_shift':
            lower_bound = param_range[0]
            upper_bound = param_range[1]
            interval = upper_bound - lower_bound
            B_approx = lambda alpha, beta: math.gamma(alpha) * beta ** (-alpha ) * (1 - (-alpha * (alpha + 1.)) / (2.0 * beta))
            B = lambda alpha, beta: (math.gamma(alpha) * math.gamma(beta)) / math.gamma(alpha + beta) if beta < 100.0 * alpha else B_approx(alpha, beta)
            fitting_function = lambda xs, A, alpha, beta, shift: (A * ((xs - lower_bound) / interval) ** (alpha - 1.0) * (1.0 - (xs - lower_bound) / interval) ** (beta - 1.0)) / B(alpha, beta) + shift
            mu_funct = lambda alpha, beta: ((alpha-1.0) / (alpha + beta - 2.0)  if alpha > 1.0 and beta > 1.0
                                                                     else 0.0 if alpha <= 1.0 and beta > 1.0
                                                                     else 1.0 if alpha > 1.0 and beta <= 1.0
                                                                     else 0.0 if abs(x_vals[max_index] - lower_bound) < abs(x_vals[max_index] - upper_bound)
                                                                     else 1.0
                                                                   )
            cumulative_prob_density = lambda alpha, beta, x: special.betainc(alpha, beta, x)
            normal_cpb = lambda sig: 0.5 * (1.0 + special.erf(sig / np.sqrt(2.0)))
            recompute_param_function = lambda params: [params[0],
                                                       mu_funct(params[1], params[2]) * interval + lower_bound,
                                                       (mu_funct(params[1], params[2]) - optimize.minimize_scalar(lambda x: np.abs(cumulative_prob_density(params[1], params[2], x) - cumulative_prob_density(params[1], params[2], mu_funct(params[1], params[2])) * normal_cpb(-1.0) * 2.0 ), bounds = [0.0, mu_funct(params[1], params[2])], method = 'bounded', tol = 10.0 ** -9.0)['x']) * interval,
                                                       (optimize.minimize_scalar(lambda x: np.abs(cumulative_prob_density(params[1], params[2], x) - (2.0 * (normal_cpb(1.0) - normal_cpb(0.0)) + cumulative_prob_density(params[1], params[2], mu_funct(params[1], params[2])) * (1.0 - 2.0 * (normal_cpb(1.0) - normal_cpb(0.0)) ))), bounds = [mu_funct(params[1], params[2]), max_n_sigma_right ], method = 'bounded', tol = 10.0 ** -9.0)['x'] - mu_funct(params[1], params[2])) * interval,
                                                       params[3]]
        #Measure initial best guess fit from data
        if guess_params is 'single_hump':
            #function is of form: A * exp (-(x-mu)**2 / (2.0 * sigma ** 2.0))
            A = np.max(y_vals) * 1.0
            mu = x_vals[np.argmax(y_vals)]
            #sigma = abs(mu - x_vals[np.argmin((np.array(y_vals) - A / math.exp(0.5)) ** 2.0)])
            sigma = (param_sampling_ranges[var_to_fit][1] - param_sampling_ranges[var_to_fit][0]) / 10.0
            guess_params = [A, mu, sigma]
        elif guess_params is 'single_hump_shift':
            #function is of form: A * exp (-(x-mu)**2 / (2.0 * sigma ** 2.0)) + shift
            A = np.max(y_vals) * 1.0
            mu = x_vals[np.argmax(y_vals)]
            #It is true that int_0^infinity (dx x e ^ (-x^2 / (2.0 * sigma ** 2.0))) = sigma ** 2.0.  So assuming the distribution is gaussian and I know the peak, this is a good first guess of the width
            sigma = np.sqrt(np.sum([(x_vals[1] - x_vals[0]) * abs(mu - x_vals[point_index]) * y_vals[point_index] / A for point_index in range(len(x_vals))]))
            power = 2.0
            shift = 0.0
            guess_params = [A, mu, sigma, power, shift]
            #guess_params = [A, mu, sigma,  shift]
        elif guess_params is 'flat_line':
            A = np.median(y_vals)
            guess_params = [A]
        elif guess_params is 'gamma_shift':
            A = np.max(y_vals) * 1.0
            lower_bound = param_range[0]
            upper_bound = param_range[1]
            interval = upper_bound - lower_bound
            mu = x_vals[np.argmax(y_vals)]
            sigma = abs(mu - x_vals[np.argmin((np.array(y_vals) - A / math.exp(0.5)) ** 2.0)])
            shift = 0.0
            mu = (mu - lower_bound) / interval
            sigma = sigma / interval
            theta = sigma ** 2.0 / mu
            k = mu/ (theta)
            guess_params = [A, k, theta, shift]
            print ('[lower_bound, interval] = ' + str([lower_bound, interval]))
        elif guess_params is 'beta_shift':
            A = np.max(y_vals) * 1.0
            mu = x_vals[np.argmax(y_vals)]
            bin_interp = interpolate.interp1d(x_vals, y_vals)
            min_x_val = np.min(x_vals)
            max_x_val = np.max(x_vals)
            mu = integrate.quad(lambda x: bin_interp(x) * x, min_x_val, max_x_val)[0] / integrate.quad(lambda x: bin_interp(x), min_x_val, max_x_val)[0]
            sigma = abs(mu - x_vals[np.argmin((np.array(y_vals) - A / math.exp(0.5)) ** 2.0)])
            sigma = np.sqrt(integrate.quad(lambda x: bin_interp(x) * x ** 2.0, min_x_val, max_x_val)[0] / integrate.quad(lambda x: bin_interp(x), min_x_val, max_x_val)[0])
            shift = 0.0
            lower_bound = param_range[0]
            upper_bound = param_range[1]
            #lower_bound = min(x_vals)
            #upper_bound = max(x_vals)
            interval = upper_bound - lower_bound
            mu = (mu - lower_bound) / interval
            sigma = sigma / interval
            print ('[interval, lower_bound, upper_bound, mu, sigma] = ' + str([interval, lower_bound, upper_bound, mu, sigma]))
            alpha = (mu ** 2.0 - mu ** 3.0 - mu * (sigma ** 2.0)) / (sigma ** 2.0)
            beta = (-1.0 + mu) * (-mu + mu ** 2.0 + (sigma ** 2.0)) / (sigma ** 2.0)
            guess_params = [A, alpha, beta, shift]
        print ('guess_params = ' + str(guess_params))
        print ('single_bounds = ' + str(single_bounds))
        if show_fit:
            plt.scatter(x_vals, y_vals)
            plt.show()
        #print ('[fitting_function, np.array(x_vals), np.array(y_vals), guess_params, np.zeros(len(x_vals)) + 1.0] = ' + str([fitting_function, np.array(x_vals), np.array(y_vals), guess_params, np.zeros(len(x_vals)) + 1.0]) )
        results = optimize.curve_fit(fitting_function, np.array(x_vals), np.array(y_vals), guess_params, np.zeros(len(x_vals)) + 1.0, maxfev = maxfev, bounds = single_bounds)
        print ('best_fit_results = ' + str(results[0]))
        recomputed_params = recompute_param_function(results[0])
        fit_results_set = fit_results_set + [recomputed_params]
        print ('recomputed_params = ' + str(recomputed_params))
        if show_fit:
            plt.scatter(x_vals, y_vals)
            plt.plot(sorted(x_vals), fitting_function(np.array(sorted(x_vals)), *(results[0])), c = 'r')
            plt.xlabel(var_to_fit)
            plt.ylabel('n steps')
            plt.show()
    print ('fit_results_set = ' + str(fit_results_set) )
    print ('fit_results_set[0][0] = ' + str(fit_results_set[0][0]))
    return fit_results_set
    #return [[0.0,np.median(x_vals)], []]


def verbose_funct_to_minimize(params, fit_funct, binned_n_visits_interp, fit_x_bin_centers, fit_y_bin_centers, deriv_of_fit):
    print ('params = ' + str(params))
    if deriv_of_fit is None:
        deriv_of_fit_to_pass = None
    else:
        deriv_of_fit_to_pass = lambda x: deriv_of_fit(x, *params)
    return -1.0 * calculateFunctionAlongCurve(lambda x: fit_funct(x, params), lambda x, y: float(binned_n_visits_interp((x, y))) + 1.0 if (y > fit_y_bin_centers[0] and y < fit_y_bin_centers[-1]) else 1.0, [fit_x_bin_centers[0], fit_x_bin_centers[-1]], curve_derivative = deriv_of_fit_to_pass)

def fitMCMCResults(results_files, vars_to_fit,
                   pops = ['MP','IM','MR'], measured_arrays = None, n_visits_array = None,
                   param_ranges_to_fit = ['buffer', 'buffer'], fit_funct = 'poly', params_for_fit = {'fit_order':3, 'deriv_of_fit':None, 'test_params':[], 'init_guess':[], 'bounds':None, 'use_fixed_point':1},
                   results_dir ='/Users/sashabrownsberger/Documents/Harvard/physics/randall/MCMCOutputs/', fixed_slope = None,
                   n_bins = [], n_ignore = 0, theta_shift = 0.0, n_fitted_points = 100, n_xy_points_to_fit = 1000, smallest_max_val_for_fit = 50.0):
    #print 'Fit info: '
    #print 'results_files = ' + str(results_files)
    #print 'vars_to_fit = ' + str(vars_to_fit)
    #print 'param_ranges_to_fit = ' + str(param_ranges_to_fit)
    #print 'fit_funct = ' + str(fit_funct)
    #print 'params_for_fit = ' + str(params_for_fit)
    #print 'results_dir = ' + str(results_dir)
    #print '[n_bins, n_ignore, theta_shift, n_fitted_points, n_xy_points_to_fit, smallest_max_val_for_fit] = ' + str([n_bins, n_ignore, theta_shift, n_fitted_points, n_xy_points_to_fit, smallest_max_val_for_fit] )
    var_cyclicity = getVarCyclicity(pops = pops)
    var_display_scalings = getVarDisplayScalings(pops = pops)
    param_display_vals = getParamDisplayVals(pops = pops)
    print ('param_display_vals = ' + str(param_display_vals) )
    bin_types = getBinTypes(pops = pops)
    x_var = vars_to_fit[0]
    y_var = vars_to_fit[1]
    if (measured_arrays is None) or (n_visits_array is None):
        measured_arrays, n_visits_array, likelihood_array = readInMCMCResults(results_files, vars_to_load = vars_to_fit, results_dir = results_dir, n_ignore = n_ignore )
    print ('[[measured_arrays[key][0], measured_arrays[key][-1]] for key in measured_arrays.keys()] = ' + str([[measured_arrays[key][0], measured_arrays[key][-1]] for key in measured_arrays.keys()]))
    print ('[n_visits_array[0],n_visits_array[-1]]  = ' + str([n_visits_array[0],n_visits_array[-1]]))
    for var in vars_to_fit:
        if var == 'theta':
            uncorrected_thetas = measured_arrays[var]
            print ('max(uncorrected_thetas) = ' + str(max(uncorrected_thetas)))
            print ('min(uncorrected_thetas) = ' + str(min(uncorrected_thetas)))
            measured_arrays[var] = [uncorrected_theta + theta_shift for uncorrected_theta in uncorrected_thetas]
        if var_cyclicity[var] != 'none':
            cyclicity = var_cyclicity[var]
            raw_measured_array = measured_arrays[var]
            measured_arrays[var] = [ ModulateCyclicParam(raw_var,cyclicity, shift = - math.pi) for raw_var in raw_measured_array ]

    vars_bin_centers = []
    vars_scalings = []
    vars_bins = []
    for i in range(len(vars_to_fit)):
        var = vars_to_fit[i]
        var_min = param_display_vals[var][0]
        var_max = param_display_vals[var][1]
        var_nbins = param_display_vals[var][2]
        var_step = (var_max - var_min) / var_nbins
        if bin_types[var] == 'lin':
            var_bins = [[var_min + var_step*i, var_min + var_step *(i+1)] for i in range(var_nbins)]
        elif bin_types[var] == 'log':
            log_spaced_borders = logList(var_min,var_max,var_nbins)
            var_bins = [ [log_spaced_borders[i],log_spaced_borders[i+1]] for i in range(len(log_spaced_borders)-1) ]

        var_bin_centers = [(var_bins[i][1] + var_bins[i][0])/ 2.0 if bin_types[var] == 'lin'
                           else 10.0 ** ((math.log10(var_bins[i][1]) + math.log10(var_bins[i][0]))/ 2.0) if bin_types[var]=='log'
                           else 0.0
                           for i in range(len(var_bins)) ]
        var_scaling = var_display_scalings[var]
        #param_display_vals[var] = param_display_vals[var] + [var_step, var_bins, var_bin_centers, var_scaling]
        vars_bin_centers = vars_bin_centers + [var_bin_centers]
        vars_scalings = vars_scalings + [var_scaling]
        vars_bins = vars_bins + [var_bins]

    total_visits = sum(n_visits_array)
    x_bins = vars_bins[0]
    y_bins = vars_bins[1]
    scaled_x_bin_centers = np.array(vars_bin_centers[0]) * vars_scalings[0]
    scaled_y_bin_centers = np.array(vars_bin_centers[1]) * vars_scalings[1]
    x_mesh, y_mesh = np.meshgrid(scaled_x_bin_centers, scaled_y_bin_centers)
    binned_n_visits = np.zeros(np.shape(x_mesh))
    for k in range(len(n_visits_array)):
        measured_x = measured_arrays[x_var][k]
        measured_y = measured_arrays[y_var][k]
        n_visits = n_visits_array[k]

        x_bin_index = -1
        y_bin_index = -1

        for l in range(len(x_bins)):
            if measured_x >=x_bins[l][0] and measured_x < x_bins[l][1]:
                x_bin_index = l
        for l in range(len(y_bins)):
            if measured_y >=y_bins[l][0] and measured_y < y_bins[l][1]:
                y_bin_index = l
        binned_n_visits[y_bin_index,x_bin_index] = binned_n_visits[y_bin_index,x_bin_index] + n_visits

    buffer_fraction = 1000
    print ('param_ranges_to_fit = ' + str(param_ranges_to_fit))
    print ('vars_scalings = ' + str(vars_scalings) )
    if param_ranges_to_fit[0] in ['all','All','ALL']:
        x_fit_index_bounds = [0, len(scaled_x_bin_centers) - 1]
    elif param_ranges_to_fit[0] in ['buff', 'bufferred', 'buffer','Buff', 'Bufferred', 'Buffer','BUFF', 'BUFFERRED', 'BUFFER']:
        x_fit_index_bounds = [len(scaled_x_bin_centers) // buffer_fraction, len(scaled_x_bin_centers) - len(scaled_x_bin_centers) // buffer_fraction]
    else:
        x_fit_index_bounds = [np.argmin(np.abs(np.array(scaled_x_bin_centers) - param_ranges_to_fit[0][0] * vars_scalings[0])), np.argmin(np.abs(np.array(scaled_x_bin_centers) - param_ranges_to_fit[0][1] * vars_scalings[0]))]

    if param_ranges_to_fit[1] in ['all','All','ALL']:
        y_fit_index_bounds = [0, len(scaled_y_bin_centers) - 1]
    elif param_ranges_to_fit[1] in ['buff', 'bufferred', 'buffer','Buff', 'Bufferred', 'Buffer','BUFF', 'BUFFERRED', 'BUFFER']:
        y_fit_index_bounds = [len(scaled_y_bin_centers) // buffer_fraction, len(scaled_y_bin_centers) - len(scaled_y_bin_centers) // buffer_fraction]
    else:
        y_fit_index_bounds = [np.argmin(np.abs(np.array(scaled_y_bin_centers) - param_ranges_to_fit[1][0] * vars_scalings[1])), np.argmin(np.abs(np.array(scaled_y_bin_centers) - param_ranges_to_fit[1][1] * vars_scalings[1]))]

    print ('x_fit_index_bounds = ' + str(x_fit_index_bounds))
    print ('y_fit_index_bounds = ' + str(y_fit_index_bounds))

    fit_x_bin_centers = scaled_x_bin_centers[x_fit_index_bounds[0]:x_fit_index_bounds[1]]
    fit_y_bin_centers = scaled_y_bin_centers[y_fit_index_bounds[0]:y_fit_index_bounds[1]]

    deriv_of_fit = params_for_fit['deriv_of_fit']
    binned_n_visits_interp = RegularGridInterpolator((scaled_x_bin_centers, scaled_y_bin_centers), binned_n_visits.transpose())
    point_in_range = lambda point: (point[0] > fit_x_bin_centers[0] and point[0] < fit_x_bin_centers[-1] and point[1] > fit_y_bin_centers[1] and point[1] < fit_y_bin_centers[-1])
    binned_n_visits_extended_interp = lambda points: [binned_n_visits_interp(point) if point_in_range(point) else 0.0 for point in points]

    if fixed_slope is None:
        fit_lines_slope = (fit_y_bin_centers[0] - fit_y_bin_centers[-1]) / (fit_x_bin_centers[-1] - fit_x_bin_centers[0])
        fit_lines_slope = -1.0 * (((param_display_vals[y_var][1] - param_display_vals[y_var][0]) * vars_scalings[1])
                                     / ((param_display_vals[x_var][1] - param_display_vals[x_var][0]) * vars_scalings[0]))
    else:
        fit_lines_slope = fixed_slope
    #print 'param_display_vals = '+ str(param_display_vals)
    #print 'fit_lines_slope = ' + str(fit_lines_slope)
    #print 'fit_x_bin_centers[-1] - fit_x_bin_centers[0] = ' + str(fit_x_bin_centers[-1] - fit_x_bin_centers[0])
    #print 'fit_y_bin_centers[0] - fit_y_bin_centers[-1] = ' + str(fit_y_bin_centers[0] - fit_y_bin_centers[-1])
    #print '(param_display_vals[x_var][1] - param_display_vals[x_var][0]) * vars_scalings[1] = ' + str((param_display_vals[x_var][1] - param_display_vals[x_var][0]) * vars_scalings[0])
    #print '(param_display_vals[y_var][1] - param_display_vals[y_var][0]) * vars_scalings[1] = ' + str((param_display_vals[y_var][1] - param_display_vals[y_var][0]) * vars_scalings[1])
    #print 'scaled_x_bin_centers[-1] - scaled_x_bin_centers[0] = ' + str(scaled_x_bin_centers[-1] - scaled_x_bin_centers[0])
    #print 'scaled_y_bin_centers[-1] - scaled_y_bin_centers[0] = ' + str(scaled_y_bin_centers[-1] - scaled_y_bin_centers[0])
    print ('fit_lines_slope = ' + str(fit_lines_slope) )
    fit_funct = lambda xs, mu, sig, A: A * np.exp(-(np.array(xs) - mu) ** 2.0 / (2.0 * sig ** 2.0))
    fit_res = []
    peak_xs = []
    peak_ys = []


    fit_box = [[fit_x_bin_centers[0], fit_x_bin_centers[-1]], [fit_y_bin_centers[0], fit_y_bin_centers[-1]]]
    print ('fit_box = ' + str(fit_box) )
    #origins are specified from top of box, and we need to account for additional x space since our sampling lines are diagonal.  We want to cover the whole box
    if abs(fit_lines_slope) > 0.0:
        fit_points_origins = [[x, fit_box[1][1]] for x in np.linspace(fit_box[0][0] - (fit_box[1][0] - fit_box[1][1]) * 1.0 / fit_lines_slope, fit_box[0][1], n_fitted_points) ]
    else:
        fit_points_origins = [[fit_box[0][1], y] for y in np.linspace(fit_box[1][0], fit_box[1][1], n_fitted_points) ]

    for m in range(n_fitted_points):
        fit_points_origin = fit_points_origins[m]
        #fit_points_origin = [fit_x_bin_centers[0] + (fit_x_bin_centers[-1] - fit_x_bin_centers[0]) * m / float((n_fitted_points - 1)),
        #                     fit_y_bin_centers[0] + (fit_y_bin_centers[-1] - fit_y_bin_centers[0]) * m / float((n_fitted_points - 1))]
        #fit_xy_points = [(fit_points_origin[0] + (fit_box[0][1] - fit_box[0][0]) * n / float((n_xy_points_to_fit - 1)),
        #                  fit_points_origin[1] + fit_lines_slope * (fit_box[0][1] - fit_box[0][0]) * n / float((n_xy_points_to_fit - 1)) )
        #                 for n in np.linspace(-n_xy_points_to_fit / 2.0, n_xy_points_to_fit / 2.0,  n_xy_points_to_fit) ]
        if abs(fit_lines_slope) > 0.0:
            fit_xy_points = [(fit_points_origin[0] + (-1.0 * (fit_box[1][1] - fit_box[1][0]) * n / float((n_xy_points_to_fit - 1))) / fit_lines_slope,
                              fit_points_origin[1] - (fit_box[1][1] - fit_box[1][0]) * n / float((n_xy_points_to_fit - 1)) )
                              for n in range(n_xy_points_to_fit) ]
        else:
            fit_xy_points = [(fit_points_origin[0] - (fit_box[0][1] - fit_box[0][0]) * n / float((n_xy_points_to_fit - 1)),
                              fit_points_origin[1] )
                              for n in range(n_xy_points_to_fit) ]
        #print ('fit_xy_points = ' + str(fit_xy_points))
    #plt.scatter([point[0] for point in fit_xy_points ], [point[1] for point in fit_xy_points ])

        try:
            n_bins_along_crossection = binned_n_visits_extended_interp(fit_xy_points)
            peak_index = np.argmax(n_bins_along_crossection)
            peak_val = float(n_bins_along_crossection[peak_index])
            peak_point = fit_xy_points[peak_index]
            width_guess_index = np.argmin((np.array(n_bins_along_crossection) - peak_val * math.exp(-0.5)) ** 2.0)
            if width_guess_index == peak_index: width_guess_index = peak_index + 1
            if peak_val < smallest_max_val_for_fit:
                print ('peak value of ' + str(peak_val) + ' less than specified threshold of ' + str(smallest_max_val_for_fit) + '.  No useful peak will be found there. ')
                fit_res  = fit_res + [None]
            else:
                init_guess = [peak_point[0], abs(peak_point[0] - fit_xy_points[width_guess_index][0]), peak_val]
                bounds = ([fit_x_bin_centers[0], 0.0, 0.0], [fit_x_bin_centers[-1], (fit_x_bin_centers[-1] - fit_x_bin_centers[0]) / 2.0, np.inf])
                new_fit_results = optimize.curve_fit(fit_funct, [point[0] for point in fit_xy_points], n_bins_along_crossection,
                                                     p0 = init_guess,bounds = bounds )
                #print ('new_fit_results[0] = ' + str(new_fit_results[0]))
                fit_res = fit_res + [new_fit_results]
                peak_xs = peak_xs + [new_fit_results[0][0]]
                peak_ys = peak_ys + [fit_points_origin[1] + (peak_xs[-1] - fit_points_origin[0]) * fit_lines_slope ]
                print ('new_peak: (x,y) = ' + str((peak_xs[-1], peak_ys[-1])) )
        except KeyboardInterrupt:
            print ('Keyboard interrupt while fitting points. ')
            return
        except:
            print ('Unable to fit ' + str(m) + ' set of points for some reason (some exception happended) . ')
            fit_res = fit_res + [None]
    print ('peak_xs = ' + str(peak_xs))
    print ('peak_ys = ' + str(peak_ys))

    return [[x / vars_scalings[0] for x in peak_xs], [y / vars_scalings[1] for y in peak_ys]]


def showMCMCResults(results_files, vars_to_disp,
                    params_to_fit = {}, fit_funct = 'poly', params_for_fit = {'fit_order':3, 'deriv_of_fit':None, 'test_params':[], 'init_guess':[], 'bounds':None, 'use_fixed_point':1},
                    pops = ['MP','IM','MR'], pop_connector = '_', pop_additional_labels = [r'$_{,MP}$', r'$_{,IM}$',r'$_{,MR}$'],
                    functions_to_display = [], show = 1, save_fig = 0,
                    results_dir = '/Users/sashabrownsberger/Documents/Harvard/physics/randall/MCMCOutputs/', extra_save_str = '', file_name = '', n_x_ticks = 5, n_y_ticks = 5,
                    n_bins = [], n_ignore = 0, theta_shift = 0.0, n_var_per_plot = 'both', fig_size_unit = 5.0, label_size_scaling = 16, tick_label_size_scaling = 13.0,
                    n_levels = 10, exp_level_powers = [-3.0, -2.0, -1.0], bin_label = r'$N_{\mathrm{bin}}$' + '\n' + r'$(10^3)$',
                    existing_subplots_array = None, return_axarr = 0, fit_line_color = 'black', n_test_lines = 100.0,
                    n_fitted_points = 100, n_xy_points_to_fit = 1000, smallest_max_val_to_be_on_curve = 50.0, fancyTex = 0,
                    y_label_coord_scalings = [-0.175, 0.5]):

    plt.rc('font', family='serif')
    #plt.tight_layout()
    if fancyTex:
        plt.rc('text', usetex=True)
    #plt.rc('xtick', labelsize='x-small')
    #plt.rc('ytick', labelsize='x-small')
    comp_archive = ComputationalArchive()
    plot_dir = comp_archive.getPlotDir()

    var_indeces = getVarIndeces()
    n_visits_index = 40
    likelihood_index = 41

    var_cyclicity = getVarCyclicity(pops = pops)
    var_display_scalings = getVarDisplayScalings(pops = pops)
    bin_types = getBinTypes(pops = pops)
    var_ticks = getVarTicks(pops = pops)
    var_ticklabels = getVarTickLabels(pops = pops)
    disp_units = getDispUnits(pops = pops)
    param_display_vals = getParamDisplayVals(pops = pops)
    display_labels = getDisplayLabels(pops = pops, pop_connector = pop_connector, pop_additional_labels = pop_additional_labels)

    vars_to_read_in = vars_to_disp

    for funct in functions_to_display:
        vars_to_read_in = vars_to_read_in + funct[1]
    vars_to_read_in = list(set(vars_to_read_in))
    vars_to_read_in = [var for var in vars_to_read_in if var in var_indeces.keys()]
    print('vars_to_disp = ' + str(vars_to_disp))
    print ('vars_to_read_in = ' + str(vars_to_read_in) )


    #n_bins_array = getNBinsArray()
    #for i in range(len(n_bins)):
    #    n_bins_array[vars_to_disp[i]] = n_bins[i]

    measured_arrays, n_visits_array, likelihood_array = readInMCMCResults(results_files, functions_to_compute = functions_to_display, vars_to_load = vars_to_read_in,
                                                                                         results_dir = results_dir, n_ignore = n_ignore )
    print ('measured_arrays.keys() = ' + str(measured_arrays.keys() ))
    pop_dependent_vars = getVarPopDependence()
    replace_vars_to_disp = []
    for i  in range(len(vars_to_disp)):
        var = vars_to_disp[i]
        print ('[var, var in pop_dependent_vars] = ' + str([var, var in pop_dependent_vars] ))
        if var in pop_dependent_vars:
            new_vars_to_disp = [var + pop_connector + pop for pop in pops]
            replace_vars_to_disp = replace_vars_to_disp +  new_vars_to_disp
        else:
            replace_vars_to_disp = replace_vars_to_disp + [var]
    vars_to_disp = replace_vars_to_disp[:]

    #print ('vars_to_disp = ' + str(vars_to_disp))
    #f2, axarr2 = plt.subplots(2,1)
    #axarr2[0].hist(measured_arrays['h_x_center'])
    #axarr2[1].hist(measured_arrays['h_z_center'])
    #plt.show()

    for funct in functions_to_display:
        var = funct[0]
        display_vals = [param_display_vals[component_var] for component_var in funct[1]]
        funct_display_vals = [funct[2](*[val[0] for val in display_vals]),
                              funct[2](*[val[2] for val in display_vals]),
                              display_vals[0][2]]
        param_display_vals[var] = funct_display_vals
        bin_types[var] = 'lin'
        var_cyclicity[var] = 'none'
        disp_units[var] = 'none'
        for component_var in funct[1]:
            if bin_types[component_var] != 'lin': bin_types[var] = bin_types[component_var]
            if var_cyclicity[component_var] != 'none': var_cyclicity[var] = var_cyclicity[component_var]
            if disp_units[component_var] != 'none': disp_units[var] = disp_units[component_var]

    print ('vars_to_disp = ' + str(vars_to_disp))
    for var in vars_to_disp:
        if var == 'theta':
            uncorrected_thetas = measured_arrays[var]
            print ('max(uncorrected_thetas = ' + str(max(uncorrected_thetas)))
            print ('min(uncorrected_thetas = ' + str(min(uncorrected_thetas)))
            measured_arrays[var] = [uncorrected_theta + theta_shift for uncorrected_theta in uncorrected_thetas]
        #print ('var_cyclicity = ' + str(var_cyclicity))
        if var_cyclicity[var] != 'none':
            cyclicity = var_cyclicity[var]
            raw_measured_array = measured_arrays[var]
            measured_arrays[var] = [ ModulateCyclicParam(raw_var,cyclicity, shift = - math.pi) for raw_var in raw_measured_array ]
            if var == 'theta':
                print ('raw_measured_array[10:30] = ')
                print (raw_measured_array[10:30])
                print ('measured_array["theta"][10:30] = ')
                print (measured_arrays[var][10:30])

    total_visits = sum(n_visits_array)

    for var in vars_to_disp:
        var_min = param_display_vals[var][0]
        var_max = param_display_vals[var][1]
        var_nbins = param_display_vals[var][2]
        var_step = (var_max - var_min) / var_nbins
        if bin_types[var] == 'lin':
            var_bins = [[var_min + var_step*i, var_min + var_step *(i+1)] for i in range(var_nbins)]
        elif bin_types[var] == 'log':
            log_spaced_borders = logList(var_min,var_max,var_nbins)
            var_bins = [ [log_spaced_borders[i],log_spaced_borders[i+1]] for i in range(len(log_spaced_borders)-1) ]

        var_bin_centers = [(var_bins[i][1] + var_bins[i][0])/ 2.0 if bin_types[var] == 'lin'
                           else 10.0 ** ((math.log10(var_bins[i][1]) + math.log10(var_bins[i][0]))/ 2.0) if bin_types[var]=='log'
                           else 0.0
                           for i in range(len(var_bins)) ]
        var_scaling = var_display_scalings[var]
        param_display_vals[var] = param_display_vals[var] + [var_step, var_bins, var_bin_centers, var_scaling]


    #plt.rcParams.update({'font.size': 6})

    if n_var_per_plot == 2 or n_var_per_plot in ['both','Both','BOTH']:

        n_x_contour_plots = len(vars_to_disp) - 1
        n_y_contour_plots = len(vars_to_disp) - 1
        if n_var_per_plot == 2:
            total_n_x_plots = n_x_contour_plots
            total_n_y_plots = n_y_contour_plots
        else:
            total_n_x_plots = n_x_contour_plots + 1
            total_n_y_plots = n_y_contour_plots + 1
        if existing_subplots_array is None:
            f, axarr = plt.subplots(total_n_x_plots, total_n_y_plots,
                                    figsize = (fig_size_unit * total_n_x_plots, fig_size_unit * total_n_y_plots),
                                    squeeze = False) #, sharex = True, sharey = True)
            #plt.tight_layout(pad = 6.0 )
            for i in range(np.shape(axarr)[0]):
                for j in range(np.shape(axarr)[1]):
                    axarr[j][i].set_xticks([])
                    axarr[j][i].set_yticks([])
                    axarr[j][i].tick_params(axis = 'x', direction = 'in')
                    axarr[j][i].tick_params(axis = 'y', direction = 'in')
            #plt.show()
        for i in range(len(vars_to_disp) - 1):
            #for j in [elem + 1 for elem in range(i+1)]:
            for j in range(i, len(vars_to_disp) - 1):
                #print 'i = ' + str(i)
                #print 'j = ' + str(j)
                if n_var_per_plot is 2:
                    axarr_x_index = i
                    #axarr_y_index = j-1
                    axarr_y_index = j
                else:
                    axarr_x_index = i + 1
                    axarr_y_index = j
                print ('[axarr_y_index, axarr_x_index] = ' + str([axarr_y_index, axarr_x_index]))
                x_var = vars_to_disp[i]
                y_var = vars_to_disp[j+1]
                x_bin_centers = param_display_vals[x_var][5]
                y_bin_centers = param_display_vals[y_var][5]

                x_scaling = param_display_vals[x_var][6]
                y_scaling = param_display_vals[y_var][6]

                scaled_x_bin_centers = np.array(x_bin_centers) * x_scaling
                scaled_y_bin_centers = np.array(y_bin_centers) * y_scaling

                x_mesh, y_mesh = np.meshgrid(scaled_x_bin_centers, scaled_y_bin_centers)

                binned_n_visits = np.zeros(np.shape(x_mesh))

                x_bins = param_display_vals[x_var][4]
                y_bins = param_display_vals[y_var][4]
                #print 'len(scaled_x_bin_centers) = ' + str(len(scaled_x_bin_centers))
                #print 'len(scaled_y_bin_centers) = ' + str(len(scaled_y_bin_centers))
                #print 'len(x_bins) = ' + str(len(x_bins))
                #print 'len(y_bins) = ' + str(len(y_bins))
                #print 'np.shape(x_mesh) = ' + str(np.shape(x_mesh))


                #print np.shape(binned_n_visits)
                for k in range(len(n_visits_array)):
                    measured_x = measured_arrays[x_var][k]
                    measured_y = measured_arrays[y_var][k]
                    n_visits = n_visits_array[k]

                    x_bin_index = -1
                    y_bin_index = -1

                    for l in range(len(x_bins)):
                        if measured_x >=x_bins[l][0] and measured_x < x_bins[l][1]:
                            x_bin_index = l
                    for l in range(len(y_bins)):
                        if measured_y >=y_bins[l][0] and measured_y < y_bins[l][1]:
                            y_bin_index = l

                    binned_n_visits[y_bin_index,x_bin_index] = binned_n_visits[y_bin_index,x_bin_index] + n_visits


                x_scaling = param_display_vals[x_var][6]
                y_scaling = param_display_vals[y_var][6]

                scaled_x_bin_centers = np.array(x_bin_centers) * x_scaling
                scaled_y_bin_centers = np.array(y_bin_centers) * y_scaling

                x_mesh, y_mesh = np.meshgrid(scaled_x_bin_centers, scaled_y_bin_centers)

                max_visits = np.max(np.abs(binned_n_visits + 1.0))
                log_levels = logList(1.0, max_visits + 1.0 , n_levels)
                exp_levels = [max_visits * np.exp(power) for power in exp_level_powers]
                levels = exp_levels

                xlabel = display_labels[x_var] + '\n' + (disp_units[x_var] if disp_units[x_var] == '' else r'$($' + disp_units[x_var] + r'$)$')
                ylabel = display_labels[y_var] + '\n' + (disp_units[y_var] if disp_units[y_var] == '' else r'$($' + disp_units[y_var] + r'$)$')
                if bin_types[x_var] == 'log':
                    axarr[axarr_y_index][axarr_x_index].set_xscale('log')
                    axarr[axarr_y_index][axarr_x_index].tick_params(top = False, bottom = False, which = 'minor')
                    if (n_var_per_plot in ['both','Both','BOTH']):
                        axarr[-1][axarr_x_index].set_xscale('log')
                        axarr[-1][axarr_x_index].tick_params(top = False, bottom = False, which = 'minor')
                if bin_types[y_var] == 'log':
                    axarr[axarr_y_index][axarr_x_index].set_yscale('log')
                    axarr[axarr_y_index][axarr_x_index].tick_params(right = False, left= False, which = 'minor')
                    if (n_var_per_plot in ['both','Both','BOTH']):
                        axarr[axarr_y_index][0].set_yscale('log')
                        axarr[axarr_y_index][0].tick_params(right = False, left= False, which = 'minor')
                CS = axarr[axarr_y_index][axarr_x_index].contour(x_mesh, y_mesh, binned_n_visits + 1.0, levels = levels, norm = LogNorm())

                if (x_var in params_to_fit and y_var in params_to_fit):
                    peak_xs, peak_ys = fitMCMCResults(results_files, [x_var, y_var], measured_arrays = measured_arrays, n_visits_array = n_visits_array, param_ranges_to_fit = [ params_to_fit[x_var],params_to_fit[y_var]],
                                                      results_dir = results_dir, n_ignore = n_ignore, theta_shift = theta_shift, n_fitted_points = n_fitted_points, n_xy_points_to_fit = n_xy_points_to_fit, smallest_max_val_for_fit = smallest_max_val_to_be_on_curve)
                    print ('peak_xs = ' + str(peak_xs))
                    print ('peak_ys = ' + str(peak_ys))
                    peak_interpolator = interpolate.interp1d( *safeSortOneListByAnother(peak_xs, [peak_xs, peak_ys]), kind = 'linear')
                    #axarr[0][1].contour(x_mesh, y_mesh, binned_n_visits_interp((x_mesh, y_mesh)) + 1.0, 20, levels = levels, norm = LogNorm())

                    axarr[axarr_y_index][axarr_x_index].scatter(np.array(peak_xs) * x_scaling, np.array(peak_ys) * y_scaling, c = 'black', s = 2.0 * (scaled_x_bin_centers[1] - scaled_x_bin_centers[0]) )
                    axarr[axarr_y_index][axarr_x_index].plot([x for x in scaled_x_bin_centers if (x > min(peak_xs) and x < max(peak_xs))], [peak_interpolator(x / x_scaling) * y_scaling for x in scaled_x_bin_centers if (x > min(peak_xs) and x < max(peak_xs))], c = fit_line_color)

                #axarr[i][j-1].colorbar(CS,format = '%.6f')
                xticks = var_ticks[x_var]
                xticklabels = var_ticklabels[x_var]
                yticks = var_ticks[y_var]
                yticklabels = var_ticklabels[y_var]
                axarr[axarr_y_index][axarr_x_index].set_xticks(xticks)
                axarr[axarr_y_index][axarr_x_index].set_yticks(yticks)
                if (axarr_y_index == len(vars_to_disp) - 2 and n_var_per_plot is 2):
                    axarr[axarr_y_index][axarr_x_index].set_xlabel(xlabel, fontsize = label_size_scaling * (fig_size_unit))
                    axarr[axarr_y_index][axarr_x_index].set_xticklabels(xticklabels, fontsize = tick_label_size_scaling * fig_size_unit, rotation = 0)
                else:
                    axarr[axarr_y_index][axarr_x_index].set_xticklabels([])
                if (axarr_x_index == 0 and n_var_per_plot is 2) :
                    axarr[axarr_y_index][axarr_x_index].set_ylabel(ylabel, fontsize = label_size_scaling * (fig_size_unit))
                    axarr[axarr_y_index][axarr_x_index].set_yticklabels(yticklabels, fontsize = tick_label_size_scaling * fig_size_unit, rotation = 90)
                else:
                    axarr[axarr_y_index][axarr_x_index].set_yticklabels([])

                if n_var_per_plot is 'both':
                    axarr[-1][axarr_x_index].set_xticks(xticks)
                    axarr[-1][axarr_x_index].set_xlabel(xlabel, fontsize = label_size_scaling * (fig_size_unit))
                    axarr[-1][axarr_x_index].set_xticklabels(xticklabels, fontsize = tick_label_size_scaling * fig_size_unit, rotation = 0)
                    axarr[axarr_y_index][0].set_yticks(yticks)
                    axarr[axarr_y_index][0].set_ylabel(ylabel, fontsize = label_size_scaling * (fig_size_unit))
                    axarr[axarr_y_index][0].yaxis.set_label_coords(y_label_coord_scalings[0]* (fig_size_unit), y_label_coord_scalings[1])
                    axarr[axarr_y_index][0].set_yticklabels(yticklabels, fontsize = tick_label_size_scaling * fig_size_unit, rotation = 90)
                else:
                    axarr[axarr_y_index][axarr_x_index].set_xticklabels([])




                #if (axarr_y_index == 0) or (axarr_y_index == 1 and n_var_per_plot in ['both','Both','BOTH']):
                #    print 'we should be adding ylabel ' + ylabel + ' to plot ' + str([[axarr_x_index],[axarr_y_index]])
                #    axarr[-1][axarr_y_index].set_ylabel(ylabel + str([axarr_x_index, axarr_y_index]), fontsize = 2.0 * (fig_size_unit) * n_y_contour_plots)
                #    axarr[-1][axarr_y_index].set_yticklabels(yticks, fontsize = 8.0 * fig_size_unit * n_y_contour_plots / (n_x_ticks), rotation = 45)
                #else:
                #    axarr[axarr_x_index][axarr_y_index].set_yticklabels([])
                axarr[axarr_y_index][axarr_x_index].set_xlim([min(scaled_x_bin_centers), max(scaled_x_bin_centers)])
                axarr[axarr_y_index][axarr_x_index].set_ylim([min(scaled_y_bin_centers), max(scaled_y_bin_centers)])
                if n_var_per_plot is 'both':
                    axarr[-1][axarr_x_index].set_xlim([min(scaled_x_bin_centers), max(scaled_x_bin_centers)])
                    axarr[axarr_y_index][0].set_ylim([min(scaled_y_bin_centers), max(scaled_y_bin_centers)])


        #for i in range(total_n_x_plots):
        #    if n_var_per_plot is 2:
        #        x_lims = axarr[-1][i].get_xlim()
        #    else:
        #        x_lims = axarr[-2][i].get_xlim()
        #    for j in range(total_n_y_plots):
        #        axarr[j][i].set_xlim(x_lims)
        #        if j != total_n_y_plots - 1:
        #            axarr[j][i].tick_params(axis='x', top = False, bottom = False, which = 'both', labelbottom = 'off')
        #        else:
        #            axarr[j][i].tick_params(axis = 'x', top = False,  which = 'both')

        #for i in range(total_n_y_plots):
        #    if n_var_per_plot is 2:
        #        y_lims = axarr[i][0].get_ylim()
        #    else:
        #        y_lims = axarr[i][1].get_ylim()
        #    for j in range(total_n_x_plots):
        #        axarr[i][j].set_ylim(y_lims)
        #        if j != 0:
        #            axarr[i][j].tick_params(axis = 'y', right = False, left = False, which = 'both', labelleft = 'off')
        #        else:
        #            axarr[i][j].tick_params(axis = 'y', right = False,  which = 'both')

        f.subplots_adjust(hspace = 0.0, wspace = 0.0)

    if n_var_per_plot == 1 or n_var_per_plot in ['both','Both','BOTH']:

        if n_var_per_plot == 1:
            n_x_plots = int( math.ceil( math.sqrt(len(vars_to_disp)) ) )
            n_y_plots = int( math.ceil( len(vars_to_disp) / float(n_x_plots) ) )
            f, axarr = plt.subplots(n_x_plots, n_y_plots, figsize = (fig_size_unit * n_x_plots, fig_size_unit * n_y_plots), squeeze = False)

        max_n_visits = 0
        for i in range(len(vars_to_disp)):
            var = vars_to_disp[i]
            xlabel = display_labels[var]
            #x_bin_centers = x_bin_centers + [param_display_vals[var][5]]
            x_bin_centers = param_display_vals[var][5]
            binned_n_visits = np.zeros((len(x_bin_centers)))
            for k in range(len(n_visits_array)):
                measured_x = measured_arrays[var][k]
                x_bins = param_display_vals[var][4]
                n_visits = n_visits_array[k]
                x_bin_index = -1
                for j in range(len(x_bins)):
                    if measured_x >=x_bins[j][0] and measured_x < x_bins[j][1]:
                        x_bin_index = j

                binned_n_visits[x_bin_index] = binned_n_visits[x_bin_index] + n_visits
            var_display_scaling = var_display_scalings[var]
            disp_x_bin_centers = np.array(x_bin_centers) * var_display_scaling
            if n_var_per_plot == 1:
               axarr[i % n_x_plots][i // n_x_plots].plot(disp_x_bin_centers, binned_n_visits)
               axarr[i % n_x_plots][i // n_x_plots].set_xlabel(xlabel, fontsize = label_size_scaling * (fig_size_unit))
            else:
                x_index = i + 1
                y_index = i - 1
                if x_index < len(vars_to_disp):
                    axarr[-1][x_index].plot(disp_x_bin_centers, binned_n_visits)
                if y_index >=  0:
                    axarr[y_index][0].plot(binned_n_visits, disp_x_bin_centers)
                max_n_visits = max(max_n_visits, max(binned_n_visits))


        print ('max_n_visits = ' + str(max_n_visits))
        if n_var_per_plot is 1:
            f.text(0.04, 0.5, 'Number of Chain Steps', va='center', rotation='vertical')

        else :
            n_visits_disp_max = max_n_visits
            axarr[-1][0].set_ylim([0.0, n_visits_disp_max])
            axarr[-1][0].set_xlim([0.0, n_visits_disp_max])
            bin_fractions = [1.0/6.0, 2.0/6.0, 3.0/6.0, 4.0/6.0, 5.0/6.0]
            bin_fractions = [1.0/6.0, 5.0/6.0 ]
            axarr[-1][0].set_xticks([int(n_visits_disp_max * i) for i in bin_fractions])
            axarr[-1][0].set_yticks([int(n_visits_disp_max * i) for i in bin_fractions])
            binned_xticks = [int(n_visits_disp_max * i ) for i in bin_fractions]
            binned_yticks = [int(n_visits_disp_max * i ) for i in bin_fractions]
            binned_xticklabels = [int(n_visits_disp_max * i / 1000.0) for i in bin_fractions]
            binned_yticklabels = [int(n_visits_disp_max * i / 1000.0) for i in bin_fractions]
            axarr[-1][0].set_xticklabels(binned_xticklabels, fontsize = tick_label_size_scaling * fig_size_unit, rotation = 0)
            axarr[-1][0].set_yticklabels(binned_yticklabels, fontsize = tick_label_size_scaling * fig_size_unit, rotation = 90)

            axarr[-1][0].set_xlabel(bin_label, fontsize = label_size_scaling * (fig_size_unit))
            axarr[-1][0].set_ylabel(bin_label, fontsize = label_size_scaling * (fig_size_unit))
            axarr[-1][0].yaxis.set_label_coords(y_label_coord_scalings[0] * (fig_size_unit), y_label_coord_scalings[1])


            for i in range(1, total_n_x_plots):
                axarr[-1][i].set_ylim([0.0, n_visits_disp_max])
                axarr[-1][i].set_yticks(binned_yticks)
                axarr[-1][i].set_yticklabels([])
                axarr[i-1][0].set_xlim([0.0, n_visits_disp_max])
                axarr[i-1][0].set_xticks(binned_xticks)
                axarr[i-1][0].set_xticklabels([])

    print ('save_fig = '  + str(save_fig) )
    if save_fig:
        if file_name is '': file_name = extra_save_str + 'MCMC_results_plots_vars_' + '_'.join([str(elem) for elem in vars_to_disp]) + '_' + str(n_var_per_plot) + '_per_plot' + '.pdf'
        print ('saving figure to ' + plot_dir + file_name)
        plt.savefig(plot_dir + file_name)

    if show:
        plt.show()
    else:
        plt.close('all')



def fitCurveToMCMCResults(results_files,  poly_order = 2, fitting_funct = 'poly', ind_var = 'rs', dep_var = 'M',
                              results_dir = '/Users/sashabrownsberger/Documents/Harvard/physics/randall/MCMCOutputs/', n_ignore = 0, n_fit_points = 200,
                              show_fit = 1, ranges_to_fit = ['all', 'all']):
    var_indeces = getVarIndeces()
    n_visits_index = 40
    likelihood_index = 41

    var_cyclicity = getVarCyclicity()
    n_bins_array = getNBinsArray()

    vars_to_read_in = [ind_var, dep_var]
    for var in vars_to_read_in:
        if not(var in var_indeces.keys()):
            print ('variable ' + str(var) + ' not a variable. Ending. ')
            return 0

    measured_arrays, n_visits_array, likelihood_array = readInMCMCResults(results_files, vars_to_load = vars_to_read_in,
                                                                          results_dir = results_dir, n_ignore = n_ignore )
    print (np.shape(n_visits_array))
    print (measured_arrays[vars_to_read_in[0]])
    x_vals, full_y_vals = binData(measured_arrays[var_to_fit], n_visits_array, n_bins = n_fit_points)
    y_vals = [full_y_vals[0][i] * full_y_vals[2][i] for i in range(len(full_y_vals[0]))]
    if ranges_to_fit[0] is 'all':
        x_vals = x_vals[:]
    else:
        x_vals = x_vals[range_to_fit[0]:range_to_fit[1]]
    if ranges_to_fit[1] is 'all':
        y_vals = y_vals[:]
    else:
        y_vals = y_vals[range_to_fit[0]:range_to_fit[1]]

    if fitting_funct in ['poly', 'polynomial', 'p']:
        fitted_curve = np.polyfit(x_vals, y_vals, poly_order)
        fitted_funct = np.poly1d(fitted_curve)
    elif fitting_funct in ['power', 'power_law']:
        fitting_funct = lambda xs, shift, power: shift + np.array(xs) ** power
        fitted_curve = optimize.curve_fit(fitting_funct, x_vals, y_vals, p0 = [0.0, 1.0], sigma = 1.0 / np.sqrt(y_vals))
        fitting_funct = lambda xs: fitting_funct(xs, *fitted_curve)
    print ('fitting_curve = ' + str(fitting_curve))


    if function_to_fit == 'gauss':
        function_to_fit = lambda x, A, mu, sigma, shift: A * np.exp(-(x - mu)**2.0 /(2.0 * sigma ** 2.0)) + shift
    #Measure initial best guess fit from data
    if guess_params is 'single_hump':
        #function is of form: A * exp (-(x-mu)**2 / (2.0 * sigma ** 2.0))
        A = np.max(y_vals) * 0.9
        mu = x_vals[np.argmax(y_vals)]
        sigma = abs(mu - x_vals[np.argmin((np.array(y_vals) - A / 2.0) ** 2.0)])
        guess_params = [A, mu, sigma]

    elif guess_params is 'single_hump_shift':
        #function is of form: A * exp (-(x-mu)**2 / (2.0 * sigma ** 2.0)) + shift
        A = np.max(y_vals) * 0.9
        mu = x_vals[np.argmax(y_vals)]
        sigma = abs(mu - x_vals[np.argmin((np.array(y_vals) - A / 2.0) ** 2.0)])
        shift = 0.0
        guess_params = [A, mu, sigma, shift]
    results = optimize.curve_fit(fitting_function, np.array(x_vals), np.array(y_vals), guess_params, np.zeros(len(x_vals)) + 1.0)
    if show_fit:
        plt.scatter(x_vals, y_vals)
        plt.plot(sorted(x_vals), fitting_function(np.array(sorted(x_vals)), *(results[0])), c = 'r')
        plt.xlabel(var_to_fit)
        plt.ylabel('n steps')
        plt.show()
    return results
