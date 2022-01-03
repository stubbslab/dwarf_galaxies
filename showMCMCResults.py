
import csv
import numpy as np
import matplotlib
import matplotlib.cm as cm
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import scipy.integrate as integrate
from logList import logList
from matplotlib.colors import LogNorm
from ModulateCyclicParam import ModulateCyclicParam
import math
import time
from ComputationalArchive import ComputationalArchive
import scipy.optimize as optimization

#function object must looks like: {[(var_name1, var_name2, ..., var_namen), string_descriptor]: function using variables }
def readInMCMCResults(results_files, functions_to_compute = [], vars_to_load = ['M', 'el', 'rs', 'phi', 'theta', 'h_x_center', 'h_y_center', 'h_z_center',
                                                     'lam', 'zeta', 'zeta', 'eps', 'a', 'b', 'd_x_center', 'd_y_center', 'd_z_center'], 
                      results_dir = '/Users/sasha/Documents/Harvard/physics/randall/MCMCOutputs/', n_ignore = 0):

    possible_operations = {'+','-','*','/'}
    var_indeces = {'M' : 1, 'el':4, 'rs': 5, 'phi': 6, 'theta':7, 'h_x_center':11, 'h_y_center':12, 'h_z_center':13,
                   'lam':15, 'zeta':16, 'eps':17, 'a': 18, 'b':19, 'd_x_center':23, 'd_y_center':24, 'd_z_center':25}
    measured_arrays = {}
    likelihood_array = []
    n_visits_array = []
    for var in vars_to_load: 
        measured_arrays[var] = []
    for function_list in functions_to_compute:
        measured_arrays[function_list[0]] = [] 

    n_visits_index = 26
    likelihood_index = 27
    
    for results_file in results_files:
        data_file = results_dir + results_file
        print 'Reading in ' + data_file 
        skip_first_line = True
        n_visits = 0
        with open(data_file) as csvfile:
            myreader = csv.reader(csvfile,delimiter=',')
            row_number = 0
            new_measured_arrays = {}
            new_likelihood_array = []
            new_n_visits_array = []
            for var in measured_arrays.keys():
                new_measured_arrays[var] = []
            for row in myreader:
                if skip_first_line:
                    skip_first_line = False
                elif n_visits < n_ignore:
                    row_number = row_number + 1
                    n_visits = n_visits + float(row[n_visits_index])
                else:
                    row_number = row_number + 1
                    for var in vars_to_load:
                        new_measured_arrays[var] = new_measured_arrays[var] + [float( row[ var_indeces[var] ] )]
                    new_likelihood_array = new_likelihood_array + [float(row[likelihood_index])]
                    new_n_visits_array = new_n_visits_array + [float(row[n_visits_index])]
            likelihood_array = likelihood_array + new_likelihood_array
            n_visits_array = n_visits_array + new_n_visits_array
            n_visits = n_visits + float(row[n_visits_index])
            for var in measured_arrays.keys():
                measured_arrays[var] = measured_arrays[var] + new_measured_arrays[var]
    for function_list in functions_to_compute:
        measured_arrays[function_list[0]] = [ function_list[2]( *[measured_arrays[var][i] for var in function_list[1]]) for i in range(len(measured_arrays[vars_to_load[0]])) ]
    

    return [measured_arrays, n_visits_array, likelihood_array ]

def measureMCMCResults(results_files, var_to_fit, function_to_fit, guess_params, functions_to_fit = [], 
                       results_dir = '/Users/sasha/Documents/Harvard/physics/randall/MCMCOutputs/', n_ignore = 0):
    var_indeces = {'M' : 1, 'el':4, 'rs': 5, 'phi': 6, 'theta':7, 'h_x_center':11, 'h_y_center':12, 'h_z_center':13, 'lam':15, 'zeta':16, 'eps':17, 'a': 18, 'b':19, 'd_x_center':23, 'd_y_center':24, 'd_z_center':25}
    n_visits_index = 26
    likelihood_index = 27
    
    var_cyclicity = {'M' : 'none', 'el':'none', 'rs':'none', 'phi':'none', 'theta':'none', 'h_x_center':'none', 'h_y_center':'none', 'h_z_center':'none', 'lam':'none', 'zeta':'none', 'eps':'none', 'a':'none', 'b':'none', 'd_x_center':'none', 'd_y_center':'none', 'd_z_center':'none' }

    vars_to_read_in = [var_to_fit]
    for funct in functions_to_display: 
        vars_to_read_in = vars_to_read_in + funct[1] 
    vars_to_read_in = list(set(vars_to_read_in))
    vars_to_read_in = [var for var in vars_to_read_in if var in var_indeces.keys()]
    print 'vars_to_read_in = ' + str(vars_to_read_in)

    measured_arrays, n_visits_array, likelihood_array = readInMCMCResults(results_files, functions_to_compute = functions_to_fit, vars_to_load = vars_to_read_in, 
                                        results_dir = results_dir, n_ignore = n_ignore )

    if function_to_fit == 'gauss':
        function_to_fit = lambda x, A, mu, sigma, shift: A * np.exp(-(x - mu)**2.0 /(2.0 * sigma ** 2.0))  
    results = optimization.curve_fit(function_to_fit, np.array(measured_arrays[var_to_fit]), np.array(n_visits_arrays), guess_params, np.zeros(len(measured_arrays[var_to_fit])) + 1.0)
    return results 
    
        
def showMCMCResults(results_files, vars_to_disp, functions_to_display = [], 
                    results_dir = '/Users/sasha/Documents/Harvard/physics/randall/MCMCOutputs/', file_name = '',
                    n_bins = [], n_ignore = 0, theta_shift = 0.0, save_fig = 0, n_var_per_plot = 1, fig_size_unit = 5.0):
    
    comp_archive = ComputationalArchive() 
    plot_dir = comp_archive.getPlotDir()
    
    var_indeces = {'M' : 1, 'el':4, 'rs': 5, 'phi': 6, 'theta':7, 'h_x_center':11, 'h_y_center':12, 'h_z_center':13, 'lam':15, 'zeta':16, 'eps':17, 'a': 18, 'b':19, 'd_x_center':23, 'd_y_center':24, 'd_z_center':25}
    n_visits_index = 26
    likelihood_index = 27
    
    var_cyclicity = {'M' : 'none', 'el':'none', 'rs':'none', 'phi':'none', 'theta':'none', 'h_x_center':'none', 'h_y_center':'none', 'h_z_center':'none', 'lam':'none', 'zeta':'none', 'eps':'none', 'a':'none', 'b':'none', 'd_x_center':'none', 'd_y_center':'none', 'd_z_center':'none' }
    disp_units = {'M' : ' (sol M)', 'el':'', 'rs':' (pc)', 'phi':'', 'theta':'', 'h_x_center':'rs', 'h_y_center':'rs', 'h_z_center':'rs', 'lam':'', 'zeta':'', 'eps':'','a':'','b':'', 'd_x_center':'rs', 'd_y_center':'rs', 'd_z_center':'rs'}
    bin_type = {'M' : ' lin', 'el':'log', 'rs':'lin', 'phi':'lin', 'theta':'lin', 'h_x_center':'lin', 'h_y_center':'lin', 'h_z_center':'lin', 'lam':'lin', 'zeta':'lin', 'eps':'lin','a':'lin','b':'lin', 'd_x_center':'lin', 'd_y_center':'lin', 'd_z_center':'lin' }

    vars_to_read_in = vars_to_disp
    for funct in functions_to_display: 
        vars_to_read_in = vars_to_read_in + funct[1] 
    vars_to_read_in = list(set(vars_to_read_in))
    vars_to_read_in = [var for var in vars_to_read_in if var in var_indeces.keys()]
    print 'vars_to_read_in = ' + str(vars_to_read_in) 
    

    n_bins_array = {'M':200, 'el':200, 'rs':500, 'phi':100, 'theta':100, 'h_x_center':200, 'h_y_center':200, 'h_z_center':200, 'lam':100, 'zeta':200, 'eps':160, 'a':100,'b':100, 'd_x_center':200, 'd_y_center':200, 'd_z_center':200}
    param_display_vals = {'M':[0.0,0.0,n_bins_array['M']],
                          'el':[0.1,10.0,n_bins_array['el']],
                          'rs':[0.0,0.0,n_bins_array['rs']],
                          'phi':[0.0,math.pi, n_bins_array['phi']],
                          'theta':[-math.pi / 2.0 , math.pi / 2.0, n_bins_array['theta']],
                          'h_x_center':[-3.0, 3.0, n_bins_array['h_x_center']],
                          'h_y_center':[-3.0, 3.0, n_bins_array['h_y_center']],
                          'h_z_center':[-3.0, 3.0, n_bins_array['h_z_center']],
                          'lam':[0.05,0.15,n_bins_array['lam']],
                          'zeta':[0.0,2.0,n_bins_array['zeta']],
                          'eps':[0.0,1.0,n_bins_array['eps']],
                          'a':[0.0,math.pi ,n_bins_array['a']],
                          'b':[-math.pi / 2.0, math.pi / 2.0, n_bins_array['b']],
                          'd_x_center':[-3.0, 3.0, n_bins_array['d_x_center']],
                          'd_y_center':[-3.0, 3.0, n_bins_array['d_y_center']],
                          'd_z_center':[-3.0, 3.0, n_bins_array['d_z_center']]
                          }
    
    measured_arrays, n_visits_array, likelihood_array = readInMCMCResults(results_files, functions_to_compute = functions_to_display, vars_to_load = vars_to_read_in, 
                                        results_dir = results_dir, n_ignore = n_ignore )
    print measured_arrays.keys() 
    #print measured_arrays
    #for var in vars_to_disp:
    #    print 'measured_arrays[var] = '
    #    print measured_arrays[var]
    #print "len(measured_arrays['phi']) = "
    #print len(measured_arrays['phi'])
    for i in range(len(n_bins)):
        n_bins_array[vars_to_disp[i]] = n_bins[i]

    for funct in functions_to_display:
        var = funct[0]
        print 'var = ' + var
        display_vals = [param_display_vals[component_var] for component_var in funct[1]]
        funct_display_vals = [funct[2](*[val[0] for val in display_vals]),
                              funct[2](*[val[2] for val in display_vals]),
                              display_vals[0][2]]
        param_display_vals[var] = funct_display_vals
        bin_type[var] = 'lin'
        var_cyclicity[var] = 'none'
        disp_units[var] = 'none' 
        for component_var in funct[1]: 
            if bin_type[component_var] != 'lin': bin_type[var] = bin_type[component_var]
            if var_cyclicity[component_var] != 'none': var_cyclicity[var] = var_cyclicity[component_var] 
            if disp_units[component_var] != 'none': disp_units[var] = disp_units[component_var] 
    
    for var in vars_to_disp:
        if var == 'theta':
            uncorrected_thetas = measured_arrays[var]
            print 'max(uncorrected_thetas = ' + str(max(uncorrected_thetas))
            print 'min(uncorrected_thetas = ' + str(min(uncorrected_thetas))
            measured_arrays[var] = [uncorrected_theta + theta_shift for uncorrected_theta in uncorrected_thetas]
        if var_cyclicity[var] != 'none':
            cyclicity = var_cyclicity[var]
            raw_measured_array = measured_arrays[var]
            measured_arrays[var] = [ ModulateCyclicParam(raw_var,cyclicity, shift = - math.pi) for raw_var in raw_measured_array ]
            if var == 'theta':
                print 'raw_measured_array[10:30] = '
                print raw_measured_array[10:30]
                print 'measured_array["theta"][10:30] = '
                print measured_arrays[var][10:30]
        
    total_visits = sum(n_visits_array)

    for var in vars_to_disp:
        if var[0:2] == 'rs':
            measured_rs_array = measured_arrays[var]
            param_display_vals[var] = [0.0,
                                       2000.0,
                                       n_bins_array[var]]
            #print 'measured_rs_array = '
            #print measured_rs_array 
            #param_display_vals[var] = [min(measured_rs_array) - 20.0,
            #                            max(measured_rs_array) + 20.0,
            #                            n_bins_array[var]]
        if var[0:1] == 'M':
            measured_M_array = measured_arrays[var]  
            param_display_vals[var] = [0.0 * 10**8,
                                       5.0 * 10**8,
                                       n_bins_array[var]]
    #print 'measured_rs_array = '
    #print measured_rs_array 
    
    #each param is min, max, nbins, step, bins, bin_centers
    for var in vars_to_disp:
        var_min = param_display_vals[var][0]
        var_max = param_display_vals[var][1]
        var_nbins = param_display_vals[var][2]
        var_step = (var_max - var_min) / var_nbins
        if bin_type[var] == 'lin':
            var_bins = [[var_min + var_step*i, var_min + var_step *(i+1)] for i in range(var_nbins)]
        elif bin_type[var] == 'log':
            log_spaced_borders = logList(var_min,var_max,var_nbins) 
            var_bins = [ [log_spaced_borders[i],log_spaced_borders[i+1]] for i in range(len(log_spaced_borders)-1) ]
        #print 'var_bins = '
        #print var_bins 
        var_bin_centers = [(var_bins[i][1] + var_bins[i][0])/ 2.0 if bin_type[var] == 'lin'
                           else 10.0 ** ((math.log10(var_bins[i][1]) + math.log10(var_bins[i][0]))/ 2.0) if bin_type[var]=='log'
                           else 0.0
                           for i in range(len(var_bins)) ]
        param_display_vals[var] = param_display_vals[var] + [var_step, var_bins, var_bin_centers]


    plt.rcParams.update({'font.size': 6})
    if  n_var_per_plot == 2:
        
        n_x_plots = len(vars_to_disp) - 1
        n_y_plots = len(vars_to_disp) - 1
        f, axarr = plt.subplots(n_x_plots, n_y_plots,
                                figsize = (fig_size_unit * n_x_plots, fig_size_unit * n_y_plots),
                                squeeze = False) #, sharex = True, sharey = True)
        print 'range(len(vars_to_disp)-1) = ' + str(range(len(vars_to_disp)-1) ) 
        for i in range(len(vars_to_disp) - 1):
            for j in [elem + 1 for elem in range(i+1)]:
                print 'i = ' + str(i)
                print 'j = ' + str(j) 
                x_var = vars_to_disp[len(vars_to_disp) - i - 2]
                y_var = vars_to_disp[len(vars_to_disp) - j]
                x_bin_centers = param_display_vals[x_var][5] 
                y_bin_centers = param_display_vals[y_var][5]
                binned_n_visits = np.zeros(( len(x_bin_centers), len(y_bin_centers) ))
                #print np.shape(binned_n_visits) 
                for k in range(len(n_visits_array)):
                    measured_x = measured_arrays[x_var][k]
                    x_bins = param_display_vals[x_var][4]
            
                    measured_y = measured_arrays[y_var][k]
                    y_bins = param_display_vals[y_var][4]

                    n_visits = n_visits_array[k]
                                   
                    x_bin_index = -1
                    y_bin_index = -1

                    for l in range(len(x_bins)):
                        if measured_x >=x_bins[l][0] and measured_x < x_bins[l][1]: 
                            x_bin_index = l
                    for l in range(len(y_bins)):
                        if measured_y >=y_bins[l][0] and measured_y < y_bins[l][1]: 
                            y_bin_index = l

                    binned_n_visits[x_bin_index,y_bin_index] = binned_n_visits[x_bin_index,y_bin_index] + n_visits
        
                #print 'np.max(binned_n_visits): ' + str(np.max(binned_n_visits))
                #print 'binned_n_visits[10:20,22:32] = '
                #print binned_n_visits[10:20,22:32]

                x_mesh, y_mesh = np.meshgrid(y_bin_centers,x_bin_centers)
        
                plt.rcParams['xtick.direction'] = 'out'
                plt.rcParams['ytick.direction'] = 'out'

                log_levels = logList(1.0,np.max(np.abs(binned_n_visits + 1.0)),30)
                #print log_levels
                #print 'np.shape(x_mesh) = '
                #print np.shape(x_mesh)
                #print 'np.shape(y_mesh) = '
                #print np.shape(y_mesh)
                #print 'np.shape(binned_n_visits): '
                #print np.shape(binned_n_visits)
                print 'x_var = ' + x_var
                print 'y_var = ' + y_var
                if bin_type[x_var] == 'log':
                    print 'setting log for x_var'
                    axarr[i][j-1].set_yscale('log')
                if bin_type[y_var] == 'log':
                    print 'setting log for y_var'
                    axarr[i][j-1].set_xscale('log')  
                CS = axarr[i][j-1].contour(x_mesh, y_mesh, binned_n_visits + 1.0, 20, levels = log_levels, norm = LogNorm())
                #axarr[i][j-1].colorbar(CS,format = '%.6f')
                if i == len(vars_to_disp) - 2:
                    print 'Printing xlabel' 
                    axarr[i][j-1].set_xlabel(y_var)
                if j-1 == 0:
                    print 'Printing ylabel'
                    axarr[i][j-1].set_ylabel(x_var)

        f.suptitle('MCMC chain outputs')
        #print 'vars_to_disp[1] = ' + str(vars_to_disp[1])
        #print 'disp_units[vars_to_disp[1]] = ' + str(disp_units[vars_to_disp[1]]) 

    if n_var_per_plot == 1:
        n_x_plots = int( math.ceil( math.sqrt(len(vars_to_disp)) ) )
        n_y_plots = int( math.ceil( len(vars_to_disp) / float(n_x_plots) ) )
        f, axarr = plt.subplots(n_x_plots, n_y_plots, figsize = (fig_size_unit * n_x_plots, fig_size_unit * n_y_plots),squeeze = False)
        #x_bin_centers = []
        #binned_n_visits = []
        for i in range(len(vars_to_disp)):
            var = vars_to_disp[i]
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
            #binned_n_visits = binned_n_visits + [individual_binned_n_visits]
            axarr[i % n_x_plots][i / n_x_plots].plot(x_bin_centers, binned_n_visits)
            axarr[i % n_x_plots][i / n_x_plots].set_xlabel(var + disp_units[var]) 

        #print 'np.shape(binned_n_visits): '
        #print np.shape(binned_n_visits)
        #print 'np.max(binned_n_visits): ' + str(np.max(binned_n_visits))
        #print 'binned_n_visits = '
        #print binned_n_visits

        plt.suptitle('MCMC chain single param results')
        #plt.xlabel(vars_to_disp[0] + disp_units[vars_to_disp[0]])
        #plt.ylabel('binned number of visits')
        f.text(0.04, 0.5, 'Number of Chain Steps', va='center', rotation='vertical')
        #print x_bin_centers
        #print binned_n_visits
        #plt.plot(x_bin_centers, binned_n_visits)
        #if bin_type[vars_to_disp[0]] == 'log':
        #    plt.xscale('log')

    print 'save_fig = '  + str(save_fig) 
    if save_fig:
        if file_name is '': file_name = 'MCMC_results_plots_vars_' + '_'.join([str(elem) for elem in vars_to_disp]) + '_' + str(n_var_per_plot) + '_per_plot' + '.png'
        print 'saving figure to ' + plot_dir + file_name
        plt.savefig(plot_dir + file_name) 
            
        
        
    plt.show()

        
