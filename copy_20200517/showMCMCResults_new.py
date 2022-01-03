import csv
import numpy as np
import matplotlib
import matplotlib.cm as cm
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import scipy.integrate as integrate
from logList import logList
from matplotlib.colors import LogNorm
import math


def showMCMCResults(results_files, vars_to_disp, results_dir = '/Users/sasha/Documents/Harvard/physics/randall/MCMCOutputs/' , n_bins = [], n_ignore = 0 ):

    var_indeces = {'rs': 5, 'phi': 6, 'eps':11, 'a': 12}
    var_cyclicity = {'rs':'none', 'phi':math.pi, 'eps':'none', 'a':math.pi}
    
    n_visits_index = 14
    likelihood_index = 15

    n_bins_array = {'rs':200, 'phi':100, 'eps':80, 'a':60}
    for i in range(len(n_bins)):
        n_bins_array[vars_to_disp[i]] = n_bins[i]
    print 'n_bins_array = '
    print n_bins_array

    measured_arrays = {}
    for var in vars_to_disp:
        measured_arrays[var] = []
    likelihood_array = []
    n_visits_array = []

    for results_file in results_files:
        data_file = results_dir + results_file

        skip_first_line = True

        with open(data_file) as csvfile:
            myreader = csv.reader(csvfile,delimiter=',')
            row_number = 0 
            for row in myreader:
                if skip_first_line:
                    skip_first_line = False
                elif row_number < n_ignore:
                    row_number = row_number + 1
                else:
                    row_number = row_number + 1
                    for var in vars_to_disp:
                        measured_arrays[var] = measured_arrays[var] + [float( row[ var_indeces[var] ] )]
                    likelihood_array = likelihood_array + [float(row[likelihood_index])]
                    n_visits_array = n_visits_array + [float(row[n_visits_index])]

    for var in vars_to_disp:
        if var_cyclicity[var] != 'none':
            raw_measured_array = measured_arrays[var]
            measured_arrays[var] = [( raw_var + min(raw_measured_array) + abs(min(raw_measured_array)) % (math.pi) ) % (math.pi) for raw_var in raw_measured_array]
            

        
    total_visits = sum(n_visits_array)
    

    
    #each param is min, max, nbins, step, bins, bin_centers
    param_display_vals = {'rs':[0.0,0.0,n_bins_array['rs']],
                          'phi':[0.0,math.pi, n_bins_array['phi']],
                          'eps':[0.0,1.0,n_bins_array['eps']],
                          'a':[0.0,math.pi,n_bins_array['a']]    }
    if 'rs' in vars_to_disp:
        measured_rs_array = measured_arrays['rs']  
        param_display_vals['rs'] = [min(measured_rs_array) - 20.0,
                                    max(measured_rs_array) + 20.0,
                                    n_bins_array['rs']]
    for var in vars_to_disp:
        var_min = param_display_vals[var][0]
        var_max = param_display_vals[var][1]
        var_nbins = param_display_vals[var][2]
        var_step = (var_max - var_min) / var_nbins
        var_bins = [[var_min + var_step*i, var_min + var_step *(i+1)] for i in range(var_nbins)]
        var_bin_centers = [(var_min + var_step*i + var_min + var_step * (i+1))/2.0 for i in range(var_nbins)]
        param_display_vals[var] = param_display_vals[var] + [var_step, var_bins, var_bin_centers]


    if len(vars_to_disp) == 2:
        binned_n_visits = np.zeros(( n_bins_array[vars_to_disp[0]], n_bins_array[vars_to_disp[1]] ))
        #print np.shape(binned_n_visits) 
        for i in range(len(n_visits_array)):
            measured_x = measured_arrays[vars_to_disp[0]][i]
            x_bins = param_display_vals[vars_to_disp[0]][4]
            
            measured_y = measured_arrays[vars_to_disp[1]][i]
            y_bins = param_display_vals[vars_to_disp[1]][4]

            n_visits = n_visits_array[i]
                                   
            x_bin_index = -1
            y_bin_index = -1

            for j in range(len(x_bins)):
                if measured_x >=x_bins[j][0] and measured_x < x_bins[j][1]: 
                    x_bin_index = j
            for j in range(len(y_bins)):
                if measured_y >=y_bins[j][0] and measured_y < y_bins[j][1]: 
                    y_bin_index = j

            binned_n_visits[x_bin_index,y_bin_index] = binned_n_visits[x_bin_index,y_bin_index] + n_visits
        
        print 'np.max(binned_n_visits): ' + str(np.max(binned_n_visits))
        #print 'binned_n_visits[10:20,22:32] = '
        #print binned_n_visits[10:20,22:32]

        x_bin_centers = param_display_vals[vars_to_disp[0]][5] 
        y_bin_centers = param_display_vals[vars_to_disp[1]][5]
        x_mesh, y_mesh = np.meshgrid(y_bin_centers,x_bin_centers)
        
        plt.rcParams['xtick.direction'] = 'out'
        plt.rcParams['ytick.direction'] = 'out'

        plt.figure(figsize=(11,9))

        log_levels = logList(1.0,np.max(np.abs(binned_n_visits + 1.0)),30)
        print log_levels
        print 'np.shape(x_mesh) = '
        print np.shape(x_mesh)
        print 'np.shape(y_mesh) = '
        print np.shape(y_mesh)
        print 'np.shape(binned_n_visits): '
        print np.shape(binned_n_visits)
        CS = plt.contour(x_mesh,y_mesh,binned_n_visits + 1.0,20, levels = log_levels, norm = LogNorm())
        plt.colorbar(CS,format = '%.6f')

        plt.title('MCMC chain')
        plt.xlabel(vars_to_disp[1])
        plt.ylabel(vars_to_disp[0])

        plt.show()
        #print 'shape x_mesh = '
        #print np.shape(x_mesh)
        #print 'shape y_mesh = '
        #print np.shape(y_mesh)
        #print x_bin_centers
        #print y_bin_centers

    if len(vars_to_disp) == 1:
        binned_n_visits = np.zeros(( n_bins_array[vars_to_disp[0]]))
        for i in range(len(n_visits_array)):
            measured_x = measured_arrays[vars_to_disp[0]][i]
            x_bins = param_display_vals[vars_to_disp[0]][4]

            n_visits = n_visits_array[i]
                                   
            x_bin_index = -1

            for j in range(len(x_bins)):
                if measured_x >=x_bins[j][0] and measured_x < x_bins[j][1]: 
                    x_bin_index = j

            binned_n_visits[x_bin_index] = binned_n_visits[x_bin_index] + n_visits        

        print 'np.shape(binned_n_visits): '
        print np.shape(binned_n_visits)
        print 'np.max(binned_n_visits): ' + str(np.max(binned_n_visits))
        print 'binned_n_visits = '
        print binned_n_visits

        x_bin_centers = param_display_vals[vars_to_disp[0]][5]
        
        plt.figure(figsize=(11,9))
        plt.title('MCMC chain')
        plt.xlabel(vars_to_disp[0])
        plt.ylabel('binned number of visits')
        plt.plot(x_bin_centers, binned_n_visits)

        plt.show()
