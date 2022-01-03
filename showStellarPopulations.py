import math 
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from logList import logList
from ObservedGalaxyStarData import ObservedGalaxyStarData

def showStellarPopulations(galaxy,pops_to_show,pop_selection_method = 'none',x_lims_arcmin = [-60,60], y_lims_arcmin = [-60,60], disp_units = 'arcmin', show_fig = 1, save_fig = 0):
    plot_save_dir = '/Users/sasha/Documents/Harvard/physics/randall/plots/'
    f, axarray = plt.subplots(2, 2,figsize = (10,10))
    for i in range(len(pops_to_show)):
        pop = pops_to_show[i]
        print 'pop = ' + pop 
        colormap = 'b' if pop == 'MP' else 'g' if pop == 'IM' else 'r' if pop == 'MR' else 'black'
        if pop == 'all' or pop == 'none':
            star_data = ObservedGalaxyStarData([galaxy,pop], pop_selection_method = 'none')
        else: 
            star_data = ObservedGalaxyStarData([galaxy,pop], pop_selection_method = pop_selection_method)
        corr_ra_one_pop = star_data.corrRa
        corr_dec_one_pop = star_data.corrDec
        if disp_units == 'arcmin':
            axarray[i / 2, i %2].scatter(corr_ra_one_pop * 60.0, corr_dec_one_pop * 60.0, s=2.0, color=colormap)
            axarray[i / 2, i % 2].set_xlim(x_lims_arcmin)
            axarray[i / 2, i % 2].set_ylim(y_lims_arcmin)            
        elif disp_units == 'degrees':
            axarray[i / 2, i % 2].scatter(corr_ra_one_pop,corr_dec_one_pop)
            axarray[i / 2, i % 2].set_xlim([x_lim / 60.0 for x_lim in x_lims_arcmin])
            axarray[i / 2, i % 2].set_ylim([y_lim / 60.0 for y_lim in y_lims_arcmin])  
        axarray[i / 2, i % 2].set_title(pop + ' stars in ' + galaxy + ' galaxy')
        axarray[i / 2, i % 2].set_xlabel('RA distance from dSph center, in arcmin')
        axarray[i / 2, i % 2].set_ylabel('Dec distance from dSph center, in arcmin')
    if show_fig:
        plt.show()
    if save_fig:
        fig_name = galaxy + '_' 
        for pop in pops_to_show:
            fig_name = fig_name + pop + '_'
        fig_name = fig_name + 'populations'
        f.savefig(plot_save_dir + fig_name + '.pdf',bbox_inches = 'tight') 
    return 0        
    

