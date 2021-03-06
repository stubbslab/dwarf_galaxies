#Do a sequence of the MCMC Series algorithms, with varying run conditions.
#You can change various parameters (present set listed in commented line right above parameter_serires assignment), and have program execute in sequence
#This script can be run from command line without having to be in python environment (bash$ python doMCMCSeriesAllPops.py)
#This method tries to be as efficient as possible, reading in the various potentialFunctionArrays only once, as that can be a time consuming part of the process.  

from runSimulMCMCForDwarfGalaxyProfilesAllPops import runSimulMCMCForDwarfGalaxyProfilesAllPops
from ArtificialParamsStorer import ArtificialParamsStorer
from readInPotentialInterpolator import readInPotentialInterpolator 

if __name__ == '__main__':
    #parameters at present are: [halo_type, disk_type, galaxy, number of iterations, population selection method, withDisk, start parameter index, outer step] 
    parameter_series = [ ['nfw','sech_disk','fornax',20,'metal_point',0,[0,1,2,],2.0],
                         ['nfw','sech_disk','fornax',20,'metal_point',0,[39,40,41],2.0]
                       ]
    saverun = 1
    show_best_fit = 0
    artificial = 0
    artificial_param_indeces = [0,1,2,3]
    if not artificial:
        artificial_param_indeces = [0]
    art_params_storer = ArtificialParamsStorer()
    artificial_storers = [art_params_storer.getStorer(index) for index in artificial_param_indeces]

    halo_index = 0
    disk_index = 1
    galaxy_index = 2
    nIter_index = 3
    pop_select_index = 4
    withDisk_index = 5
    start_indeces_index = 6
    step_index = 7
    unique_halo_types = []
    unique_disk_types = []
    #Create all of the potential files used in the MCMC algorithms.
    # That way, we don't need to reread in a potential file every time it is used.
    for parameter_set in parameter_series:
        if parameter_set[halo_index] not in unique_halo_types:
            print 'assigning halo type ' + str(parameter_set[halo_index]) + ' to unique_halo_types'
            unique_halo_types = unique_halo_types + [parameter_set[halo_index]]
        if parameter_set[disk_index] not in unique_disk_types:
            print 'assigning disk type ' + str(parameter_set[disk_index]) + ' to unique_disk_types'
            unique_disk_types = unique_disk_types + [parameter_set[disk_index]]
    if artificial:
        for storer in artificial_storers:
            art_el = storer.el
            art_lam = storer.lam
            if art_el not in unique_els:
                print 'assigning el ' + str(art_el) + ' to unique_els'
                unique_els = unique_els + [art_el]
            if parameter_set[lam_index] not in unique_lams:
                print 'assigning lam ' + str(art_lam) + ' to unique_lams'            
                unique_lams = unique_lams + [art_lam]
    #print 'unique_els = '
    #print unique_els
    #print 'unique_lams = '
    #print unique_lams
    halo_funct_array = {}
    disk_funct_array = {}
    for unique_halo_type in unique_halo_types:
        print 'Loading potential interpolator for halo type ' + unique_halo_type
        halo_funct_array[unique_halo_type] = readInPotentialInterpolator(unique_halo_type)
    for unique_disk_type in unique_disk_types:
        print 'Loading potential interpolator for disk type ' + unique_disk_type
        disk_funct_array[unique_disk_type] = readInPotentialInterpolator(unique_disk_type)

    #Now actually do each MCMC run.  This can take a while.  
    for i in range(len(parameter_series)):
        for j in range(len(artificial_param_indeces)): 
            print 'Starting MCMC series ' + str(i) +', ' + str(j)
            if artificial: 
                runSimulMCMCForDwarfGalaxyProfilesAllPops(
                    parameter_series[i][el_index], parameter_series[i][lam_index], parameter_series[i][galaxy_index],
                    halo_funct = halo_funct_array[parameter_series[i][el_index]], disk_funct = disk_funct_array[parameter_series[i][lam_index]],
                    withDisk=parameter_series[i][withDisk_index], nIterations = parameter_series[i][nIter_index], saverun = saverun, show_best_fit = show_best_fit,
                    pop_selection_method = parameter_series[i][pop_select_index], start_param_index = parameter_series[i][start_index_index], outer_step = parameter_series[i][step_index],
                    use_artificial_data = artificial, generation_params_storer_index = artificial_param_indeces[j], generation_disk_funct = disk_funct_array[artificial_storers[j].lam], generation_halo_funct = halo_funct_array[artificial_storers[j].el])
            else:
                runSimulMCMCForDwarfGalaxyProfilesAllPops(
                    parameter_series[i][galaxy_index],
                    halo_funct = halo_funct_array[parameter_series[i][halo_index]], disk_funct = disk_funct_array[parameter_series[i][disk_index]],
                    withDisk=parameter_series[i][withDisk_index], nIterations = parameter_series[i][nIter_index], saverun = saverun, show_best_fit = show_best_fit,
                    pop_selection_method = parameter_series[i][pop_select_index], start_param_indeces = parameter_series[i][start_indeces_index], outer_step = parameter_series[i][step_index],
                    use_artificial_data = artificial, generation_params_storer_index = 'None', generation_disk_funct = 0, generation_halo_funct = 0,
                    h_type = parameter_series[i][halo_index], d_type = parameter_series[i][disk_index])
