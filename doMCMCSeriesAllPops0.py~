#Do a sequence of the MCMC Series algorithms, with varying run conditions.
#You can change various parameters (present set listed in commented line right above parameter_serires assignment), and have program execute in sequence
#This script can be run from command line without having to be in python environment (bash$ python doMCMCSeriesAllPops.py)
#This method tries to be as efficient as possible, reading in the various potentialFunctionArrays only once, as that can be a time consuming part of the process.  

from runMCMCForDwarfGalaxyProfilesAllPops import runMCMCForDwarfGalaxyProfilesAllPops
from ArtificialParamsStorer import ArtificialParamsStorer
from readInPotentialInterpolator import readInPotentialInterpolator 

if __name__ == '__main__':
    #parameters at present are: [halo_type, disk_type, galaxy, number of iterations, population selection method, withDisk, start parameter index, outer step, halo edge on, disk edge on, fixed disk angle, fixed_params ] 
    parameter_series = [ ['cored','sech_disk','fornax',10,'metal_point',0,0,2.0,0,0,None,{}],
                         ['burkert','sech_disk','fornax',10,'metal_point',0,1,2.0,0,0,None,{}],
                         ['nfw','sech_disk','fornax',10,'metal_point',0,2,2.0,0,0,None,{}],
                         #['nfw','sech_disk','fornax',10000,'metal_point',0,3,2.0,0,0,None,{}],
                         #['nfw','sech_disk','fornax',10000,'metal_point',0,4,2.0,0,0,None,{}],
                         #['nfw','sech_disk','fornax',10000,'metal_point',0,5,2.0,0,0,None,{}],
                         #['nfw','sech_disk','fornax',10000,'metal_point',0,6,2.0,0,0,None,{}],
                         #['nfw','sech_disk','fornax',10000,'metal_point',0,7,2.0,0,0,None,{}],
                         #['nfw','sech_disk','fornax',10000,'metal_point',0,8,2.0,0,0,None,{}],
                         #['nfw','sech_disk','fornax',10000,'metal_point',0,9,2.0,0,0,None,{}],
                         #['nfw','sech_disk','fornax',10000,'metal_point',0,10,2.0,0,0,None,{}],
                         #['nfw','sech_disk','fornax',10000,'metal_point',0,11,2.0,0,0,None,{}],
                         #['nfw','sech_disk','fornax',10000,'metal_point',0,12,2.0,0,0,None,{}],
                         #['nfw','sech_disk','fornax',10000,'metal_point',0,13,2.0,0,0,None,{}],
                         #['nfw','sech_disk','fornax',10000,'metal_point',0,14,2.0,0,0,None,{}],
                         #['nfw','sech_disk','fornax',10000,'metal_point',0,15,2.0,0,0,None,{}],
                         #['nfw','sech_disk','fornax',10000,'metal_point',0,16,2.0,0,0,None,{}],
                         #['nfw','sech_disk','fornax',10000,'metal_point',0,17,2.0,0,0,None,{}],
                         #['nfw','sech_disk','fornax',10000,'metal_point',0,18,2.0,0,0,None,{}],
                         #['nfw','sech_disk','fornax',10000,'metal_point',0,19,2.0,0,0,None,{}],
                         #['nfw','sech_disk','fornax',10000,'metal_point',0,20,2.0,0,0,None,{}],
                         #['nfw','sech_disk','fornax',10000,'metal_point',0,21,2.0,0,0,None,{}],
                         #['nfw','sech_disk','fornax',10000,'metal_point',0,22,2.0,0,0,None,{}],
                         #['nfw','sech_disk','fornax',10000,'metal_point',0,23,2.0,0,0,None,{}],
                         #['nfw','sech_disk','fornax',10000,'metal_point',0,24,2.0,0,0,None,{}],
                         #['nfw','sech_disk','fornax',10000,'metal_point',0,25,2.0,0,0,None,{}],
                         #['nfw','sech_disk','fornax',10000,'metal_point',0,26,2.0,0,0,None,{}],
                         #['nfw','sech_disk','fornax',10000,'metal_point',0,27,2.0,0,0,None,{}],
                         #['nfw','sech_disk','fornax',10000,'metal_point',0,28,2.0,0,0,None,{}],
                         #['nfw','sech_disk','fornax',10000,'metal_point',0,29,2.0,0,0,None,{}],
                         #['nfw','sech_disk','fornax',10000,'metal_point',0,30,2.0,0,0,None,{}],
                         #['nfw','sech_disk','fornax',10000,'metal_point',0,31,2.0,0,0,None,{}],
                         #['nfw','sech_disk','fornax',10000,'metal_point',0,32,2.0,0,0,None,{}],
                         #['nfw','sech_disk','fornax',10000,'metal_point',0,33,2.0,0,0,None,{}],
                         #['nfw','sech_disk','fornax',10000,'metal_point',0,34,2.0,0,0,None,{}],
                         #['nfw','sech_disk','fornax',10000,'metal_point',0,35,2.0,0,0,None,{}],
                         #['nfw','sech_disk','fornax',10000,'metal_point',0,36,2.0,0,0,None,{}],
                         #['nfw','sech_disk','fornax',10000,'metal_point',0,37,2.0,0,0,None,{}],
                         #['nfw','sech_disk','fornax',10000,'metal_point',0,38,2.0,0,0,None,{}],
                         #['nfw','sech_disk','fornax',10000,'metal_point',0,39,2.0,0,0,None,{}],
                         #['nfw','sech_disk','fornax',10000,'metal_point',0,40,2.0,0,0,None,{}],
                         #['nfw','sech_disk','fornax',10000,'metal_point',0,41,2.0,0,0,None,{}]
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
    start_index_index = 6
    step_index = 7
    h_edge_on_index = 8
    d_edge_on_index = 9
    fixed_disk_angle_index = 10
    fixed_params_index = 11
    
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
                runMCMCForDwarfGalaxyProfilesAllPops(
                    parameter_series[i][el_index], parameter_series[i][lam_index], parameter_series[i][galaxy_index],
                    halo_funct = halo_funct_array[parameter_series[i][el_index]], disk_funct = disk_funct_array[parameter_series[i][lam_index]],
                    withDisk=parameter_series[i][withDisk_index], nIterations = parameter_series[i][nIter_index], saverun = saverun, show_best_fit = show_best_fit,
                    pop_selection_method = parameter_series[i][pop_select_index], start_param_index = parameter_series[i][start_index_index], outer_step = parameter_series[i][step_index],
                    use_artificial_data = artificial, generation_params_storer_index = artificial_param_indeces[j], generation_disk_funct = disk_funct_array[artificial_storers[j].lam], generation_halo_funct = halo_funct_array[artificial_storers[j].el], fixed_params = parameter_series[i][fixed_params_index])
            else:
                runMCMCForDwarfGalaxyProfilesAllPops(
                    parameter_series[i][galaxy_index],
                    halo_funct = halo_funct_array[parameter_series[i][halo_index]], disk_funct = disk_funct_array[parameter_series[i][disk_index]],
                    withDisk=parameter_series[i][withDisk_index], nIterations = parameter_series[i][nIter_index], saverun = saverun, show_best_fit = show_best_fit,
                    pop_selection_method = parameter_series[i][pop_select_index], start_param_index = parameter_series[i][start_index_index], outer_step = parameter_series[i][step_index],
                    use_artificial_data = artificial, generation_params_storer_index = 'None', generation_disk_funct = 0, generation_halo_funct = 0,
                    h_type = parameter_series[i][halo_index], d_type = parameter_series[i][disk_index], halo_edge_on = parameter_series[i][h_edge_on_index], disk_edge_on = parameter_series[i][d_edge_on_index], fixed_disk_angle = parameter_series[i][fixed_disk_angle_index], fixed_params = parameter_series[i][fixed_params_index])
