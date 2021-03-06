#Do a sequence of the MCMC Series algorithms, with varying run conditions.
#You can change various parameters (present set listed in commented line right above parameter_serires assignment), and have program execute in sequence
#This script can be run from command line without having to be in python environment (bash$ python doMCMCSeriesAllPops.py)
#This method tries to be as efficient as possible, reading in the various potentialFunctionArrays only once, as that can be a time consuming part of the process.  

from runMCMCForDwarfGalaxyProfilesAllPops import runMCMCForDwarfGalaxyProfilesAllPops
from PotentialFunctionArray import PotentialFunctionArray
from PotentialArchive import PotentialArchive
from ArtificialParamsStorer import ArtificialParamsStorer
import time

if __name__ == '__main__':
    start = time.time()
    #parameters at present are: [el, lam, galaxy, number of iterations, population selection method, start parameter index] 
    parameter_series = [ #[0.2,0.1,'fornax',20,'metal_point',1,0],
                         #[0.2,0.1,'fornax',20,'metal_point',1,1],
                         #[0.2,0.1,'fornax',20,'metal_point',1,2],
                         #[0.2,0.1,'fornax',20,'metal_point',1,3],
                         #[0.2,0.1,'fornax',2000,'metal_point',1,4],
                         #[0.2,0.1,'fornax',2000,'metal_point',1,5],
                         #[0.2,0.1,'fornax',2000,'metal_point',1,6],
                         #[0.2,0.1,'fornax',2000,'metal_point',1,7],
                         #[0.2,0.1,'fornax',2000,'metal_point',1,8],
                         #[0.2,0.1,'fornax',2000,'metal_point',1,9],
                         #[0.2,0.1,'fornax',2000,'metal_point',1,10],
                         #[0.2,0.1,'fornax',2000,'metal_point',1,11]
                         #[0.2,0.1,'fornax',2000,'metal_point',1,12],
                         #[0.2,0.1,'fornax',2000,'metal_point',1,13],
                         #[0.2,0.1,'fornax',2000,'metal_point',1,14],
                         #[0.2,0.1,'fornax',2000,'metal_point',1,15],
                         #[0.2,0.1,'fornax',2000,'metal_point',1,16],
                         #[0.2,0.1,'fornax',2000,'metal_point',1,17],
                         #[0.2,0.1,'fornax',2000,'metal_point',1,18],
                         #[0.2,0.1,'fornax',2000,'metal_point',1,19],
                         #[0.2,0.1,'fornax',2000,'metal_point',1,20],
                         #[0.2,0.1,'fornax',2000,'metal_point',1,21],
                         #[0.2,0.1,'fornax',2000,'metal_point',1,22],
                         #[0.2,0.1,'fornax',2000,'metal_point',1,23],
                         #[0.2,0.1,'fornax',2000,'metal_point',1,24],
                         #[0.2,0.1,'fornax',2000,'metal_point',1,25],
                         #[0.2,0.1,'fornax',2000,'metal_point',1,26],
                         #[0.2,0.1,'fornax',2000,'metal_point',1,27],
                         #[0.2,0.1,'fornax',2000,'metal_point',1,28],
                         #[0.2,0.1,'fornax',2000,'metal_point',1,29],
                         #[0.2,0.1,'fornax',2000,'metal_point',1,30],
                         #[0.2,0.1,'fornax',2000,'metal_point',1,31],
                         #[0.2,0.1,'fornax',2000,'metal_point',1,32],
                         #[0.2,0.1,'fornax',2000,'metal_point',1,33],
                         #[0.2,0.1,'fornax',2000,'metal_point',1,34],
                         #[0.2,0.1,'fornax',2000,'metal_point',1,35],
                         #[0.2,0.1,'fornax',2000,'metal_point',1,36],
                         #[0.2,0.1,'fornax',2000,'metal_point',1,37],
                         #[0.2,0.1,'fornax',2000,'metal_point',1,38],
                         #[0.2,0.1,'fornax',2000,'metal_point',1,39],
                         #[0.2,0.1,'fornax',2000,'metal_point',1,40],
                         #[0.2,0.1,'fornax',2000,'metal_point',1,41],
                         #[0.2,0.1,'fornax',2000,'metal_point',1,42],
                         #[0.2,0.1,'fornax',2000,'metal_point',1,43],
                         #[0.2,0.1,'fornax',2000,'metal_point',1,44],
                         #[0.2,0.1,'fornax',2000,'metal_point',1,45],
                         #[0.2,0.1,'fornax',2000,'metal_point',1,46],
                         [0.2,0.1,'fornax',20,'metal_point',1,47]
                       ]
    saverun = 1
    show_best_fit = 0
    artificial = 1
    artificial_param_indeces = [4,5,6,7]
    if not artificial:
        artificial_param_indeces = ['None']
    art_params_storer = ArtificialParamsStorer()
    artificial_storers = [art_params_storer.getStorer(index) for index in artificial_param_indeces]

    el_index = 0
    lam_index = 1
    galaxy_index = 2
    nIter_index = 3
    pop_select_index = 4
    withDisk_index = 5
    start_index_index = 6
    pot_archive = PotentialArchive()
    unique_els = []
    unique_lams = []
    #Create all of the potential files used in the MCMC algorithms.
    # That way, we don't need to reread in a potential file every time it is used.
    for parameter_set in parameter_series:
        if parameter_set[el_index] not in unique_els:
            print 'assigning el ' + str(parameter_set[el_index]) + ' to unique_els'
            unique_els = unique_els + [parameter_set[el_index]]
        if parameter_set[lam_index] not in unique_lams:
            print 'assigning lam ' + str(parameter_set[lam_index]) + ' to unique_lams'
            unique_lams = unique_lams + [parameter_set[lam_index]]
    for storer in artificial_storers:
        art_el = storer.el
        art_lam = storer.lam
        if art_el not in unique_els:
            print 'assigning el ' + str(art_el) + ' to unique_els'
            unique_els = unique_els + [art_el]
        if parameter_set[lam_index] not in unique_lams:
            print 'assigning lam ' + str(art_lam) + ' to unique_lams'            
            unique_lams = unique_lams + [art_lam]
    print 'unique_els = '
    print unique_els
    print 'unique_lams = '
    print unique_lams
    halo_funct_array = {}
    disk_funct_array = {}
    for unique_el in unique_els:
        halo_file = pot_archive.getHaloFile(unique_el)
        print 'Defining halo potential file ' + halo_file
        halo_funct_array[unique_el] = PotentialFunctionArray(halo_file)
    for unique_lam in unique_lams:
        disk_file = pot_archive.getDiskFile(unique_lam)
        print 'Defining disk potential file ' + disk_file
        disk_funct_array[unique_lam] = PotentialFunctionArray(disk_file)

    #Now actually do each MCMC run.  This can take a while.
    end = time.time()
    print 'Took ' + str(end -start) + 's to get to MCMC chains. '
    for i in range(len(parameter_series)):
        for j in range(len(artificial_param_indeces)): 
            print 'Starting MCMC series ' + str(i) +', ' + str(j) 
            runMCMCForDwarfGalaxyProfilesAllPops(
                parameter_series[i][el_index], parameter_series[i][lam_index], parameter_series[i][galaxy_index],
                halo_funct = halo_funct_array[parameter_series[i][el_index]], disk_funct = disk_funct_array[parameter_series[i][lam_index]],
                withDisk=parameter_series[i][withDisk_index], nIterations = parameter_series[i][nIter_index], saverun = saverun, show_best_fit = show_best_fit,
                pop_selection_method = parameter_series[i][pop_select_index], start_param_index = parameter_series[i][start_index_index],
                use_artificial_data = artificial, generation_params_storer_index = artificial_param_indeces[j], generation_disk_funct = disk_funct_array[artificial_storers[j].lam], generation_halo_funct = halo_funct_array[artificial_storers[j].el])
