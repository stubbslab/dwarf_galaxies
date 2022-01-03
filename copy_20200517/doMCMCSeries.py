from runMCMCForDwarfGalaxyProfiles import runMCMCForDwarfGalaxyProfiles
from PotentialFunctionArray import PotentialFunctionArray
from PotentialArchive import PotentialArchive
from AstronomicalParameterArchive import AstronomicalParameterArchive 

if __name__ == '__main__':
    pot_archive = PotentialArchive() 
    parameter_series = [ #[0.2,0.1,['fornax','MR'],100,'metal_point',0],
                         #[0.3,0.1,['fornax','MR'],100,'metal_point',1],
                         #[0.5,0.1,['fornax','MR'],100,'metal_point',2],
                         #[0.2,0.2,['fornax','MR'],100,'metal_point',3],
                         #[0.3,0.2,['fornax','MR'],100,'metal_point',4], 
                         [0.3,0.2,['fornax','MR'],0,'metal_point',5]  ]
    unique_els = []
    unique_lams = []
    for parameter_set in parameter_series:
        if parameter_set[0] not in unique_els:
            unique_els = unique_els + [parameter_set[0]]
        if parameter_set[1] not in unique_lams:
            unique_lams = unique_lams + [parameter_set[1]]
    #print 'unique_els = '
    #print unique_els
    #print 'unique_lams = '
    #print unique_lams
    disk_funct_array = {}
    halo_funct_array = {}
    #pot_dir = '/Users/sasha/Documents/Harvard/physics/randall/potential_tables/'
    for unique_el in unique_els:
        halo_file = pot_archive.getHaloFile(unique_el)
        print 'Defining halo potential file ' + halo_file
        halo_funct_array[unique_el] = PotentialFunctionArray(halo_file)
    for unique_lam in unique_lams:
        disk_file = pot_archive.getDiskFile(unique_lam)
        print 'Defining disk potential file ' + disk_file
        disk_funct_array[unique_lam] = PotentialFunctionArray(disk_file)
    

    #disk_file = 'disk_pot_lambda_1_5_master_fixed.csv'
    #print 'reading in disk potential file from: ' + pot_dir + disk_file 
    #disk_funct = PotentialFunctionArray(pot_dir + disk_file)
    for i in range(len(parameter_series)):
        runMCMCForDwarfGalaxyProfiles(parameter_series[i][0], parameter_series[i][1], parameter_series[i][2], halo_funct = halo_funct_array[parameter_series[i][0]], disk_funct = disk_funct_array[parameter_series[i][1]], withDisk=0, nIterations = parameter_series[i][3], saverun = 0, show_best_fit = 0, pop_selection_method = parameter_series[i][4], start_param_index = parameter_series[i][5])
