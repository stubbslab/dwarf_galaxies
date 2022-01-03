#Do a sequence of the MCMC Series algorithms, with varying run conditions.
#You can change various parameters (present set listed in commented line right above parameter_serires assignment), and have program execute in sequence
#This script can be run from command line without having to be in python environment (bash$ python doMCMCSeriesAllPops.py)
#This method tries to be as efficient as possible, reading in the various potentialFunctionArrays only once, as that can be a time consuming part of the process.

from runMCMCForDwarfGalaxyProfilesAllPops import runMCMCForDwarfGalaxyProfilesAllPops
from ArtificialParamsStorer import ArtificialParamsStorer
from readInPotentialInterpolator import readInPotentialInterpolator
import AstronomicalParameterArchive as apa

if __name__ == '__main__':
    extra_save_str = 'Elliptical_Oblate_WD_'
    astro_arch = apa.AstronomicalParameterArchive()
    deg_to_rad = astro_arch.getDegToRad()
    #parameters at present are: [halo_type, disk_type, galaxy, withDisk, number of iterations, population selection method, start parameter index, mask_N_R, mask_N_z, max_sky_angle, target_n_sky_pixels, n_los_bins, include morph,  include_kinem, apply_obs_mask, n_save_iters, halo_edge_on, disk_edge_on, fixed_disk_angle, fixed_halo_angle, fixed_params ]
    parameter_series = [ ['nfw','sech_disk','fornax',1,15000,'metal_point',0,100,100,None,500,21,1,1,1,101,1,1, 178.0 * deg_to_rad, 2.390, 'oblate', {}, {'lam':[0.1, 0.1], 'el':[0.427, 0.427], }],
                         #['nfw','sech_disk','fornax',1,15000,'metal_point',1,100,100,None,500,21,1,1,1,101,1,1, 178.0 * deg_to_rad, 2.390, 'oblate', {}, {'lam':[0.1, 0.1], 'el':[0.427, 0.427], }],
                         #['nfw','sech_disk','fornax',1,15000,'metal_point',2,100,100,None,500,21,1,1,1,101,1,1, 178.0 * deg_to_rad, 2.390, 'oblate', {}, {'lam':[0.1, 0.1], 'el':[0.427, 0.427], }],
                         #['nfw','sech_disk','fornax',1,15000,'metal_point',3,100,100,None,500,21,1,1,1,101,1,1, 178.0 * deg_to_rad, 2.390, 'oblate', {}, {'lam':[0.1, 0.1], 'el':[0.427, 0.427], }],
                         #['nfw','sech_disk','fornax',1,15000,'metal_point',4,100,100,None,500,21,1,1,1,101,1,1, 178.0 * deg_to_rad, 2.390, 'oblate', {}, {'lam':[0.1, 0.1], 'el':[0.427, 0.427], }],
                         #['nfw','sech_disk','fornax',1,15000,'metal_point',5,100,100,None,500,21,1,1,1,101,1,1, 178.0 * deg_to_rad, 2.390, 'oblate', {}, {'lam':[0.1, 0.1], 'el':[0.427, 0.427], }],
                         #['nfw','sech_disk','fornax',1,15000,'metal_point',6,100,100,None,500,21,1,1,1,101,1,1, 178.0 * deg_to_rad, 2.390, 'oblate', {}, {'lam':[0.1, 0.1], 'el':[0.427, 0.427], }],
                         #['nfw','sech_disk','fornax',1,15000,'metal_point',7,100,100,None,500,21,1,1,1,101,1,1, 178.0 * deg_to_rad, 2.390, 'oblate', {}, {'lam':[0.1, 0.1], 'el':[0.427, 0.427], }],
                         #['nfw','sech_disk','fornax',1,15000,'metal_point',8,100,100,None,500,21,1,1,1,101,1,1, 178.0 * deg_to_rad, 2.390, 'oblate', {}, {'lam':[0.1, 0.1], 'el':[0.427, 0.427], }],
                         #['nfw','sech_disk','fornax',1,15000,'metal_point',9,100,100,None,500,21,1,1,1,101,1,1, 178.0 * deg_to_rad, 2.390, 'oblate', {}, {'lam':[0.1, 0.1], 'el':[0.427, 0.427], }],
                         #['nfw','sech_disk','fornax',1,15000,'metal_point',10,100,100,None,500,21,1,1,1,101,1,1, 178.0 * deg_to_rad, 2.390, 'oblate', {}, {'lam':[0.1, 0.1], 'el':[0.427, 0.427], }],
                         #['nfw','sech_disk','fornax',1,15000,'metal_point',11,100,100,None,500,21,1,1,1,101,1,1, 178.0 * deg_to_rad, 2.390, 'oblate', {}, {'lam':[0.1, 0.1], 'el':[0.427, 0.427], }],
                         #['nfw','sech_disk','fornax',1,15000,'metal_point',12,100,100,None,500,21,1,1,1,101,1,1, 178.0 * deg_to_rad, 2.390, 'oblate', {}, {'lam':[0.1, 0.1], 'el':[0.427, 0.427], }],
                         #['nfw','sech_disk','fornax',1,15000,'metal_point',13,100,100,None,500,21,1,1,1,101,1,1, 178.0 * deg_to_rad, 2.390, 'oblate', {}, {'lam':[0.1, 0.1], 'el':[0.427, 0.427], }],
                         #['nfw','sech_disk','fornax',1,15000,'metal_point',14,100,100,None,500,21,1,1,1,101,1,1, 178.0 * deg_to_rad, 2.390, 'oblate', {}, {'lam':[0.1, 0.1], 'el':[0.427, 0.427], }],
                         #['nfw','sech_disk','fornax',1,15000,'metal_point',15,100,100,None,500,21,1,1,1,101,1,1, 178.0 * deg_to_rad, 2.390, 'oblate', {}, {'lam':[0.1, 0.1], 'el':[0.427, 0.427], }],
                         #['nfw','sech_disk','fornax',1,15000,'metal_point',16,100,100,None,500,21,1,1,1,101,1,1, 178.0 * deg_to_rad, 2.390, 'oblate', {}, {'lam':[0.1, 0.1], 'el':[0.427, 0.427], }],
                         #['nfw','sech_disk','fornax',1,15000,'metal_point',17,100,100,None,500,21,1,1,1,101,1,1, 178.0 * deg_to_rad, 2.390, 'oblate', {}, {'lam':[0.1, 0.1], 'el':[0.427, 0.427], }],
                         #['nfw','sech_disk','fornax',1,15000,'metal_point',18,100,100,None,500,21,1,1,1,101,1,1, 178.0 * deg_to_rad, 2.390, 'oblate', {}, {'lam':[0.1, 0.1], 'el':[0.427, 0.427], }],
                         #['nfw','sech_disk','fornax',1,15000,'metal_point',19,100,100,None,500,21,1,1,1,101,1,1, 178.0 * deg_to_rad, 2.390, 'oblate', {}, {'lam':[0.1, 0.1], 'el':[0.427, 0.427], }],
                         #['nfw','sech_disk','fornax',1,15000,'metal_point',20,100,100,None,500,21,1,1,1,101,1,1, 178.0 * deg_to_rad, 2.390, 'oblate', {}, {'lam':[0.1, 0.1], 'el':[0.427, 0.427], }],
                         #['nfw','sech_disk','fornax',1,15000,'metal_point',21,100,100,None,500,21,1,1,1,101,1,1, 178.0 * deg_to_rad, 2.390, 'oblate', {}, {'lam':[0.1, 0.1], 'el':[0.427, 0.427], }],
                         #['nfw','sech_disk','fornax',1,15000,'metal_point',22,100,100,None,500,21,1,1,1,101,1,1, 178.0 * deg_to_rad, 2.390, 'oblate', {}, {'lam':[0.1, 0.1], 'el':[0.427, 0.427], }],
                         #['nfw','sech_disk','fornax',1,15000,'metal_point',23,100,100,None,500,21,1,1,1,101,1,1, 178.0 * deg_to_rad, 2.390, 'oblate', {}, {'lam':[0.1, 0.1], 'el':[0.427, 0.427], }],
                         #['nfw','sech_disk','fornax',1,15000,'metal_point',24,100,100,None,500,21,1,1,1,101,1,1, 178.0 * deg_to_rad, 2.390, 'oblate', {}, {'lam':[0.1, 0.1], 'el':[0.427, 0.427], }],
                         #['nfw','sech_disk','fornax',1,15000,'metal_point',25,100,100,None,500,21,1,1,1,101,1,1, 178.0 * deg_to_rad, 2.390, 'oblate', {}, {'lam':[0.1, 0.1], 'el':[0.427, 0.427], }],
                         #['nfw','sech_disk','fornax',1,15000,'metal_point',26,100,100,None,500,21,1,1,1,101,1,1, 178.0 * deg_to_rad, 2.390, 'oblate', {}, {'lam':[0.1, 0.1], 'el':[0.427, 0.427], }],
                         #['nfw','sech_disk','fornax',1,15000,'metal_point',27,100,100,None,500,21,1,1,1,101,1,1, 178.0 * deg_to_rad, 2.390, 'oblate', {}, {'lam':[0.1, 0.1], 'el':[0.427, 0.427], }],
                         #['nfw','sech_disk','fornax',1,15000,'metal_point',28,100,100,None,500,21,1,1,1,101,1,1, 178.0 * deg_to_rad, 2.390, 'oblate', {}, {'lam':[0.1, 0.1], 'el':[0.427, 0.427], }],
                         #['nfw','sech_disk','fornax',1,15000,'metal_point',29,100,100,None,500,21,1,1,1,101,1,1, 178.0 * deg_to_rad, 2.390, 'oblate', {}, {'lam':[0.1, 0.1], 'el':[0.427, 0.427], }],
                         #['nfw','sech_disk','fornax',1,15000,'metal_point',30,100,100,None,500,21,1,1,1,101,1,1, 178.0 * deg_to_rad, 2.390, 'oblate', {}, {'lam':[0.1, 0.1], 'el':[0.427, 0.427], }],
                         #['nfw','sech_disk','fornax',1,15000,'metal_point',31,100,100,None,500,21,1,1,1,101,1,1, 178.0 * deg_to_rad, 2.390, 'oblate', {}, {'lam':[0.1, 0.1], 'el':[0.427, 0.427], }],
                         #['nfw','sech_disk','fornax',1,15000,'metal_point',32,100,100,None,500,21,1,1,1,101,1,1, 178.0 * deg_to_rad, 2.390, 'oblate', {}, {'lam':[0.1, 0.1], 'el':[0.427, 0.427], }],
                         #['nfw','sech_disk','fornax',1,15000,'metal_point',33,100,100,None,500,21,1,1,1,101,1,1, 178.0 * deg_to_rad, 2.390, 'oblate', {}, {'lam':[0.1, 0.1], 'el':[0.427, 0.427], }],
                         #['nfw','sech_disk','fornax',1,15000,'metal_point',34,100,100,None,500,21,1,1,1,101,1,1, 178.0 * deg_to_rad, 2.390, 'oblate', {}, {'lam':[0.1, 0.1], 'el':[0.427, 0.427], }],
                         #['nfw','sech_disk','fornax',1,15000,'metal_point',35,100,100,None,500,21,1,1,1,101,1,1, 178.0 * deg_to_rad, 2.390, 'oblate', {}, {'lam':[0.1, 0.1], 'el':[0.427, 0.427], }],
                         #['nfw','sech_disk','fornax',1,15000,'metal_point',36,100,100,None,500,21,1,1,1,101,1,1, 178.0 * deg_to_rad, 2.390, 'oblate', {}, {'lam':[0.1, 0.1], 'el':[0.427, 0.427], }],
                         #['nfw','sech_disk','fornax',1,15000,'metal_point',37,100,100,None,500,21,1,1,1,101,1,1, 178.0 * deg_to_rad, 2.390, 'oblate', {}, {'lam':[0.1, 0.1], 'el':[0.427, 0.427], }],
                         #['nfw','sech_disk','fornax',1,15000,'metal_point',38,100,100,None,500,21,1,1,1,101,1,1, 178.0 * deg_to_rad, 2.390, 'oblate', {}, {'lam':[0.1, 0.1], 'el':[0.427, 0.427], }],
                         #['nfw','sech_disk','fornax',1,15000,'metal_point',39,100,100,None,500,21,1,1,1,101,1,1, 178.0 * deg_to_rad, 2.390, 'oblate', {}, {'lam':[0.1, 0.1], 'el':[0.427, 0.427], }],
                         #['nfw','sech_disk','fornax',1,15000,'metal_point',40,100,100,None,500,21,1,1,1,101,1,1, 178.0 * deg_to_rad, 2.390, 'oblate', {}, {'lam':[0.1, 0.1], 'el':[0.427, 0.427], }],
                         #['nfw','sech_disk','fornax',1,15000,'metal_point',41,100,100,None,500,21,1,1,1,101,1,1, 178.0 * deg_to_rad, 2.390, 'oblate', {}, {'lam':[0.1, 0.1], 'el':[0.427, 0.427], }]
                       ]
    saverun = 1
    show_best_fit = 0
    artificial = 0
    artificial_param_indeces = [0]
    if not artificial:
        artificial_param_indeces = [0]
    art_params_storer = ArtificialParamsStorer()
    #artificial_storers = [art_params_storer.getStorer(index) for index in artificial_param_indeces]


    halo_index = 0
    disk_index = halo_index + 1
    galaxy_index = disk_index + 1
    withDisk_index = galaxy_index + 1
    nIter_index = withDisk_index + 1
    pop_select_index = nIter_index + 1
    start_index_index = pop_select_index + 1
    mask_R_index = start_index_index + 1
    mask_z_index = mask_R_index + 1
    max_sky_index = mask_z_index + 1
    target_n_sky_index = max_sky_index + 1
    n_los_bins_index = target_n_sky_index + 1
    include_morph_index = n_los_bins_index + 1
    include_kinem_index = include_morph_index + 1
    obs_mask_index = include_kinem_index + 1
    n_save_iters_index = obs_mask_index + 1
    halo_edge_on_index = n_save_iters_index + 1
    disk_edge_on_index = halo_edge_on_index + 1
    fixed_disk_angle_index = disk_edge_on_index + 1
    fixed_halo_angle_index = fixed_disk_angle_index + 1
    prol_or_obl_index = fixed_halo_angle_index + 1
    fixed_params_index = prol_or_obl_index + 1
    specified_params_index = fixed_params_index + 1

    #Now actually do each MCMC run.  This can take a while.
    for i in range(len(parameter_series)):
        parameter_set = parameter_series[i]
        halo_funct = readInPotentialInterpolator(parameter_set[halo_index])
        disk_funct = readInPotentialInterpolator(parameter_set[disk_index])
        for j in range(len(artificial_param_indeces)):
            print ('Starting MCMC series ' + str(i) +', ' + str(j))
            if artificial:
                runMCMCForDwarfGalaxyProfilesAllPops(
                    parameter_set[el_index], parameter_set[lam_index], parameter_set[galaxy_index],
                    halo_funct = halo_funct_array[parameter_set[el_index]], disk_funct = disk_funct_array[parameter_set[lam_index]],
                    withDisk=parameter_set[withDisk_index], nIterations = parameter_set[nIter_index], saverun = saverun, show_best_fit = show_best_fit,
                    pop_selection_method = parameter_set[pop_select_index], start_param_index = parameter_set[start_index_index], outer_step = parameter_set[step_index],
                    use_artificial_data = artificial, generation_params_storer_index = artificial_param_indeces[j], generation_disk_funct = disk_funct_array[artificial_storers[j].lam], generation_halo_funct = halo_funct_array[artificial_storers[j].el], fixed_params = parameter_series[i][fixed_params_index])
            else:
                runMCMCForDwarfGalaxyProfilesAllPops(
                    parameter_set[galaxy_index],
                    h_type = parameter_set[halo_index], d_type = parameter_set[disk_index],
                    withDisk=parameter_set[withDisk_index], nIterations = parameter_set[nIter_index], saverun = 1, show_best_fit = 0, start_param_index = parameter_set[start_index_index],
                    mask_R_bins = parameter_set[mask_R_index], mask_z_bins = parameter_set[mask_z_index], max_sky_angle_deg = parameter_set[max_sky_index], target_n_sky_pixels = parameter_set[target_n_sky_index],
                    pop_selection_method = parameter_set[pop_select_index], n_los_bins =  parameter_set[n_los_bins_index], #n_los_vel_bins = n_los_vel_bins, max_v_los = max_v_los,
                    halo_edge_on = parameter_set[halo_edge_on_index], disk_edge_on = parameter_set[disk_edge_on_index], halo_funct = halo_funct, disk_funct = disk_funct,
                    apply_observation_mask = parameter_set[obs_mask_index], include_morph_prob = parameter_set[include_morph_index], include_kinem_prob = parameter_set[include_kinem_index],
                    fixed_params = parameter_set[fixed_params_index], specified_param_ranges = parameter_set[specified_params_index],
                    fixed_disk_angle = parameter_series[i][fixed_disk_angle_index], fixed_halo_angle = parameter_series[i][fixed_halo_angle_index],
                    n_iters_to_save = parameter_set[n_save_iters_index], extra_save_str = extra_save_str , prolate_or_oblate = parameter_series[i][prol_or_obl_index])
