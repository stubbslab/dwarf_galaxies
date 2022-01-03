

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
from compRZLos import compRZLos
from GalaxyMask import GalaxyMask 
import math
import numpy as np 

if __name__ == '__main__':
    #parameters at present are: [el, lam, galaxy, number of iterations, population selection method, start parameter index]
    print 'Starting computation of likelihood series. '
    astro_archive = AstronomicalParameterArchive ()
    compute_params_archive = ComputationalArchive() 
    dSph_archive = DwarfGalDataArchive()
    
    saverun = 1
    save_array = []
    results_number = '3'
    gamma = astro_archive.getGamma() 
    
    halo_index = 0
    el_index = 1
    M_index = 2
    rs_index = 3
    phi_index = 4
    theta_index = 5
    c_index = 6
    halo_sym_axis_index = 7
    halo_center_index = 8
    disk_index = 9
    lam_index = 10
    zeta_index = 11
    eps_index = 12
    a_index = 13
    b_index = 14
    disk_sym_axis_index = 15
    disk_center_index = 16
    gal_index = 17
    pops_index = 18
    pop_select_index = 19
    step_index = 20 #step for determining density of profile sampling 
    withDisk_index = 21
    apply_observation_mask_index = 22 
    
    #                     halo       el     M                       rs         phi    theta c     halo_ax, halo_center,         disk          lam   zeta  eps     a             b              disk_ax disk_center,       gal       pops             pop_select     step, WD, apply_observation_mask 
    #Gaussian-on-one-variable peaks
    parameter_series = [ ['nfw',     0.319, 1.57 * 10**8,           644.4,     2.162, 0.0,  5.0, 'none',  [-0.203,0.0,-0.106], 'sech_disk',  0.10, 1.0,   0.0,   -math.pi*0.0, math.pi * 0.0, 'none', [0.0,0.0,0.0],     'fornax', ['MR','IM','MP'], 'metal_point', 1.0,  0, 0],
                         ['nfw',     2.680, 1.57 * 10**8,           308.4,     0.592, 0.0,  5.0, 'none',  [-0.421,0.0,-0.221], 'sech_disk',  0.10, 1.0,   0.0,   -math.pi*0.0, math.pi * 0.0, 'none', [0.0,0.0,0.0],     'fornax', ['MR','IM','MP'], 'metal_point', 1.0,  0, 0],
                         ['cored',   0.369, 1.57 * 10**8,           545.8,     2.161, 0.0,  5.0, 'none',  [-0.245,0.0,-0.147], 'sech_disk',  0.10, 1.0,   0.1,   -math.pi*0.0, math.pi * 0.0, 'none', [0.0,0.0,0.0],     'fornax', ['MR','IM','MP'], 'metal_point', 1.0,  0, 0],
                         ['cored',   2.544, 1.57 * 10**8,           270.2,     0.588, 0.0,  5.0, 'none',  [-0.492,0.0,-0.301], 'sech_disk',  0.10, 1.0,   0.0,   -math.pi*0.0, math.pi * 0.0, 'none', [0.0,0.0,0.0],     'fornax', ['MR','IM','MP'], 'metal_point', 1.0,  0, 0],
                         ['burkert', 0.392, 1.57 * 10**8,           1174.9,    2.176, 0.0,  5.0, 'none',  [-0.114,0.0,-0.077], 'sech_disk',  0.10, 1.0,   0.1,   -math.pi*0.0, math.pi * 0.0, 'none', [0.0,0.0,0.0],     'fornax', ['MR','IM','MP'], 'metal_point', 1.0,  0, 0],
                         ['burkert', 2.299, 1.57 * 10**8,           610.1,     0.604, 0.0,  5.0, 'none',  [-0.224,0.0,-0.143], 'sech_disk',  0.10, 1.0,   0.0,   -math.pi*0.0, math.pi * 0.0, 'none', [0.0,0.0,0.0],     'fornax', ['MR','IM','MP'], 'metal_point', 1.0,  0, 0],

                         ['nfw',     0.429, 9788.656 * 10.0 ** 6.0, 7282.797,  2.343, 0.0,  None, 'none', [-70.294, 0.0, -143.808], 'sech_disk',  0.10, 1.0,   0.0,   -math.pi*0.0, math.pi * 0.0, 'none', [0.0,0.0,0.0],     'fornax', ['MR','IM','MP'], 'metal_point', 1.0,  0, 0],
                         ['nfw',     2.160, 9773.098 * 10.0 ** 6.0, 3936.0,    0.592, 0.0,  None, 'none',  [-69.961, 0.0, -144.074], 'sech_disk',  0.10, 1.0,   0.0,   -math.pi*0.0, math.pi * 0.0, 'none', [0.0,0.0,0.0],     'fornax', ['MR','IM','MP'], 'metal_point', 1.0,  0, 0],
                         ['cored',   0.467, 3556.0 * 10.0 ** 6.0,   1904.382,  2.161, 0.0,  None, 'none',  [-71.922, 0.0, -144.965], 'sech_disk',  0.10, 1.0,   0.1,   -math.pi*0.0, math.pi * 0.0, 'none', [0.0,0.0,0.0],     'fornax', ['MR','IM','MP'], 'metal_point', 1.0,  0, 0],
                         ['cored',   2.039, 3055.992 * 10.0 ** 6.0, 994.021,   0.588, 0.0,  None, 'none',  [-73.764, 0.0, -142.656], 'sech_disk',  0.10, 1.0,   0.0,   -math.pi*0.0, math.pi * 0.0, 'none', [0.0,0.0,0.0],     'fornax', ['MR','IM','MP'], 'metal_point', 1.0,  0, 0],
                         ['burkert', 0.466, 473.006 * 10.0 ** 6.0 , 1434.028,  2.176, 0.0,  None, 'none',  [-72.094, 0.0, -144.010], 'sech_disk',  0.10, 1.0,   0.1,   -math.pi*0.0, math.pi * 0.0, 'none', [0.0,0.0,0.0],     'fornax', ['MR','IM','MP'], 'metal_point', 1.0,  0, 0],
                         ['burkert', 2.018, 409.973 * 10.0 ** 6.0 , 751.177,   0.604, 0.0,  None, 'none',  [-70.108, 0.0, -144.233], 'sech_disk',  0.10, 1.0,   0.0,   -math.pi*0.0, math.pi * 0.0, 'none', [0.0,0.0,0.0],     'fornax', ['MR','IM','MP'], 'metal_point', 1.0,  0, 0],

                         ['nfw',     0.429, 9788.656 * 10.0 ** 6.0, 7282.797,  2.343, 0.0,  None, 'none', [-70.294, 0.0, -143.808], 'sech_disk',  0.10, 1.0,   0.0,   -math.pi*0.0, math.pi * 0.0, 'none', [0.0,0.0,0.0],     'fornax', ['MR','IM','MP'], 'metal_rigid', 2.0,  0, 1],
                         ['nfw',     2.160, 9773.098 * 10.0 ** 6.0, 3936.0,    0.592, 0.0,  None, 'none',  [-69.961, 0.0, -144.074], 'sech_disk',  0.10, 1.0,   0.0,   -math.pi*0.0, math.pi * 0.0, 'none', [0.0,0.0,0.0],     'fornax', ['MR','IM','MP'], 'metal_rigid', 2.0,  0, 1],
                         ['cored',   0.467, 3556.0 * 10.0 ** 6.0,   1904.382,  2.161, 0.0,  None, 'none',  [-71.922, 0.0, -144.965], 'sech_disk',  0.10, 1.0,   0.1,   -math.pi*0.0, math.pi * 0.0, 'none', [0.0,0.0,0.0],     'fornax', ['MR','IM','MP'], 'metal_rigid', 2.0,  0, 1],
                         ['cored',   2.039, 3055.992 * 10.0 ** 6.0, 994.021,   0.588, 0.0,  None, 'none',  [-73.764, 0.0, -142.656], 'sech_disk',  0.10, 1.0,   0.0,   -math.pi*0.0, math.pi * 0.0, 'none', [0.0,0.0,0.0],     'fornax', ['MR','IM','MP'], 'metal_rigid', 2.0,  0, 1],
                         ['burkert', 0.466, 473.006 * 10.0 ** 6.0 , 1434.028,  2.176, 0.0,  None, 'none',  [-72.094, 0.0, -144.010], 'sech_disk',  0.10, 1.0,   0.1,   -math.pi*0.0, math.pi * 0.0, 'none', [0.0,0.0,0.0],     'fornax', ['MR','IM','MP'], 'metal_rigid', 2.0,  0, 1],
                         ['burkert', 2.018, 409.973 * 10.0 ** 6.0 , 751.177,   0.604, 0.0,  None, 'none',  [-70.108, 0.0, -144.233], 'sech_disk',  0.10, 1.0,   0.0,   -math.pi*0.0, math.pi * 0.0, 'none', [0.0,0.0,0.0],     'fornax', ['MR','IM','MP'], 'metal_rigid', 2.0,  0, 1],

                         #dist,                     M,                        omegaphi,                  sigsqr ,                  el,                       rs ,                      phi ,                     theta ,                   h_xhat,                   h_yhat,                   h_zhat,                    h_x_center,               h_y_center,              h_z_center,               c,                       lam,                     zeta,                    eps,                     a,                       b,                        d_xhat,                  d_yhat,                  d_zhat,                  d_x_center,              d_y_center,              d_z_center,              nVisits,                 logLikelihood
                         #1.470000000000000000e+05, 6.805955054625613213e+09, -3.024615280204695296e-17, 1.561654949446706269e+02, 4.705847116190558510e-01, 2.602688826492321368e+03, 2.315298980873465418e+00, 0.000000000000000000e+00, 7.354249979691316330e-01, 0.000000000000000000e+00, -6.776061336514765943e-01, -8.247976085986836381e+01,0.000000000000000000e+00,-1.533811367329362554e+02,1.982036072616902089e+01,1.000000000000000056e-01,1.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,-1.570796326794896558e+00,0.000000000000000000e+00,0.000000000000000000e+00,1.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,0.000000000000000000e+00,2.000000000000000000e+00,2.142793232362578419e+04
                         #['nfw',     0.429, 9788.656 * 10.0 ** 6.0, 7282.797,  2.343, 0.0,  None, 'none', [-70.294, 0.0, -143.808], 'sech_disk',  0.10, 1.0,   0.0,   -math.pi*0.0, math.pi * 0.0, 'none', [0.0,0.0,0.0],     'fornax', ['MR','IM','MP'], 'metal_rigid', 1.0,  0, 1],
                         #['nfw',     2.160, 9773.098 * 10.0 ** 6.0, 3936.0,    0.592, 0.0,  None, 'none',  [-69.961, 0.0, -144.074], 'sech_disk',  0.10, 1.0,   0.0,   -math.pi*0.0, math.pi * 0.0, 'none', [0.0,0.0,0.0],     'fornax', ['MR','IM','MP'], 'metal_rigid', 1.0,  0, 1],
                         ['cored',   4.705847116190558510e-01, 6.805955054625613213 * 10.0 ** 9.0,   2.602688826492321368e+03,  2.315298980873465418e+00, 0.0,  None, 'none',  [-8.247976085986836381e+01, 0.0, -1.533811367329362554e+02], 'sech_disk',  0.10, 1.0,   0.1,   -math.pi*0.0, math.pi * 0.0, 'none', [0.0,0.0,0.0],     'fornax', ['MR','IM','MP'], 'metal_rigid', 2.0,  0, 1],
                         #['cored',   2.039, 3055.992 * 10.0 ** 6.0, 994.021,   0.588, 0.0,  None, 'none',  [-73.764, 0.0, -142.656], 'sech_disk',  0.10, 1.0,   0.0,   -math.pi*0.0, math.pi * 0.0, 'none', [0.0,0.0,0.0],     'fornax', ['MR','IM','MP'], 'metal_rigid', 1.0,  0, 1],
                         #['burkert', 0.466, 473.006 * 10.0 ** 6.0 , 1434.028,  2.176, 0.0,  None, 'none',  [-72.094, 0.0, -144.010], 'sech_disk',  0.10, 1.0,   0.1,   -math.pi*0.0, math.pi * 0.0, 'none', [0.0,0.0,0.0],     'fornax', ['MR','IM','MP'], 'metal_rigid', 1.0,  0, 1],
                         #['burkert', 2.018, 409.973 * 10.0 ** 6.0 , 751.177,   0.604, 0.0,  None, 'none',  [-70.108, 0.0, -144.233], 'sech_disk',  0.10, 1.0,   0.0,   -math.pi*0.0, math.pi * 0.0, 'none', [0.0,0.0,0.0],     'fornax', ['MR','IM','MP'], 'metal_rigid', 1.0,  0, 1],
                         
                         
                         ]


    pot_archive = PotentialArchive()
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
    halo_funct_array = {}
    disk_funct_array = {}
    for unique_halo_type in unique_halo_types:
        print 'Loading potential interpolator for halo type ' + unique_halo_type
        halo_funct_array[unique_halo_type] = readInPotentialInterpolator(unique_halo_type)
    for unique_disk_type in unique_disk_types:
        print 'Loading potential interpolator for disk type ' + unique_disk_type
        disk_funct_array[unique_disk_type] = readInPotentialInterpolator(unique_disk_type)

    #arcmin_limits_z =50.0
    likelihood_array = []
    save_array = []
    #Now actually compute the likelihood for each of the specified point. 
    for parameter_set in parameter_series:
        print 'Computing likelihood for parameter_set '
        print parameter_set
        parameter_storers = []
        surface_profiles = []
        
        halo_type = parameter_set[halo_index]
        galaxy = parameter_set[gal_index]
        withDisk = parameter_set[withDisk_index]
        el = parameter_set[el_index]
        M = parameter_set[M_index]
        rs = parameter_set[rs_index]
        dist = dSph_archive.getDistanceFromSun([galaxy,'dummy_pop_variable'])
        phi = parameter_set[phi_index]
        theta = parameter_set[theta_index]
        halo_sym_axis = parameter_set[halo_sym_axis_index]
        halo_center = parameter_set[halo_center_index]
        disk_type = parameter_set[disk_index]
        c = parameter_set[c_index]
        lam = parameter_set[lam_index]
        zeta = parameter_set[zeta_index]
        eps = parameter_set[eps_index]
        a = parameter_set[a_index]
        b = parameter_set[b_index]
        disk_sym_axis = parameter_set[disk_sym_axis_index]
        disk_center = parameter_set[disk_center_index]
        pops = parameter_set[pops_index]
        if len(pops) == 0:
            pops = dSph_archive.getPopulations(galaxy)
        pop_selection_method = parameter_set[pop_select_index]
        outer_step = parameter_set[step_index]
        apply_observation_mask = parameter_set[apply_observation_mask_index]

        #print ('[halo_type, galaxy, withDisk, el, M, rs, dist, phi, theta, halo_sym_axis, halo_center, disk_type, c, lam, zeta, eps, a, b, disk_sym_axis, disk_center, pops, pop_selection_method, outer_step, apply_observation_mask]')
        #print ([halo_type, galaxy, withDisk, el, M, rs, dist, phi, theta, halo_sym_axis, halo_center, disk_type, c, lam, zeta, eps, a, b, disk_sym_axis, disk_center, pops, pop_selection_method, outer_step, apply_observation_mask])
        #halo_file = pot_archive.getHaloFile(el)
        #disk_file = pot_archive.getDiskFile(lam)
        
        #print 'pops = '
        #print pops

        if not withDisk:
            zeta = 1.0
            eps = 0.0 

        #print 'max, min R in arcminutes are: ' + str(max(R_outer)) + ', ' + str(min(R_outer))
        #print 'max, min z in arcminutes are: ' + str(max(z_outer)) + ', ' + str(min(z_outer))
        storer = 0
        log_likelihood = 0.0

        arcmin_limits_R,arcmin_limits_z = dSph_archive.getObservationBounds([galaxy,'dummy_var'],return_unit = 'arcmin')
        #arcmin_limits_R = [-60.0,60.0]
        #arcmin_limits_z = [-60.0,60.0]
        arcmin_limits_los = [-50.0,50.0] 
        #R,z,los_bins = compRZLos(zeta, zeta * outer_step, outer_step, arcmin_limits_R, arcmin_limits_z, arcmin_limits_los, rs, dist)
        R, z, los_bins = compRZLos(1.0, 1.0 * outer_step, outer_step, arcmin_limits_R, arcmin_limits_z, arcmin_limits_los, rs, dist)
        Rmesh_degrees, zmesh_degrees = np.meshgrid(R * (rs / dist) * (180.0 / math.pi), z * (rs / dist) * (180.0 / math.pi))
        if apply_observation_mask:
            observation_mask = GalaxyMask(Rmesh_degrees, zmesh_degrees, galaxy, mask_types = ['n_vel_meas']).final_mask
        else:
            observation_mask = np.zeros(np.shape(Rmesh_degrees)) + 1.0 
        for pop in pops:
            population = [galaxy,pop]
            star_data = ObservedGalaxyStarData(population, pop_selection_method = pop_selection_method)
            #print 'for population ' + pop + ', len(star_data.corrRa) ' + str(len(star_data.corrRa))
            storer =  DwarfGalaxyParametersStorer(population, withDisk, sigsqr = star_data.sigSqr, 
                                                  el = el, M = M, rs = rs, phi = phi, theta = theta, halo_sym_axis = halo_sym_axis, c = c, halo_center = halo_center,
                                                  lam = lam, zeta = zeta, eps = eps, a = a, b = b, disk_sym_axis = disk_sym_axis, disk_center = disk_center )
            #storer.printContents()
            #print 'storer.sigsqr = ' + str(storer.sigsqr) + ' for population ' + pop
            #'storer.printContents() = '
            #print storer.printContents()

            parameter_storers = parameter_storers + [storer]
            surface_profile = SurfaceBrightnessProfile(R, z, los_bins, gamma, storer,
                                                       disk_interpolating_function = disk_funct_array[disk_type],
                                                       halo_interpolating_function = halo_funct_array[halo_type],
                                                       observation_mask = observation_mask )
            #if pop == 'MR':
            #        #print 'pop == ' + pop 
            #        #print 'surface_profile.Rlikeaxis  = '
            #        #print surface_profile.Rlikeaxis
            #        #print 'surface_profile.zlikeaxis  = '
            #        #print surface_profile.zlikeaxis    
            #        brightness_mesh = surface_profile.surfaceBrightness
            #        #print 'surface_profile.surfaceBrightness = '
            #        #print surface_profile.surfaceBrightness
            surface_profiles = surface_profiles + [surface_profile]
            log_likelihood = log_likelihood + surface_profile.sumLogSurfaceBrightness(star_data.proj_x, star_data.proj_y)

            #plt.scatter(np.array(star_data.corrRa) * 60.0 ,np.array(star_data.corrDec) * 60.0)
        pop_string = ' '.join(pops) 
        halo_sym_string = ''
        if halo_sym_axis is 'none':
            halo_sym_string = 'none'
        else:
            halo_sym_string = ' '.join([str(elem) for elem in halo_sym_axis])
        halo_center_string = ''
        if halo_center is 'none':
            halo_center_string = 'none'
        else:
            halo_center_string = ' '.join([str(elem) for elem in halo_center])
        disk_sym_string = ''
        if disk_sym_axis is 'none':
            disk_sym_string = 'none'
        else:
            disk_sym_string = ' '.join([str(elem) for elem in disk_sym_axis])
        disk_center_string = ''
        if disk_center is 'none':
            disk_center_string = 'none'
        else:
            disk_center_string = ' '.join([str(elem) for elem in disk_center])
        save_array = save_array + [[storer.dist, storer.omegaphi, storer.sigsqr, storer.c]
                                    + parameter_set[0:halo_sym_axis_index] + [halo_sym_string]
                                    + parameter_set[halo_sym_axis_index+1:halo_center_index] + [halo_center_string]
                                    + parameter_set[halo_center_index+1:disk_sym_axis_index] + [disk_sym_string] 
                                    + parameter_set[disk_sym_axis_index+1:disk_center_index] + [disk_center_string]
                                    + parameter_set[disk_center_index+1:pops_index] + [pop_string] + parameter_set[pops_index+1:] + [log_likelihood]]
        likelihood_array = likelihood_array + [log_likelihood]


    #print likelihood_array
    #save_array = [parameter_series[i] + [likelihood_array[i]] for i in range(len(likelihood_array)) ]
    #print save_array 
    #print np.array(save_array)

    if saverun:
        save_dir = compute_params_archive.getSpotProbabilityDir()
        print 'save_dir = ' + save_dir 
        file_name = 'best_fit_probabilities_number_' + str(results_number) + '.csv'
        print 'Saving series to ' + save_dir + file_name
        header = 'dist, omega_phi, sigsqr, c, halo, el, M, rs , phi , theta , halo_sym_axis, halo_center, disk, lam, zeta, eps, a, b, disk_sym_axis, disk_center, galaxy, pops, pop_select, step, withDisk, obs_mask, logLikelihood' 
        np.savetxt(save_dir + file_name, np.array(save_array), delimiter = ',',fmt = '%10s',header = header) 
         
