

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
    print ('Starting computation of likelihood series. ') 
    astro_archive = AstronomicalParameterArchive ()
    compute_params_archive = ComputationalArchive() 
    dSph_archive = DwarfGalDataArchive()
    r_half_light = 791.0
    M_star = 10.0 ** 7.39
    
    saverun = 1
    save_array = []
    results_number = 'recalibrating'
    gamma = astro_archive.getGamma() 

    #                     halo       disk          gal       pops               pop_select      step,  WD, apply_observation_mask,  M                rs                      halo_sym_axis                  halo_center                        Rd                     el     eps         disk_sym_axis                disk_center                        lam 
    
    halo_index = 0
    disk_index = halo_index + 1
    gal_index = disk_index + 1
    pops_index = gal_index + 1
    pop_select_index = pops_index + 1
    step_index = pop_select_index + 1 #step for determining density of profile sampling
    withDisk_index = step_index + 1
    apply_observation_mask_index = withDisk_index + 1
    M_index = apply_observation_mask_index + 1
    rs_index = M_index + 1 
    halo_sym_axis_index = rs_index + 1
    halo_center_index = halo_sym_axis_index + 1
    Rd_index = halo_center_index + 1
    el_index = Rd_index + 1
    eps_index = el_index + 1
    disk_sym_axis_index = eps_index + 1
    disk_center_index = disk_sym_axis_index + 1
    lam_index = disk_center_index + 1 
    
    
    
    #                     halo       disk          gal       pops               pop_select      step,  WD, apply_observation_mask,  M                rs                      halo_sym_axis                                                                            halo_center                                                                     Rd                     el     eps         disk_sym_axis                disk_center                        lam 
    #Gaussian-on-one-variable peaks
    parameter_series = [ ['nfw',     'sech_disk',  'fornax', ['MR','IM','MP'],  'metal_rigid',  2.0,   1,  1,                       8.587455914154990196e+09, 3.829752462557314175e+03, [7.009092642998508982e-01,0.000000000000000000e+00,7.132504491541815650e-01], [-1.388017101740483383e+02,0.000000000000000000e+00,-1.400154202388922613e+02], 4.511197670577145602e-01 * 8.587455914154990196e+09, 2.160000000000000142e+00,1.409738201437771114e-02, [7.660512571018597283e-01,0.000000000000000000e+00,-6.427794890105472669e-01], [4.426018416512122826e+02,0.000000000000000000e+00,3.924504606951299479e+01], 0.1],
                         ['nfw',     'sech_disk',  'fornax', ['MR','IM','MP'],  'metal_rigid',  2.0,   1,  1,                      8.587455914154990196e+09, 3.829752462557314175e+03, [7.009092642998508982e-01,0.000000000000000000e+00,7.132504491541815650e-01], [-1.388017101740483383e+02,0.000000000000000000e+00,-1.400154202388922613e+02], 4.511197670577145602e-01 * 8.587455914154990196e+09, 2.160000000000000142e+00,1.409738201437771114e-02, [7.660512571018597283e-01,0.000000000000000000e+00,-6.427794890105472669e-01], [4.426018416512122826e+02,0.000000000000000000e+00,3.924504606951299479e+01], 0.1],
                         ['cored',   'sech_disk',  'fornax', ['MR','IM','MP'],  'metal_rigid',  2.0,   1,  1,                       2.023644847128327847e+09, 3.011400279999999839e+03, [6.114354707349896056e-01, 0.0, -7.912942974185279699e-01], [-1.276999999999999886e+02, 0.0, 8.750000000000000000e+01], 1.0 * 2.023644847128327847e+09, 1.749999999999999889e-01, 0.0, [0.0, 0.0, 1.0], [0.0, 0.0, 0.0], 0.1],
                         ['cored',   'sech_disk',  'fornax', ['MR','IM','MP'],  'metal_rigid',  2.0,   1,  1,                       2.023644847128327847e+09, 3.011400279999999839e+03, [6.114354707349896056e-01, 0.0, -7.912942974185279699e-01], [-1.276999999999999886e+02, 0.0, 8.750000000000000000e+01],  1.0 * 2.023644847128327847e+09, 5.714285714285714413e+00, 0.0, [0.0, 0.0, 1.0], [0.0, 0.0, 0.0], 0.1],
                         ['burkert', 'sech_disk',  'fornax', ['MR','IM','MP'],  'metal_rigid',  2.0,   1,  1,                       1.190855390975193739e+09 , 5.566678319999999985e+03, [5.453209358290151965e-01, 0.0, 8.382273420418634435e-01], [-1.276999999999999886e+02, 0.0, -2.791999999999999886e+02], 1.0 * 1.190855390975193739e+09, 1.185817621249851822e-01, 0.0, [0.0, 0.0, 1.0], [0.0, 0.0, 0.0], 0.1],
                         ['burkert', 'sech_disk',  'fornax', ['MR','IM','MP'],  'metal_rigid',  2.0,   1,  1,                       2.023644847128327847e+09 , 3.011400279999999839e+03, [6.114354707349896056e-01, 0.0, -7.912942974185279699e-01], [-1.276999999999999886e+02, 0.0, 8.750000000000000000e+01],  1.0 * 2.023644847128327847e+09, 5.714285714285714413e+00, 0.0, [0.0, 0.0, 1.0], [0.0, 0.0, 0.0], 0.1],
                         ['nfw',     'sech_disk',  'fornax', ['MR','IM','MP'],  'metal_rigid',  2.0,   1,  1,                       8.587455914154990196e+09, 3.829752462557314175e+03, [7.009092642998508982e-01,0.000000000000000000e+00,7.132504491541815650e-01], [-1.388017101740483383e+02,0.000000000000000000e+00,-1.400154202388922613e+02], 4.511197670577145602e-01 * 8.587455914154990196e+09, 2.160000000000000142e+00,1.409738201437771114e-02, [7.660512571018597283e-01,0.000000000000000000e+00,-6.427794890105472669e-01], [4.426018416512122826e+02,0.000000000000000000e+00,3.924504606951299479e+01], 0.1],
                         ['nfw',     'sech_disk',  'fornax', ['MR','IM','MP'],  'metal_rigid',  2.0,   1,  1,                       9394654968.17372,       5192.67178574411,              [0.716910607650482, 0.0, -0.697165102854564], [-101.713746912193, 0.0, -138.145020724826], 6842.1842448351, 0.429, 0.39814545233736, [0.729615424741271, 0.0, -0.683857684010067], [148.713370317747, 0.0, -115.143124951233], 0.1],
                         ['burkert',     'sech_disk',  'fornax', ['MR','IM','MP'],  'metal_rigid',  2.0,   1,  1,                   3600565968.35128,       4940.5199353749,              [0.886270223331628, 0.0, 0.463168534375668], [-247.497640981089, 0.0, -121.336054888426], 2847.86282096028, 0.467, 0.496269371569995, [0.666939789196733, 0.0, -0.745111614180195], [7.60594665982927, 0.0, -76.9250271886467], 0.1],
                         ['burkert',     'sech_disk',  'fornax', ['MR','IM','MP'],  'metal_rigid',  2.0,   1,  1,                   2.45471e+7,       5.74458932e+02,              [0.886270223331628, 0.0, 0.463168534375668], [-38.924277, 0.0, -188.091593], 2847.86282096028, 0.467, 0.496269371569995, [0.666939789196733, 0.0, -0.745111614180195], [7.60594665982927, 0.0, -76.9250271886467], 0.1],
                         ['nfw',     'sech_disk',  'fornax', ['MR','IM','MP'],  'metal_rigid',  2.0,   0,  1,                       9.800e+09,       7.280e+03,              [0.7181262977631888, 0.0, -0.6959127965923143], [-70.0, 0.0, -143.0], 3.940e+03, 0.427, 0.0, [0.0, 0.0, 1.0], [0.0, 0.0, 0.0], 0.1],
                         ['nfw',     'sech_disk',  'fornax', ['MR','IM','MP'],  'metal_rigid',  2.0,   0,  1,                       9.800e+09,       3.940e+03,              [0.6971651028545646, 0.0, 0.7169106076504828], [-70.0, 0.0, -144.0], 3.940e+03, 2.16, 0.0, [0.0, 0.0, 1.0], [0.0, 0.0, 0.0], 0.1],
                         ['cored',     'sech_disk',  'fornax', ['MR','IM','MP'],  'metal_rigid',  2.0,   0,  1,                     3.560e+09,       1.900e+03,              [0.7229671459095686, 0.0, -0.690882411076858], [-72.0, 0.0, -144.0], 3.940e+03, 0.466, 0.0, [0.0, 0.0, 1.0], [0.0, 0.0, 0.0], 0.1],
                         ['cored',     'sech_disk',  'fornax', ['MR','IM','MP'],  'metal_rigid',  2.0,   0,  1,                     3.060e+09,       0.990e+03,              [0.6921431738704067, 0.0, 0.7217602280983623], [-72.0, 0.0, -144.0], 3.940e+03, 2.04, 0.0, [0.0, 0.0, 1.0], [0.0, 0.0, 0.0], 0.1],
                         ['burkert',     'sech_disk',  'fornax', ['MR','IM','MP'],  'metal_rigid',  2.0,   0,  1,                   0.470e+09,       1.430e+03,              [0.7229671459095686, 0.0, -0.690882411076858], [-72.0, 0.0, -144.0], 3.940e+03, 0.466, 0.0, [0.0, 0.0, 1.0], [0.0, 0.0, 0.0], 0.1],
                         ['burkert',     'sech_disk',  'fornax', ['MR','IM','MP'],  'metal_rigid',  2.0,   0,  1,                   0.410e+09,       0.750e+03,              [0.693401828275813, 0.0, 0.7205511116803305], [-71.0, 0.0, -144.0], 3.940e+03, 2.02, 0.0, [0.0, 0.0, 1.0], [0.0, 0.0, 0.0], 0.1],
                         ['nfw',     'sech_disk',  'fornax', ['MR','IM','MP'],  'metal_rigid',  2.0,   1,  1,                       9.700e+09,       7.100e+03,              [0.7181262977631888, 0.0, -0.6959127965923143], [-77.0, 0.0, -138.0], 6.000e+03, 0.427, 0.030, [0.0348995, 0, -0.999391], [-6.0, 0.0, -370.0], 0.1],
                         ['nfw',     'sech_disk',  'fornax', ['MR','IM','MP'],  'metal_rigid',  2.0,   1,  1,                       9.700e+09,       3.900e+03,              [0.6971651028545646, 0.0, 0.7169106076504828], [-76.0, 0.0, -137.0], 6.000e+03, 2.16, 0.039, [0.0348995, 0, -0.999391], [-80.0, 0.0, -350.0], 0.1],
                         ['cored',     'sech_disk',  'fornax', ['MR','IM','MP'],  'metal_rigid',  2.0,   1,  1,                     3.600e+09,       1.900e+03,              [0.7229671459095686, 0.0, -0.690882411076858], [-73.0, 0.0, -144.0], 6.800e+03, 0.466, 0.001, [0.0348995, 0, -0.999391], [-120.0, 0.0, 220.0], 0.1],
                         ['cored',     'sech_disk',  'fornax', ['MR','IM','MP'],  'metal_rigid',  2.0,   1,  1,                     2.870e+09,       1.000e+03,              [0.6921431738704067, 0.0, 0.7217602280983623], [-74.0, 0.0, -146.0], 1.290e+03, 2.04, 0.001, [0.0348995, 0, -0.999391], [-86.0, 0.0, -150.0], 0.1],
                         ['burkert',     'sech_disk',  'fornax', ['MR','IM','MP'],  'metal_rigid',  2.0,   1,  1,                   0.390e+09,       1.290e+03,              [0.7229671459095686, 0.0, -0.690882411076858], [-73.0, 0.0, -147.0], 8.100e+03, 0.466, 0.0041, [0.0348995, 0, -0.999391], [-180.0, 0.0, 0.0], 0.1],
                         ['burkert',     'sech_disk',  'fornax', ['MR','IM','MP'],  'metal_rigid',  2.0,   1,  1,                   0.520e+09,       0.840e+03,              [0.693401828275813, 0.0, 0.7205511116803305], [-72.0, 0.0, -145.0], 8.300e+03, 2.02, 0.0040, [0.0348995, 0, -0.999391], [-170.0, 0.0, -140.0], 0.1],
                         ['nfw',     'sech_disk',  'fornax', ['MR','IM','MP'],  'metal_rigid',  2.0,   1,  1,                       6.500e+09,       5.070e+03,              [0.7181, 0, -0.6959], [-78.0, 0.0, -145.0], 6.400e+03, 0.427, 0.00, [0.0348995, 0, -0.999391], [90.0, 0.0, 30.0], 0.1],
                         ['nfw',     'sech_disk',  'fornax', ['MR','IM','MP'],  'metal_rigid',  2.0,   1,  1,                       6.800e+09,       2.800e+03,              [0.6972, 0, 0.7169], [-73.0, 0.0, -146.0], 6.800e+03, 2.16, 0.039, [0.0348995, 0, -0.999391], [-140.0, 0.0, 15.0], 0.1],
                         ['cored',     'sech_disk',  'fornax', ['MR','IM','MP'],  'metal_rigid',  2.0,   1,  1,                     3.400e+09,       1.800e+03,              [0.7230, 0, -0.6909], [-72.0, 0.0, -145.0], 8.200e+03, 0.466, 0.00, [0.0348995, 0, -0.999391], [240.0, 0.0, 118.0], 0.1],
                         ['cored',     'sech_disk',  'fornax', ['MR','IM','MP'],  'metal_rigid',  2.0,   1,  1,                     3.200e+09,       0.920e+03,              [0.6921, 0, 0.7218], [-72.0, 0.0, -145.0], 1.044e+03, 2.04, 0.072, [0.0348995, 0, -0.999391], [-90.0, 0.0, -145.0], 0.1],
                         ['burkert',     'sech_disk',  'fornax', ['MR','IM','MP'],  'metal_rigid',  2.0,   1,  1,                   0.720e+09,       1.500e+03,              [0.7230, 0, -0.6909], [-72.0, 0.0, -144.0], 7.900e+03, 0.466, 0.25, [0.0348995, 0, -0.999391], [-280.0, 0.0, 64.0], 0.1],
                         ['burkert',     'sech_disk',  'fornax', ['MR','IM','MP'],  'metal_rigid',  2.0,   1,  1,                   0.580e+09,       0.770e+03,              [0.6934, 0, 0.7206], [-71.0, 0.0, -144.0], 7.900e+03, 2.02, 0.0, [0.0348995, 0, -0.999391], [-87.0, 0.0, -133.0], 0.1],
                         ['nfw',     'sech_disk',  'fornax', ['MR','IM','MP'],  'metal_rigid',  2.0,   1,  1,                       2.6420206199837527e+09,   3.1636456422890765e+03,       [0.7181262977631888,0.0,-0.6959127965923143], [-72.0, 0.0, -138.0], 3.1636456422890765e+03 * 1.8892760681235155, 0.427, 0.021, [0.0348995, 0, -0.999391], [-5.7, 0.0, -367.0], 0.1]
                         
                         ]


    pot_archive = PotentialArchive()
    unique_halo_types = []
    unique_disk_types = []
    #Create all of the potential files used in the MCMC algorithms.
    # That way, we don't need to reread in a potential file every time it is used.
    for parameter_set in parameter_series:
        if parameter_set[halo_index] not in unique_halo_types:
            print ('assigning halo type ' + str(parameter_set[halo_index]) + ' to unique_halo_types')
            unique_halo_types = unique_halo_types + [parameter_set[halo_index]]
        if parameter_set[disk_index] not in unique_disk_types:
            print ('assigning disk type ' + str(parameter_set[disk_index]) + ' to unique_disk_types')
            unique_disk_types = unique_disk_types + [parameter_set[disk_index]]
    halo_funct_array = {}
    disk_funct_array = {}
    for unique_halo_type in unique_halo_types:
        print ('Loading potential interpolator for halo type ' + unique_halo_type)
        halo_funct_array[unique_halo_type] = readInPotentialInterpolator(unique_halo_type)
    for unique_disk_type in unique_disk_types:
        print ('Loading potential interpolator for disk type ' + unique_disk_type)
        disk_funct_array[unique_disk_type] = readInPotentialInterpolator(unique_disk_type)

    #arcmin_limits_z =50.0
    likelihood_array = []
    save_array = []
    #Now actually compute the likelihood for each of the specified point. 
    for parameter_set in parameter_series:
        print ('Computing likelihood for parameter_set: ' + str(parameter_set))
        parameter_storers = []
        surface_profiles = []

        
        halo_type = parameter_set[halo_index]
        galaxy = parameter_set[gal_index]
        withDisk = parameter_set[withDisk_index]
        el = parameter_set[el_index]
        M = parameter_set[M_index]
        rs = parameter_set[rs_index]
        dist = dSph_archive.getDistanceFromSun([galaxy,'dummy_pop_variable'])
        #phi = parameter_set[phi_index]
        #theta = parameter_set[theta_index]
        halo_sym_axis = parameter_set[halo_sym_axis_index]
        halo_center = parameter_set[halo_center_index]
        disk_type = parameter_set[disk_index]
        #c = parameter_set[c_index]
        lam = parameter_set[lam_index]
        Rd = parameter_set[Rd_index]
        eps = parameter_set[eps_index]
        #a = parameter_set[a_index]
        #b = parameter_set[b_index]
        disk_sym_axis = parameter_set[disk_sym_axis_index]
        disk_center = parameter_set[disk_center_index]
        pops = parameter_set[pops_index]
        if len(pops) == 0:
            pops = dSph_archive.getPopulations(galaxy)
        pop_selection_method = parameter_set[pop_select_index]
        outer_step = parameter_set[step_index]
        apply_observation_mask = parameter_set[apply_observation_mask_index]
        #halo_file = pot_archive.getHaloFile(el)
        #disk_file = pot_archive.getDiskFile(lam)

        print ('Rd = ' + str(parameter_set[12])) 
        
        #print 'pops = '
        #print pops

        if not withDisk:
            zeta = 1.0
            eps = 0.0
        else:
            zeta = Rd / rs 

        #print 'max, min R in arcminutes are: ' + str(max(R_outer)) + ', ' + str(min(R_outer))
        #print 'max, min z in arcminutes are: ' + str(max(z_outer)) + ', ' + str(min(z_outer))
        storer = 0
        log_likelihood = 0.0

        arcmin_limits_R,arcmin_limits_z = dSph_archive.getObservationBounds([galaxy,'dummy_var'],return_unit = 'arcmin')
        #arcmin_limits_R = [-60.0,60.0]
        #arcmin_limits_z = [-60.0,60.0]
        arcmin_limits_los = [-50.0,50.0] 
        #R,z,los_bins = compRZLos(zeta, zeta * outer_step, outer_step, arcmin_limits_R, arcmin_limits_z, arcmin_limits_los, rs, dist)
        R,z,los_bins = compRZLos(1.0, 1.0 * outer_step, outer_step, arcmin_limits_R, arcmin_limits_z, arcmin_limits_los, rs, dist)
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
                                                  el = el, M = M, rs = rs, halo_sym_axis = halo_sym_axis, halo_center = halo_center,
                                                  lam = lam, zeta = zeta, eps = eps, disk_sym_axis = disk_sym_axis, disk_center = disk_center, stellar_data = star_data )
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
        #                     halo       disk          gal       pops               pop_select      step,  WD, apply_observation_mask,  dist,  omegaphi, sigsqr, c, M                rs                      halo_sym_axis                  halo_center                        Rd                     el     eps         disk_sym_axis                disk_center                        lam
        save_array = save_array + [[halo_type, disk_type, galaxy, '/'.join(pops), pop_selection_method, outer_step, withDisk, apply_observation_mask]
                                    + [storer.dist, storer.omegaphi, storer.sigsqr, storer.c]
                                    + [storer.M, storer.rs]
                                    + halo_sym_axis
                                    + [storer.phi, storer.theta]
                                    + halo_center
                                    + [Rd, el, eps]
                                    + disk_sym_axis
                                    + [storer.a, storer.b]
                                    + disk_center
                                    + [lam]
                                    + [log_likelihood] 
                                   ] 
                                   
                                    #+ parameter_set[0:halo_sym_axis_index] + [halo_sym_string]
                                    #+ parameter_set[halo_sym_axis_index+1:halo_center_index] + [halo_center_string]
                                    #+ parameter_set[halo_center_index+1:disk_sym_axis_index] + [disk_sym_string] 
                                    #+ parameter_set[disk_sym_axis_index+1:disk_center_index] + [disk_center_string]
                                    #+ parameter_set[disk_center_index+1:pops_index] + [pop_string] + parameter_set[pops_index+1:] + [log_likelihood]]
        likelihood_array = likelihood_array + [log_likelihood]


    #print likelihood_array
    #save_array = [parameter_series[i] + [likelihood_array[i]] for i in range(len(likelihood_array)) ]
    #print save_array 
    #print np.array(save_array)

    if saverun:
        save_dir = compute_params_archive.getSpotProbabilityDir()
        print ('save_dir = ' + save_dir )
        file_name = 'best_fit_probabilities_number_' + str(results_number) + '.csv'
        print ('Saving series to ' + save_dir + file_name) 
        header = 'halo, disk, gal, pops, pop_select, step, WD, obs_mask, dist, omega_phi, sigsqr, c, M, rs, hsym_x, hsym_y, hsym_z, phi, theta, hc_x, hc_y, hc_z, Rd, el, eps, hsym_x, hsym_y, hsym_z, a, b, hc_x, hc_y, hc_z, lam, logLikelihood' 
        np.savetxt(save_dir + file_name, np.array(save_array), delimiter = ',',fmt = '%10s',header = header) 
         
