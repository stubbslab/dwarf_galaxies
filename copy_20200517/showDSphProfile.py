import math
import numpy as np
from AstronomicalParameterArchive import AstronomicalParameterArchive
from DwarfGalDataArchive import DwarfGalDataArchive
from DwarfGalaxyParametersStorer import DwarfGalaxyParametersStorer
from PotentialArchive import PotentialArchive
from ObservedGalaxyStarData import ObservedGalaxyStarData
from PotentialFunctionArray import PotentialFunctionArray
from SurfaceBrightnessProfile import SurfaceBrightnessProfile
from VisitedDsphMCMCPoint import VisitedDsphMCMCPoint
from distributeRandomStars import distributeRandomStars
from logList import logList
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from compRZLos import compRZLos
from ComputationalArchive import ComputationalArchive
from GalaxyMask import GalaxyMask
from readInPotentialInterpolator import readInPotentialInterpolator


#def showDSphProfiles(populations, pop_selection_methods = ['none'],
#                    apply_crowding_masking = [0], apply_observation_masking = [1], masking_fields = [[]], show_contours = [1], show_stars = [1], fields_of_interest = [[]],
#                    withDisk = [0], disk_funct = [0], halo_funct = [0], return_prof = [0],
#                    dist = ['none'], M = ['none'], vphi = ['none'], sigsqr = ['none'], figname = ['unnamed_dSph_profile.png'], save = [0], show = [1], colorbar = [0],
#                    el=[0.3], rs = [350], phi = [math.pi*0.0], theta = [math.pi * 0.0], halo_sym_axis = ['none'], c = [5.0], halo_center = [[0.0, 0.0, 0.0]], halo_edge_on = [0],
#                    lam=[0.1], zeta = [0.2], eps = [0.05], a = [math.pi*0.0], b = [math.pi*0.0, disk_sym_axis = 'none', disk_center = [0.0, 0.0, 0.0],  disk_edge_on = 0,
#                    x_lims_arcmin = [-40.3,35.05], y_lims_arcmin = [-48.2,34.05], outer_step_arcmin = 1.0, inner_step_arcmin = 1.0, inner_lim_factor = 0.0, disp_units = 'arcmin',
#                    artificial_data = 0, artificial_parameters_storer = 'None',art_disk_funct = 0, art_halo_funct = 0, plot_title = None, figsize = (10,10))

def plotAllFornaxPopulations(param_dicts, disk_functs, halo_functs, withDisks, titles,
                             pop_selection_method = 'metal_point', x_lims_arcmin = [-50.0, 50.0], y_lims_arcmin =  [-50.0, 50.0], n_xticks = 8, n_yticks = 8):
    pops = ['MP','IM','MR']
    f, axarr = plt.subplots(3, len(param_dicts), sharex = True, sharey = True, figsize = [2.0 * len(param_dicts), 2.0 * len(pops)  ])
    plt.subplots_adjust(wspace = 0.0, hspace = 0.0)
    x_ticks = np.around(np.linspace(x_lims_arcmin[0], x_lims_arcmin[1], n_xticks+1)[0:n_xticks] + (x_lims_arcmin[1] - x_lims_arcmin[0]) / (2.0 * (n_xticks + 1)), 1)
    y_ticks = np.around(np.linspace(y_lims_arcmin[0], y_lims_arcmin[1], n_yticks+1)[0:n_yticks] + (y_lims_arcmin[1] - y_lims_arcmin[0]) / (2.0 * (n_yticks + 1)), 1)
    print ('x_ticks = ' + str(x_ticks))
    print ('y_ticks = ' + str(y_ticks))
    for i in range(len(pops)):
        pop = pops[i]
        for j in range(len(param_dicts)):
            params = param_dicts[j]
            disk_funct = disk_functs[j]
            halo_funct = halo_functs[j]
            withDisk = withDisks[j]
            showDSphProfiles(['fornax',pop], pop_selection_method = pop_selection_method, axarr = axarr, axarr_indeces = [i, j],
                            apply_crowding_masking = 0, apply_observation_masking = 1, masking_fields = [], show_contours = 1, show_stars = 1, fields_of_interest = [],
                            withDisk = withDisk, disk_funct = disk_funct, halo_funct = halo_funct, return_prof = 0,
                            dist = params['dist'], M = params['M'], vphi = params['vphi'], sigsqr = 'none', save = 0, show = 0, colorbar = 0,
                            el=params['el'], rs = params['rs'], phi = params['phi'], theta = params['theta'],
                            halo_sym_axis = 'none', c = params['c'], halo_center = [params['hx'], params['hy'], params['hz']], halo_edge_on = 0,
                            lam=params['lam'], zeta = params['zeta'], eps = params['eps'], a = params['a'], b = params['b'],
                            disk_sym_axis = 'none', disk_center = [params['dx'], params['dy'], params['dz']],  disk_edge_on = 0,
                            x_lims_arcmin = x_lims_arcmin, y_lims_arcmin = y_lims_arcmin, outer_step_arcmin = 2.0, inner_step_arcmin = 1.0, inner_lim_factor = 0.0, disp_units = 'arcmin',
                            artificial_data = 0, artificial_parameters_storer = 'None',art_disk_funct = 0, art_halo_funct = 0, plot_title = None, figsize = (10,10))
            if j == 0:
                axarr[i,j].set_ylabel(pop + ' population \n Dist. from dSph center \n parallel to Dec (amin)')
                axarr[i,j].set_yticks(y_ticks)
                axarr[i,j].set_yticklabels(np.around(y_ticks, 1), rotation = 45)
            #else:
            #    axarr[i,j].set_yticks(y_ticks)
            #    axarr[i,j].set_yticklabels([])
            if i == len(pops) -1:
                axarr[i,j].set_xlabel('Dist. from dSph center \n parallel to RA (amin)')
                axarr[i,j].set_xticks(x_ticks)
                axarr[i,j].set_xticklabels(np.around(x_ticks, 1) , rotation = 45)
            else:
                axarr[i,j].set_xticks(x_ticks)
                axarr[i,j].set_xticklabels([])
            if i == 0:
                axarr[i,j].set_title(titles[j])

    plt.show()

def showDSphProfiles(galaxy, populations,
                     pop_selection_method = 'metal_rigid', disk_type = 'sech_disk', halo_type = 'nfw', plot_titles = ['MR', 'IM', 'MP'], master_title = '',
                     axarr = None, axarr_indeces = None,
                     apply_crowding_masking = 0, apply_observation_mask = 1, masking_fields = [], show_contours = 1, show_stars = 1, fields_of_interest = [],
                     withDisk = 0, disk_funct = None, halo_funct = None, return_prof = 0,
                     dist = 'none', M = 'none', vphi = 'none', sigsqr_RR = 'none',
                     dispersion_bs = [1.0, 1.0, 1.0], gamma0s = [0.0, 0.0, 0.0],  beta_exps = [1.0, 1.0, 1.0], beta_R1s = [1.0, 1.0, 1.0], beta_R2s = [1.0, 1.0, 1.0],
                     save = 0, show = 1, colorbar = 0, outer_step = 0.5,
                     el=0.3, rs = 350.0, phi = math.pi*0.0, theta = math.pi * 0.0, halo_sym_axis = 'none', c = 'none', halo_center = [0.0, 0.0, 0.0], halo_edge_on = 0,
                     lam=0.1, zeta = 0.2, eps = 0.05, a = math.pi*0.0, b = math.pi*0.0, disk_sym_axis = 'none', disk_center = [0.0, 0.0, 0.0],  disk_edge_on = 0,
                     x_lims_arcmin = [-50.0, 50.0], y_lims_arcmin =  [-50.0, 50.0], outer_step_arcmin = 1.0, inner_step_arcmin = 1.0, inner_lim_factor = 0.0, disp_units = 'arcmin',
                     artificial_data = 0, artificial_parameters_storer = 'None',art_disk_funct = 0, art_halo_funct = 0, figsize = (15,5)):


    dSph_archive = DwarfGalDataArchive ()
    potential_archive = PotentialArchive()
    computational_params_archive = ComputationalArchive()
    matplotlib.rcParams['xtick.direction'] = 'out'
    matplotlib.rcParams['ytick.direction'] = 'out'
    #disk_file = potential_archive.getDiskFile(lam)
    #halo_file = potential_archive.getHaloFile(el)

    astro_archive = AstronomicalParameterArchive ()

    gamma = astro_archive.getGamma()
    deg_to_rad = astro_archive.getDegToRad()
    #(param_indeces_to_vary, param_ranges, gauss_sampling_width) =  start_dsph_parameter_storer.getMCMCInfo()
    #print 'start_dsph_parameter_storer.dist = ' + str(start_dsph_parameter_storer.dist)

    arcmin_limits_R,arcmin_limits_z = dSph_archive.getObservationBounds([galaxy,'dummy_var'],return_unit = 'arcmin')
    #arcmin_limits_R = [-60.0, 60.0]
    #arcmin_limits_z = [-60.0, 60.0]
    arcmin_limits_los = [-50.0,50.0]

    if not disk_funct:
        print ('Loading potential interpolator for disk type ' + disk_type)
        disk_funct= readInPotentialInterpolator(disk_type)

    if not halo_funct:
        print ('Loading potential interpolator for halo type ' + halo_type)
        halo_funct= readInPotentialInterpolator(halo_type)

    if axarr is None:
        fig, axarr = plt.subplots(1, len(populations), figsize=figsize, squeeze = False)
        axarr_indeces = [[0,0], [0,1], [0,2]]

    for pop_index in range(len(populations)):
        population = populations[pop_index]
        dispersion_b = dispersion_bs[pop_index]
        gamma0 = gamma0s[pop_index]
        beta_exp = beta_exps[pop_index]
        beta_R1 = beta_R1s[pop_index]
        beta_R2 = beta_R2s[pop_index]
        print ('population = ' + str(population))
        local_axarr_indeces = axarr_indeces[pop_index]
        star_data = ObservedGalaxyStarData([galaxy,population], pop_selection_method = pop_selection_method)
        fields = star_data.Field


        dsph_parameter_storer = DwarfGalaxyParametersStorer([galaxy, population],withDisk=withDisk, stellar_data = star_data,
                                                        dist = dist, M = M,
                                                        vphi = vphi, sigsqr_RR = sigsqr_RR, dispersion_b = dispersion_b, gamma0 = gamma0, beta_exp = beta_exp, beta_R1 = beta_R1, beta_R2 = beta_R2,
                                                        el=el, rs = rs, phi = phi, theta = theta, halo_sym_axis = halo_sym_axis, c = c, halo_center = halo_center, halo_edge_on = halo_edge_on,
                                                        lam=lam, zeta = zeta, eps = eps, a = a, b = b, disk_sym_axis = disk_sym_axis, disk_center = disk_center, disk_edge_on = disk_edge_on  )
        print ('dsph_parameter_storer.printContents() = ')
        dsph_parameter_storer.printContents()

        dist = dsph_parameter_storer.dist
        #print 'dist = ' + str(dist)
        rs = dsph_parameter_storer.rs
        R,z,los_bins = compRZLos(1.0, 1.0 * outer_step, outer_step, arcmin_limits_R, arcmin_limits_z, arcmin_limits_los, rs, dist)
        print ('arcmin_limits_R = ' + str(arcmin_limits_R))
        Rmesh_degrees, zmesh_degrees = np.meshgrid(R * (rs / dist) * (180.0 / math.pi), z * (rs / dist) * (180.0 / math.pi))
        Rmesh_arcmin, zmesh_arcmin = [Rmesh_degrees * 60.0, zmesh_degrees * 60.0]

        if apply_observation_mask:
            observation_mask = GalaxyMask(Rmesh_degrees, zmesh_degrees, galaxy, mask_types = ['n_vel_meas']).final_mask
        else:
            observation_mask = np.zeros(np.shape(Rmesh_degrees)) + 1.0

        dsph_brightness_profile  = SurfaceBrightnessProfile(R, z, los_bins, gamma, dsph_parameter_storer,
                                                           disk_interpolating_function = disk_funct,
                                                           halo_interpolating_function = halo_funct,
                                                           observation_mask = observation_mask )
        print ('dsph_brightness_profile.onSkyInterpolator([[0.0, 0.1, 0.0, 0.1], [0.0, 0.0, 0.1, 0.1]] ) = ' +str(dsph_brightness_profile.onSkyInterpolator( np.dstack([[0.0, 0.1, 0.0, 0.1], [0.0, 0.0, 0.1, 0.1]]) [0]  )))
        Rmesh_scale_radii, zmesh_scale_radii = np.meshgrid(R, z)
        surface_prob_mesh = dsph_brightness_profile.onSkyInterpolator((Rmesh_scale_radii, zmesh_scale_radii) )
        print ('surface_prob_mesh[0:5, 0:5] = ' + str(surface_prob_mesh[0:5, 0:5] ))
        colormap='b' if population == 'MP' else 'g' if population == 'IM' else 'r'
        if pop_selection_method == 'none':
            colormap = 'black'

        max_prob = np.max(np.abs(surface_prob_mesh))
        min_prob = np.min(np.abs(surface_prob_mesh))
        if min_prob <= 0.0:
            #print 'np.unique(np.abs(surface_prob_mesh).flatten()) = ' + str(np.unique(np.abs(surface_prob_mesh).flatten()))
            min_prob = np.partition(np.unique(np.abs(surface_prob_mesh).flatten()), 1)[1]
            #min_prob = min(10**-8,max_prob * 10**-8)
        #print 'min_prob = '
        #print min_prob
        #print 'max_prob = '
        #print max_prob
        nlevels = 20
        log_levels1=logList(min_prob,100.0,nlevels / 4)
        log_levels2=logList(max_prob * 0.01,max_prob,3*nlevels/4)
        #print 'log_levels1 = '
        #print log_levels1
        #print 'log_levels2 = '
        #print log_levels2
        #print type(log_levels1)

        #log_levels = np.union1d(log_levels1, log_levels2)
        log_levels = logList(max_prob * 10.0 ** (-3.0), max_prob, nlevels)
        test_levels = np.arange(min_prob,max_prob + (max_prob-min_prob)/(2.0 * 20.0) ,(max_prob-min_prob)/(20.0) )
        #print 'surface_prob_mesh = '
        # print surface_prob_mesh
        print ('local_axarr_indeces = ' + str(local_axarr_indeces ))
        if show_contours:
            if apply_observation_mask:
                masking_population = population
            else:
                masking_population = None
            CS=axarr[local_axarr_indeces[0], local_axarr_indeces[1]].contour(Rmesh_arcmin, zmesh_arcmin, surface_prob_mesh, norm=LogNorm(),levels=log_levels)
            if colorbar:
                CB_contour=axarr[local_axarr_indeces[0], local_axarr_indeces[1]].colorbar(shrink=0.8, extend='both', format='%.2e')

        if artificial_data:
            print ('Nstars = ' + str(len(star_data.corrRa)))
            if artificial_parameters_storer == 'None':
                artificial_profile = dsph_brightness_profile
            else:
                 art_disk_file = potential_archive.getDiskFile(artificial_parameters_storer.lam)
                 art_halo_file = potential_archive.getHaloFile(artificial_parameters_storer.el)
                 artificial_parameters_storer.sigsqr = sigsqr
                 print ('artificial_parameters_storer.sigsqr = ' + str(artificial_parameters_storer.sigsqr) )
                 artificial_profile = SurfaceBrightnessProfile(R,z,los_bins,gamma,artificial_parameters_storer,art_disk_file,art_halo_file, disk_interpolating_function = art_disk_funct, halo_interpolating_function = art_halo_funct)
            corr_ra_one_pop_in_scale_radii,corr_dec_one_pop_in_scale_radii = distributeRandomStars(len(star_data.corrRa),artificial_profile,400,400,withDisk)
            corr_ra_one_pop = [RA * rs  / dist * 1 / deg_to_rad for RA in corr_ra_one_pop_in_scale_radii]
            corr_dec_one_pop = [Dec * rs  / dist * 1 / deg_to_rad for Dec in corr_dec_one_pop_in_scale_radii]

        else:
            corr_ra_one_pop = star_data.corrRa
            corr_dec_one_pop = star_data.corrDec
        if len(fields_of_interest) > 0:
            corr_ra_one_pop_untrimmed = corr_ra_one_pop[:]
            corr_dec_one_pop_untrimmed = corr_dec_one_pop[:]
            corr_ra_one_pop = []
            corr_dec_one_pop = []
            for i in range(len(fields)):
                field = fields[i]
                if field in fields_of_interest:
                    corr_ra_one_pop = corr_ra_one_pop + [corr_ra_one_pop_untrimmed[i]]
                    corr_dec_one_pop = corr_dec_one_pop + [corr_dec_one_pop_untrimmed[i]]


        if show_stars:

            if disp_units == 'arcmin':
                 print ('local_axarr_indeces = ' + str(local_axarr_indeces) )
                 axarr[local_axarr_indeces[0], local_axarr_indeces[1]].scatter(np.array(corr_ra_one_pop) * 60.0,np.array(corr_dec_one_pop) * 60.0 ,s=4.,color=colormap)
            elif disp_units == 'degrees':
                 axarr[axarr_indeces[0], axarr_indeces[1]].scatter(corr_ra_one_pop, corr_dec_one_pop, s=4., color=colormap)
        #axarr[axarr_indeces[0], axarr_indeces[1]].scatter(corr_ra_one_pop, corr_dec_one_pop)

        if disp_units == 'arcmin':
            y_lims = y_lims_arcmin
            x_lims = x_lims_arcmin
        elif disp_units == 'degrees':
            x_lims = [x_lim / 60.0 for x_lim in x_lims_arcmin]
            y_lims = [y_lim / 60.0 for y_lim in y_lims_arcmin]
        print ('[x_lims, y_lims] =' + str([x_lims, y_lims]))

        axarr[local_axarr_indeces[0], local_axarr_indeces[1]].set_ylims = y_lims #ylimits from Walker plot
        axarr[local_axarr_indeces[0], local_axarr_indeces[1]].set_xlims = x_lims #xlimits from Walker plot
        axarr[local_axarr_indeces[0], local_axarr_indeces[1]].set_xlabel('Arcmin from Center Parallel to RA')
        axarr[local_axarr_indeces[0], local_axarr_indeces[1]].set_ylabel('Arcmin from Center Parallel to Dec')
        plot_title = plot_titles[pop_index]
        axarr[local_axarr_indeces[0], local_axarr_indeces[1]].set_title(plot_title)

    matplotlib.rcParams.update({'font.size': 12})
    #axarr[axarr_indeces[0], axarr_indeces[1]].set_xlabel('ra sep')
    #axarr[axarr_indeces[0], axarr_indeces[1]].set_ylabel('dec sep')
    if plot_title is None:
        plot_title = 'Fornax Stars Position on Sky'
    #axarr[axarr_indeces[0], axarr_indeces[1]].set_title(plot_title)
    pltFileDir = computational_params_archive.getPlotDir()
    #Uncomment the following line if you want to save the file

    #axarr[axarr_indeces[0], axarr_indeces[1]].axis('scaled')
    #axarr[local_axarr_indeces[0], local_axarr_indeces[1]].axis(x_lims + y_lims) #combined limits from Walker plot
    plt.suptitle(master_title)
    if show:
        fig.tight_layout(rect=[0, 0.03, 1.0, 0.95])
        plt.show()
    if save:
        axarr[local_axarr_indeces[0], local_axarr_indeces[1]].savefig(pltFileDir + fig_name)
    if return_prof:
        if show_contours:
            return dsph_brightness_profile
        else:
            print ('No contour profile generated, so no profile to return. ')
    else:
        return 1
