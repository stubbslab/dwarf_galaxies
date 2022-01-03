import math
import numpy as np
from AstronomicalParameterArchive import AstonomicalParameterArchive
from DwarfGalDataArchive import DwarfGalDataArchive
from ObservedGalaxyStarData import ObservedGalaxyStarData
from distributeRandomStars import distributeRandomStars
from SurfaceBrightnessProfile import SurfaceBrightnessProfile
from PotentialFunctionArray import PotentialFunctionArray
from compRZLos import compRZLos

def runMCMCForArtificialDwarfGalaxies (gal_to_replicate, generation_params_storer, generation_Rbin_walls, generatio_zbin_walls, withDisk = 0, pop_selection_method = 'None', MCMC_el = 'None', MCMC_lam = 'None' , generation_disk_funct = 0, generation_halo_funct = 0, disk_funct = 0, halo_funct = 0, start_param_index = 0, arcmin_limits_R = 60.0, arcmin_limits_z = 60.0, outer_step = 1.0):
    astro_archive = AstronomicalParameterArchiv()
    dSph_archive = DwarfGalDataArchiv()
    potetial_archive = PotentialArchive()

    deg_to_rad = astro_archive.getDegToRad()
    gamma = astro_archive.getGamma() 

    generation_disk_file = potential_archive.getDiskFile(generation_params_storer.lam)
    generation_halo_file = potential_archive.getHalofile(generation_params_storer.el)

    if not generation_disk_funct:
        generation_disk_funct = PotentialFunctionArray(generation_disk_file)
    if not generation_halo_funct:
        generation_halo_funct = PotentialFunctionArray(generation_halo_file) 

    start_params = start_param_storer.getStartParameters(start_param_index, withDisk) 
    
    if MCMC_el == 'None':
        MCMC_el = generation_params_storer.el
    if MCMC_lam == 'None':
        MCMC_lam = generation_params_storer.lam

    pops = dSph_archive.getPopulations(gal_to_replicate)
    print 'pops = '
    print pops

    gen_R,gen_z,gen_los = compRZLos(generation_params_storer.zeta, inner_step = outer_step *  generation_param_storer.zeta,
                                    outer_step, arcmin_limits_R, arcmin_limits_z, generation_param_storer.rs, generation_param_storer.dist) 
    real_star_data = []
    popGalPairs = []
    start_dsph_parameter_storers  = []
    generation_params_storers = []
    for i in range(len(pops)):
        pop = pops[i]
        popGalPairs = popGalPairs + [[galaxy,pop]]
        real_star_data = real_star_data + [ObservedGalaxyStarData([galaxy,pop], pop_selection_method = pop_selection_method)]
        start_dsph_parameter_storers = start_dsph_parameter_storers + [
            DwarfGalaxyParametersStorer([galaxy,pop],withDisk=withDisk,el=MCMC_el, lam=MCMC_lam,
                                        rs = start_params[0], phi = start_params[1], theta = start_params[2],
                                        zeta = start_params[3], eps = start_params[4], b= start_params[5], sigsqr = real_star_data[i].sigSqr)
             ]
        generation_params_storers = generation_params_storers + [generation_params_storer]
        generation_params_storers[i].sigsqr = real_star_data[i].sigsqr
        generation_surface_profiles = generation_surface_profiles + [SurfaceBrightnessProfile(gen_R, gen_z, gen_los, gamma, generation_params_storers[i],
                                                                                              generation_disk_file, generation_halo_file,
                                                                                              disk_interpolating_function = generation_disk_funct,
                                                                                              halo_interpolating_function = generation_halo_funct)]
    artificial_star_data = []
    
    for i in range(len(pops)):
        artificial_star_RAs_in_scale_radii, artificial_star_Decs_in_scale_radii = distributeRandomStars(len(real_star_data[i].corr_Ra), generation_surface_profiles[i], generation_Rbin_walls, generation_zbin_walls)
        artificial_star_RAs_in_deg = [RA * generation_params_storers[i].rs / generation_params_storers[i].dist * 1 / deg_to_rad for RA in artificial_star_RAs_in_scale_radii]
        artificial_star_Decs_in_deg = [Dec * generation_params_storers[i].rs / generation_params_storers[i].dist * 1 / deg_to_rad for Dec in artificial_star_Decs_in_scale_radii]
        artificial_star_data = artificial_star_data + [artificial_star_RAs_in_deg,artificial_star_Decs_in_deg]

        

    #So now we have artificial star RAs and Decs to take the place of our real RAs and Decs.  Everything else should follow from here identically to how it did before.  
