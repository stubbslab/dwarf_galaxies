import math 
import numpy as np
from AstronomicalParameterArchive import AstronomicalParameterArchive
from DwarfGalaxyParametersStorer import DwarfGalaxyParametersStorer
from PotentialArchive import PotentialArchive
from ObservedGalaxyStarData import ObservedGalaxyStarData
from PotentialFunctionArray import PotentialFunctionArray
from SurfaceBrightnessProfile import SurfaceBrightnessProfile
from VisitedDsphMCMCPoint import VisitedDsphMCMCPoint
from logList import logList 
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm


def showDSphProfiles(populations, show_stars = 1, withDisk = 0, disk_funct = 0, halo_funct = 0, pop_selection_method='none',
                    dist = 'none', M = 'none', vphi = 'none', sigsqr_RR = 'none', dispersion_b = 0.0, gamma0 = 0.0, gammainf = 0.0, beta_exp = 1.0, beta_R0 = 1.0, 
                    el=0.3, rs = 350, phi = math.pi*0.5, theta = 0.0, c = 5.0, 
                    lam=0.2, zeta = 0.2, eps = 0.05, a = math.pi*0.8, b = math.pi*0.5):

    e = math.sqrt(1-(1-el)**2) # the eccentricity of the halo (this is the e that appears in my expression for ...) 
    #f = math.log(1+c) - c/(1+c)
    astro_archive = AstronomicalParameterArchive ()
    potential_archive = PotentialArchive()
    gamma = astro_archive.getGamma()
    

    
    disk_file = potential_archive.getDiskFile(lam)
    halo_file = potential_archive.getHaloFile(el)
    if not disk_funct: 
        print ('Creating the disk_funct from file: ' + disk_file)
        disk_funct = PotentialFunctionArray(disk_file)
    if not halo_funct:
        print ('Creating the halo_funct from file: ' + halo_file)
        halo_funct = PotentialFunctionArray(halo_file)
    print ('Halo and disk interpolating functions created. ')

    arcmin_limits_R = 60.0
    arcmin_limits_z = 60.0
    R_inner=np.arange(0,arcmin_limits_R*0.05 + 0.1,0.1)
    z_inner=np.arange(0,arcmin_limits_z*0.05 + 0.1,0.1)
    R_outer=np.arange(0,arcmin_limits_R + 0.5,0.5)
    z_outer=np.arange(0,arcmin_limits_z + 0.5,0.5)
    R_arcmin = np.unique(np.concatenate((-R_inner,-R_outer,R_inner,R_outer)))
    z_arcmin = np.unique(np.concatenate((-z_inner,-z_outer,z_inner,z_outer)))
    print ('max, min R in arcminutes are: ' + str(max(R_outer)) + ', ' + str(min(R_outer)))
    print ('max, min z in arcminutes are: ' + str(max(z_outer)) + ', ' + str(min(z_outer)))
    print ('dist = ' + str(dist))
    R=R_arcmin*1/60.0*1/180.0*math.pi*dist/rs
    z=z_arcmin*1/60.0*1/180.0*math.pi*dist/rs
    
    los_bins = np.unique(np.around(np.concatenate((np.arange(max(R)*-1.0,max(R)*1.0 + 0.25,0.25,dtype=np.float),np.arange(max(R)*-0.5*0.1,max(R)*0.5*0.1+0.05,0.05,dtype=np.float))),5))
    #print 'los_bins = ' 
    #print los_bins
    #print ' R = '
    #print R
    #print 'z = '
    #print z 
    
    
    fig=plt.figure(1,figsize=(11,9))
    for population in populations: 
        star_data = ObservedGalaxyStarData(population, pop_selection_method = pop_selection_method)
        print ('a 0 = ' + str(a))
        print ('b 0 = ' + str(b) )
        dsph_parameter_storer = DwarfGalaxyParametersStorer(population,withDisk=withDisk,
                                                            dist = dist, M = M, vphi = vphi, sigsqr_RR = 'none', dispersion_b = 0.0, gamma0 = 0.0, gammainf = 0.0, beta_exp = 1.0, beta_R0 = 1.0, 
                                                            el=el, rs = rs, phi = phi, theta = theta, c = c, 
                                                            lam=lam, zeta = zeta, eps = eps, a = a, b = b)
        dsph_brightness_profile = SurfaceBrightnessProfile(R, z, los_bins, gamma, dsph_parameter_storer,
                                                                  disk_file, halo_file, 
                                                                  disk_interpolating_function = disk_funct,
                                                                  halo_interpolating_function = halo_funct       ) 



        c = dsph_parameter_storer.c
        dist = dsph_parameter_storer.dist
        rs = dsph_parameter_storer.rs 
        pop_of_interest = population[1] # can be 'MP', 'IM', 'MR'
        colormap='b' if pop_of_interest == 'MP' else 'g' if pop_of_interest == 'IM' else 'r'

        corr_ra_one_pop = star_data.corrRa
        corr_dec_one_pop = star_data.corrDec
        print ('NStars for population ' + population[0] + ' ' + population[1] + ' is ' + str(len(corr_dec_one_pop))())
    
        if show_stars:
            plt.scatter(corr_ra_one_pop * 60.0,corr_dec_one_pop * 60.0 ,s=4.,color=colormap)
        
    zmesh,Rmesh=np.meshgrid(dsph_brightness_profile.zlikeaxis,dsph_brightness_profile.Rlikeaxis)
    zmesh = zmesh*180.0/math.pi*rs/dist*60.0
    Rmesh = Rmesh*180.0/math.pi*rs/dist*60.0
    #print 'a 2 = ' + str(dsph_brightness_profile.a)
    #print 'b 2 = ' + str(dsph_brightness_profile.b) 
    brightness_mesh = dsph_brightness_profile.surfaceBrightness
    log_levels=logList(np.min(np.abs(brightness_mesh)),np.max(np.abs(brightness_mesh)),20)
    CS=plt.contour(Rmesh,zmesh,brightness_mesh,norm=LogNorm(),levels=log_levels)
    CB_contour=plt.colorbar(shrink=0.8,extend='both',format='%.2e')

    matplotlib.rcParams['xtick.direction'] = 'out'
    matplotlib.rcParams['ytick.direction'] = 'out'

    plt.ylim = (-20,20) #ylimits from Walker plot
    plt.xlim = (20,-20) #xlimits from Walker plot
    plt.xlabel('ra sep')
    plt.ylabel('dec sep')
    plt.title('Fornax Stars Position on Sky')
    pltFileDir="/Users/sasha/Documents/Harvard/physics/randall/dwarfGalaxyProject/plots/"
    #Uncomment the following line if you want to save the file
    #plt.savefig(pltFileDir + gal_of_interest + '_' + pop_of_interest + '_metal_cuts_' + str((metallicity_cuts[gal_of_interest])[0]) + '_' + str((metallicity_cuts[gal_of_interest])[1]) + '.png')
    #plt.savefig(pltFileDir + gal_of_interest + '_' + pop_of_interest + '_full_pop_division' + '.png')
    plt.show()
