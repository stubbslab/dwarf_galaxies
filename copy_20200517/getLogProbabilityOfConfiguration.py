#Write MCMC type procedure for comparing stellar distributions 
#Due to computational constraints, assumes a fixed value of
#  halo ellipticity and disk height-to-radius scale.
#Once I get this working, perhaps we can generalize

from ObservedGalaxyStarData import ObservedGalaxyStarData
from PotentialFunctionArray import PotentialFunctionArray
#from generateGalaxyDensityProfile import GalaxyDensityProfile
from PotentialArchive import PotentialArchive 
from SurfaceBrightnessProfile import SurfaceBrightnessProfile 
from DwarfGalaxyParametersStorer import DwarfGalaxyParametersStorer
#from DwarfGalaxyParametersStorer import getMCMCInfo
from AstronomicalParameterArchive import AstronomicalParameterArchive
from VisitedDsphMCMCPoint import VisitedDsphMCMCPoint
from logList import logList 
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import copy 
import random
import math 
import numpy as np 

def getLogProbabilityOfConfiguration(population, withDisk = 0, disk_funct = 0, halo_funct = 0, pop_selection_method='none',
                    dist = 'none', M = 'none', vphi = 'none', sigsqr = 'none',
                    el=0.3, rs = 350, phi = math.pi*0.5, theta = 0.0, c = 5.0, 
                    lam=0.2, zeta = 0.2, eps = 0.05, a = math.pi*0.8, b = math.pi*0.5):
    e = math.sqrt(1-(1-el)**2) # the eccentricity of the halo (this is the e that appears in my expression for ...) 
    f = math.log(1+c) - c/(1+c)
    astro_archive = AstronomicalParameterArchive ()
    star_data = ObservedGalaxyStarData(population, pop_selection_method = pop_selection_method)
    dsph_parameter_storer = DwarfGalaxyParametersStorer(population, withDisk, 
                    dist = dist, M = M, vphi = vphi, sigsqr = sigsqr,
                    el=el, rs = rs, phi = phi, theta = theta, c = c, 
                    lam=lam, zeta = zeta, eps = eps, a = a, b = b)
    e = math.sqrt(1-(1-el)**2) # the eccentricity of the halo (this is the e that appears
    potential_archive = PotentialArchive()
    disk_file = potential_archive.getDiskFile(lam)
    halo_file = potential_archive.getHaloFile(el)
    if not disk_funct: 
        print 'Creating the disk_funct from file: ' + disk_file
        disk_funct = PotentialFunctionArray(disk_file)
    if not halo_funct:
        print 'Creating the halo_funct from file: ' + halo_file
        halo_funct = PotentialFunctionArray(halo_file)
    print 'Halo and disk interpolating functions created. '

    c = dsph_parameter_storer.c
    dist = dsph_parameter_storer.dist
    rs = dsph_parameter_storer.rs 
    gamma = astro_archive.getGamma()
    #(param_indeces_to_vary, param_ranges, gauss_sampling_width) =  start_dsph_parameter_storer.getMCMCInfo()
    #print 'start_dsph_parameter_storer.dist = ' + str(start_dsph_parameter_storer.dist)
    
    arcmin_limits_R = 60.0
    arcmin_limits_z = 60.0
    R_inner=np.arange(0,arcmin_limits_R*0.1 + 0.5,0.5)
    z_inner=np.arange(0,arcmin_limits_z*0.1 + 0.5,0.5)
    R_outer=np.arange(0,arcmin_limits_R + 2.0,2.0)
    z_outer=np.arange(0,arcmin_limits_z + 2.0,2.0)
    R_arcmin = np.unique(np.concatenate((-R_inner,-R_outer,R_inner,R_outer)))
    z_arcmin = np.unique(np.concatenate((-z_inner,-z_outer,z_inner,z_outer)))
    print 'max, min R in arcminutes are: ' + str(max(R_outer)) + ', ' + str(min(R_outer))
    print 'max, min z in arcminutes are: ' + str(max(z_outer)) + ', ' + str(min(z_outer))
    R=R_arcmin*1/60.0*1/180.0*math.pi*dist/rs
    z=z_arcmin*1/60.0*1/180.0*math.pi*dist/rs
    
    los_bins = np.unique(np.around(np.concatenate((np.arange(max(R)*-1.0,max(R)*1.0 + 0.25,0.25,dtype=np.float),np.arange(max(R)*-0.5*0.1,max(R)*0.5*0.1+0.05,0.05,dtype=np.float))),5))
    #print 'los_bins = ' 
    #print los_bins
    #print ' R = '
    #print R
    #print 'z = '
    #print z 
    
    dsph_brightness_profile = SurfaceBrightnessProfile(R, z, los_bins, gamma, dsph_parameter_storer,
                                                                  disk_file, halo_file, 
                                                                  disk_interpolating_function = disk_funct,
                                                                  halo_interpolating_function = halo_funct       )
    print 'The total log likelihood of the configuration is ' + str(dsph_brightness_profile.sumLogSurfaceBrightness(star_data.corrRa,star_data.corrDec))
    
    return dsph_brightness_profile.starProbabilityValues(star_data.corrRa,star_data.corrDec)
    
