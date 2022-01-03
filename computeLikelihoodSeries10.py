#Do a sequence of the MCMC Series algorithms, with varying run conditions.
#You can change various parameters (present set listed in commented line right above parameter_serires assignment), and have program execute in sequence
#This script can be run from command line without having to be in python environment (bash$ python doMCMCSeriesAllPops.py)
#This method tries to be as efficient as possible, reading in the various potentialFunctionArrays only once, as that can be a time consuming part of the process.  

from runMCMCForDwarfGalaxyProfilesAllPops import runMCMCForDwarfGalaxyProfilesAllPops
from PotentialFunctionArray import PotentialFunctionArray
from PotentialArchive import PotentialArchive
from SurfaceBrightnessProfile import SurfaceBrightnessProfile
from DwarfGalaxyParametersStorer import DwarfGalaxyParametersStorer
from AstronomicalParameterArchive import AstronomicalParameterArchive
from DwarfGalDataArchive import DwarfGalDataArchive
from ObservedGalaxyStarData import ObservedGalaxyStarData 
import math
import numpy as np 

if __name__ == '__main__':
    #parameters at present are: [el, lam, galaxy, number of iterations, population selection method, start parameter index]
    print 'Starting. '
    astro_archive = AstronomicalParameterArchive ()
    dSph_archive = DwarfGalDataArchive()
    
    saverun = 1
    save_array = []
    results_number = '10'
    gamma = astro_archive.getGamma() 
    
    el_index = 0
    rs_index = 1
    phi_index = 2
    theta_index = 3
    c_index = 4
    lam_index = 5
    zeta_index = 6
    eps_index = 7
    a_index = 8
    b_index = 9
    gal_index = 10
    pop_select_index = 11
    withDisk_index = 12
    
    #                     el    rs   phi   thet   c    lam  zeta eps     a             b               gal      pop_select     WD
    parameter_series = [ [0.2 , 258, 0.65, -0.05, 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.2 , 263, 0.65, -0.05, 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.2 , 268, 0.65, -0.05, 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.2 , 273, 0.65, -0.05, 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.2 , 278, 0.65, -0.05, 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.2 , 258, 0.7 , -0.05, 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.2 , 263, 0.7 , -0.05, 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.2 , 268, 0.7 , -0.05, 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.2 , 273, 0.7 , -0.05, 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.2 , 278, 0.7 , -0.05, 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.2 , 258, 0.75, -0.05, 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.2 , 263, 0.75, -0.05, 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.2 , 268, 0.75, -0.05, 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.2 , 273, 0.75, -0.05, 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.2 , 278, 0.75, -0.05, 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.2 , 258, 0.65, 0.0  , 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.2 , 263, 0.65, 0.0  , 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.2 , 268, 0.65, 0.0  , 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.2 , 273, 0.65, 0.0  , 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.2 , 278, 0.65, 0.0  , 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.2 , 258, 0.7 , 0.0  , 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.2 , 263, 0.7 , 0.0  , 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.2 , 268, 0.7 , 0.0  , 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.2 , 273, 0.7 , 0.0  , 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.2 , 278, 0.7 , 0.0  , 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.2 , 258, 0.75, 0.0  , 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.2 , 263, 0.75, 0.0  , 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.2 , 268, 0.75, 0.0  , 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.2 , 273, 0.75, 0.0  , 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.2 , 278, 0.75, 0.0  , 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.2 , 258, 0.65, 0.05 , 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.2 , 263, 0.65, 0.05 , 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.2 , 268, 0.65, 0.05 , 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.2 , 273, 0.65, 0.05 , 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.2 , 278, 0.65, 0.05 , 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.2 , 258, 0.7 , 0.05 , 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.2 , 263, 0.7 , 0.05 , 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.2 , 268, 0.7 , 0.05 , 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.2 , 273, 0.7 , 0.05 , 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.2 , 278, 0.7 , 0.05 , 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.2 , 258, 0.75, 0.05 , 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.2 , 263, 0.75, 0.05 , 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.2 , 268, 0.75, 0.05 , 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.2 , 273, 0.75, 0.05 , 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.2 , 278, 0.75, 0.05 , 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.25, 275, 0.65, -0.05, 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.25, 280, 0.65, -0.05, 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.25, 285, 0.65, -0.05, 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.25, 290, 0.65, -0.05, 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.25, 295, 0.65, -0.05, 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.25, 275, 0.7 , -0.05, 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.25, 280, 0.7 , -0.05, 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.25, 285, 0.7 , -0.05, 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.25, 290, 0.7 , -0.05, 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.25, 295, 0.7 , -0.05, 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.25, 275, 0.75, -0.05, 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.25, 280, 0.75, -0.05, 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.25, 285, 0.75, -0.05, 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.25, 290, 0.75, -0.05, 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.25, 295, 0.75, -0.05, 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.25, 275, 0.65, 0.0  , 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.25, 280, 0.65, 0.0  , 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.25, 295, 0.65, 0.0  , 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.25, 290, 0.65, 0.0  , 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.25, 295, 0.65, 0.0  , 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.25, 275, 0.7 , 0.0  , 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.25, 280, 0.7 , 0.0  , 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.25, 285, 0.7 , 0.0  , 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.25, 290, 0.7 , 0.0  , 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.25, 295, 0.7 , 0.0  , 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.25, 275, 0.75, 0.0  , 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.25, 280, 0.75, 0.0  , 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.25, 285, 0.75, 0.0  , 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.25, 290, 0.75, 0.0  , 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.25, 295, 0.75, 0.0  , 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.25, 275, 0.65, 0.05 , 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.25, 280, 0.65, 0.05 , 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.25, 285, 0.65, 0.05 , 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.25, 290, 0.65, 0.05 , 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.25, 295, 0.65, 0.05 , 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.25, 275, 0.7 , 0.05 , 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.25, 280, 0.7 , 0.05 , 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.25, 285, 0.7 , 0.05 , 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.25, 290, 0.7 , 0.05 , 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.25, 295, 0.7 , 0.05 , 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.25, 275, 0.75, 0.05 , 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.25, 280, 0.75, 0.05 , 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.25, 285, 0.75, 0.05 , 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.25, 290, 0.75, 0.05 , 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.25, 295, 0.75, 0.05 , 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.3 , 283, 0.65, -0.05, 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.3 , 288, 0.65, -0.05, 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.3 , 293, 0.65, -0.05, 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.3 , 298, 0.65, -0.05, 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.3 , 303, 0.65, -0.05, 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.3 , 283, 0.7 , -0.05, 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.3 , 288, 0.7 , -0.05, 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.3 , 293, 0.7 , -0.05, 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.3 , 298, 0.7 , -0.05, 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.3 , 303, 0.7 , -0.05, 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.3 , 283, 0.75, -0.05, 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.3 , 288, 0.75, -0.05, 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.3 , 293, 0.75, -0.05, 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.3 , 298, 0.75, -0.05, 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.3 , 303, 0.75, -0.05, 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.3 , 283, 0.65, 0.0  , 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.3 , 288, 0.65, 0.0  , 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.3 , 293, 0.65, 0.0  , 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.3 , 298, 0.65, 0.0  , 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.3 , 303, 0.65, 0.0  , 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.3 , 283, 0.7 , 0.0  , 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.3 , 288, 0.7 , 0.0  , 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.3 , 293, 0.7 , 0.0  , 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.3 , 298, 0.7 , 0.0  , 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.3 , 303, 0.7 , 0.0  , 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.3 , 283, 0.75, 0.0  , 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.3 , 288, 0.75, 0.0  , 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.3 , 293, 0.75, 0.0  , 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.3 , 298, 0.75, 0.0  , 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.3 , 303, 0.75, 0.0  , 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.3 , 283, 0.65, 0.05 , 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.3 , 288, 0.65, 0.05 , 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.3 , 293, 0.65, 0.05 , 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.3 , 298, 0.65, 0.05 , 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.3 , 303, 0.65, 0.05 , 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.3 , 283, 0.7 , 0.05 , 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.3 , 288, 0.7 , 0.05 , 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.3 , 293, 0.7 , 0.05 , 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.3 , 298, 0.7 , 0.05 , 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.3 , 303, 0.7 , 0.05 , 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.3 , 283, 0.75, 0.05 , 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.3 , 288, 0.75, 0.05 , 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.3 , 293, 0.75, 0.05 , 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.3 , 298, 0.75, 0.05 , 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.3 , 303, 0.75, 0.05 , 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.4 , 298, 0.65, -0.05, 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.4 , 303, 0.65, -0.05, 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.4 , 308, 0.65, -0.05, 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.4 , 313, 0.65, -0.05, 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.4 , 318, 0.65, -0.05, 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.4 , 298, 0.7 , -0.05, 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.4 , 303, 0.7 , -0.05, 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.4 , 308, 0.7 , -0.05, 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.4 , 313, 0.7 , -0.05, 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.4 , 318, 0.7 , -0.05, 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.4 , 298, 0.75, -0.05, 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.4 , 303, 0.75, -0.05, 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.4 , 308, 0.75, -0.05, 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.4 , 313, 0.75, -0.05, 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.4 , 318, 0.75, -0.05, 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.4 , 298, 0.65, 0.0  , 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.4 , 303, 0.65, 0.0  , 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.4 , 308, 0.65, 0.0  , 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.4 , 313, 0.65, 0.0  , 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.4 , 318, 0.65, 0.0  , 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.4 , 298, 0.7 , 0.0  , 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.4 , 303, 0.7 , 0.0  , 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.4 , 308, 0.7 , 0.0  , 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.4 , 313, 0.7 , 0.0  , 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.4 , 318, 0.7 , 0.0  , 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.4 , 298, 0.75, 0.0  , 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.4 , 303, 0.75, 0.0  , 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.4 , 308, 0.75, 0.0  , 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.4 , 313, 0.75, 0.0  , 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.4 , 318, 0.75, 0.0  , 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.4 , 298, 0.65, 0.05 , 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.4 , 303, 0.65, 0.05 , 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.4 , 308, 0.65, 0.05 , 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.4 , 313, 0.65, 0.05 , 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.4 , 318, 0.65, 0.05 , 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.4 , 298, 0.7 , 0.05 , 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.4 , 303, 0.7 , 0.05 , 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.4 , 308, 0.7 , 0.05 , 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.4 , 313, 0.7 , 0.05 , 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.4 , 318, 0.7 , 0.05 , 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.4 , 298, 0.75, 0.05 , 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.4 , 303, 0.75, 0.05 , 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.4 , 308, 0.75, 0.05 , 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.4 , 313, 0.75, 0.05 , 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.4 , 318, 0.75, 0.05 , 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.5 , 270, 0.65, -0.05, 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.5 , 275, 0.65, -0.05, 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.5 , 280, 0.65, -0.05, 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.5 , 285, 0.65, -0.05, 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.5 , 290, 0.65, -0.05, 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.5 , 270, 0.7 , -0.05, 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.5 , 275, 0.7 , -0.05, 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.5 , 280, 0.7 , -0.05, 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.5 , 285, 0.7 , -0.05, 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.5 , 290, 0.7 , -0.05, 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.5 , 270, 0.75, -0.05, 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.5 , 275, 0.75, -0.05, 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.5 , 280, 0.75, -0.05, 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.5 , 285, 0.75, -0.05, 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.5 , 290, 0.75, -0.05, 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.5 , 270, 0.65, 0.0  , 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.5 , 275, 0.65, 0.0  , 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.5 , 280, 0.65, 0.0  , 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.5 , 285, 0.65, 0.0  , 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.5 , 290, 0.65, 0.0  , 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.5 , 270, 0.7 , 0.0  , 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.5 , 275, 0.7 , 0.0  , 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.5 , 280, 0.7 , 0.0  , 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.5 , 285, 0.7 , 0.0  , 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.5 , 290, 0.7 , 0.0  , 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.5 , 270, 0.75, 0.0  , 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.5 , 275, 0.75, 0.0  , 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.5 , 280, 0.75, 0.0  , 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.5 , 285, 0.75, 0.0  , 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.5 , 290, 0.75, 0.0  , 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.5 , 270, 0.65, 0.05 , 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.5 , 275, 0.65, 0.05 , 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.5 , 280, 0.65, 0.05 , 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.5 , 285, 0.65, 0.05 , 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.5 , 290, 0.65, 0.05 , 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.5 , 270, 0.7 , 0.05 , 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.5 , 275, 0.7 , 0.05 , 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.5 , 280, 0.7 , 0.05 , 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.5 , 285, 0.7 , 0.05 , 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.5 , 290, 0.7 , 0.05 , 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.5 , 270, 0.75, 0.05 , 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.5 , 275, 0.75, 0.05 , 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.5 , 280, 0.75, 0.05 , 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.5 , 285, 0.75, 0.05 , 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0],
                         [0.5 , 290, 0.75, 0.05 , 5.0, 0.2, 1.0,  0.0,   -math.pi*0.0, math.pi * 0.5, 'fornax', 'metal_point', 0] 
                         ]


    pot_archive = PotentialArchive()
    unique_els = []
    unique_lams = []
    #Create all of the potential files used in the MCMC algorithms.
    # That way, we don't need to reread in a potential file every time it is used.  
    for parameter_set in parameter_series:
        if parameter_set[el_index] not in unique_els:
            unique_els = unique_els + [parameter_set[el_index]]
        if parameter_set[lam_index] not in unique_lams:
            unique_lams = unique_lams + [parameter_set[lam_index]]
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


    arcmin_limits_R = 60.0
    arcmin_limits_z = 60.0
    R_inner=np.arange(0,arcmin_limits_R*0.1 + 0.05,0.05)
    z_inner=np.arange(0,arcmin_limits_z*0.1 + 0.05,0.05)
    R_outer=np.arange(0,arcmin_limits_R + 1.0,1.0)
    z_outer=np.arange(0,arcmin_limits_z + 1.0,1.0)
    R_arcmin = np.unique(np.around(np.concatenate( (-R_outer,-R_inner,R_inner,R_outer) ),5))
    z_arcmin = np.unique(np.around(np.concatenate( (-z_outer,-z_inner,z_inner,z_outer) ),5))
    R_arcmin = np.sort(R_arcmin)
    z_arcmin = np.sort(z_arcmin)
    likelihood_array = []
    save_array = []
    #Now actually compute the likelihood for each of the specified point. 
    for parameter_set in parameter_series:
        print 'Computing likelihood for parameter_set '
        print parameter_set
        parameter_storers = []
        surface_profiles = []
        
        galaxy = parameter_set[gal_index]
        withDisk = parameter_set[withDisk_index]
        el = parameter_set[el_index]
        rs = parameter_set[rs_index]
        phi = parameter_set[phi_index]
        theta = parameter_set[theta_index]
        c = parameter_set[c_index]
        lam = parameter_set[lam_index]
        zeta = parameter_set[zeta_index]
        eps = parameter_set[eps_index]
        a = parameter_set[a_index]
        b = parameter_set[b_index]
        pop_selection_method = parameter_set[pop_select_index]
        halo_file = pot_archive.getHaloFile(el)
        disk_file = pot_archive.getDiskFile(lam)
        pops = dSph_archive.getPopulations(galaxy)
        print 'pops = '
        print pops

        #print 'max, min R in arcminutes are: ' + str(max(R_outer)) + ', ' + str(min(R_outer))
        #print 'max, min z in arcminutes are: ' + str(max(z_outer)) + ', ' + str(min(z_outer))
        storer = 0
        log_likelihood = 0.0
        for pop in pops:
            population = [galaxy,pop]
            star_data = ObservedGalaxyStarData(population, pop_selection_method = pop_selection_method)
            print 'for population ' + pop + ', len(star_data.corrRa) ' + str(len(star_data.corrRa))
            storer =  DwarfGalaxyParametersStorer(population, withDisk, sigsqr = star_data.sigSqr, 
                                                  el = el, rs = rs, phi = phi, theta = theta, c = c,
                                                  lam = lam, zeta = zeta, eps = eps, a = a, b = b)
            print 'storer.sigsqr = ' + str(storer.sigsqr) + ' for population ' + pop 

            dist = storer.dist
            R=R_arcmin*1/60.0*1/180.0*math.pi*dist/rs
            z=z_arcmin*1/60.0*1/180.0*math.pi*dist/rs
            los_bins_outer = np.arange(max(R)*-1.0,max(R)*1.0 + 0.25,0.25,dtype=np.float)
            los_bins_inner = np.arange(max(R)*-0.1,max(R)*0.1+0.05,0.05,dtype=np.float)
            los_bins = np.unique(np.around(np.concatenate( (los_bins_outer,los_bins_inner) ),5))

            parameter_storers = parameter_storers + [storer]
            surface_profile = SurfaceBrightnessProfile(R, z, los_bins, gamma, storer,
                                                       disk_file, halo_file,
                                                       disk_interpolating_function = disk_funct_array[lam],
                                                       halo_interpolating_function = halo_funct_array[el] )
            if pop == 'MR':
                    print 'pop == ' + pop 
                    print 'surface_profile.Rlikeaxis  = '
                    print surface_profile.Rlikeaxis
                    print 'surface_profile.zlikeaxis  = '
                    print surface_profile.zlikeaxis    
                    brightness_mesh = surface_profile.surfaceBrightness
                    print 'surface_profile.surfaceBrightness = '
                    print surface_profile.surfaceBrightness
            surface_profiles = surface_profiles + [surface_profile]
            log_likelihood = log_likelihood + surface_profile.sumLogSurfaceBrightness(star_data.corrRa,star_data.corrDec)
        save_array = save_array + [[storer.dist, storer.M, storer.vphi, storer.sigsqr] + parameter_set + [log_likelihood]]
        likelihood_array = likelihood_array + [log_likelihood]


    #print likelihood_array
    #save_array = [parameter_series[i] + [likelihood_array[i]] for i in range(len(likelihood_array)) ]
    #print save_array 
    #print np.array(save_array)

    if saverun:
        save_dir = '/Users/sasha/Documents/Harvard/physics/randall/probabilityTables/'
        file_name = 'spot_computed_probabilities_number_' + str(results_number) + '.csv'
        print 'Saving series to ' + save_dir + file_name
        header = 'dist, M, vphi, sigsqr ,el, rs , phi , theta , c, lam, zeta, eps, a, b, galaxy, pop_select, withDisk, logLikelihood' 
        np.savetxt(save_dir + file_name, np.array(save_array), delimiter = ',',fmt = '%10s',header = header) 
         
