#This file generates a distribution of random stars, modeled after some real distribution.
# The purpose is to see how likely it is that a given structure (eg, flattening), could have arisen by sheer chance.

import math
import numpy as np
import random 

def generateRandomStarDistribution(ref_pop_RA, ref_pop_Dec, prob_funct_method = 'simple', prob_funct_shape = 'spherical'):
    ref_pop_radial_dist = [ math.sqrt( ref_pop_RA[i] ** 2 + ref_pop_Dec[i] **2 ) for i in range(len(ref_pop_RA)) ]
    N_stars = len(ref_pop_RA)
    
    if prob_funct_method == 'simple':
        prob_density_funct = lambda 
    if prob_funct_method == 'pop_weighted':
        #Do a polynomial fit to radial distribution of actual population radial distribution, then divide by R 
