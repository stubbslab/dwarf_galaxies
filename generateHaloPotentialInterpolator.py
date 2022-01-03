from scipy.interpolate import RegularGridInterpolator
import numpy as np
import math
from PotentialFunctionArray import PotentialFunctionArray

def generatePotentialInterpolator(Rlikeaxis, zlikeaxis, varying_parameter, potential_functions):

    np.array()
    
    potential_meshes = [function.surfaceProbDensity for function in potential_functions]

    potential_meshes = np.array(potential_meshes) 
    
    interpolator = RegularGridInterpolator((Rlikeaxis, zlikeaxis, varying_parameter), potential_functions.transpose())

    return interpolator 
