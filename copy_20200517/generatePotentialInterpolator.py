from scipy.interpolate import RegularGridInterpolator
import numpy as np
import math
from PotentialFunctionArray import PotentialFunctionArray

def generatePotentialInterpolator(varying_parameter, potential_functions):
    
    potential_meshes = [function.potential_mesh.transpose() for function in potential_functions]

    potential_meshes = np.array(potential_meshes)

    Rlikeaxis = potential_functions[0].Rlikeaxis
    zlikeaxis = potential_functions[0].zlikeaxis
    
    interpolator = RegularGridInterpolator((Rlikeaxis, zlikeaxis, varying_parameter), potential_meshes.transpose())

    return interpolator 
