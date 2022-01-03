from scipy.interpolate import RegularGridInterpolator
import numpy as np
import math
from PotentialFunctionArray import PotentialFunctionArray
from ComputationalArchive import ComputationalArchive
from PotentialArchive import PotentialArchive 

class FullPotentialInterpolator:

    def getPotentialInterpolator(self):
        return self.potential_interpolator 

    def generatePotentialInterpolator(self, Rlikeaxis, zlikeaxis, varying_parameter, potential_functions):
    
        #potential_meshes = [function.potential_mesh.transpose() for function in potential_functions]

        potential_meshes = np.array(potential_meshes) 
    
        interpolator = RegularGridInterpolator((Rlikeaxis, zlikeaxis, varying_parameter), potential_meshes.transpose())

        return interpolator


    def __init__(self, potential_type, generate_new_interpolator = 'no', save_interpolator = 'no', interpolator_name = None, params_to_interpolate = None):
        potential_type = potential_type.lower()
        generate_new_interpolator = generate_new_interpolator.lower()
        
        comp_array = ComputationalArchive()
        pot_archive = PotentialArchive() 
        potential_dir = comp_array.getPotentialDir()
        interpolator_suffix = '_full_interpolator.npy'
        potential_info = {'nfw': [potential_dir + 'nfw'+interpolator_suffix],
                          'cored': [potential_dir + 'cored'+interpolator_suffix],
                          'sech_disk': [potential_dir + 'sech_disk'+interpolator_suffix]}
        
        if interpolator_name:
            potential_info[potential_type][0] = interpolator_name 
        if 'generate_new_interpolator' == 'no':
            self.potential_interpolator = np.load(potential_info[potential_type][0]).item()
        else:
            potential_files = [ pot_archive.getGeneralFile(param, potential_type) for param in params_to_interpolate ]
            potential_functions = [ PotentialFunctionArray(file) for file in potential_files ]
            self.potential_interpolator = self.generatePotentialInterpolator (potential_functions[0].Rlikeaxis,potential_functions[0].zlikeaxis,params_to_interpolate, potential_functions)
