from PotentialArchive import PotentialArchive
import numpy as np
from PotentialFunctionArray import PotentialFunctionArray
from ComputationalArchive import ComputationalArchive 
from scipy.interpolate import RegularGridInterpolator


class PreGeneratedPotentialInterpolator:

    def saveToFile(self,file_name = None):
        computational_params_archive = ComputationalArchive()
        potential_dir = computational_params_archive.getPotentialDir()
        if file_name is None:
            file_name = self.potential_type + '_potential_interpolator.npy'
        print (file_name )
        print ('file_name = ' + file_name)
        print ('potential_dir = '+ potential_dir )
        np.save(potential_dir + file_name, self)

        return 1.0


    def generatePotentialInterpolator(self, varying_parameter, potential_functions):
    
        potential_meshes = [function.potential_mesh.transpose() for function in potential_functions]

        potential_meshes = np.array(potential_meshes)

        Rlikeaxis = potential_functions[0].Rlikeaxis
        zlikeaxis = potential_functions[0].zlikeaxis
    
        interpolator = RegularGridInterpolator((Rlikeaxis, zlikeaxis, varying_parameter), potential_meshes.transpose())

        return interpolator 

    def __init__(self, potential_type, shape_params, potential_function_arrays = None):
        pot_archive = PotentialArchive()
        if potential_function_arrays == None:
            potential_function_arrays = [PotentialFunctionArray( pot_archive.getGeneralFile( shape_param, potential_type ) ) for shape_param in shape_params]

        if len(shape_params) != len(potential_function_arrays):
            print ('length of shape parameters array not same as length of potential function array! Truncating longest...') 
            if len (shape_params) > len(potential_function_arrays): shape_params = shape_params[0:len(potential_function_arrays)]
            if len (shape_params) < len(potential_function_arrays): potential_function_array = potential_function_array[0:len(shape_params)]

        self.potential_type = potential_type 
        self.shape_params = shape_params 
        self.Rlikeaxis = potential_function_arrays[0].Rlikeaxis
        self.Rlike_range = (np.min(self.Rlikeaxis), np.max(self.Rlikeaxis))
        self.zlikeaxis = potential_function_arrays[0].zlikeaxis
        self.zlike_range = (np.min(self.zlikeaxis), np.max(self.zlikeaxis))
        for pot_function in potential_function_arrays:
            if pot_function.Rlikeaxis.all() != self.Rlikeaxis.all(): print ('Not all R-like axes match!')
            if pot_function.zlikeaxis.all() != self.zlikeaxis.all(): print ('Not all z-like axes match!')

        self.interpolator = self.generatePotentialInterpolator(shape_params, potential_function_arrays) 
