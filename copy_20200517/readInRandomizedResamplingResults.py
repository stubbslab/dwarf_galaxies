import numpy as np
import math
from DirectoryArchive import DirectoryArchive 
from ComputationalArchive import ComputationalArchive 

def readInRandomizedResamplingResults(file_name = 'best_fit_probabilities_number_bootstrapping.csv', n_dsph_configurations = 6):
    compute_params_archive = ComputationalArchive() 
    load_dir = compute_params_archive.getSpotProbabilityDir()

    loaded_randomized_vals = (np.genfromtxt(load_dir + file_name, delimiter=',').transpose()).tolist()

    return loaded_randomized_vals

    
    
