import numpy as np
from PotentialArchive import PotentialArchive 

def readInPotentialInterpolator(pot_type):
    pot_archive = PotentialArchive()
    pot_file = pot_archive.getPreGenFile(pot_type)
    print ('Loading file '+ pot_file ) 
    return np.load(pot_file, encoding = 'latin1', allow_pickle=True).item() 
