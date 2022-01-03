import numpy as np
import math
from DwarfGalaxyParametersStorer import DwarfGalaxyParametersStorer 

class ArtificialParamsStorer:

    def getStorer(self, index):
        artificial_params = self.artificialParameters[index]
        return DwarfGalaxyParametersStorer(artificial_params[self.indeces['population']],
                                           artificial_params[self.indeces['withDisk']],
                                           el = artificial_params[self.indeces['el']],
                                           lam = artificial_params[self.indeces['lam']],
                                           rs = artificial_params[self.indeces['rs']],
                                           phi = artificial_params[self.indeces['phi']],
                                           theta = artificial_params[self.indeces['theta']],
                                           zeta = artificial_params[self.indeces['zeta']],
                                           eps = artificial_params[self.indeces['eps']],
                                           a = artificial_params[self.indeces['a']],
                                           b = artificial_params[self.indeces['b']]
                                           ) 


    def __init__(self):
        #population, withDisk, el, lam, rs, phi, theta, zeta, eps, a, b, WD
        self.indeces = {'population' : 0, 'withDisk' : 1, 'el' : 2, 'lam' : 3, 'rs' : 4, 'phi' : 5, 'theta' : 6, 'zeta' : 7, 'eps' : 8, 'a' : 9, 'b' : 10, 'disk_type' : 11, 'halo_type' : 12}
        self.artificialParameters = [[ ['fornax','MR'], 0, 0.2, 0.15, 500, math.pi * 0.0 , math.pi * 0.0 , 1.0, 0.0 , math.pi * 0.0, 0, 'sech_disk', 'nfw'],
                                     [ ['fornax','MR'], 0, 0.2, 0.15, 500, math.pi * 0.25, math.pi * 0.0 , 1.0, 0.0 , math.pi * 0.0, 0, 'sech_disk', 'cored'],
                                     [ ['fornax','MR'], 0, 0.2, 0.15, 500, math.pi * 0.0 , math.pi * 0.25, 1.0, 0.0 , math.pi * 0.0, 0, 'sech_disk', 'nfw'],
                                     [ ['fornax','MR'], 0, 0.2, 0.15, 250, math.pi * 0.25, math.pi * 0.0 , 1.0, 0.0 , math.pi * 0.0, 0, 'sech_disk', 'cored'],
                                     [ ['fornax','MR'], 1, 0.2, 0.1, 268, 0.7           , math.pi * 0.0 , 0.2, 0.0 , math.pi * 0.0, 0, 'sech_disk', 'nfw'],        
                                     [ ['fornax','MR'], 1, 0.2, 0.1, 268, 0.7           , math.pi * 0.0 , 0.2, 0.03, math.pi * 0.0, 0, 'sech_disk', 'cored'],
                                     [ ['fornax','MR'], 1, 0.2, 0.1, 268, 0.7           , math.pi * 0.0 , 0.2, 0.05, math.pi * 0.0, 0, 'sech_disk', 'nfw'],
                                     [ ['fornax','MR'], 1, 0.2, 0.1, 268, 0.7           , math.pi * 0.0 , 0.2, 0.1 , math.pi * 0.0, 0, 'sech_disk', 'cored'],
                                     [ ['fornax','MR'], 1, 0.2, 0.1, 268, 0.7           , math.pi * 0.0 , 0.2, 0.15, math.pi * 0.0, 0, 'sech_disk', 'nfw']
                                    ]
