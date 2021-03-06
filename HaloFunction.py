from PotentialArchive import PotentialArchive
from PotentialFunctionArray import PotentialFunctionArray
import math
import numpy as np

class HaloFunction:

    def getHFunction(self):
        return self.HFunction

    def getCutScaling(self):
        return self.cut_scaling
    
    def getCutOffset(self):
        return self.cut_offset 

    def __init__(type, c, el, halo_interpolating_function = None):
        self.c = c
        self.el = el
        self.e = math.sqrt(1 - el**2) 
        type = type.lower()
        pot_archive = PotentialArchive()
        if halo_interpolating_function:
            self.HFunction = halo_interpolating_function
        else:
            if type == 'nfw':
                pot_file = pot_archive.getHaloFile(el)
            elif type == 'cored':
                pot_file = pot_archive.getCHaloFile(el) 
            self.HFunction = PotentialFunctionArray(pot_file)
        if type == 'nfw':
            self.cut_scaling = math.log(1.0+self.c) - self.c / (1.0+self.c)   # f(c)
            self.cut_offset = (-self.c)/(self.c+1.0) * (math.arcsin(self.e)) / self.e 
        elif type == 'cored':
            self.cut_scaling = math.log(1.0+self.c) - (3.0 * self.c**2 + 2.0 * self.c) / (2.0*(1+self.c)**2)    # h(c)
            self.cut_offset = (-2.0 * self.c - 1.0)/(2.0 * (self.c+1)**2) * (math.arcsin(self.e)) / self.e 
                
    
