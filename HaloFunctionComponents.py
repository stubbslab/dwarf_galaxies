from PotentialArchive import PotentialArchive
from PotentialFunctionArray import PotentialFunctionArray
import math
import numpy as np

class HaloFunctionComponents:

    def getHFunction(self):
        return self.HFunction

    def getOverallScaling(self):
        return self.overall_scale_factor 
    
    def getCutOffset(self):
        return self.cut_offset 

    def __init__(self,halo_type, c, el, halo_interpolating_function = None):
        self.c = c
        self.el = el
        halo_type = halo_type.lower()
        pot_archive = PotentialArchive()
        if halo_interpolating_function:
            self.HFunction = halo_interpolating_function
        else:
            if halo_type == 'nfw':
                pot_file = pot_archive.getHaloFile(e)
            elif halo_type == 'cored':
                pot_file = pot_archive.getCHaloFile(e)
            elif halo_type == 'burkert':
                pot_file = pot_archive.getCHaloFile(e) 
            self.HFunction = PotentialFunctionArray(pot_file)
        if halo_type == 'nfw':
            self.cut_scale_factor = math.log(1.0+self.c) - self.c / (1.0+self.c)   # f(c)
            self.overall_scale_factor = 1.0 / (2.0 * self.cut_scale_factor) 
            #self.cut_offset = (-1.0)/(self.c+1.0) * (math.asin(self.e)) / self.e
            self.cut_offset = 0.0
        elif halo_type == 'cored':
            self.cut_scale_factor = math.log(1.0+self.c) - (3.0 * self.c**2 + 2.0 * self.c) / (2.0*(1.0+self.c)**2.0)    # h(c)
            self.overall_scale_factor = 1.0 / (4.0 * self.cut_scale_factor) 
            #self.cut_offset = (-2.0 * self.c - 1.0)/(2.0 * (self.c+1)**2) * (math.asin(self.e)) / self.e
            self.cut_offset = 0.0
        elif halo_type == 'burkert':
            self.cut_scale_factor = 1.0 / 4.0 * ( 2.0 * math.log(1.0 + c) + math.log(1.0 + c ** 2.0) - 2.0 * math.atan(c) )
            self.overall_scale_factor= 1.0 / (8.0 * self.cut_scale_factor) 
            self.cut_offset = 0.0
        else:
            self.overall_scale_factor = 1.0 
                
    
