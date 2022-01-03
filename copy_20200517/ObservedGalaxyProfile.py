#Here we define a python object that basically includes all of the information on an dSph observation

from DwarfGalDataArchive import DwarfGalDataArchive
from readInGalaxyData import readInGalaxyData
from DwarfGalDataArchive import DwarfGalDataArchive 

class ObservedGalaxyProfile:
    def __init__(self,population,pop_selection_method='none'):
        dataArchive = DwarfGalDataArchive()
        #self.dataFile = dataArchive.getFile(population)
        #self.dist = dataArchive.getDistance(population)
        self.dist = data_archive.getDistanceFromSun(population) #distance from sun  in pc {"carina":105000,"fornax":147000,"sculptor":86000,"sextans":86000}
        self.M = data_archive.getTotalMass(population) # total mass in M_sun {"carina":1.28*10**8,"fornax":1.28*10**8,"sculptor":1.28*10**8,"sextans":1.28*10**8} 
        self.gal_center_tuple = data_archive.getRADec(population) 
        (self.dist, self.M) = ExtraGalDataArchive.readInExtraGalData(population, pop_selection_method)       
        ( self.Target,        self.Date,       self.RA,      self.DEC,    self.Vmag,  self.VImag,
          self.Vhel,          self.VhelE,      self.SigFe,   self.SigFeE, self.SigMg, self.SigMgE,
          self.Pm,            self.RFe,        self.RMg,     self.RFeE,   self.RMgE,  self.V_Intensities, 
          self.SigMg_corr,    self.corrRa,     self.corrDec
          ) = readInGalaxyData(population,pop_selection_method)

        
        #print 'self.dist = distance from sun = ' + str(self.dist)
        
        
     