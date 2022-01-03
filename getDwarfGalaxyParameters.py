#Defines all the parameters need to do the dwarf galaxy simulations, in the following list:
# [dist, M, vphi, el, rs, phi, theta ,c, e, f, lam, zeta, eps, a, b]
import math 
from DwarfGalDataArchive import DwarfGalDataArchive



def getAllParameters(self):
    return self.all_params
    
    
#I need to choose the parameters we want to vary, depending on whether we are working with the disk + halo or just disk models
def getMCMCInfo(self):
    param_indeces_to_vary = []
    param_ranges = []
    gauss_sampling_width = []
    if self.withDisk: 
        param_indeces_to_vary = [4,5,9,10,11] #vary rs, phi, zeta, eps, a
        param_ranges = [[100,1000], [-math.pi,math.pi], [0.0,0.5],,[0.0,0.5],[-math.pi,math.pi]]
        gauss_sampling_width = [100, math.pi * 0.2, 0.1,0.1, math.pi * 0.2, 0.0]
    else:
        param_indeces_to_vary = [4,5] #vary rs, phi 
        param_ranges = [[100,1000],[-math.pi,math.pi]]
        gauss_sampling_width = [100, math.pi * 0.2]
    
    return [param_indeces_to_vary, param_ranges, gauss_sampling_width]
        
        
Class DwarfGalaxyParameters: 

    def __init__(self, population, withDisk=0, dist = 'none', M = 'none', vphi = 0.0,
                  el=0.3, rs = 500, phi = math.pi*0.0, theta = 0.0, c = 5.0,
                  lam=0.2, zeta = 0.2, eps = 0.1, a = math.pi*0.1, b = math.pi*0.5):
        
        data_archive = DwarfGalDataArchive() 
        
        self.population = population
        self.withDisk = withDisk
        #Parameters for which we have observational data 
        if dist == 'none':
            self.dist  = data_archive.getDistanceFromSun(population) #distance from sun  in pc {"carina":105000,"fornax":147000,"sculptor":86000,"sextans":86000}
        else:
            self.dist = dist 
        if M == 'none'
            self.M = data_archive.getTotalMass(population) # total mass in M_sun {"carina":1.28*10**8,"fornax":1.28*10**8,"sculptor":1.28*10**8,"sextans":1.28*10**8} 
        else:
            self.M = M 
        if vphi == 'none'
            self.vphi = data_archive.getRotationalVelocity(population) # rotational_velocity, in km/s {"carina":0.0,"fornax":0.0,"sculptor":0.0,"sextans":0.0} 
        else:
            self.vphi = vphi 
          
        #Parameters we must specify about the halo
        self.el    = el #ellipticity of halo
        self.rs    = rs  #scale radius in pc 
        self.phi   = phi # angle of halo axis of symmety off of on sky north
        self.theta = theta # rotation angle of halo rotated into sky 
        self.c     = c # ratio of rs to rvir.  Since we're plotting in r/rs and rmax=rvir, this is the distance past which we strongly suspect our model to be innappropriate 
        self.e     = math.sqrt(1-(1-el)**2) # the eccentricity of the halo (this is the e that appears in my expression for ...) 
        self.f     = math.log(1+c) - c/(1+c) #function accounting for cutoff of potential at r/rs = c
       
        #Parameters we must specify for the disk 
        self.lam   = lam # scale height to radius ratio of disk
        self.zeta  = zeta #ratio of disk and halo scale radii 
        self.eps   = eps #self interacting DM fraction 
        self.a     = a # angle of disk axis with respect to halo axis 
        self.b     = b # angle of disk rated around central axis 
        
    if withDisk:
        self.all_params = [self.dist,self.M, self.vphi,
                             self.el, self.rs, self.phi, self.theta ,self.c,
                             self.lam, self.zeta, self.eps, self.a, self.b]
    else:
        self.all_params = [self.dist,self.M, self.vphi,
                             self.el, self.rs, self.phi, self.theta ,self.c,
                             0.0, 0.0, 0.0, 0.0, 0.0]
    return parameter_set 
def getParamatersToVary(withDisk):
    if withDisk:
        return [2,3,4,5,7]
    else:
        return [4,5]