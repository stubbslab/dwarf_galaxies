#Defines the DwarfGalaxyParametersStorer class.
#This class is basically a repository fo the parameters needed to do the dwarf galaxy simulations.
#We define it so that there is a standard way to pass parameters between classes without having to keep track of parameter ordering.
#This class ALSO keeps track of the parameters that we intend to vary in the MCMC simulation, as well as the parameters of the simulation.
#Which parameters are varied is based on whether we are doing the simulation with disk or without disk.  
#Defines all the parameters needed to do the dwarf galaxy simulations, in the following list:
# [dist, M, vphi, el, rs, phi, theta ,c, e, f, lam, zeta, eps, a, b]
#This is a file with hard coded parameters.   
import math 
from DwarfGalDataArchive import DwarfGalDataArchive

        
class DwarfGalaxyParametersStorer: 

    # Returns the array of all the parameters stored in this class, in the standard order
    def getParameterSet(self):
        return self.all_params_ordered

    #Returns the MCMC info, which consists of:
    # the indeces of the standard array to vary
    # the allowed parameter range, for eac of the parameters to vary
    # the width of the gaussian jump taken in each step
    def getMCMCInfo(self):
        param_indeces_to_vary = []
        param_ranges = []
        gauss_sampling_width = []
        if self.withDisk: 
            params_to_vary = ['rs','phi','eps','a']#set params to vary in case of disk + halo
            param_indeces_to_vary = [self.param_order_array.index(param) for param in params_to_vary] #set params to vary in case of disk + halo
            param_ranges = [[50,5000], [-1000*math.pi,1000*math.pi], [0.0,0.5],[-1000*math.pi,1000*math.pi]]
            gauss_sampling_width = [10, math.pi * 0.03, 0.005, math.pi * 0.1]
        else:
            params_to_vary = ['rs','phi']#set params to vary in case of disk + halo
            param_indeces_to_vary = [self.param_order_array.index(param) for param in params_to_vary] #set params to vary in case of just halo
            param_ranges = [[50,5000],[-1000*math.pi,1000*math.pi]]
            gauss_sampling_width = [10, math.pi * 0.03]
        
        return [param_indeces_to_vary, param_ranges, gauss_sampling_width]
        

    def __init__(self, population, withDisk, dist = 'none', M = 'none', vphi = 'none', sigsqr = 'none',
                  el=0.3, rs = 500, phi = math.pi*0.0, theta = 0.0, c = 5.0, 
                  lam=0.2, zeta = 0.2, eps = 0.05, a = math.pi*0.8, b = math.pi*0.5):
        
        data_archive = DwarfGalDataArchive() 
        
        self.population = population
        self.withDisk = withDisk
        #Parameters for which we have observational data 
        if dist == 'none':
            #print "dist == 'none' so reassigning to be: "
            self.dist  = data_archive.getDistanceFromSun(population) #distance from sun  in pc {"carina":105000,"fornax":147000,"sculptor":86000,"sextans":86000}
            #print self.dist
        else:
            #print "dist != 'none', so assigning dist to be " + str(dist)
            self.dist = dist 
        if M == 'none':
            self.M = data_archive.getTotalMass(population) # total mass in M_sun {"carina":1.28*10**8,"fornax":1.28*10**8,"sculptor":1.28*10**8,"sextans":1.28*10**8} 
        else:
            self.M = M 
        if vphi == 'none':
            self.vphi = data_archive.getRotationalVelocity(population) # rotational_velocity, in km/s {"carina":0.0,"fornax":0.0,"sculptor":0.0,"sextans":0.0} 
        else:
            self.vphi = vphi 
        if sigsqr == 'none':
            self.sigsqr = data_archive.getVelocityDispersion(population) 
        else:
            self.sigsqr = sigsqr
        
        #Parameters we must specify about the halo
        self.el    = el #ellipticity of halo
        self.rs    = rs  #scale radius in pc 
        self.phi   = phi # angle of halo axis of symmety off of on sky north
        self.theta = theta # rotation angle of halo rotated into sky 
        self.c     = c # ratio of rs to rvir.  Since we're plotting in r/rs and rmax=rvir, this is the distance past which we strongly suspect our model to be innappropriate 
        self.e     = math.sqrt(1-(1-el)**2) # the eccentricity of the halo (this is the e that appears in my expression for ...) 
        self.f = math.log(1+c) - c/(1+c) #function accounting for cutoff of potential at r/rs = c
        
        #Parameters we must specify for the disk 
        self.lam   = lam # scale height to radius ratio of disk
        self.zeta = zeta #ratio of disk and halo scale radii 
        self.a     = a # angle of disk axis with respect to halo axis 
        self.b     = b # angle of disk rated around central axis 
        if withDisk:
            self.eps   = eps #self interacting DM fraction          
        else:
            self.eps   = 0.0 #self interacting DM fraction 
        
        self.param_order_array = ['dist','M', 'vphi' , 'sigsqr' ,
                  'el', 'rs' , 'phi' , 'theta' , 'c', 
                  'lam', 'zeta', 'eps', 'a', 'b' ] #The order of the parameters in all_params_ordered array
        #Can get the index of a parameter in the all_params_ordered array using .index('name of parameter')
        self.all_params_ordered = [self.dist, self.M, self.vphi , self.sigsqr ,
                  self.el, self.rs , self.phi , self.theta , self.c, 
                  self.lam, self.zeta, self.eps, self.a, self.b ]
        
        
