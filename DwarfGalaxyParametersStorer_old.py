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
            #params_to_vary = ['M', 'rs', 'phi', 'theta', 'el', 'zeta', 'eps', 'a', 'b', 'lam']#set params to vary in case of disk + halo
            params_to_vary = ['rs', 'h_axis_dist', 'el', 'zeta', 'd_axis_dist', 'eps',  'lam']#set params to vary in case of disk + halo
            param_indeces_to_vary = [self.param_order_array.index(param) for param in params_to_vary] #set params to vary in case of disk + halo
            #param_ranges = [[0.0,1000.0 * 10**7], [200.0,5000.0], [-1000.0 * math.pi, 1000.0 * math.pi], [-1000.0 * math.pi, 1000.0 * math.pi], [0.1, 10.0],
            #                [0.00001,1.0],        [0.0,1.0],      [-1000.0 * math.pi, 1000.0 * math.pi], [-1000.0 * math.pi, 1000.0 * math.pi], [0.1, 1.0]   ]
            param_ranges = [[200.0,5000.0], [- float('inf'), float('inf')], [0.1, 10.0],
                            [0.00001,1.0],  [- float('inf'), float('inf')], [0.0,1.0],  [0.1, 1.0]   ]
            gauss_sampling_width = [10.0, math.pi * 0.03, math.pi * 0.03, 0.01,
                                    0.01,  0.01, math.pi * 0.03, math.pi * 0.03, 0.01   ]
        else:
            #params_to_vary = ['M', 'rs', 'phi', 'theta', 'el']#set params to vary in case of just halo
            params_to_vary = ['rs', 'phi', 'theta', 'el']#set params to vary in case of just halo
            param_indeces_to_vary = [self.param_order_array.index(param) for param in params_to_vary] #set params to vary in case of just halo
            #param_ranges = [[0.0,1000 * 10**7], [200,5000], [-1000*math.pi,1000*math.pi], [-1000*math.pi,1000*math.pi], [0.1,10.0]]
            param_ranges = [[200,5000], [-1000*math.pi,1000*math.pi], [-1000*math.pi,1000*math.pi], [0.1,10.0]]
            #gauss_sampling_width = [10**6, 10, math.pi * 0.03, math.pi*0.03, 0.01]
            gauss_sampling_width = [10, math.pi * 0.03, math.pi*0.03, 0.01]
        
        return [param_indeces_to_vary, param_ranges, gauss_sampling_width]

    def printContents(self):
        print 'Printing full contents of storer. '
        print 'self.population = '
        print self.population
        print 'self.withDisk = ' + str(self.withDisk)
        print 'self.dist = ' + str(self.dist)
        print 'self.M = ' + str(self.M)
        print 'self.vphi = ' + str(self.vphi)
        print 'self.sigsqr = ' + str(self.sigsqr)
        print 'self.el = ' + str(self.el)
        print 'self.rs = ' + str(self.rs)
        print 'self.phi = ' + str(self.phi)
        print 'self.halo_sym_axis = [' + ','.join(self.halo_sym_axis) + ']'  
        print 'self.theta = ' + str(self.theta)
        print 'self.c = ' + str(self.c)
        print 'self.lam = ' + str(self.lam)
        print 'self.a = ' + str(self.a)
        print 'self.b = ' + str(self.b)
        print 'self.disk_sym_axis = [' + ','.join(self.disk_sym_axis) + ']'  
        print 'self.eps = ' + str(self.eps)
        print 'self.zeta = ' + str(self.zeta) 

    def __init__(self, population, withDisk, dist = 'none', M = 'none', vphi = 'none', sigsqr = 'none',
                  el=0.3, rs = 500, phi = 'none', theta = 'none', halo_sym_axis = [0.0,0.0,1.0], c = 5.0, 
                  lam=0.2, zeta = 1.0, eps = 0.05, a = 'none', b = 'none', disk_sym_axis = [0.0,0.0,1.0],
                  fixed_disk_angle = 0):
        
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
        #phi is angle of halo axis of symmetry off of sky north. Theta is rotation angle along halo axis of symmetry
        # halo_sym_axis is unit vector defining orientatin of disk symmetry axis in sky coordinates. 
        if (phi == 'none' and theta == 'none'): 
            self.theta = math.arctan( halo_sym_axis[1] / halo_sym_axis[0] ) 
            self.phi = math.arccos(halo_sym_axis[2] / math.sqrt(sum([component ** 2.0 for component in halo_sym_axis])) 
        else:
            self.phi = phi
            self.theta = theta
            halo_sym_axis = [math.cos(self.theta) * math.sin(self.phi), math.sin(self.theta) * math.sin(self.phi), math.cos(self.phi)]
            
        self.c     = c # ratio of rs to rvir.  Since we're plotting in r/rs and rmax=rvir, this is the distance past which we strongly suspect our model to be innappropriate 
        self.f = math.log(1+c) - c/(1+c) #function accounting for cutoff of potential at r/rs = c
        
        #Parameters we must specify for the disk
        self.lam   = lam # scale height to radius ratio of disk
        # a and b are defined for the disk as phi and theta are defined for the halo 
        # disk_sym_axis is unit vector defining orientation of disk symmetry axis in sky coordinates. 
        if (a == 'none' and b == 'none'): 
            self.b = math.arctan(disk_sym_axis[1] / disk_sym_axis[0] ) 
            self.a = math.arccos(disk_sym_axis[2] / math.sqrt(sum([component ** 2.0 for component in disk_sym_axis])) ) 
        else:
            self.a = a
            self.b = b
            disk_sym_axis = [math.cos(self.b) * math.sin(self.a), math.sin(self.b) * math.sin(self.a), math.cos(self.a)] 
        if withDisk:
            self.eps   = eps #self interacting DM fraction
            self.zeta = zeta #ratio of disk and halo scale radii
        else:
            self.eps   = 0.0 #self interacting DM fraction
            self.zeta = 1.0
        
        self.param_order_array = ['dist','M', 'vphi' , 'sigsqr' ,
                  'el', 'rs' , 'phi' , 'theta' , 'c', 
                  'lam', 'zeta', 'eps', 'a', 'b' ] #The order of the parameters in all_params_ordered array
        #Can get the index of a parameter in the all_params_ordered array using .index('name of parameter')
        self.all_params_ordered = [self.dist, self.M, self.vphi , self.sigsqr ,
                  self.el, self.rs , self.phi , self.theta , self.c, 
                  self.lam, self.zeta, self.eps, self.a, self.b ]
        
        
