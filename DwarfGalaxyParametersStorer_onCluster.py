#Defines the DwarfGalaxyParametersStorer class.
#This class is basically a repository fo the parameters needed to do the dwarf galaxy simulations.
#We define it so that there is a standard way to pass parameters between classes without having to keep track of parameter ordering.
#This class ALSO keeps track of the parameters that we intend to vary in the MCMC simulation, as well as the parameters of the simulation.
#Which parameters are varied is based on whether we are doing the simulation with disk or without disk.
#Defines all the parameters needed to do the dwarf galaxy simulations, in the following list:
# [dist, M, vphi, el, rs, phi, theta ,c, e, f, lam, zeta, eps, a, b]
#This is a file with hard coded parameters.
import math
import numpy as np
from DwarfGalDataArchive import DwarfGalDataArchive
import samplingFunctions
from samplingFunctions import constrainedGauss
from samplingFunctions import constrainedLogGauss
from samplingFunctions import varyUnitVector
from samplingFunctions import varyUpperHemisphereVector
from samplingFunctions import elSurrogateGauss
import random
from BackgroundMCMCInformationStorer import BackgroundStorer
from AstronomicalParameterArchive import AstronomicalParameterArchive
from CosmologicalParameterArchive import CosmologicalParameterArchive
import scipy.special as special

#Calculate virial radius, as described at https://en.wikipedia.org/wiki/Virial_mass#Virial_radius
# overdensity is Delta_c
def ComputeVirialRadius(mass_in_MMsun, scale_radius_pc, ellipticity, overdensity):
    astro_arch = AstronomicalParameterArchive()
    cosmo_arch = CosmologicalParameterArchive ()
    parsec_to_m = astro_arch.getParsecToM()
    solar_mass_in_kg = astro_arch.getSolarMass()
    gravitational_constant = astro_arch.getGravitationalConstant()
    hubble_constant = cosmo_arch.getH0('s')[0] #Get hubble constant, in seconds
    unit_scaling = (solar_mass_in_kg * gravitational_constant) / (hubble_constant ** 2.0 * (parsec_to_m) ** 3.0)
    #print 'unit_scaling = ' + str(unit_scaling)
    scaled_virial_rad = ((mass_in_MMsun) / (ellipticity * scale_radius_pc ** 3.0 * overdensity ) * (2.0 * unit_scaling)) ** (1.0 / 3.0) #scaled by rs
    return scaled_virial_rad


class DwarfGalaxyParametersStorer:

    def defineBetaIntFunctionForVelocityDispersion(self):
        gamma0, beta_exp, beta_R1, beta_R2 = [self.gamma0, self.beta_exp, self.beta_R1, self.beta_R2]
        beta0 = 2.0 * gamma0/ (1.0 + gamma0)
        unscaled_max_value = np.exp(-(beta_exp/2)) * (beta_R1/beta_R2 * np.sqrt(beta_exp / 2.0) ) ** (beta_exp)
        beta_scaling = beta0 / unscaled_max_value
        print ('beta_scaling = ' + str(beta_scaling))
        beta_funct = lambda Rs, beta_scaling = beta_scaling, beta_exp = beta_exp, beta_R1 = beta_R1, beta_R2 = beta_R2: beta_scaling * 0.5 * (beta_R1 / beta_R2) ** beta_exp * special.gammaincc(beta_exp/2.0, (Rs/beta_R1) ** 2.0) / special.gamma(beta_exp/2.0)
        return beta_funct

    # Returns the array of all the parameters stored in this class, in the standard order
    def getParameterSet(self):
        return self.all_params_ordered

    #Returns the MCMC info, which consists of:
    # the indeces of the standard array to vary
    # the allowed parameter range, for eac of the parameters to vary
    # the width of the gaussian jump taken in each step
    def getMCMCInfo(self):
        param_indeces_to_vary = []
        gauss_sampling_width = []
        if self.withDisk:
            params_to_vary = ['M', 'rs', 'halo_sym_axis', 'halo_center', 'sigsqr_RR', 'omega_phi', 'el', 'Rd', 'disk_sym_axis', 'disk_center', 'eps',  'lam'] #set params to vary in case of disk + halo
            #params_to_vary = ['rs', 'halo_sym_axis', 'halo_center', 'el', 'zeta', 'disk_sym_axis', 'disk_center', 'eps',  'lam'] #set params to vary in case of disk + halo
            param_indeces_to_vary = [self.param_order_array.index(param) for param in params_to_vary] #set params to vary in case of disk + halo
            #param_ranges = [[0.0,1000.0 * 10**7], [200.0,5000.0], [-1000.0 * math.pi, 1000.0 * math.pi], [-1000.0 * math.pi, 1000.0 * math.pi], [0.1, 10.0],
            #                [0.00001,1.0],        [0.0,1.0],      [-1000.0 * math.pi, 1000.0 * math.pi], [-1000.0 * math.pi, 1000.0 * math.pi], [0.1, 1.0]   ]

            param_functions = [lambda scaler, jump_multiplier: constrainedGauss (scaler, self.param_step_sizes['M'] * jump_multiplier, self.param_ranges['M'] ), #M
                               lambda scaler, jump_multiplier: constrainedGauss (scaler, self.param_step_sizes['rs'] * jump_multiplier, self.param_ranges['rs']  ), #rs
                               lambda array,  jump_multiplier: varyUnitVector (array, abs(random.gauss(0.0, 0.05 * jump_multiplier))) if self.halo_edge_on == 0
                                                               else varyUnitVector (array, abs(random.gauss(0.0, 0.05 * jump_multiplier)), fixed_values = {1: 0.0}), # halo_sym_axis
                               lambda vector, jump_multiplier: [constrainedGauss(vector[0], self.param_step_sizes['h_x_center'] * jump_multiplier, self.param_ranges['h_x_center'] ),
                                                                0.0,
                                                                constrainedGauss(vector[2], self.param_step_sizes['h_z_center'] * jump_multiplier, self.param_ranges['h_z_center'])], #halo_center
                               lambda scaler, jump_multiplier: constrainedGauss(scaler, self.param_step_sizes['sigsqr_RR'] * jump_multiplier, self.param_ranges['sigsqr_RR'] ), #sigsqr_RR
                               lambda scaler, jump_multiplier: constrainedGauss(scaler, self.param_step_sizes['omega_phi'] * jump_multiplier, self.param_ranges['omega_phi'] ), #omega_phi
                               lambda scaler, jump_multiplier: elSurrogateGauss(scaler, self.param_step_sizes['el'] * jump_multiplier, self.param_ranges['el'] ), #el
                               lambda scaler, jump_multiplier: constrainedGauss(scaler, self.param_step_sizes['Rd'] * jump_multiplier, self.param_ranges['Rd'] ), #zeta
                               lambda array,  jump_multiplier: varyUnitVector (array, abs(random.gauss(0.0, 0.05 * jump_multiplier))) if self.disk_edge_on == 0
                                                               else varyUnitVector (array, abs(random.gauss(0.0, 0.05 * jump_multiplier)), fixed_values = {1: 0.0}), # disk_sym_axis
                               lambda vector, jump_multiplier: [constrainedGauss(vector[0], self.param_step_sizes['d_x_center'] * jump_multiplier, self.param_ranges['d_x_center']),
                                                                0.0,
                                                                constrainedGauss(vector[2], self.param_step_sizes['d_z_center'] * jump_multiplier, self.param_ranges['d_z_center'])], #disk_center
                               lambda scaler, jump_multiplier: constrainedGauss(scaler, self.param_step_sizes['eps'] * jump_multiplier, self.param_ranges['eps'] ), #eps
                               lambda scaler, jump_multiplier: constrainedGauss(scaler, self.param_step_sizes['lam'] * jump_multiplier, self.param_ranges['lam'] ) #lam
                               ]
            #param_ranges = [[200.0,5000.0], [- float('inf'), float('inf')], [0.1, 10.0],
            #                [0.00001,1.0],  [- float('inf'), float('inf')], [0.0,1.0],  [0.1, 1.0]   ]
            #gauss_sampling_width = [10.0, math.pi * 0.03, math.pi * 0.03, 0.01,
            #                        0.01,  0.01, math.pi * 0.03, math.pi * 0.03, 0.01   ]
        else:
            params_to_vary = ['M', 'rs', 'halo_sym_axis', 'halo_center', 'sigsqr_RR', 'omega_phi', 'el']#set params to vary in case of just halo
            #params_to_vary = ['rs', 'halo_sym_axis', 'halo_center', 'el']#set params to vary in case of just halo
            param_indeces_to_vary = [self.param_order_array.index(param) for param in params_to_vary] #set params to vary in case of just halo
            #param_ranges = [[0.0,1000 * 10**7], [200,5000], [-1000*math.pi,1000*math.pi], [-1000*math.pi,1000*math.pi], [0.1,10.0]]
            param_functions = [lambda scaler, jump_multiplier: constrainedGauss (scaler, self.param_step_sizes['M'] * jump_multiplier, self.param_ranges['M'] ), #M
                               lambda scaler, jump_multiplier: constrainedGauss (scaler, self.param_step_sizes['rs'] * jump_multiplier, self.param_ranges['rs']  ), #rs
                               lambda array,  jump_multiplier: varyUnitVector (array, abs(random.gauss(0.0, 0.05 * jump_multiplier))) if self.halo_edge_on == 0
                                                               else varyUnitVector (array, abs(random.gauss(0.0, 0.05 * jump_multiplier)), fixed_values = {1: 0.0}), # halo_sym_axis
                               lambda vector, jump_multiplier: [constrainedGauss(vector[0], self.param_step_sizes['h_x_center'] * jump_multiplier, self.param_ranges['h_x_center'] ),
                                                                0.0,
                                                                constrainedGauss(vector[2], self.param_step_sizes['h_z_center'] * jump_multiplier, self.param_ranges['h_z_center'])], #halo_center
                               lambda scaler, jump_multiplier: constrainedGauss(scaler, self.param_step_sizes['sigsqr_RR'] * jump_multiplier, self.param_ranges['sigsqr_RR'] ), #sigsqr_RR
                               lambda scaler, jump_multiplier: constrainedGauss(scaler, self.param_step_sizes['omega_phi'] * jump_multiplier, self.param_ranges['omega_phi'] ), #omega_phi
                               lambda scaler, jump_multiplier: elSurrogateGauss(scaler, self.param_step_sizes['el'] * jump_multiplier, self.param_ranges['el'] ), #el
                               ]
            #param_ranges = [[200,5000], [-1000*math.pi,1000*math.pi], [-1000*math.pi,1000*math.pi], [0.1,10.0]]
            #gauss_sampling_width = [10**6, 10, math.pi * 0.03, math.pi*0.03, 0.01]
            #gauss_sampling_width = [10, math.pi * 0.03, math.pi*0.03, 0.01]

        return [param_indeces_to_vary, param_functions]

    def printContents(self):
        print ('Printing full contents of storer. ')
        print ('self.population = ' + str(self.population))
        print ('self.withDisk = ' + str(self.withDisk))
        print ('self.dist = ' + str(self.dist))
        print ('self.M = ' + str(self.M))
        print ('self.omega_phi = ' + str(self.omega_phi))
        print ('self.sigsqr_RR = ' + str(self.sigsqr_RR))
        print ('self.dispersion_b = ' + str(self.dispersion_b))
        print ('self.gamma0 = ' + str(self.gamma0))
        print ('self.beta_exp = ' + str(self.beta_exp))
        print ('self.beta_R1 = ' + str(self.beta_R1))
        print ('self.beta_R2 = ' + str(self.beta_R2))
        #print ('self.gamma_turn2 = ' + str(self.gamma_turn2))
        print ('self.el = ' + str(self.el))
        print ('self.rs = ' + str(self.rs))
        print ('self.phi = ' + str(self.phi))
        print ('self.halo_sym_axis = [' + ','.join([str(axis_comp) for axis_comp in self.halo_sym_axis]) + ']')
        print ('self.halo_center = [' + ','.join([str(center_comp) for center_comp in self.halo_center]) + ']')
        print ('self.theta = ' + str(self.theta))
        print ('self.overdensity = ' + str(self.overdensity) )
        print ('self.c = ' + str(self.c))
        print ('self.lam = ' + str(self.lam))
        print ('self.a = ' + str(self.a))
        print ('self.b = ' + str(self.b))
        print ('self.disk_sym_axis = [' + ','.join([str(axis_comp) for axis_comp in self.disk_sym_axis]) + ']')
        print ('self.disk_center = [' + ','.join([str(center_comp) for center_comp in self.disk_center]) + ']')
        print ('self.eps = ' + str(self.eps))
        print ('self.Rd = ' + str(self.Rd))
        print ('self.halo_edge_on = ' + str(self.halo_edge_on))
        print ('self.disk_edge_on = ' + str(self.disk_edge_on))

    def recomputeCutoffRadius(self):
        self.c = ComputeVirialRadius (self.M, self.rs, self.el, self.overdensity)  #CALCULATE c from M, rs, Q, Delta

    #We assume beta = beta0 + (betainf - beta0) * 1.0 / (1.0 + (R / beta_R0) ** beta_exp)
    #beta = 2.0 * gamma / (1.0 + gamma)
    #Give your answers in terms of gamma0 and gammainf (they can be between -1 and 1)
    #Now a new function: beta = beta0 * exp(-(R/beta_R1) ** 2) / ((R/beta_R2) ** beta_exp) (has beta-> 0 at R = 0 and R = inf; otherwise, things seem to diverge)
    # You now give only gamma0 between -1 and 1

    def __init__(self, population, withDisk, dist = 'none', M = 'none', omega_phi = 0.0, sigsqr_RR = 100.0,
                  dispersion_b = 1.0, gamma0 = 0.0, beta_exp = 1.0, beta_R1 = 1.0, beta_R2 = 1.0,
                  el=0.3, rs = 500.0, phi = 'none', theta = 'none', halo_sym_axis = [0.0,0.0,1.0], halo_center = [0.0, 0.0, 0.0], c = 5.0,
                  lam=0.15, Rd = 500.0, eps = 0.05, a = 'none', b = 'none', disk_sym_axis = [0.0,0.0,1.0], disk_center = [0.0, 0.0, 0.0],
                  fixed_positive_index = 0, halo_edge_on = 0, disk_edge_on = 0, fixed_disk_angle = None, fixed_halo_angle = None, specified_param_ranges = {}, stellar_data = None,
                  overdensity = 200.0):

        BackgroundInfo = BackgroundStorer()
        data_archive = DwarfGalDataArchive()

        self.overdensity = overdensity
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
        self.dispersion_b = dispersion_b
        self.gamma0 = gamma0
        self.beta_exp = beta_exp
        self.beta_R1 = beta_R1
        self.beta_R2 = beta_R2
        #self.beta_int_dispersion_funct = self.defineBetaIntFunctionForVelocityDispersion()

        #Parameters we must specify about the halo
        self.el    = el #ellipticity of halo
        self.rs    = rs  #scale radius in pc
        #self.halo_center = [ component / self.rs for component in halo_center ]
        self.halo_center = halo_center
        #phi is angle of halo axis of symmetry off of sky north. Theta is rotation angle along halo axis of symmetry
        # halo_sym_axis is unit vector defining orientatin of disk symmetry axis in sky coordinates.
        self.halo_edge_on = halo_edge_on
        if (halo_sym_axis is 'none' or halo_sym_axis is None):
            if halo_edge_on is 0:
                self.theta = 0.0
            else:
                self.theta = theta
            if fixed_halo_angle is None:
                self.phi = phi
            else:
                self.phi = fixed_halo_angle
            self.halo_sym_axis = [math.cos(self.theta) * math.sin(self.phi), math.sin(self.theta) * math.sin(self.phi), math.cos(self.phi)]
        else:
            if halo_edge_on is 1: halo_sym_axis = [halo_sym_axis[0], 0.0, halo_sym_axis[2]]
            if not(fixed_halo_angle is None):
                halo_x_z_mag = np.sqrt(halo_sym_axis[0] ** 2.0 + halo_sym_axis[2] ** 2.0)
                halo_sym_axis[0] = math.sin(fixed_halo_angle) * halo_x_z_mag
                halo_sym_axis[2] = math.cos(fixed_halo_angle) * halo_x_z_mag
            if (fixed_positive_index != -1 and halo_sym_axis[fixed_positive_index] < 0.0):
                halo_sym_axis = [component * -1.0 for component in halo_sym_axis]
            halo_sym_axis = (np.array(halo_sym_axis) / math.sqrt(sum([elem ** 2.0 for elem in halo_sym_axis]))).tolist()
            self.halo_sym_axis = halo_sym_axis
            #print 'self.halo_sym_axis = ' + str(self.halo_sym_axis)
            if (halo_sym_axis[0] > 0.0):
                #print 'here 1'
                self.theta = math.atan( halo_sym_axis[1] / halo_sym_axis[0] )
            elif halo_sym_axis[0] == 0.0:
                if halo_sym_axis[1] > 0.0:
                    #print 'here 2'
                    self.theta = math.pi / 2.0
                else:
                    #print 'here 3'
                    self.theta = - math.pi / 2.0
            else:
                if halo_sym_axis[1] < 0.0:
                    #print 'here 4'
                    self.theta = -math.pi  + math.atan( halo_sym_axis[1] / halo_sym_axis[0] )
                else:
                    #print 'here 5'
                    self.theta = math.pi + math.atan( halo_sym_axis[1] / halo_sym_axis[0] )
            self.phi = math.acos(halo_sym_axis[2] / math.sqrt(sum([component ** 2.0 for component in halo_sym_axis])) )
        #print 'self.theta = ' + str(self.theta)
        #print 'self.phi = ' + str(self.phi)

        self.c = c # ratio of rs to rvir.  Since we're plotting in r/rs and rmax=rvir, this is the distance past which we strongly suspect our model to be innappropriate
        self.c = ComputeVirialRadius (self.M, self.rs, self.el, self.overdensity)  #CALCULATE c from M, rs, Q, Delta
        self.f = math.log(1+self.c) - self.c/(1+self.c) #function accounting for cutoff of potential at r/rs = c

        #Parameters we must specify for the disk
        self.lam   = lam # scale height to radius ratio of disk
        #self.disk_center = [ component / self.rs for component in disk_center ]
        self.disk_center = disk_center
        # a and b are defined for the disk as phi and theta are defined for the halo
        # disk_sym_axis is unit vector defining orientation of disk symmetry axis in sky coordinates.
        self.disk_edge_on = disk_edge_on
        self.fixed_disk_angle =  fixed_disk_angle
        if (disk_sym_axis is 'none' or disk_sym_axis is None):
            if disk_edge_on is 0:
                self.b = 0.0
            else:
                self.b = b
            if fixed_disk_angle is None:
                self.a = a
            else:
                self.a = fixed_disk_angle
            self.disk_sym_axis = [math.cos(self.b) * math.sin(self.a), math.sin(self.b) * math.sin(self.a), math.cos(self.a)]
        else:
            if (fixed_positive_index != -1 and disk_sym_axis[fixed_positive_index] < 0.0):
                disk_sym_axis = [component * -1.0 for component in disk_sym_axis]
            if disk_edge_on is 1: disk_sym_axis = [disk_sym_axis[0], 0.0, disk_sym_axis[2]]
            self.disk_sym_axis = disk_sym_axis
            if (disk_sym_axis[0] > 0.0):
                self.b = math.atan( disk_sym_axis[1] / disk_sym_axis[0] )
            elif disk_sym_axis[0] == 0.0:
                if disk_sym_axis[1] > 0.0: self.b = math.pi / 2.0
                else: self.b = - math.pi / 2.0
            else:
                if disk_sym_axis[1] < 0.0:
                    self.b = -math.pi  + math.atan( disk_sym_axis[1] / disk_sym_axis[0] )
                else:
                    self.b = math.pi + math.atan( disk_sym_axis[1] / disk_sym_axis[0] )
            if fixed_disk_angle is None:
                self.a = math.acos(disk_sym_axis[2] / math.sqrt(sum([component ** 2.0 for component in disk_sym_axis])) )
            else:
                self.a = fixed_disk_angle

        if withDisk:
            self.eps   = eps #self interacting DM fraction
            self.Rd = Rd
        else:
            self.eps   = 0.0 #self interacting DM fraction
            self.Rd = self.rs
        self.zeta = self.Rd/ self.rs#ratio of disk and halo scale radii

        self.omega_phi = omega_phi

        if sigsqr_RR == 'none':
            if stellar_data is None:
                self.sigsqr_RR = data_archive.getVelocityDispersion(population)
            else:
                self.sigsqr_RR = stellar_data.measureDispersionSignal( population, self.dist, self.halo_center, self.halo_sym_axis, self.dispersion_b,
                                                                       self.gamma0, self.beta_exp, self.beta_R1, self.beta_R2 )
        else:
            self.sigsqr_RR = sigsqr_RR

        self.param_order_array = ['dist','M', 'omega_phi' , 'sigsqr_RR' ,
                                  'dispersion_b', 'gamma0', 'beta_exp', 'beta_R1','beta_R2',
                                  'el', 'rs' , 'phi' , 'theta' , 'halo_sym_axis', 'halo_center', 'c',
                                  'lam', 'Rd', 'eps', 'a', 'b', 'disk_sym_axis', 'disk_center' ] #The order of the parameters in all_params_ordered array
        #Can get the index of a parameter in the all_params_ordered array using .index('name of parameter')
        self.all_params_ordered = [self.dist, self.M, self.omega_phi , self.sigsqr_RR,
                                   self.dispersion_b, self.gamma0, self.beta_exp, self.beta_R1,self.beta_R2,
                                   self.el, self.rs , self.phi , self.theta , self.halo_sym_axis, self.halo_center, self.c,
                                   self.lam, self.Rd, self.eps, self.a, self.b, self.disk_sym_axis, self.disk_center ]
        #The range over which the parameters are allowed to vary
        self.param_ranges = BackgroundInfo.param_ranges
        self.param_step_sizes = BackgroundInfo.param_step_sizes
        #User allowed to override.  Note a range of [a, a] means that the parameter will always have the value a.
        for param in specified_param_ranges.keys():
            self.param_ranges[param] = specified_param_ranges[param]
