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
from SphericalBackgroundMCMCInformationStorer import BackgroundStorer
from AstronomicalParameterArchive import AstronomicalParameterArchive
from CosmologicalParameterArchive import CosmologicalParameterArchive
import scipy.special as special

#Calculate virial radius, as described at https://en.wikipedia.org/wiki/Virial_mass#Virial_radius
# overdensity is Delta_c
def ComputeVirialRadius(mass_in_MMsun, scale_radius_pc, overdensity):
    astro_arch = AstronomicalParameterArchive()
    cosmo_arch = CosmologicalParameterArchive ()
    parsec_to_m = astro_arch.getParsecToM()
    solar_mass_in_kg = astro_arch.getSolarMass()
    gravitational_constant = astro_arch.getGravitationalConstant()
    hubble_constant = cosmo_arch.getH0('s')[0] #Get hubble constant, in seconds
    unit_scaling = (solar_mass_in_kg * gravitational_constant) / (hubble_constant ** 2.0 * (parsec_to_m) ** 3.0)
    #print 'unit_scaling = ' + str(unit_scaling)
    scaled_virial_rad = ((mass_in_MMsun) / ( scale_radius_pc ** 3.0 * overdensity ) * (2.0 * unit_scaling)) ** (1.0 / 3.0) #scaled by rs
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
        params_to_vary = ['M', 'rs', 'halo_center',
                          'sigsqr_rr_0', 'sigsqr_rr_inf_to_0_rat', 'r_sigsqr_rr0', 'alpha_sigsqr_rr', 'gamma_for_beta_inf', 'r_beta0', 'alpha_beta' ] #set params to vary in case of disk + halo
        #params_to_vary = ['rs', 'halo_sym_axis', 'halo_center', 'el', 'zeta', 'disk_sym_axis', 'disk_center', 'eps',  'lam'] #set params to vary in case of disk + halo
        param_indeces_to_vary = [self.param_order_array.index(param) for param in params_to_vary] #set params to vary in case of disk + halo
        #param_ranges = [[0.0,1000.0 * 10**7], [200.0,5000.0], [-1000.0 * math.pi, 1000.0 * math.pi], [-1000.0 * math.pi, 1000.0 * math.pi], [0.1, 10.0],
        #                [0.00001,1.0],        [0.0,1.0],      [-1000.0 * math.pi, 1000.0 * math.pi], [-1000.0 * math.pi, 1000.0 * math.pi], [0.1, 1.0]   ]

        param_functions = [lambda scaler, jump_multiplier: constrainedGauss (scaler, self.param_step_sizes['M'] * jump_multiplier, self.param_ranges['M'] ), #M
                              lambda scaler, jump_multiplier: constrainedGauss (scaler, self.param_step_sizes['rs'] * jump_multiplier, self.param_ranges['rs']  ), #rs
                              lambda vector, jump_multiplier: [constrainedGauss(vector[0], self.param_step_sizes['h_x_center'] * jump_multiplier, self.param_ranges['h_x_center'] ),
                                                                0.0,
                                                                constrainedGauss(vector[2], self.param_step_sizes['h_z_center'] * jump_multiplier, self.param_ranges['h_z_center'])], #halo_center
                               lambda scaler, jump_multiplier: constrainedGauss(scaler, self.param_step_sizes['sigsqr_rr_0'] * jump_multiplier, self.param_ranges['sigsqr_rr_0'] ), #sigsqr_rr_0_MR
                               lambda scaler, jump_multiplier: constrainedGauss(scaler, self.param_step_sizes['sigsqr_rr_inf_to_0_rat'] * jump_multiplier, self.param_ranges['sigsqr_rr_inf_to_0_rat'] ), #sigsqr_rr_inf_MR
                               lambda scaler, jump_multiplier: constrainedGauss(scaler, self.param_step_sizes['r_sigsqr_rr0'] * jump_multiplier, self.param_ranges['r_sigsqr_rr0'] ), #r_sigsqr_rr0_MR
                               lambda scaler, jump_multiplier: constrainedGauss(scaler, self.param_step_sizes['alpha_sigsqr_rr'] * jump_multiplier, self.param_ranges['alpha_sigsqr_rr'] ), #alpha_sigsqr_rr_MR
                               lambda scaler, jump_multiplier: constrainedGauss(scaler, self.param_step_sizes['gamma_for_beta_inf'] * jump_multiplier, self.param_ranges['gamma_for_beta_inf'] ), #gamma_for_beta_inf_MR
                               lambda scaler, jump_multiplier: constrainedGauss(scaler, self.param_step_sizes['r_beta0'] * jump_multiplier, self.param_ranges['r_beta0'] ), #r_beta0_MR
                               lambda scaler, jump_multiplier: constrainedGauss(scaler, self.param_step_sizes['alpha_beta'] * jump_multiplier, self.param_ranges['alpha_beta'] ),  #alpha_beta_MR ),  #alpha_beta_MP
                               ]

        return [param_indeces_to_vary, param_functions]

    def printContents(self):
        print ('Printing full contents of storer. ')
        print ('self.population = ' + str(self.population))
        print ('self.dist = ' + str(self.dist))
        print ('self.M = ' + str(self.M))
        print ('self.sigsqr_rr_0 = ' + str(self.sigsqr_rr_0))
        print ('self.sigsqr_rr_inf_to_0_rat = ' + str(self.sigsqr_rr_inf_to_0_rat))
        print ('self.sigsqr_rr_inf = ' + str(self.sigsqr_rr_inf))
        print ('self.r_sigsqr_rr0 = ' + str(self.r_sigsqr_rr0 ))
        print ('self.alpha_sigsqr_rr = ' + str(self.alpha_sigsqr_rr))
        print ('self.gamma_for_beta_inf = ' + str(self.gamma_for_beta_inf))
        print ('self.r_beta0 = ' + str(self.r_beta0))
        print ('self.alpha_beta = ' + str(self.alpha_beta ))
        print ('self.rs = ' + str(self.rs))
        print ('self.halo_center = [' + ','.join([str(center_comp) for center_comp in self.halo_center]) + ']')
        print ('self.overdensity = ' + str(self.overdensity) )
        print ('self.c = ' + str(self.c))

    def recomputeCutoffRadius(self):
        self.c = ComputeVirialRadius (self.M, self.rs, self.overdensity)  #CALCULATE c from M, rs, Q, Delta

    #We assume beta = beta0 + (betainf - beta0) * 1.0 / (1.0 + (R / beta_R0) ** beta_exp)
    #beta = 2.0 * gamma / (1.0 + gamma)
    #Give your answers in terms of gamma0 and gammainf (they can be between -1 and 1)
    #Now a new function: beta = beta0 * exp(-(R/beta_R1) ** 2) / ((R/beta_R2) ** beta_exp) (has beta-> 0 at R = 0 and R = inf; otherwise, things seem to diverge)
    # You now give only gamma0 between -1 and 1
    def __init__(self, population, dist = 'none', M = 'none',
                  sigsqr_rr_0 = 100, sigsqr_rr_inf_to_0_rat = 1, r_sigsqr_rr0 = 1.0, alpha_sigsqr_rr = 2.0, gamma_for_beta_inf = 0.0, r_beta0 = 1.0, alpha_beta = 1.0, #dispersion_rr_params = [100.0, 100.0, 1.0, 2.0], dispersion_beta_params = [0.0, 1.0, 1.0],
                  rs = 500.0, halo_center = [0.0, 0.0, 0.0], c = 5.0,
                  specified_param_ranges = {}, stellar_data = None, overdensity = 200.0):

        BackgroundInfo = BackgroundStorer()
        data_archive = DwarfGalDataArchive()

        dispersion_rr_params = [sigsqr_rr_0, sigsqr_rr_inf_to_0_rat, r_sigsqr_rr0, alpha_sigsqr_rr]
        dispersion_beta_params = [gamma_for_beta_inf, r_beta0, alpha_beta]

        self.overdensity = overdensity
        self.population = population
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
        self.dispersion_rr_params = dispersion_rr_params[:]
        self.dispersion_beta_params = dispersion_beta_params[:]
        self.sigsqr_rr_0, self.sigsqr_rr_inf_to_0_rat, self.r_sigsqr_rr0, self.alpha_sigsqr_rr = dispersion_rr_params
        self.sigsqr_rr_inf = self.sigsqr_rr_0 * self.sigsqr_rr_inf_to_0_rat
        self.gamma_for_beta_inf, self.r_beta0, self.alpha_beta = dispersion_beta_params[:]
        #self.beta_int_dispersion_funct = self.defineBetaIntFunctionForVelocityDispersion()

        #Parameters we must specify about the halo
        self.rs    = rs  #scale radius in pc
        #self.halo_center = [ component / self.rs for component in halo_center ]
        self.halo_center = halo_center
        #phi is angle of halo axis of symmetry off of sky north. Theta is rotation angle along halo axis of symmetry
        # halo_sym_axis is unit vector defining orientatin of disk symmetry axis in sky coordinates.


        self.c = c # ratio of rs to rvir.  Since we're plotting in r/rs and rmax=rvir, this is the distance past which we strongly suspect our model to be innappropriate
        self.c = ComputeVirialRadius (self.M, self.rs, self.overdensity)  #CALCULATE c from M, rs, Q, Delta

        self.param_order_array = ['dist','M',
                                  'sigsqr_rr_0', 'sigsqr_rr_inf_to_0_rat', 'r_sigsqr_rr0', 'alpha_sigsqr_rr', 'gamma_for_beta_inf', 'r_beta0', 'alpha_beta',
                                  'rs' , 'halo_center', 'c',  ] #The order of the parameters in all_params_ordered array
        #Can get the index of a parameter in the all_params_ordered array using .index('name of parameter')
        self.all_params_ordered = [self.dist, self.M,
                                   self.sigsqr_rr_0, self.sigsqr_rr_inf_to_0_rat, self.r_sigsqr_rr0, self.alpha_sigsqr_rr,
                                   self.gamma_for_beta_inf, self.r_beta0, self.alpha_beta,
                                   self.rs, self.halo_center, self.c,  ]
        #The range over which the parameters are allowed to vary
        self.param_ranges = BackgroundInfo.param_ranges
        self.param_step_sizes = BackgroundInfo.param_step_sizes
        #User allowed to override.  Note a range of [a, a] means that the parameter will always have the value a.
        for param in specified_param_ranges.keys():
            self.param_ranges[param] = specified_param_ranges[param]
