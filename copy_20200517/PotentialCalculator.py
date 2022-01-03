from HaloFunctionComponents import HaloFunctionComponents 
import numpy as np
import scipy.integrate as integrate 


class CylindricalPotentialCaculator:

    #Assume R, z passed as points on a mesh 
    def getDPhiDZ(self, R, z):
        el = self.el
        c = self.c 
        mu_with_tau = lambda localR, localz, tau: np.sqrt(localR ** 2.0 / (1.0 + tau) + localz ** 2.0 / (el ** 2.0 + tau)) 
        raw_tauc = self.tauc_funct(R,z)
        tauc = np.where(raw_tauc > 0.0, raw_tauc, 0.0)
        raw_dtauc_dz = self.dtauc_dz_funct(R, z)
        dtauc_dz = np.where(tauc > 0.0, raw_dtauc_dz, 0.0)

        if len(np.shape(R)) == 0: 
            mass_int_bounds = [mu_with_tau(R, z, tauc) , c ]
            mass_int = self.mass_int_funct (mass_int_bounds[1]) - self.mass_int_funct (mass_int_bounds[0])
            term2 = integrate.quad (lambda tau: ((tau  + 1) ** 2.0 + (tau + el ** 2.0)) ** -0.5 * z / (el ** 2.0 + tau) * self.mass_funct(mu_with_tau(R, z, tau)), tauc, np.inf )[0]
        else:
            mass_int = np.zeros(np.shape(R)) 
            term2 = np.zeros(np.shape(R))
            for i in range(np.shape(term2)[0]):
                for j in range(np.shape(term2)[1]):
                    ind_R = R[i][j]
                    ind_z = z[i][j]
                    ind_tauc = tauc[i][j]
                    ind_mu_with_tau = mu_with_tau(ind_R, ind_z, tauc)[i][j]
                    mass_int_bounds = [ind_mu_with_tau, c]
                    mass_int[i][j] = self.mass_int_funct (mass_int_bounds[1]) - self.mass_int_funct (mass_int_bounds[0])
                    int_res =    integrate.quad (lambda tau: ((tau  + 1) ** 2.0 + (tau + el ** 2.0)) ** -0.5 * ind_z / (el ** 2.0 + tau) * self.mass_funct(mu_with_tau(ind_R, ind_z, tau)), ind_tauc, np.inf )[0]
                    ind_term2 = int_res 
                    term2[i][j] = ind_term2 

        term1 = np.where(R ** 2.0 + z ** 2.0 > 0.0, dtauc_dz * ((tauc  + 1) ** 2.0 + (tauc + el ** 2.0)) ** -0.5 * mass_int, 0.0)
        return -2 * np.pi * self.el * (-term1 - term2)

    #Assume R, z passed as points on a mesh 
    def getDPhiDR(self, R, z):
        el = self.el
        c = self.c 
        mu_with_tau = lambda localR, localz, tau: np.sqrt(localR ** 2.0 / (1.0 + tau) + localz ** 2.0 / (el ** 2.0 + tau))
        raw_tauc = self.tauc_funct(R,z)
        tauc = np.where(raw_tauc > 0.0, raw_tauc, 0.0)
        raw_dtauc_dR = self.dtauc_dR_funct(R, z)
        dtauc_dR = np.where(tauc > 0.0, raw_dtauc_dR, 0.0)

        if len(np.shape(R)) == 0: 
            mass_int_bounds = [mu_with_tau(R, z, tauc) , c ]
            mass_int = self.mass_int_funct (mass_int_bounds[1]) - self.mass_int_funct (mass_int_bounds[0])
            term2 = integrate.quad (lambda tau: ((tau  + 1) ** 2.0 + (tau + el ** 2.0)) ** -0.5 * z / (el ** 2.0 + tau) * self.mass_funct(mu_with_tau(R, z, tau)), tauc, np.inf )[0]
        else:
            mass_int = np.zeros(np.shape(R)) 
            term2 = np.zeros(np.shape(R))
            for i in range(np.shape(term2)[0]):
                for j in range(np.shape(term2)[1]):
                    ind_R = R[i][j]
                    ind_z = z[i][j]
                    ind_tauc = tauc[i][j]
                    ind_mu_with_tau = mu_with_tau(ind_R, ind_z, tauc)[i][j]
                    mass_int_bounds = [ind_mu_with_tau, c]
                    mass_int[i][j] = self.mass_int_funct (mass_int_bounds[1]) - self.mass_int_funct (mass_int_bounds[0])
                    int_res =    integrate.quad (lambda tau: ((tau  + 1) ** 2.0 + (tau + el ** 2.0)) ** -0.5 * ind_R / (1.0 + tau) * self.mass_funct(mu_with_tau(ind_R, ind_z, tau)), ind_tauc, np.inf )[0]
                    ind_term2 = int_res 
                    term2[i][j] = ind_term2 

        term1 = np.where(R ** 2.0 + z ** 2.0 > 0.0, dtauc_dR * ((tauc  + 1) ** 2.0 + (tauc + el ** 2.0)) ** -0.5 * mass_int, 0.0) 
        return -2 * np.pi * self.el * (-term1 - term2)

    def getDeltaR(self, R, z): 

        return 1 

    #Note: this gives properties of tau up to a scaling by (G * r_s) * (M / r_s ** 3.0) 
    def __init__(self, el, c, halo_type):

        self.el = el
        self.c = c 
        self.halo_type = halo_type 
        hfc = HaloFunctionComponents(halo_type, self.c, self.el, halo_interpolating_function = 'dummy')
        halo_scaling = hfc.getOverallScaling() 
        #self.tauc_funct = lambda R, z: np.where((R ** 2.0 + z ** 2.0) > 0.0, - (1 + el ** 2.0) / 2.0 + 0.5 * ((1 - el) ** 2.0 + 4.0 * c ** 2.0 / (R ** 2.0 + z ** 2.0)) ** 0.5, np.inf )
        #self.tauc_funct = lambda R, z: np.where((R ** 2.0 + z ** 2.0) > 0.0, -0.5 * ((1 + el ** 2.0) - (R ** 2.0 + z ** 2.0) / (c ** 2.0)) + 0.5 * np.sqrt( ((1.0 + el ** 2.0) - (R ** 2.0 + z ** 2.0) / (c ** 2.0) ) ** 2.0 - 4.0 * (el ** 2.0 + (R ** 2.0 * el ** 2.0 + z ** 2.0) / (c ** 2.0))), 0.0)
        self.tauc_funct = lambda R, z:  -0.5 * ((1 + el ** 2.0) - (R ** 2.0 + z ** 2.0) / (c ** 2.0)) + 0.5 * np.sqrt( ((1.0 + el ** 2.0) - (R ** 2.0 + z ** 2.0) / (c ** 2.0) ) ** 2.0 - 4.0 * (el ** 2.0 - (R ** 2.0 * el ** 2.0 + z ** 2.0) / (c ** 2.0))) 
        #self.dtauc_dz_funct = lambda R, z: np.where((R ** 2.0 + z ** 2.0) > 0.0, (z / c ** 2.0) * (1.0 + 1.0 / np.sqrt((1.0 + el ** 2.0 - (R ** 2.0 + z ** 2.0) / (c ** 2.0) ) ** 2.0 - 4.0 * (el ** 2.0 - (R ** 2.0 * el ** 2.0 + z ** 2.0 ) / c ** 2.0 ) ) * z / c ** 2.0 * (2 - (1+ el ** 2.0 - (R ** 2.0 + z ** 2.0) / c ** 2.0 )) ), np.inf )
        self.dtauc_dz_funct = lambda R, z: (z / c ** 2.0) * (1.0 + 1.0 / np.sqrt((1.0 + el ** 2.0 - (R ** 2.0 + z ** 2.0) / (c ** 2.0) ) ** 2.0 - 4.0 * (el ** 2.0 - (R ** 2.0 * el ** 2.0 + z ** 2.0 ) / c ** 2.0 ) ) * z / c ** 2.0 * (2 - (1+ el ** 2.0 - (R ** 2.0 + z ** 2.0) / c ** 2.0 )) )
        #self.dtauc_dR_funct = lambda R, z: np.where((R ** 2.0 + z ** 2.0) > 0.0, (z / c ** 2.0) * (1.0 + 1.0 / np.sqrt((1.0 + el ** 2.0 - (R ** 2.0 + z ** 2.0) / (c ** 2.0) ) ** 2.0 - 4.0 * (el ** 2.0 - (R ** 2.0 * el ** 2.0 + z ** 2.0 ) / c ** 2.0 ) ) * z / c ** 2.0 * (2 * el ** 2.0 - (1+ el ** 2.0 - (R ** 2.0 + z ** 2.0) / c ** 2.0 )) ), np.inf )
        self.dtauc_dR_funct = lambda R, z: (z / c ** 2.0) * (1.0 + 1.0 / np.sqrt((1.0 + el ** 2.0 - (R ** 2.0 + z ** 2.0) / (c ** 2.0) ) ** 2.0 - 4.0 * (el ** 2.0 - (R ** 2.0 * el ** 2.0 + z ** 2.0 ) / c ** 2.0 ) ) * z / c ** 2.0 * (2 * el ** 2.0 - (1+ el ** 2.0 - (R ** 2.0 + z ** 2.0) / c ** 2.0 )) ) 

        mu_funct = lambda R, z: np.sqrt(R ** 2.0 + (z / self.el) ** 2.0) 

        #mass_funct_int is: int d\mu \mu \rho(\mu), a term that shows up a lot 
        if halo_type in ['NFW', 'nfw','Nfw']: #\rho(\mu') = M / (4 \pi Q rs ** 3.0) * 1.0 / f(c) * 1.0 / mu(R, z) * 1.0 / (1.0 + mu(R,z)) ** 2.0
            min_mu = 0.0001 
            self.mass_funct = lambda mu: np.where(mu > min_mu, halo_scaling / (4.0 * np.pi * self.el) * mu ** (-1.0) * (1.0 + mu) ** (-2.0), halo_scaling / (4.0 * np.pi * self.el) * min_mu ** (-1.0) * (1.0 + min_mu) ** (-2.0)) 
            self.mass_int_funct = lambda mu: halo_scaling / (4.0 * np.pi * self.el) * (-1.0 / (1.0 + mu) ** 2.0) 

        elif halo_type in ['acored', 'Acored','ACORED']:
            self.mass_funct = lambda mu: halo_scaling / (4.0 * np.pi * self.el) * (1.0 + mu) ** (-3.0)
            self.mass_int_funct = lambda mu: halo_scaling/ (4.0 * np.pi * self.el) * (-(2.0 * mu + 1)/ (2.0 * (1.0 + mu) ** 2.0)) 
        
        elif halo_type in ['burkert','Burkert', 'BURKERT']:
            self.mass_funct = lambda mu: halo_scaling / (4.0 * np.pi * self.el) * (1.0 + mu) ** (-1.0) * (1.0 + mu ** 2.0) ** (-1.0)
            self.mass_int_funct = lambda mu: halo_scaling  / (4.0 * np.pi * self.el ) * mu ** 3.0 / 3.0 + mu ** 2.0 / 2.0

        elif halo_type in ['test', 'Test', 'TEST', 'uniform', 'Uniform', 'UNIFORM', 'constant','Constant', 'CONSTANT']:
            halo_scaling = 1.0 
            self.mass_funct = lambda mu: halo_scaling / (4.0 * np.pi * self.el)
            self.mass_int_funct = lambda mu: halo_scaling / (4.0 * np.pi * self.el) * mu
        else:
            print ('Halo type: ' + halo_type + 'not recognized!. ' ) 

        
        
        
        

        
