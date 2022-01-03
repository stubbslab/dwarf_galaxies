#Defines a class that relates the characteristic values for disks and halos to the files that contain their computed values.
#The halo files are characterized by an ellipticity. 
#The disk files are characterized by a value of lambda (ratio of scale radius to radius height).
#Will need to be updated as new files with different values are created.  

from ComputationalArchive import ComputationalArchive 

class PotentialArchive:

    def getGeneralFile(self,param,keyword):
        keyword = keyword.lower()
        if keyword == 'nfw':
            return self.getHaloFile(param)
        elif keyword == 'cored':
            return self.getCHaloFile(param)
        elif keyword == 'burkert':
            return self.getBHaloFile(param) 
        elif keyword == 'sech_disk':
            return self.getDiskFile(param)
        elif keyword in ['exp_sphere', 'exp_sph', 'exponential_sphere', 'exponential_exp']:
            return self.getExpSphereFile()
        else:
            print ("keyword '" + keyword + "' not recognized. ") 
            return 0 

    def getDiskFile(self):
        #print 'lam = ' + str(lam) 
        return self.expSpherePotentialFileLibrary[1] 
        
    def getDiskFile(self, lam):
        #print 'lam = ' + str(lam) 
        return self.diskPotentialFileLibrary[lam]
        
    def getHaloFile(self, el):
        #print 'el = ' + str(el) 
        return self.haloPotentialFileLibrary[el]

    def getCHaloFile(self, el):
        #print 'el = ' + str(el)
        return self.coredHaloPotentialFileLibrary[el]
    
    def getBHaloFile(self, el):
        #print 'el = ' + str(el)
        return self.burkertHaloPotentialFileLibrary[el]

    def getPreGenFile(self, pot_type):
        pot_type = pot_type.lower() 
        return self.pre_generated_potential_interpolators[pot_type]

    def __init__(self):
        computational_params_archive = ComputationalArchive()
        potential_dir = computational_params_archive.getPotentialDir() 
        self.diskPotentialFileLibrary = {     #maps lambda values to disk files 
                                         0.2: potential_dir + 'disk_pot_lambda_1_5_master_fixed.csv',
                                         0.15: potential_dir + 'sech_disk_potential_lam_3_20_malpha_800_mbeta_800.csv',
                                         0.1: potential_dir + 'sech_disk_potential_lam_1_10_malpha_800_mbeta_800.csv',
                                         0.05: potential_dir + 'sech_disk_potential_lam_1_20_malpha_800_mbeta_800.csv'}
        self.haloPotentialFileLibrary = {     #maps ellipticity values to halo files
                                         0.05: potential_dir + 'halo_potential_el_1_20_c_5_rmax_8c.csv',
                                         0.1: potential_dir + 'halo_potential_el_1_10_c_5_rmax_8c.csv',
                                         0.15: potential_dir + 'halo_potential_el_3_20_rmax_40.csv',
                                         0.2: potential_dir +  'halo_potential_el_2_10_c_5_rmax_8c.csv', 
                                         0.25: potential_dir + 'halo_potential_el_1_4_rmax_40.csv',
                                         0.3: potential_dir + 'halo_potential_el_3_10_c_5_rmax_8c.csv',  #'halo_potential_el_3_10_c_5_small.csv',  #
                                         0.35: potential_dir + 'halo_potential_el_7_20_rmax_40.csv',
                                         0.4: potential_dir + 'halo_potential_el_4_10_c_5_rmax_8c.csv',
                                         0.45: potential_dir + 'halo_potential_el_9_20_rmax_40.csv', 
                                         0.5: potential_dir +  'halo_potential_el_5_10_c_5_rmax_8c.csv',
                                         0.6: potential_dir +  'halo_potential_el_6_10_c_5_rmax_8c.csv',
                                         0.7: potential_dir +  'halo_potential_el_7_10_c_5_rmax_8c.csv',
                                         0.8: potential_dir +  'halo_potential_el_8_10_c_5_rmax_8c.csv',
                                         0.9: potential_dir +  'halo_potential_el_9_10_c_5_rmax_8c.csv',
                                         0.99: potential_dir + 'halo_potential_el_99_100_rmax_40.csv',
                                         1.0: potential_dir +  'halo_potential_el_10_10_c_5_rmax_8c.csv',
                                         1.1: potential_dir + 'halo_potential_el_11_10_c_5_rmax_8c.csv',
                                         1.11: potential_dir + 'halo_potential_el_10_9_c_5_rmax_8c.csv',
                                         1.2: potential_dir + 'halo_potential_el_12_10_c_5_rmax_8c.csv',
                                         1.25: potential_dir + 'halo_potential_el_10_8_c_5_rmax_8c.csv',
                                         1.43: potential_dir + 'halo_potential_el_10_7_c_5_rmax_8c.csv',
                                         1.67: potential_dir + 'halo_potential_el_10_6_c_5_rmax_8c.csv',
                                         2.0: potential_dir + 'halo_potential_el_10_5_c_5_rmax_8c.csv',
                                         2.5: potential_dir + 'halo_potential_el_10_4_c_5_rmax_8c.csv',
                                         3.33: potential_dir + 'halo_potential_el_10_3_c_5_rmax_8c.csv',
                                         5.0: potential_dir + 'halo_potential_el_10_2_c_5_rmax_8c.csv',
                                         10.0: potential_dir + 'halo_potential_el_10_1_c_5_rmax_8c.csv',
                                         20.0: potential_dir + 'halo_potential_el_20_1_c_5_rmax_8c.csv'} 
                                         # 'test' : 'halo_potential_el_3_10_c_5_small.csv'} #halo_potential_e_1_2_c_5.csv'}
        self.coredHaloPotentialFileLibrary = { #maps ellipticity values to cored halo files
                                              #maps ellipticity values to halo files

                                         0.05: potential_dir + 'cored_halo_potential_el_1_20_c_5_rmax_8c.csv',
                                         0.1: potential_dir + 'cored_halo_potential_el_1_10_c_5_rmax_8c.csv',
                                         0.2: potential_dir +  'cored_halo_potential_el_2_10_c_5_rmax_8c.csv', 
                                         0.3: potential_dir + 'cored_halo_potential_el_3_10_c_5_rmax_8c.csv',  #'cored_halo_potential_el_3_10_c_5_small.csv',  #
                                         0.4: potential_dir + 'cored_halo_potential_el_4_10_c_5_rmax_8c.csv',
                                         0.5: potential_dir +  'cored_halo_potential_el_5_10_c_5_rmax_8c.csv',
                                         0.6: potential_dir +  'cored_halo_potential_el_6_10_c_5_rmax_8c.csv',
                                         0.7: potential_dir +  'cored_halo_potential_el_7_10_c_5_rmax_8c.csv',
                                         0.8: potential_dir +  'cored_halo_potential_el_8_10_c_5_rmax_8c.csv',
                                         0.9: potential_dir +  'cored_halo_potential_el_9_10_c_5_rmax_8c.csv',
                                         1.0: potential_dir +  'cored_halo_potential_el_10_10_c_5_rmax_8c.csv',
                                         1.1: potential_dir + 'cored_halo_potential_el_11_10_c_5_rmax_8c.csv',
                                         1.11: potential_dir + 'cored_halo_potential_el_10_9_c_5_rmax_8c.csv',
                                         1.25: potential_dir + 'cored_halo_potential_el_10_8_c_5_rmax_8c.csv',
                                         1.43: potential_dir + 'cored_halo_potential_el_10_7_c_5_rmax_8c.csv',
                                         1.67: potential_dir + 'cored_halo_potential_el_10_6_c_5_rmax_8c.csv',
                                         2.0: potential_dir + 'cored_halo_potential_el_10_5_c_5_rmax_8c.csv',
                                         2.5: potential_dir + 'cored_halo_potential_el_10_4_c_5_rmax_8c.csv',
                                         3.33: potential_dir + 'cored_halo_potential_el_10_3_c_5_rmax_8c.csv',
                                         5.0: potential_dir + 'cored_halo_potential_el_10_2_c_5_rmax_8c.csv',
                                         10.0: potential_dir + 'cored_halo_potential_el_10_1_c_5_rmax_8c.csv',
                                         20.0: potential_dir + 'cored_halo_potential_el_10_1_c_5_rmax_8c.csv'}
        self.burkertHaloPotentialFileLibrary = { #maps ellipticity values to burkert halo files
                                              #maps ellipticity values to halo files

                                         0.05: potential_dir + 'burkert_halo_potential_el_1_20_c_5_rmax_8c.csv',
                                         0.1: potential_dir + 'burkert_halo_potential_el_1_10_c_5_rmax_8c.csv',
                                         0.2: potential_dir +  'burkert_halo_potential_el_2_10_c_5_rmax_8c.csv', 
                                         0.3: potential_dir + 'burkert_halo_potential_el_3_10_c_5_rmax_8c.csv',  #'burkert_halo_potential_el_3_10_c_5_small.csv',  #
                                         0.4: potential_dir + 'burkert_halo_potential_el_4_10_c_5_rmax_8c.csv',
                                         0.5: potential_dir +  'burkert_halo_potential_el_5_10_c_5_rmax_8c.csv',
                                         0.6: potential_dir +  'burkert_halo_potential_el_6_10_c_5_rmax_8c.csv',
                                         0.7: potential_dir +  'burkert_halo_potential_el_7_10_c_5_rmax_8c.csv',
                                         0.8: potential_dir +  'burkert_halo_potential_el_8_10_c_5_rmax_8c.csv',
                                         0.9: potential_dir +  'burkert_halo_potential_el_9_10_c_5_rmax_8c.csv',
                                         1.0: potential_dir +  'burkert_halo_potential_el_10_10_c_5_rmax_8c.csv',
                                         1.11: potential_dir + 'burkert_halo_potential_el_10_9_c_5_rmax_8c.csv',
                                         1.25: potential_dir + 'burkert_halo_potential_el_10_8_c_5_rmax_8c.csv',
                                         1.43: potential_dir + 'burkert_halo_potential_el_10_7_c_5_rmax_8c.csv',
                                         1.67: potential_dir + 'burkert_halo_potential_el_10_6_c_5_rmax_8c.csv',
                                         2.0: potential_dir + 'burkert_halo_potential_el_10_5_c_5_rmax_8c.csv',
                                         2.5: potential_dir + 'burkert_halo_potential_el_10_4_c_5_rmax_8c.csv',
                                         3.33: potential_dir + 'burkert_halo_potential_el_10_3_c_5_rmax_8c.csv',
                                         5.0: potential_dir + 'burkert_halo_potential_el_10_2_c_5_rmax_8c.csv',
                                         10.0: potential_dir + 'burkert_halo_potential_el_10_1_c_5_rmax_8c.csv',
                                         20.0: potential_dir + 'burkert_halo_potential_el_10_1_c_5_rmax_8c.csv'}
        self.pre_generated_potential_interpolators = { #maps file type to pre-generated potential interpolators
                                                   'nfw' : potential_dir + 'nfw_potential_interpolator.npy',
                                                   'cored' : potential_dir + 'cored_potential_interpolator.npy',
                                                   'burkert' : potential_dir + 'burkert_potential_interpolator.npy', 
                                                   'sech_disk' : potential_dir + 'sech_disk_potential_interpolator.npy',
                                                   'exp_sphere': potential_dir + 'exponential_sphere_pot.npy' }
