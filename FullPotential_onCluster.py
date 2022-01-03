#Define FullPotential object,
# Basically combines the disk and ellipticity computed gravitational potentials to create an interpolator.
# Then gives you, for a given coordinate R and z, the value of a potential at a point (calculated by interpolation).
# Note that this class returns values based on TRUE R and z coordinates (in scale radii).
# It relies on using the scale parameters (zeta, lam, el) to appropriately convert R and z to the R and z -like coordinates used in the potential files.
from PotentialFunctionArray import PotentialFunctionArray
from HaloFunctionComponents import HaloFunctionComponents
from readInPotentialInterpolator import readInPotentialInterpolator
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm
from logList import logList
import time
import math
from PreGeneratedPotentialInterpolator import PreGeneratedPotentialInterpolator

class FullPotential:

    def full_value_single_point(self, R_halo, z_halo, R_disk, z_disk):
        #for reasons of speed, we'll just return the full potential, without worrying about whether it is valid in a region.
        #start = time.time()
        R_to_disk_scale = 1.0/(self.zeta)
        #z_to_disk_scale = 1/(self.zeta*self.lam)
        z_to_disk_scale = 1.0/(self.zeta)
        R_to_halo_scale = 1.0
        z_to_halo_scale = 1.0
        disk_coefficient = 1.0/(self.zeta)
        total_scaling = (-1.0)*(self.gamma*self.M)/(self.rs)
        #end = time.time()
        #print ('Took ' + str(end - start) + 's for first part.')
        #start = time.time()
        disk_part = (self.eps * self.potential_value_single_point(R_disk, z_disk, disk_coefficient, R_to_disk_scale, z_to_disk_scale, self.lam, self.ffunction))
        #end = time.time()
        #print ('Took ' + str(end - start) + 's for second part.')
        start = time.time()
        halo_part = (1-self.eps) * (self.halo_function.getCutOffset() + self.potential_value_single_point(R_halo, z_halo, self.halo_coefficient, R_to_halo_scale, z_to_halo_scale, self.el, self.gammafunction))
        #end = time.time()
        #print ('Took ' + str(end - start) + 's for third part.')
        #disk_coefficient = 1.0

        return (total_scaling * (disk_part + halo_part ))

    def full_value(self,R_halo,z_halo,R_disk,z_disk):

        #for reasons of speed, we'll just return the full potential, without worrying about whether it is valid in a region.
        R_to_disk_scale = 1.0/(self.zeta)
        #z_to_disk_scale = 1/(self.zeta*self.lam)
        z_to_disk_scale = 1.0/(self.zeta)
        R_to_halo_scale = 1.0
        z_to_halo_scale = 1.0
        disk_coefficient = 1.0/(self.zeta)
        #disk_coefficient = 1.0

        return ((-1.0)*(self.gamma*self.M)/(self.rs)
                * (self.eps * self.potential_value(R_disk, z_disk, disk_coefficient, R_to_disk_scale, z_to_disk_scale, self.lam, self.ffunction) +
                 (1-self.eps) * (self.halo_function.getCutOffset() + self.potential_value(R_halo, z_halo, self.halo_coefficient, R_to_halo_scale, z_to_halo_scale, self.el, self.gammafunction)))
                 )


    #should handle R and z at single point
    def potential_value_single_point(self, R, z, overall_coefficient, R_scale, z_scale, shape_parameter, potential_function):
       R_axis = abs(R) * R_scale
       z_axis = abs(z) * z_scale
       valid = (R_axis <= potential_function.Rlike_range[1]) * (R_axis >= potential_function.Rlike_range[0]) * (z_axis <= potential_function.zlike_range[1]) * (z_axis >= potential_function.zlike_range[0])
       if valid:
           potential_val = potential_function.interpolator((R_axis, z_axis, shape_parameter))
       else:
           potential_val = 1.0 / np.sqrt(R ** 2.0 + z ** 2.0)

       return potential_val

    #should handle R and z being meshgrids
    def potential_value(self, R, z, overall_coefficient, R_scale, z_scale, shape_parameter, potential_function, takeAbs = 1):
        Rmesh = 0.0
        zmesh = 0.0
        if takeAbs:
            Rmesh = np.abs(R)
            zmesh = np.abs(z)
        else:
            Rmesh = np.array(R)
            zmesh = np.array(z)
        R_axis_mesh=Rmesh*R_scale
        z_axis_mesh=zmesh*z_scale
        if (np.any(np.isnan(R))):
            print ('!!!!! R contains np.nan !!!!!')
        if (np.any(np.isnan(z))):
            print ('!!!!! z contains np.nan !!!!!')

        where_valid_new = (R_axis_mesh <= potential_function.Rlike_range[1]) * (R_axis_mesh >= potential_function.Rlike_range[0])
        where_valid_new = where_valid_new * (z_axis_mesh <= potential_function.zlike_range[1]) * (z_axis_mesh >= potential_function.zlike_range[0])

        shape_param_array = np.zeros(np.shape(R_axis_mesh)) + shape_parameter

        if (np.any(np.isnan(where_valid_new))):
            print ('!!!!! where_valid_new contains np.nan !!!!!')
        if (np.any(np.isnan(R_axis_mesh))):
            print ('!!!!! R_axis_mesh contains np.nan !!!!!')
        if (np.any(np.isnan(z_axis_mesh))):
            print ('!!!!! z_axis_mesh contains np.nan !!!!!')
        potential_mesh_new = potential_function.interpolator((R_axis_mesh * where_valid_new, z_axis_mesh * where_valid_new, shape_param_array))

        external_array = np.zeros(np.shape(potential_mesh_new))
        external_array[where_valid_new == 0] = np.nan_to_num(1/np.sqrt(Rmesh[where_valid_new == 0]**2 + zmesh[where_valid_new == 0]**2))

        potential_mesh_new = overall_coefficient * potential_mesh_new * where_valid_new + external_array

        return potential_mesh_new

    def showPotential(self,Rhalo,zhalo,Rdisk,zdisk):
        print ('Plotting potential...')
        Rmesh_halo,zmesh_halo = np.meshgrid(Rhalo,zhalo)
        Rmesh_disk,zmesh_disk = np.meshgrid(Rdisk,zdisk)
        potential_mesh = self.full_value(np.abs(Rmesh_halo),np.abs(zmesh_halo),np.abs(Rmesh_disk),np.abs(zmesh_disk))
        #log_levels = logList(min_prob,max_prob,nlevels)
        CS=plt.contour(Rmesh_halo,zmesh_halo,potential_mesh,20)
        CB_contour=plt.colorbar(shrink=0.8,extend='both',format='%.2e')
        plt.show()

        return 0

    def __init__(self, gamma, M, rs, eps, c, lam, zeta, el, disk_potential_file = None, halo_potential_file = None, disk_interpolating_function=None, halo_interpolating_function=None, halo_type = 'NFW', disk_type = 'sech_disk'):
        #print halo_interpolating_function
        if disk_interpolating_function:
            self.ffunction = disk_interpolating_function
        else:
            #self.ffunction=PotentialFunctionArray(disk_potential_file)
            self.ffunction = readInPotentialInterpolator(disk_type)
        if halo_interpolating_function:
            self.gammafunction = halo_interpolating_function
        else:
            #self.gammafunction=PotentialFunctionArray(halo_potential_file)
            self.gammafunction = PotentialFunctionArray(halo_type)
        self.eps=eps
        #self.cut_corr=cutoff_correction
        self.gamma=gamma
        self.M=M
        self.rs=rs
        self.zeta=zeta
        self.lam=lam
        self.c=c
        self.Mhalo=(1-self.eps)*M
        self.Mdisk=self.eps*M
        self.el = el
        #print self.gammafunction
        #print halo_type
        #print self.c
        self.halo_function = HaloFunctionComponents(halo_type,self.c,self.el,halo_interpolating_function = self.gammafunction )
        self.halo_coefficient = self.halo_function.getOverallScaling()
