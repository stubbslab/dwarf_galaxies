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
    def full_value(self,R_halo,z_halo,R_disk,z_disk):
        #for reasons of speed, we'll just return the full potential, without worrying about whether it is valid in a region.  
        R_to_disk_scale = 1.0/(self.zeta)
        #z_to_disk_scale = 1/(self.zeta*self.lam)
        z_to_disk_scale = 1.0/(self.zeta) 
        R_to_halo_scale = 1.0
        z_to_halo_scale = 1.0
        disk_coefficient = 1.0/(self.zeta)
        #disk_coefficient = 1.0
        halo_coefficient = self.halo_function.getOverallScaling()
        if (np.any(np.isnan(R_halo))):
            print ('!!!!! R_halo contains np.nan !!!!!') 
        if (np.any(np.isnan(z_halo))):
            print ('!!!!! z_halo contains np.nan !!!!!') 
        if (np.any(np.isnan(R_disk))):
            print ('!!!!! R_disk contains np.nan !!!!!') 
        if (np.any(np.isnan(z_disk))):
            print ('!!!!! z_disk contains np.nan !!!!!') 
        return ((-1.0)*(self.gamma*self.M)/(self.rs) 
                * (self.eps * self.potential_value(R_disk, z_disk, disk_coefficient, R_to_disk_scale, z_to_disk_scale, self.lam, self.ffunction) + 
                 (1-self.eps) * (self.halo_function.getCutOffset() + self.potential_value(R_halo,z_halo,halo_coefficient,R_to_halo_scale,z_to_halo_scale, self.el, self.gammafunction))) 
                 )
    #should handle R and z being meshgrids
    def potential_value(self, R, z, overall_coefficient, R_scale, z_scale, shape_parameter, potential_function, takeAbs = 1):
        #print (Ibreak)
        #print 'R_scale = ' + str(R_scale)
        #print 'z_scale = ' + str(z_scale) 
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
        
        
        #print 'np.unique(R_axis_mesh) = '
        #print np.unique(R_axis_mesh)
        #print 'np.unique(z_axis_mesh) = '
        #print np.unique(z_axis_mesh) 
        potential_mesh=np.zeros(np.shape(R_axis_mesh),float)
        potential_zero=0
        potential_nonzero=0
        #if not np.shape(R_axis_mesh): #This checks if the R and z arrays have only one value; ie, we don't want to try to iterate
        #    if (R_axis_mesh <= potential_function.Rlike_range[1] and R_axis_mesh >= potential_function.Rlike_range[0]
        #        and
        #        z_axis_mesh <= potential_function.zlike_range[1] and z_axis_mesh >= potential_function.zlike_range[0] ):
        #        potential_mesh_new=potential_function.interpolating_function([R_axis_mesh,z_axis_mesh])
        #    else:
        #        potential_mesh_new=0
        #    return (-1)*(self.gamma*self.M)/(self.rs)*(self.eps/self.zeta)*disk_mesh_new
        #print 'np.shape(R_axis_mesh) = '
        #print np.shape(R_axis_mesh)
        #print 'np.shape(z_axis_mesh) = '
        #print np.shape(z_axis_mesh)
        #print 'potential_function.Rlike_range = '
        #print potential_function.Rlike_range
        #print 'potential_function.zlike_range = '
        #print potential_function.zlike_range
        where_valid_new = (R_axis_mesh <= potential_function.Rlike_range[1]) * (R_axis_mesh >= potential_function.Rlike_range[0])
        where_valid_new = where_valid_new * (z_axis_mesh <= potential_function.zlike_range[1]) * (z_axis_mesh >= potential_function.zlike_range[0])
        #where_valid_new = np.zeros(np.shape(where_valid_new)) + 1.0 
        #print 'where_valid_new = '
        #print where_valid_new
        
        #print 'coordinate_array = '
        #print coordinate_array 
        #coordinate_array=np.dstack((R_axis_mesh * where_valid_new, z_axis_mesh * where_valid_new))
        #print 'np.shape(coordinate_array) = '
        #print np.shape(coordinate_array) 
        #full_mesh = np.zeros(list(np.shape(coordinate_array)[0:3]) + [3]) + shape_parameter
        #full_mesh[:,:,:,0:2] = coordinate_array

        shape_param_array = np.zeros(np.shape(R_axis_mesh)) + shape_parameter

        #print '(R_axis_mesh * where_valid_new).tolist() = '
        #print (R_axis_mesh * where_valid_new).tolist()
        #print '(z_axis_mesh * where_valid_new).tolist() = '
        #print (z_axis_mesh * where_valid_new).tolist()
        #flat_R = (R_axis_mesh * where_valid_new).flatten().tolist()
        #flat_z = (z_axis_mesh * where_valid_new).flatten().tolist()
        #for i in range(len(R_axis_mesh.flatten().tolist())):
        #    print 'flat_R[i] = ' + str(flat_R[i])
        #    print 'flat_z[i] = ' + str(flat_z[i])
        #    print potential_function.interpolator((flat_R[i], flat_z[i], shape_param_array[0]))
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
        #potential_mesh_new = overall_coefficient * potential_mesh_new * where_valid_new + (1 - where_valid_new) * np.nan_to_num(1/np.sqrt(Rmesh**2 + zmesh**2))

        #print 'np.shape(potential_mesh_new) = '
        #print np.shape(potential_mesh_new)
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
