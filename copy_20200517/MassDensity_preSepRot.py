from FullPotential import FullPotential 
import numpy as np 
import time 
import math

class MassDensity:
    
    #Here, we rotate the mass distribution by arbitrary angles and project it onto the sky
    #  (the normal spherical angles describing how the rotated the observer sees the symmetry axis of the 
    #  galaxy) and then projecting that rotated shape onto the 2d viewing surface of the ski.  The
    #  x-axis is defined as the observer's (pre-rotation) los.
    #We've replaced R with x and y so we now have a full set of 3d Cartesian coordiantes.
    #z is the viewer's vertical sky coordinate
    #x is the distance from the center along the viewing direction
    #y is the viewer's horizontal sky coordinate.
    #So the final plot will be one of mass density as a function of z and y; x is the projection axis.
    #There are also 3 coordinates that define the configuration of the potential on the sky.
    #theta is the angle rotating the x axis into the z axis of the halo
    #a and b, the two other angles are defined as the difference between the disk coordinates and the halo coordinates
    #a is the angle rotating x into z
    #b is the angle rotating x into y 
    #sum ('integrate') along line of sight range for a given angle.  Use stored R_range to determine number and size of integration steps.  
    
    def generate_star_prof(self,theta,sky_project_prob_array,star_number):
        #sky_project_prob should be the result of self.integrate_los(theta,integration_rect_centers) called already
        #that leaves us with an (unnormalized) array of numbers 
        # whose relative sizes indicate the likelihood of finding a star at that point on the sky
        # (rho, remember, is really a probability density, so probability per unit volue
        #  We then just integrated along los, so that turned it into a probability density per unit area).
        Rmesh,zmesh=np.meshgrid(self.get_R_range(),self.get_z_range())
        summed_prob_array = sky_project_prob_array/np.sum(sky_project_prob_array)
        stars=np.zeros(np.shape(sky_project_prob_array))
        #print 'pre_summed_array'
        #print summed_prob_array
        prev_prob=0
        for  i in range(np.shape(summed_prob_array)[0]):
            for j in range(np.shape(summed_prob_array)[1]): 
                summed_prob_array[i][j] = prev_prob + summed_prob_array[i][j]
                prev_prob=summed_prob_array[i][j]
        #print 'post_summed_array'
        #print summed_prob_array
    
        for star in range(star_number):
            placement=random.random()
            #print 'placement = ' + str(placement)
            index = np.argmin(np.abs(placement-summed_prob_array))
            #print 'index = ' + str(index)
            i = index/(np.shape(summed_prob_array)[1])
            j = index%(np.shape(summed_prob_array)[1])
            stars[i][j] = stars[i][j] + 1
        #print 'stars = '
        #print stars
        return stars
            
                     
        
        #Now we randomly drop a star onto the field at a weighted random location
    
    def get_density(self,R_halo,z_halo,R_disk,z_disk):
        return ((self.C*R_halo**((self.vphi**2)/self.sigSqr))*np.exp(-self.potential.full_value(R_halo,z_halo,R_disk,z_disk)/self.sigSqr))

    #Turn a 3-dimensional dwarf galaxy distribution into a 2-dimensional distribution by summing a series of cross sections perpendicular to the viewing line of site.
    #The integration_rect_centers tell you the positions along los where you wish to measure a cross section, measured in scale radii from the center of the galaxy.
    #I think I should change this to be the bin borders 
    #Each cross section is weighted by the portion of the viewing axis that lies between the bounds of the two bins.

    #|   |   |   |   | | | |   |   |   |
    #  x   x   x   x  x x x  x   x   x 
    def integrate_los(self,phi,theta,a,b,integration_bin_borders): #,integration_rect_centers):
        integration_bin_borders.sort()
        full_start=time.time()
        integrated_array=np.zeros((len(self.R_values),len(self.z_values)))
        #next_x=integration_rect_centers[1]
        #To have a well defined calculation for first step, have to effectively add an extra point
        #prev_x=integration_rect_centers[0] - (next_x - integration_rect_centers[0])
        #Note that we have one more bin border than the number of bins, so the for loop goes up to the length of the array minus 1
        for i in range(len(integration_bin_borders)-1):
            start_los = time.time()
            #Calculate a cross section at the center of the next two bin borders 
            x = (integration_bin_borders[i + 1] + integration_bin_borders[i]) / 2.0  #integration_rect_centers[i]
            #print 'adding x = ' + str(x)
            start_int_array=time.time()
            #print 'self.get_rotated_cross_section(phi,theta,a,b,' + str(x) + ')[' + str(len(self.R_values)/2) + ',' + str(len(self.z_values/2)) + '] = ' + str(self.get_rotated_cross_section(phi,theta,a,b,x)[len(self.R_values)/2,len(self.z_values)/2]) 
            integrated_array=integrated_array+self.get_rotated_cross_section(phi,theta,a,b,x)*(integration_bin_borders[i + 1] - integration_bin_borders[i])  
            end_int_array=time.time()
            #print "Took " + str(end_int_array-start_int_array) + " to get new cross section."
            #prev_x=x
            #if i < len(integration_rect_centers)-2:
            #    next_x=integration_rect_centers[i+2]
            #else: 
            #    next_x = next_x + (next_x - x)
            end_los = time.time()
            #print "Took " + str(end_los-start_los) + " to compute."
        full_end=time.time()
        #print "Took " + str(full_end-full_start) + " to do full computation. "
        return integrated_array
    
    #We want to mass density on some plane that intersects our distribution at some arbitrary angle and depth
    #The basic way we do the rotations is to define 6, 2d meshgrids.
    #Each of these corresponds to the value of the natural x, y, and z coordinates of the DM disk and halos on the 2d cross-section of interest.  
    def get_rotated_cross_section(self,phi,theta,a,b,los_depth):
        #proj_array=np.zeros((len(self.R_values),len(self.z_values)))+0.5
        #offset=np.shape(proj_array)[1]/2
        #x_sky=los_depth
        
        #ymesh_sky,zmesh_sky = np.meshgrid(self.R_values,self.z_values)
        
        #xmesh_halo=math.cos(theta) * x_sky - math.sin(theta) * zmesh_sky
        #ymesh_halo=ymesh_sky
        #zmesh_halo=math.sin(theta) * x_sky + math.cos(theta) * zmesh_sky
        
        #Now we must add on the extra phi rotation, which just affects angle of viewing box on sky
        #ymesh_halo_pre_phi = ymesh_halo 
        #ymesh_halo = math.cos(phi) * ymesh_halo + math.sin(phi) * zmesh_halo
        #zmesh_halo = - math.sin(phi) * ymesh_halo_pre_phi + math.cos(phi) * zmesh_halo

        #x goes right on sky, z goes up on sky, y goes along los AWAY from observer
        y_sky = los_depth # the distance PAST gal center, alogn los
        xmesh_sky,zmesh_sky = np.meshgrid(self.R_values,self.z_values)

        #theta rotates x into y (unfortunately reverse of typical convention, but I fear it's coded in deep)
        xmesh_halo = math.cos(theta) * xmesh_sky + math.sin(theta) * y_sky
        ymesh_halo = -math.sin(theta) * xmesh_sky + math.cos(theta) * y_sky
        zmesh_halo = zmesh_sky

        #Now perform phi rotation, which rotates z into x (up on sky into right on sky). Again, reverse of usual convention.
        
        xmesh_halo_pre_phi = xmesh_halo
        xmesh_halo = xmesh_halo_pre_phi * math.cos(phi) - zmesh_halo * math.sin(phi)
        ymesh_halo = ymesh_halo
        zmesh_halo = xmesh_halo_pre_phi * math.sin(phi) + zmesh_halo * math.cos(phi) 
        
        

        xmesh_disk= (math.cos(b) * xmesh_halo + math.sin(b) * ymesh_halo) * math.cos(a) - math.sin(a) * zmesh_halo
        ymesh_disk= -math.sin(b) * xmesh_halo + math.cos(b) * ymesh_halo
        zmesh_disk= (math.cos(b) * xmesh_halo + math.sin(b) * ymesh_halo) * math.sin(a) + math.cos(a) * zmesh_halo
         
        proj_array_where_valid=np.sqrt((ymesh_halo**2+xmesh_halo**2) + zmesh_halo**2) <= self.c
        proj_array_correction = proj_array_where_valid + (np.sqrt((ymesh_halo**2+xmesh_halo**2) + zmesh_halo**2) >= self.c) * math.e**(self.vSup * (1 - np.sqrt((ymesh_halo**2+xmesh_halo**2) + zmesh_halo**2) / self.c) ) 
       
        proj_array=self.get_density(np.sqrt(xmesh_halo**2+ymesh_halo**2),np.abs(zmesh_halo),np.sqrt(xmesh_disk**2+ymesh_disk**2),np.abs(zmesh_disk))

        #print 'Applying correction.' 
        #return np.transpose(proj_array * proj_array_where_valid + proj_array_correction )
        #return np.transpose(proj_array * proj_array_correction) # if you impose exponential cut off past some limit
        return np.transpose(proj_array) # if you don't want to impose any exponential cut off

    #Computes a non-rotated cross section 
    def compute_cross_section(self,x_depth):
        zmesh,Rmesh=np.meshgrid(self.z_values,self.R_values)
        where_valid=np.sqrt(Rmesh**2+x_depth**2+zmesh**2) <= c
        return (self.C*np.exp(-self.potential.full_value(np.sqrt(Rmesh**2+x_depth**2),abs(zmesh))/self.sigSqr)) * where_valid
    
    def get_R_range(self):
        return self.R_values
    
    def get_z_range(self):
        return self.z_values

    #In defining the MassDensity, we only need to define the key pieces that are needed to specify the parameters.
    #Angles aren't specified here, so that we can examine the same mass density from different angles, using the above functions.  
    def __init__(self,R,z,gamma,M,rs,eps,f,c,lam,zeta,e,vphi,disk_file,halo_file,C,sigSqr,virial_suppression_factor = 10, disk_interpolating_function=None,halo_interpolating_function=None):
        #Using our Jeans analysis, we know that there is a 1-to-1 correspondences between potential value and mass density value. 
        self.potential=FullPotential(gamma,M,rs,eps,f,c,lam,zeta,e,disk_file,halo_file,disk_interpolating_function,halo_interpolating_function)
        self.C=C
        self.c=c
        self.vphi=vphi
        self.sigSqr=sigSqr
        #R and z specify the ON SKY coordinates at which we wish to know the value of the mass density profile.  
        self.z_values=z 
        self.R_values=R
        self.R_range=[np.min(R),np.max(R)]
        self.z_range=[np.min(z),np.max(z)]
        self.vSup = virial_suppression_factor
