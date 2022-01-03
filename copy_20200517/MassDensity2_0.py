from FullPotential import FullPotential 
import numpy as np 
import time 
import math
from AstronomicalParameterArchive import AstronomicalParameterArchive 

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
        term1 = (self.C*R_halo**((self.vphi**2)/self.sigsqr))
        exponent = -self.potential.full_value(R_halo,z_halo,R_disk,z_disk)/self.sigsqr
        term2 = np.exp(exponent)
        
        return term1 * term2

    #Turn a 3-dimensional dwarf galaxy distribution into a 2-dimensional distribution by summing a series of cross sections perpendicular to the viewing line of site.
    #The integration_rect_centers tell you the positions along los where you wish to measure a cross section, measured in scale radii from the center of the galaxy.
    #I think I should change this to be the bin borders 
    #Each cross section is weighted by the portion of the viewing axis that lies between the bounds of the two bins.

    #|   |   |   |   | | | |   |   |   |
    #  x   x   x   x  x x x  x   x   x 
    def integrate_los(self, halo_sym_axis, disk_sym_axis, integration_bin_centers, apply_tidal_constraint = 0, los_bin_step = 1.0): #,integration_rect_centers):
        if apply_tidal_constraint:
            integration_bin_limits = self.getLosTidalLimits(halo_sym_axis, disk_sym_axis, tide_source = 'MW',sampling_density = los_bin_step)
            integration_bin_borders = np.arange(integration_bin_limits[0], integration_bin_limits[1] + los_bin_step / 2.0, los_bin_step)
        else:
            n_centers = len(integration_bin_centers) 
            integration_bin_borders = ( [integration_bin_centers[0] - (integration_bin_centers[1] - integration_bin_centers[0]) / 2.0]
                                        + [(integration_bin_centers[i+1] + integration_bin_centers[i]) / 2.0 for i in range(n_centers - 1)]
                                        + [integration_bin_centers[n_centers-1] + (integration_bin_centers[n_centers-1] - integration_bin_centers[n_centers - 2]) /2.0 ] )
            
        integration_bin_borders.sort()
        integrated_array=np.zeros((len(self.R_values),len(self.z_values)))
        #next_x=integration_rect_centers[1]
        #To have a well defined calculation for first step, have to effectively add an extra point
        #prev_x=integration_rect_centers[0] - (next_x - integration_rect_centers[0])
        #Note that we have one more bin border than the number of bins, so the for loop goes up to the length of the array minus 1
        integrated_array = sum([ self.get_rotated_cross_section(halo_sym_axis, disk_sym_axis, integration_bin_centers[i])
                                  * (integration_bin_borders[i + 1] - integration_bin_borders[i]) for i in range(len(integration_bin_centers))])
        #for i in range(len(integration_bin_borders)-1):
        #    #Calculate a cross section at the center of the next two bin borders 
        #    x = integration_bin_centers[i]
        #    #print 'adding x = ' + str(x)
        #    #print 'self.get_rotated_cross_section(phi,theta,a,b,' + str(x) + ')[' + str(len(self.R_values)/2) + ',' + str(len(self.z_values/2)) + '] = ' + str(self.get_rotated_cross_section(phi,theta,a,b,x)[len(self.R_values)/2,len(self.z_values)/2])
        #    #print 'self.get_rotated_cross_section(phi,theta,a,b,x)*(integration_bin_borders[i + 1] - integration_bin_borders[i]) = '
        #    #print self.get_rotated_cross_section(phi,theta,a,b,x)*(integration_bin_borders[i + 1] - integration_bin_borders[i])
        #    integrated_array=integrated_array+self.get_rotated_cross_section(halo_sym_axis, disk_sym_axis, x)*(integration_bin_borders[i + 1] - integration_bin_borders[i])  
        #    end_int_array=time.time()
        #    #print "Took " + str(end_int_array-start_int_array) + " to get new cross section."
        #    #prev_x=x
        #    #if i < len(integration_rect_centers)-2:
        #    #    next_x=integration_rect_centers[i+2]
        #    #else: 
        #    #    next_x = next_x + (next_x - x)
        #    #print "Took " + str(end_los-start_los) + " to compute."
        #print "Took " + str(full_end-full_start) + " to do full computation. "
        return integrated_array 
    
    #We want to mass density on some plane that intersects our distribution at some arbitrary angle and depth
    #The basic way we do the rotations is to define 6, 2d meshgrids.
    #Each of these corresponds to the value of the natural x, y, and z coordinates of the DM disk and halos on the 2d cross-section of interest.  
    def get_rotated_cross_section(self, halo_sym_axis, disk_sym_axis, los_depth):

        #x goes right on sky, z goes up on sky, y goes along los AWAY from observer
        y_sky = los_depth # the distance PAST gal center, alogn los
        xmesh_sky,zmesh_sky = np.meshgrid(self.R_values,self.z_values)

        #We have a vector (x_sky,y_sky,z_sky) in sky coordinates and have the unit vector (x_zp_hat,y_zp_hat,z_zp_hat)
        # of the symmetry (z) axes in the sky coordinates.  So, z_sym = (x_sky,y_sky,z_sky) dot (x_zp_hat,y_zp_hat,z_zp_hat)
        # (just the component of the sky vector along the symmetry z axis)
        zmesh_halo = xmesh_sky * halo_sym_axis[0] + y_sky * halo_sym_axis[1] + zmesh_sky * halo_sym_axis[2]
        zmesh_disk = xmesh_sky * disk_sym_axis[0] + y_sky * disk_sym_axis[1] + zmesh_sky * disk_sym_axis[2]

        #The R_sym known since R_sky ** 2 + z_sky ** 2 = R_sym ** 2 + z_sym ** 2
        # => R_sym ** 2 = R_sky ** 2 + z_sky **2 - z_sym ** 2
        Rmesh_halo = np.sqrt( xmesh_sky * xmesh_sky + y_sky * y_sky + zmesh_sky * zmesh_sky - zmesh_halo * zmesh_halo )
        Rmesh_disk = np.sqrt( xmesh_sky * xmesh_sky + y_sky * y_sky + zmesh_sky * zmesh_sky - zmesh_disk * zmesh_disk )

        #theta rotates x into y (unfortunately reverse of typical convention, but I fear it's coded in deep)
        #xmesh_halo = math.cos(theta) * xmesh_sky + math.sin(theta) * y_sky
        #ymesh_halo = -math.sin(theta) * xmesh_sky + math.cos(theta) * y_sky
        #zmesh_halo = zmesh_sky

        #Now perform phi rotation, which rotates z into x (up on sky into right on sky). Again, reverse of usual convention 
        #xmesh_halo_pre_phi = xmesh_halo
        #xmesh_halo = xmesh_halo_pre_phi * math.cos(phi) - zmesh_halo * math.sin(phi)
        #ymesh_halo = ymesh_halo
        #zmesh_halo = xmesh_halo_pre_phi * math.sin(phi) + zmesh_halo * math.cos(phi) 
        
        #Now I want to rotate disk and halo INDEPENDENTLY.  This is different from my old convention 
        #b rotates x into y for disk
        #xmesh_disk = math.cos(b) * xmesh_sky + math.sin(b) * y_sky
        #ymesh_disk = -math.sin(b) * xmesh_sky + math.cos(b) * y_sky
        #zmesh_disk = zmesh_sky

        #Now perform a rotation, which rotates z into x for disk
        #xmesh_disk_pre_a = xmesh_disk
        #xmesh_disk = xmesh_disk_pre_a * math.cos(a) - zmesh_disk * math.sin(a)
        #ymesh_disk = ymesh_disk
        #zmesh_disk = xmesh_disk_pre_a * math.sin(a) + zmesh_disk * math.cos(a)
         
        #proj_array_where_valid=np.sqrt((ymesh_halo**2+xmesh_halo**2) + zmesh_halo**2) <= self.c
        #proj_array_where_valid should now be determined by TIDAL RADIUS; ie, where potential of Milky Way MATCHES my computed potential
        #phi_MW = -self.gamma * self.M_MW / ( self.dist + los_sky )#Milky Way potential at dSph
        #proj_array_where_valid = FullPotential() > 
        #proj_array_correction = proj_array_where_valid + (np.sqrt((ymesh_halo**2+xmesh_halo**2) + zmesh_halo**2) >= self.c) * math.e**(self.vSup * (1 - np.sqrt((ymesh_halo**2+xmesh_halo**2) + zmesh_halo**2) / self.c) ) 

        #R_halo = np.sqrt(xmesh_halo**2 + ymesh_halo ** 2)
        #R_disk = np.sqrt(xmesh_disk**2 + ymesh_disk ** 2)
        #print 'R_halo = '
        #print R_halo
        #print 'R_disk = '
        #print R_disk 
        proj_array=self.get_density(Rmesh_halo,zmesh_halo,Rmesh_disk,zmesh_disk)

        #NOTE WE RETURN THE TRANSPOSE SINCE THAT IS THE THING THAT WE CAN FEED INTO AN INTERPOLATOR
        #print 'Applying correction.' 
        #return np.transpose(proj_array * proj_array_where_valid + proj_array_correction )
        #return np.transpose(proj_array * proj_array_correction) # if you impose exponential cut off past some limit
        return np.transpose(proj_array) # if you don't want to impose any exponential cut off

    #Computes a non-rotated cross section 
    def compute_cross_section(self,x_depth):
        zmesh,Rmesh=np.meshgrid(self.z_values,self.R_values)
        where_valid=np.sqrt(Rmesh**2+x_depth**2+zmesh**2) <= self.c
        return (self.C*np.exp(-self.potential.full_value(np.sqrt(Rmesh**2+x_depth**2),abs(zmesh))/self.sigsqr)) * where_valid

    def getLosTidalLimits(self,theta, phi, a, b, tide_source = 'MW',sampling_density = 1.0,default_limits = [5.0,5.0]):
        start = time.time()
        time_limit = 10.0
        tidal_force = lambda los: 0.0
        #print 'tide_source = ' + str(tide_source) 
        if tide_source == 'MW':
            tidal_force = lambda los: - self.gamma * self.M_MW / ((self.dist + los * self.rs) ** 2)
            tidal_potential = lambda los: - self.gamma * self.M_MW / ((self.dist + los * self.rs))
            
        R_halo = lambda los: abs(los * math.sqrt((math.sin(theta) * math.cos(phi))**2  +  math.cos(theta)**2  )) 
        z_halo = lambda los: abs(los * math.sin(theta) * math.sin(phi))
        R_disk = lambda los: abs(los * math.sqrt((math.sin(b) * math.cos(a))**2  +  math.cos(b)**2  ))
        z_disk = lambda los: abs(los * math.sin(b) * math.sin(a))
        tidal_limits = [0.0,0.0]
        test_limits = [0.0, -sampling_density]
        while tidal_limits[0] == 0.0:
            dwarf_force = ( -self.potential.full_value(R_halo(test_limits[1]),z_halo(test_limits[1]),R_disk(test_limits[1]),z_disk(test_limits[1]))
                            +self.potential.full_value(R_halo(test_limits[0]),z_halo(test_limits[0]),R_disk(test_limits[0]),z_disk(test_limits[1])) ) / (self.rs*(test_limits[1] - test_limits[0]))
            print 'test_limits = (' + str(test_limits[0]) +', '+ str(test_limits[1]) + ')' 
            print 'dwarf_force = ' + str(dwarf_force)
            print 'tidal_force = ' + str(tidal_force(test_limits[1]))
            if tidal_force(test_limits[1]) + dwarf_force <= 0.0:
                print 'Found tidal limit: ' + str(test_limits[0]) 
                tidal_limits[0] = test_limits[0]
            else:
                if time.time() - start > time_limit / 2.0:
                    print 'Ran out of time finding front tidal limit. Using default'
                    tidal_limits[0] = 5.0
                else:
                    test_limits = [limit - sampling_density for limit in test_limits]
        limiting_potential = tidal_potential (tidal_limits[0]) + self.potential.full_value(R_halo(tidal_limits[0]),z_halo(tidal_limits[0]),R_disk(tidal_limits[0]),z_disk(tidal_limits[1]))
        print 'Found limiting_potential = '+ str(limiting_potential) 
        test_limit = sampling_density
        while tidal_limits[1] == 0.0:
            #if time.time() - start > 2.0:
                #print 'Ran out of time. '
                #break 
            print 'Testing potential at test_limit = ' + str(test_limit) + ' which is ' + str(tidal_potential (test_limit) + self.potential.full_value(R_halo(test_limit),z_halo(test_limit),R_disk(test_limit),z_disk(test_limit))) 
            if  tidal_potential (test_limit) + self.potential.full_value(R_halo(test_limit),z_halo(test_limit),R_disk(test_limit),z_disk(test_limit)) >= limiting_potential:
                tidal_limits[1] = test_limit
            else:
                if time.time() - start > time_limit:
                    print 'Ran out of time finding back tidal limit. Using default.'
                    tidal_limits[1] = 5.0
                else:
                    test_limit = test_limit + sampling_density
                
        return tidal_limits 
                            
    
    def get_R_range(self):
        return self.R_values
    
    def get_z_range(self):
        return self.z_values

    #In defining the MassDensity, we only need to define the key pieces that are needed to specify the parameters.
    #Angles aren't specified here, so that we can examine the same mass density from different angles, using the above functions.  
    def __init__(self, R, z, gamma, M, rs, eps, c, lam, zeta, el, vphi, C, sigsqr, disk_file=None, halo_file=None, dist = 147000.0, virial_suppression_factor = 10, disk_interpolating_function=None, halo_interpolating_function=None, halo_type = 'nfw', disk_type = 'sech_disk'):
        #Using our Jeans analysis, we know that there is a 1-to-1 correspondences between potential value and mass density value.
        astro_archive = AstronomicalParameterArchive()
        self.M_MW = astro_archive.getMWMass()
        self.gamma = gamma
        self.potential=FullPotential(gamma, M, rs, eps, c, lam, zeta, el, disk_potential_file=disk_file, halo_potential_file=halo_file, disk_interpolating_function=disk_interpolating_function, halo_interpolating_function=halo_interpolating_function, halo_type = halo_type, disk_type = disk_type)
        self.C=C
        self.c=c
        self.rs = rs
        self.vphi=vphi
        self.el = el
        self.sigsqr=sigsqr
        #R and z specify the ON SKY coordinates at which we wish to know the value of the mass density profile.  
        self.z_values=z 
        self.R_values=R
        self.R_range=[np.min(R),np.max(R)]
        self.z_range=[np.min(z),np.max(z)]
        self.vSup = virial_suppression_factor
        self.dist = dist
