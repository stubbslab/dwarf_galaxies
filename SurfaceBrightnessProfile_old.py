#Define the SurfaceBrightnessProfile class.
#This takes in a MassDensity class, integrates it along a los, and then creates an interpolator for the observed surface brightness.  

from MassDensity import MassDensity 
import numpy as np
from scipy.interpolate import RegularGridInterpolator 
from DwarfGalaxyParametersStorer import DwarfGalaxyParametersStorer
from AstronomicalParameterArchive import AstronomicalParameterArchive

class SurfaceBrightnessProfile:

    #Sum the log of the value of the dSph surface brightness profile at every observed star position.
    #Must have corrRa and corrDec given in degrees from center of galaxy.  
    def sumLogSurfaceBrightness(self,corrRa,corrDec):
        astro_archive = AstronomicalParameterArchive()
        deg_to_rad = astro_archive.getDegToRad()
        corr_Ra_in_scale_radii = corrRa*deg_to_rad*self.dist/self.rs
        corr_Dec_in_scale_radii = corrDec*deg_to_rad*self.dist/self.rs
        #print 'np.max(corr_Ra_in_scale_radii) = ' + str(np.max(corr_Ra_in_scale_radii))
        #print 'np.min(corr_Ra_in_scale_radii) = ' + str(np.min(corr_Ra_in_scale_radii))
        #print 'np.max(corr_Dec_in_scale_radii) = ' + str(np.max(corr_Dec_in_scale_radii))
        #print 'np.min(corr_Dec_in_scale_radii) = ' + str(np.min(corr_Dec_in_scale_radii))
        log_likelihoods = np.log(self.onSkyInterpolator(np.dstack((corr_Ra_in_scale_radii,corr_Dec_in_scale_radii))[0]))
        return np.sum(log_likelihoods)

    #Sum the value of the dSph surface brightness profile at every observed star position
    #Must have corrRa and corrDec given in degrees from center of galaxy
    #Note that this quantity CANNOT serve as a good proxy for the likelihood of an array, since we're adding together values where a true probability would demand multiplication. 
    def sumSurfaceBrightness(self,corrRa,corrDec):
        #print np.dstack(np.meshgrid(corrRa,corrDec))
        #print np.shape(np.dstack(np.meshgrid(corrRa,corrDec)))
        #print np.meshgrid(corrRa,corrDec)
        astro_archive = AstronomicalParameterArchive()
        deg_to_rad = astro_archive.getDegToRad()
        corr_Ra_in_scale_radii = corrRa*deg_to_rad*self.dist/self.rs 
        corr_Dec_in_scale_radii = corrDec*deg_to_rad*self.dist/self.rs
        return np.sum(self.onSkyInterpolator(np.dstack((corr_Ra_in_scale_radii,corr_Dec_in_scale_radii))[0]))

    #Multiply the value of the dSph surface brightness profile at every observed star position
    #Must have corrRa and corrDec given in degrees from center of galaxy 
    #Note that this quantity CAN serve as a proxy for probability, but this value can quickly grow too small for python to effectively handle well. 
    def multiplySurfaceBrightness(self,corrRA,corrDec):
        astro_archive = AstronomicalParameterArchive()
        deg_to_rad = astro_archive.getDegToRad()
        corr_Ra_in_scale_radii = corrRa*deg_to_rad*self.dist/self.rs 
        corr_Dec_in_scale_radii = corrDec*deg_to_rad*self.dist/self.rs
        return np.prod(self.onSkyInterpolator(np.dstack((corr_Ra_in_scale_radii,corr_Dec_in_scale_radii))[0]))

    #Return the value of the dSph surface brightness profile at every observed star position
    #Must have corrRa and corrDec given in degrees from center of galaxy  
    def starProbabilityValues(self,corrRa,corrDec):
        astro_archive = AstronomicalParameterArchive()
        deg_to_rad = astro_archive.getDegToRad()
        corr_Ra_in_scale_radii = corrRa*deg_to_rad*self.dist/self.rs 
        corr_Dec_in_scale_radii = corrDec*deg_to_rad*self.dist/self.rs
        #print 'np.dstack(np.meshgrid(corr_Ra_in_scale_radii,corr_Dec_in_scale_radii))'
        #print np.dstack(np.meshgrid(corr_Ra_in_scale_radii,corr_Dec_in_scale_radii))
        #print 'np.shape(np.dstack((corr_Ra_in_scale_radii,corr_Dec_in_scale_radii))[0])'
        #print np.shape(np.dstack((corr_Ra_in_scale_radii,corr_Dec_in_scale_radii))[0])
        return self.onSkyInterpolator(np.dstack((corr_Ra_in_scale_radii,corr_Dec_in_scale_radii))[0])

    #The old method I was using the match up the stars on the sky.
    # Obviosly wrong, since it was taking the inner product of the RA and Dec lists, rather than just matching them up elementwise.  
    #def starProbabilityValuesOld(self,corrRa,corrDec ):
    #    astro_archive = AstronomicalParameterArchive()
    #    deg_to_rad = astro_archive.getDegToRad()
    #    corr_Ra_in_scale_radii = corrRa*deg_to_rad*self.dist/self.rs 
    #    corr_Dec_in_scale_radii = corrDec*deg_to_rad*self.dist/self.rs
    #    #print 'np.dstack(np.meshgrid(corr_Ra_in_scale_radii,corr_Dec_in_scale_radii))'
    #    #print np.dstack(np.meshgrid(corr_Ra_in_scale_radii,corr_Dec_in_scale_radii))
    #    print 'np.shape(np.dstack(np.meshgrid(corr_Ra_in_scale_radii,corr_Dec_in_scale_radii)))'
    #    print np.shape(np.dstack(np.meshgrid(corr_Ra_in_scale_radii,corr_Dec_in_scale_radii)))
    #    return self.onSkyInterpolator(np.dstack(np.meshgrid(corr_Ra_in_scale_radii,corr_Dec_in_scale_radii)))

    #I now also need to account for the fact that 


    #A surface brightness profile is defined by an underlying massDensity, along with some angles that allow one to integrate the density along the los.
    #We then create an interpolator over the grid generated after integration along los
    # (note this is different than the graviational potential interpolator used in computing the mass density profile). 
    def __init__(self,R,z,los_bins,gamma,parameter_storer,disk_file,halo_file,C=1.0,disk_interpolating_function=None,halo_interpolating_function=None):
        
        self.Rlikeaxis = R
        self.zlikeaxis = z
        self.los_bins = los_bins 
        self.gamma = gamma
        self.C = C
        self.dist = parameter_storer.dist 
        self.rs = parameter_storer.rs 
        self.el = parameter_storer.el
        self.lam = parameter_storer.lam
        self.eps = parameter_storer.eps
        self.zeta = parameter_storer.zeta
        self.rs = parameter_storer.rs
        self.phi = parameter_storer.phi
        self.theta = parameter_storer.theta
        self.a = parameter_storer.a 
        self.b = parameter_storer.b
        self.M = parameter_storer.M
        self.f = parameter_storer.f
        self.c = parameter_storer.c
        self.e = parameter_storer.e
        self.vphi = parameter_storer.vphi
        self.sigsqr = parameter_storer.sigsqr 

        #print 'self.Rlikaxis = '
        #print self.Rlikeaxis
        #print 'self.zlikaxis = '
        #print self.zlikeaxis
        #print 'self.gamma = ' + str(self.gamma) 
        #print 'self.M = ' + str(self.M)
        #print 'self.rs = ' + str(self.rs) 
        #print 'self.eps = ' + str(self.eps) 
        #print 'self.f = ' + str(self.f) 
        #print 'self.c = ' + str(self.c) 
        #print 'self.lam = ' + str(self.lam)
        #print 'self.zeta = ' + str(self.zeta) 
        #print 'self.e = ' + str(self.e)
        #print 'self.vphi = ' + str(self.vphi) 
        #print 'disk_fil = ' + str(disk_file) 
        #print 'halo_file = ' + str(halo_file) 
        #print 'self.C = ' + str(self.C) 
        #print 'self.sigsqr = ' + str(self.sigsqr)
        #print 'self.zeta = ' + str(self.zeta)
        #print 'self.eps = ' + str(self.eps)
        #print 'self.dist = ' + str(self.dist) 
        
        self.massDensity = MassDensity(self.Rlikeaxis,self.zlikeaxis,self.gamma,self.M,self.rs,self.eps,self.f,self.c,self.lam,self.zeta,self.e,self.vphi,disk_file,halo_file,self.C,self.sigsqr,disk_interpolating_function=disk_interpolating_function,halo_interpolating_function=halo_interpolating_function)
        self.surfaceBrightness = self.massDensity.integrate_los(self.phi,self.theta,self.a,self.b,self.los_bins)
        #Normalize surface brightness to unity?  Now is this truly a normalization, since I'm only summing over points I've measured?  I don't think so...
        #I should normalize somehow.
        #Since this is a SURFACE PROBABILITY DENSITY, I need to multiply by the surface area of the portion of the sky that each point describes.  Yes?
        # Here, we normalize the surface density to get probability density per square scale radius. That normalization is convenient for ...
        # But for comparison of two profiles on the sky, I need to normalize by angular size on sky
        # (that way, I'm guarenteed to have each differential area element the same).  
        #Rmesh, zmesh = np.meshgrid(self.Rlikeaxis,self.zlikeaxis)
        RAmesh,Decmesh = np.meshgrid(self.Rlikeaxis * self.rs / self.dist, self.zlikeaxis * self.rs / self.dist)
        RAmesh_scaleRadii,Decmesh_scaleRadii = np.meshgrid(self.Rlikeaxis, self.zlikeaxis)
        #print 'RAmesh[10:20,10:20]'
        #print RAmesh[10:20,10:20]
        #print 'Decmesh[10:20,10:20]'
        #print Decmesh[10:20,10:20]
        normalization_constant = 0.0
        surface_prob = self.surfaceBrightness
        for i in range(np.shape(surface_prob)[0]):
            for j in range(np.shape(surface_prob)[1]):
                prob_at_point = surface_prob[i][j]
                #print 'surface_prob[i][j] = ' + str(surface_prob[i][j]) 
                if i == 0:
                    prob_at_point = prob_at_point * (Decmesh_scaleRadii[i+1][j]-Decmesh_scaleRadii[i][j])
                elif i == np.shape(surface_prob)[0]-1:
                    prob_at_point = prob_at_point * (Decmesh_scaleRadii[i][j]-Decmesh_scaleRadii[i-1][j])
                else:
                    prob_at_point = prob_at_point * (Decmesh_scaleRadii[i+1][j]-Decmesh_scaleRadii[i-1][j]) / 2.0
                if j == 0:
                    prob_at_point = prob_at_point * (RAmesh_scaleRadii[i][j+1]-RAmesh_scaleRadii[i][j])
                elif j == np.shape(surface_prob)[1]-1:
                    prob_at_point = prob_at_point * (RAmesh_scaleRadii[i][j]-RAmesh_scaleRadii[i][j-1])
                else:
                    prob_at_point = prob_at_point * (RAmesh_scaleRadii[i][j+1]-RAmesh_scaleRadii[i][j-1]) / 2.0
                normalization_constant = normalization_constant + prob_at_point
        #print 'Surface profile normalization_constant = ' + str(normalization_constant) 
        self.surface_prob = surface_prob 
        #print 'self.surfaceBrightness[10:20,10:20] before normalization: '
        #print self.surfaceBrightness[10:20,10:20]
        self.surface_prob = self.surface_prob / normalization_constant
        self.surfaceBrightness = self.surfaceBrightness / normalization_constant 
        self.surfaceProbDensity_perSqrScaleRadius = self.surfaceBrightness
        self.surfaceProbDensity = self.surfaceProbDensity_perSqrScaleRadius * (self.dist / self.rs) ** 2 
        #print 'self.surfaceProbDensity[10:20,10:20] after normalization: '
        #print self.surfaceProbDensity[10:20,10:20]
        #print 'self.surfaceProbDensity[10:20,10:20] = '
        #print self.surfaceProbDensity[10:20,10:20] 
                    
        #self.surfaceBrightness = self.surfaceBrightness / np.sum(self.surfaceBrightness)
        #print 'self.Rlikeaxis = '
        #print self.Rlikeaxis
        #print 'self.zlikeaxis = '
        #print self.zlikeaxis
        #print 'self.surfaceBrightness = '
        #print self.surfaceBrightness 
        #self.onSkyBrightnessInterpolator = RegularGridInterpolator((self.Rlikeaxis,self.zlikeaxis),self.surfaceBrightness,method="linear")

        #This gives an interpolator for surface probability density per square radian (not per square degree or scale radii, but square radian)
        self.scaleRadiusInterpolator = RegularGridInterpolator((self.Rlikeaxis,self.zlikeaxis),self.surfaceProbDensity_perSqrScaleRadius,method = 'linear') 
        self.onSkyInterpolator = RegularGridInterpolator((self.Rlikeaxis,self.zlikeaxis), self.surfaceProbDensity, method = 'linear')

        #print 'self.onSkyBrightnessInterpolator((-1.27760528,0.60904147)) = ' + str(self.onSkyBrightnessInterpolator((-1.27760528,0.60904147)))
        #print 'self.onSkyInterpolator((-1.27760528,0.60904147)) = ' + str(self.onSkyInterpolator((-1.27760528,0.60904147)))
        
        
        
