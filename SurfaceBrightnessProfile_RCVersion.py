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
        
        self.massDensity = MassDensity(self.Rlikeaxis,self.zlikeaxis,self.gamma,self.M,self.rs,self.eps,self.f,self.c,self.lam,self.zeta,self.e,self.vphi,disk_file,halo_file,self.C,self.sigsqr,disk_interpolating_function=disk_interpolating_function,halo_interpolating_function=halo_interpolating_function)
        self.surfaceBrightness = self.massDensity.integrate_los(self.phi,self.theta,self.a,self.b,self.los_bins)
        #Normalize surface brightness to unity?  Now is this truly a normalization, since I'm only summing over points I've measured?  I don't think so...
        #I should normalize somehow.
        #Since this is a SURFACE PROBABILITY DENSITY, I need to multiply by the surface area of the portion of the sky that each point describes.  Yes?
        # Now note that the R and z here describe the physical distance from the center, in scale radii.  So we need to rescale back to, say, radii on sky.  
        #Rmesh, zmesh = np.meshgrid(self.Rlikeaxis,self.zlikeaxis)
        RAmesh_scaleRadii,Decmesh_scaleRadii = np.meshgrid(self.Rlikeaxis * self.rs / self.dist, self.zlikeaxis * self.rs / self.dist)
        RAmesh,Decmesh = np.meshgrid(self.Rlikeaxis, self.zlikeaxis)
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
                    prob_at_point = prob_at_point * (Decmesh[i+1][j]-Decmesh[i][j])
                elif i == np.shape(surface_prob)[0]-1:
                    prob_at_point = prob_at_point * (Decmesh[i][j]-Decmesh[i-1][j])
                else:
                    prob_at_point = prob_at_point * (Decmesh[i+1][j]-Decmesh[i-1][j]) / 2.0
                if j == 0:
                    prob_at_point = prob_at_point * (RAmesh[i][j+1]-RAmesh[i][j])
                elif j == np.shape(surface_prob)[1]-1:
                    prob_at_point = prob_at_point * (RAmesh[i][j]-RAmesh[i][j-1])
                else:
                    prob_at_point = prob_at_point * (RAmesh[i][j+1]-RAmesh[i][j-1]) / 2.0
                normalization_constant = normalization_constant + prob_at_point
        #print 'self.surfaceBrightness[100:111,100:111] before normalization: '
        #print self.surfaceBrightness[100:111,100:111]
        self.surfaceBrightness = self.surfaceBrightness / normalization_constant
        #print 'self.surfaceBrightness[100:111,100:111] after normalization: '
        #print self.surfaceBrightness[100:111,100:111]
                    
        #self.surfaceBrightness = self.surfaceBrightness / np.sum(self.surfaceBrightness) 
        self.onSkyInterpolator = RegularGridInterpolator((self.Rlikeaxis,self.zlikeaxis),self.surfaceBrightness,method="linear")
        
        
