#Here, I will write a simple galaxy creation routine for various values of DDM.  
#For Fornax, we extend it to allow for multiple stellar populations (only thing that changes is sigSqr)

from MassDensity import MassDensity
from FullPotential import FullPotential
from PotentialFunctionArray import PotentialFunctionArray
from logList import logList
import numpy as np
import random
import csv
import math
import scipy as sci
import scipy.integrate as integrate
from scipy import ndimage
import os
import matplotlib
import matplotlib.cm as cm
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import time
import scipy.optimize as optimize
from numpy import inf
from operator import itemgetter
from scipy.interpolate import RegularGridInterpolator 
from matplotlib.colors import LogNorm
from matplotlib import colors, ticker
#%matplotlib notebook
print 'test'
#One will need to generate two potential function arrays: one for the disk, one for the halo.
# Do this using the potentialFunctionArray object, one for the elliptical halo and one for the disk.
# The results of those output should be fed in as the disk_funct and halo_funct variables below.  
# Make sure that the zeta and el parameters match the files you read in. 


def generateGalaxyDensityProfile(lam, zeta, eps, el, theta, a, b, disk_file, halo_file, disk_funct, halo_funct):
    #initialize hard coded parameters

    galaxy = 'fornax' #name of galaxy, for documentation purposes
    gal_dist = 147000 #distance from sun to galaxy, in pc
    
    #velocity dispersion squared, in (km/s)^2, for each population.  At present, still assuming constant
    sigSqr=10.0**2  
    
    #really, G,M,rs are only important as (G M)/(rs sigsqr)
    gamma=4.307 * 10**(-3) #constant scaling between (G M)/ (rs sigsqr) when M, rs, sigsqr are expressed in particular units
    M=1.28*10**8 #Mass of ALL dynamical matter, in terms of solar masses (for our purposes, set to DM)
    M_star=20*10**6
    
    #scaling constants between probability density and luminosity contribution for each population
    #equal to typicaly luminosity of the stars, multiplied by abundance, then normalized
    tot_stars=562.0
    typ_mag = (21.2+21.6)/2.0
    
    #scaling constant of star mass densities.  
    C=10**(typ_mag/-2.5)*tot_stars
    
    disk_frac = 0.3 #fraction of young stars in disk
    mass_to_other = 3 #relative typical mass of the young stars vs old stars
    eps_stars=(M_star/M) * tot_stars/tot_stars * disk_frac * mass_to_other #fraction of total mass in stellar disk purely due to stars
    
    print 'the epsilon due to the stars alone is roughly ' + str(eps_stars)
    
    #really, vphi is only important as (vphi^2)/(sigsqr)
    vphi=0
    #a=1
    #eps=0.05 #eps_stars #fraction of DM in DM disk
    
    #potential_dir="/media/sf_Sasha/Documents/Harvard/physics/randall/dwarfGalaxyProject/sech_potential_computations/"
    #disk_file=potential_dir+'Potential_1to10.csv'
    potential_dir="/media/sf_Sasha/Documents/Harvard/physics/randall/dwarfGalaxyProject/potential_tables/"
    disk_file=potential_dir+'disk_pot_lambda_1_10_master3.csv'
    halo_file=potential_dir+'halo_potential_el_3_10_rmax_20.csv'
    #el=0.3 #ellipticity of the halo 
    e=math.sqrt(1-(1-el)**2) #eccentricity of the halo
    #lam=0.1 #ratio of zd to Rd: lam=zd/Rd
    #zeta=0.1 # ratio of rs to Rd: zeta=Rd/rs
    #Q=1 # paramater defining how elliptical a dwarf galaxy is 
    c=5  # ratio of rs to rvir.  Since we're plotting in r/rs and rmax=rvir, this is the maximum value for our plots
    f= math.log(1+c) - c/(1+c) #function accounting for cutoff of potential at r/rs = c
    print 'Parameters set.'
    
    #Now we need to estimate rs from other paramaters, using a known mass, Mk inside some known radius rk
    #This is still done assuming a spherical distribution.  
    #  This is becuase doing this sort of thing analytically for an elliptic halo is challenging.
    #  My hope is that it still gives us the right rough estimation, but more thought on this wouldn't be a bad idea. 
    rk = 300 #radius inside of which total mass is known, in parsecs
    Mk = 0.7 * 10**7 #total mass inside sphere of radius rk, in solar masses
    #Define function for finding eta_k = rk/rs, breaking into two parts for easier writing.  
    part1 = lambda eta_k: (1-eps)/f * (np.log(1 + eta_k) - eta_k/(1 + eta_k))
    part2 = lambda eta_k: eps * np.tanh(eta_k/(2 * zeta * lam))*(1-np.exp(-eta_k/zeta)*(1 + eta_k/zeta))
    our_funct = lambda eta_k: M/Mk*(part1(eta_k) + part2(eta_k)) - 1 
    eta_k = optimize.ridder(our_funct,0,100)
    rs = rk/eta_k
    print 'Halo scale radius is computed to be ' + str(rs) + ' parsecs. '
    
    print str(gamma*M/(sigSqr * rs)) + ' is our exponential scaling.'
    
    #Let's only go out to 60 arcmin (the size typically shown in the observational data)
    arcmin_limits = 60.0
    num_pnts_R = 400
    num_pnts_z = 500
    
    R_inner=np.arange(0,arcmin_limits*0.05 + 0.01,0.01)
    z_inner=np.arange(0,arcmin_limits*0.05 + 0.01,0.01)
    R_outer=np.arange(0,arcmin_limits + 0.25,0.25)
    z_outer=np.arange(0,arcmin_limits + 0.25,0.25)
    R = np.unique(np.concatenate((-R_inner,-R_outer,R_inner,R_outer)))*1/60.0*1/180.0*math.pi*gal_dist/rs
    z = np.unique(np.concatenate((-z_inner,-z_outer,z_inner,z_outer)))*1/60.0*1/180.0*math.pi*gal_dist/rs
    
    mass_prof=MassDensity(R,z,gamma,M,rs,eps,f,c,lam,zeta,e,vphi,disk_file,halo_file,C,sigSqr,disk_funct,halo_funct)
    
    #Now we can try adding up los.  We can also check to see how far we differ from rotational symmetry for a spherically 
    #  symmetric potential.  
    #theta=math.pi*0.0
    #a=math.pi*0.25
    #b=math.pi*0.5
    x_integration_pnts=np.unique(np.concatenate((np.arange(c*-1.0,c*1.002,0.1,dtype=np.float),np.arange(c*-0.5*0.1,c*0.5*0.1+0.005,0.005,dtype=np.float))))
    #print 'integration points for rectangle integration method = ' + str(x_integration_pnts)
    los_total=mass_prof.integrate_los(theta,a,b,x_integration_pnts)
    
    zmesh,Rmesh=np.meshgrid(mass_prof.get_z_range(),mass_prof.get_R_range())
    zmesh = zmesh*180.0/math.pi*rs/gal_dist*60.0
    Rmesh = Rmesh*180.0/math.pi*rs/gal_dist*60.0
    #print mass_prof_MR.get_z_range()
    #print mass_prof_MR.get_R_range()
    
    matplotlib.rcParams['xtick.direction'] = 'out'
    matplotlib.rcParams['ytick.direction'] = 'out'
    
    plt.figure(figsize=(9,7))
    
    #im=plt.imshow(cross_sect,interpolation='spline16',origin='lower',cmap=cm.cool,extent=(-5.4,5.5,-2.5,2.5)) #not working, for some reaso
    
    log_levels=logList(np.min(np.abs(los_profs[pop_of_interest])),np.max(np.abs(los_profs[pop_of_interest])),20)
    CS=plt.contour(Rmesh,zmesh,los_profs[pop_of_interest],norm=LogNorm(),levels=log_levels)
    CB_contour=plt.colorbar(shrink=0.8,extend='both',format='%.2f')
    
    #CB_gradient=plt.colorbar(im,shrink=0.8, orientation='horizontal')
    
    l,b,w,h = plt.gca().get_position().bounds
    ll,bb,ww,hh = CB_contour.ax.get_position().bounds
    
    CB_contour.ax.set_position ([ll,b+0.1*h,ww,h*0.8])
    
    plt.title('Total Probability Density Profile for Having a Star')
    plt.xlabel('R')
    plt.ylabel('z')
    
    pltFileDir="/media/sf_Sasha/Documents/Harvard/physics/randall/dwarfGalaxyProject/dSphDataFiles/plots/"
    #plt.savefig(pltFileDir + galaxy + '_' + 'pop_of_interest'  + '_probability_isophotes' + '.png')
    
    plt.show()
    
    header='galaxy = ' + galaxy + ' zeta='+str(zeta)+",c="+str(c)+",lambda="+str(lam)+",theta="+str(theta)+",epsilon="+str(eps)+",el="+str(el)+",a="+str(a)+",b="+str(b)
    full_file_dir = '/media/sf_Sasha/Documents/Harvard/physics/randall/dwarfGalaxyProject/simulatedProfileFiles/'
    file_name = 'DGP_' + galaxy + '_' + '_zeta_' + str(zeta) + '_lambda_' + str(lam) + '_epsilon_'+ str(eps) + '_el_' + str(el) + '_theta_' + str(theta) + '_a_' + str(a) + '_b_' + str(b) + '.txt'
    np.savetxt(full_file_dir + file_name,los_total,delimiter=',')#,header=header)
    file_name_R_coords = 'DGP_RAs_' + galaxy + '_zeta_' + str(zeta) + '_lambda_' + str(lam) + '_epsilon_'+ str(eps) + '_el_' + str(el) + '_theta_' + str(theta) + '_a_' + str(a) + '_b_' + str(b) +'.txt'
    file_name_z_coords = 'DGP_Dec_' + galaxy + '_zeta_' + str(zeta) + '_lambda_' + str(lam) + '_epsilon_'+ str(eps) + '_el_' + str(el) + '_theta_' + str(theta) + '_a_' + str(a) + '_b_' + str(b) +'.txt'
    np.savetxt(full_file_dir + file_name_R_coords,Rmesh,delimiter=',')
    np.savetxt(full_file_dir + file_name_z_coords,zmesh,delimiter=',')
    
print 'done'