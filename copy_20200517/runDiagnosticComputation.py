#Want to run a code that checks to make sure I am correctly computing the probability
#I need to define my own star positions at points where I know the values of the potential

from PotentialArchive import PotentialArchive
from PotentialFunctionArray import PotentialFunctionArray
from DwarfGalaxyParametersStorer import DwarfGalaxyParametersStorer
from DwarfGalDataArchive import DwarfGalDataArchive 
from SurfaceBrightnessProfile import SurfaceBrightnessProfile
from AstronomicalParameterArchive import AstronomicalParameterArchive
from compRZLos import compRZLos
from MassDensity import MassDensity
from FullPotential import FullPotential
import math
import numpy as np

def runDiagnosticComputation(el=0.2,  lam=0.1, rs=500,  phi=0.0, theta=0.0, eps=0.0, zeta=1.0, a=0.0,   b=0.0,   galaxy='fornax',  pop='MR',  pop_selection_method='none', withDisk=0, outer_step=5.0,disk_interpolating_function = None,halo_interpolating_function = None): 
    #parameters at present are: [el,  lam, rs,  phi, theta, eps, zeta, a,   b,   galaxy,  pop,  pop_selection_method, withDisk, outer_step]
    
    astro_archive = AstronomicalParameterArchive()
    gamma = astro_archive.getGamma()
    pot_archive = PotentialArchive()
    dSph_archive = DwarfGalDataArchive()

    halo_file = pot_archive.getHaloFile(el)
    if not halo_interpolating_function:
        halo_interpolating_function = PotentialFunctionArray(halo_file)
    disk_file = pot_archive.getDiskFile(lam) 
    if not disk_interpolating_function:
        disk_interpolating_function = PotentialFunctionArray(disk_file)

    storer = DwarfGalaxyParametersStorer([galaxy,pop],withDisk,el = el, rs = rs, phi = phi, theta = theta,  lam = lam, zeta = zeta, eps = eps, a = a, b = b)
    f = storer.f
    c = storer.c
    dist = storer.dist
    M = storer.M
    e = storer.e
    sigsqr = storer.sigsqr

    
    arcmin_limits_R = [-90.0,90.0]
    arcmin_limits_z = [-60.0,60.0]
    arcmin_limits_los = [-60.0,60.0]
    R, z, los = compRZLos(zeta, zeta * outer_step, outer_step, arcmin_limits_R, arcmin_limits_z, arcmin_limits_los, rs, dist)
    print 'R = '
    print R
    print 'z = '
    print z
    print 'los = '
    print los

    Rmesh,zmesh = np.meshgrid(R,z)

    test_pot = FullPotential(gamma,M,rs,eps,f,c,lam,zeta,e,disk_file,halo_file,disk_interpolating_function = disk_interpolating_function, halo_interpolating_function = halo_interpolating_function)

    #print 'test_pot.full_value = '
    #print test_pot.full_value(Rmesh,zmesh,Rmesh,zmesh) 


        
    
    test_mass_profile = MassDensity(R,z,gamma,M,rs,eps,f,c,lam,zeta,e,0.0,disk_file,halo_file,1.0,sigsqr,dist = dist, disk_interpolating_function = disk_interpolating_function, halo_interpolating_function = halo_interpolating_function)
    print 'sigsqr = ' + str(sigsqr)
    #print 'phi = ' +str(phi)
    #print 'theta = ' +str(theta)
    #print 'a = ' +str(a)
    #print 'b = ' +str(b)

    #print 'Rmesh = '
    #print Rmesh
    #print 'zmesh = '
    #print zmesh 
    
    #print 'test_mass_profile.get_rotated_cross_section(phi,theta,a,b,0.0) = '
    #print test_mass_profile.get_rotated_cross_section(phi,theta,a,b,0.0)

    #print 'np.transpose(np.exp  (-1 *  test_pot.full_value(Rmesh,zmesh,Rmesh,zmesh) / sigsqr )) = '
    #print np.transpose(np.exp  (-1 *  test_pot.full_value(Rmesh,zmesh,Rmesh,zmesh) / sigsqr ) )

    x_centers = los
    x_borders = [los[0] - (los[1] - los[0]) / 2.0] + [(los[j+1] + los[j]) / 2.0 for j in range(len(los) -1) ] + [los[len(los)-1] + (los[len(los)-1] - los[len(los)-2]) / 2.0]
    x_bin_widths = [x_borders[j+1] - x_borders[j] for j in range(len(x_centers))]
    for i in range(len(x_centers)):
        #print 'x_centers[i] = ' + str(x_centers[i]) 
        x_bin_width = x_bin_widths[i]
        #print 'test_mass_profile.get_rotated_cross_section(phi,theta,a,b,x_centers[i]) = '
        #print test_mass_profile.get_rotated_cross_section(phi,theta,a,b,x_centers[i])
        #print 'bin_width = ' + str(bin_width)
    print 'test_mass_profile.integrate_los(phi,theta,a,b,los) = '
    print test_mass_profile.integrate_los(phi,theta,a,b,los)
    print 'x_borders = '
    print x_borders
    print 'x_centers = '
    print x_centers
    print 'x_bin_widths = '
    print x_bin_widths
    normalized_mass_profile = test_mass_profile.integrate_los(phi,theta,a,b,los)
    R_borders = [R[0] - (R[1] - R[0]) / 2.0] + [(R[j+1] + R[j]) / 2.0 for j in range(len(R) -1) ] + [R[len(R)-1] + (R[len(R)-1] - R[len(R)-2]) / 2.0]
    z_borders = [z[0] - (z[1] - z[0]) / 2.0] + [(z[j+1] + z[j]) / 2.0 for j in range(len(z) -1) ] + [z[len(z)-1] + (z[len(z)-1] - z[len(z)-2]) / 2.0]
    R_bin_widths = [R_borders[j+1] - R_borders[j] for j in range(len(R))]
    z_bin_widths = [z_borders[j+1] - z_borders[j] for j in range(len(z))]
    area_array = np.array([[R_bin * z_bin for z_bin in z_bin_widths] for R_bin in R_bin_widths] )
    print 'rs /dist = ' + str(rs/dist) 
    area_array = area_array * (rs / dist) ** 2.0
    print 'area_array = '
    print area_array 
    normalized_mass_profile = normalized_mass_profile / np.sum(normalized_mass_profile * area_array)
    print 'normalized_mass_profile = '
    print normalized_mass_profile 
    

    #observation_mask = dSph_archive.getObservationMask([galaxy,pop],R * rs / dist * 180.0 / math.pi, z * rs / dist * 180.0 / math.pi)
    observation_mask = None 
    
    test_surface_profile = SurfaceBrightnessProfile(R,z,los, gamma, storer,disk_file,halo_file,disk_interpolating_function = disk_interpolating_function, halo_interpolating_function = halo_interpolating_function, observation_mask = observation_mask )
    test_interpolator = test_surface_profile.onSkyInterpolator

    print 'test_interpolator((0.0,0.0)) = ' + str( test_interpolator((0.0,0.0)) )

    for r in R:
        for Z in z:
            print 'test_interpolator((' + str(r) + ',' + str(Z) + ')) = ' + str( test_interpolator((r,Z)) )

    test_RA_scaleRadii = [0.0,R[0],R[1],R[2]]
    test_Dec_scaleRadii = [0.0,z[0],z[1],z[2]]

    for i in range(len(test_RA_scaleRadii)):
        print 'test_surface_profile.sumLogSurfaceBrightness(' + str(test_RA_scaleRadii[i] * rs/dist * 180 / math.pi) + ',' + str(test_Dec_scaleRadii[i] * rs/dist * 180 / math.pi) + ') = '
        print test_surface_profile.sumLogSurfaceBrightness(test_RA_scaleRadii[i] * rs/dist * 180 / math.pi, test_Dec_scaleRadii[i] * rs/dist * 180 / math.pi)

    print 'test_surface_profile.sumLogSurfaceBrightness(test_RA_scaleRadii, test_Dec_scaleRadii) = '
    print test_surface_profile.sumLogSurfaceBrightness(np.array([RA * rs/dist * 180 / math.pi for RA in test_RA_scaleRadii]), np.array([Dec * rs/dist * 180 / math.pi for Dec in test_Dec_scaleRadii]))

    print 'Done. '

    


