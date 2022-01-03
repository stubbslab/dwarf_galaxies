import scipy.integrate as integrate
from matplotlib.colors import LogNorm
from logList import logList
import csv
import math 
import numpy as np
from smoothedStarDist import smoothedStarDist 
from DwarfGalDataArchive import DwarfGalDataArchive
from AstronomicalParameterArchive import AstronomicalParameterArchive

def readInBoylanArtificialStarData(halo_number, viewer_phi = 0.0, viewer_theta = 0.0, dist_from_sun = None, arcmin_limits_R = None, arcmin_limits_z = None, inclusion_range = None):
    data_archive = DwarfGalDataArchive()
    astro_params = AstronomicalParameterArchive()
    deg_to_rad = astro_params.getDegToRad() #math.pi/180.0)
    data_file = data_archive.getFile(['Boylan' + str(halo_number), 'dummy'])
    print (data_file)
    #dist_from_sun = data_archive.getDistanceFromSun(population) #{"carina":105000,"fornax":147000,"sculptor":86000,"sextans":86000} in pc
    #total_mass = data_archive.getTotalMass(population) #{"carina":1.28*10**8,"fornax":1.28*10**8,"sculptor":1.28*10**8,"sextans":1.28*10**8} in M_sun 
    #astro_params = AstronomicalParameterArchive() 
    
    #Read the data into various variables 
    xs = []
    x_index = 0
    ys = []
    y_index = x_index + 1
    zs = []
    z_index = y_index + 1
    vxs = []
    vx_index = z_index + 1
    vys = []
    vy_index = vx_index + 1 
    vzs = []
    vz_index = vy_index + 1

    #file = open(data_file, 'rb')
    start_row = 1
    start_column = 0
    with open(data_file) as csvfile:
        
        myreader = csv.reader(csvfile, delimiter = ' ')
        
        row_num = 0
        new_xs = []
        new_ys = []
        new_zs = []
        new_vxs = []
        new_vys = []
        new_vzs = []
        for row in myreader:
            if row_num >= start_row:
                row = row[start_column:]
                if row_num % 1000 == 0:
                    print ('Reading in row ' + str(row_num)) 
                    xs = xs + new_xs[:]
                    ys = ys + new_ys[:]
                    zs = zs + new_zs[:]
                    vxs = vxs + new_vxs[:]
                    vys = vys + new_vys[:]
                    vzs = vzs + new_vzs[:]
                    new_xs = [float(row[x_index]) * 10 ** 3.0] #File has in kpc; I want in pc 
                    new_ys = [float(row[y_index]) * 10 ** 3.0] #File has in kpc; I want in pc 
                    new_zs = [float(row[z_index]) * 10 ** 3.0] #File has in kpc; I want in pc 
                    new_vxs = [float(row[vx_index])]
                    new_vys = [float(row[vy_index])]
                    new_vzs = [float(row[vz_index])]
                else:
                    new_xs = new_xs + [float(row[x_index]) * 10 ** 3.0] #File has in kpc; I want in pc 
                    new_ys = new_ys + [float(row[y_index]) * 10 ** 3.0] #File has in kpc; I want in pc 
                    new_zs = new_zs + [float(row[z_index]) * 10 ** 3.0] #File has in kpc; I want in pc 
                    new_vxs = new_vxs + [float(row[vx_index])]
                    new_vys = new_vys + [float(row[vy_index])]
                    new_vzs = new_vzs + [float(row[vz_index])]
            row_num = row_num + 1
        xs = xs + new_xs[:]
        ys = ys + new_ys[:]
        zs = zs + new_zs[:]
        vxs = vxs + new_vxs[:]
        vys = vys + new_vys[:]
        vzs = vzs + new_vzs[:]
    if inclusion_range is None:
        inclusion_range = [0, len(xs)]
    print ('inclusion_range = ' + str(inclusion_range) )
    xs = xs[min(inclusion_range[0], len(xs) - 1): min(inclusion_range[1], len(xs))]
    ys = ys[min(inclusion_range[0], len(ys) - 1): min(inclusion_range[1], len(ys))]
    zs = zs[min(inclusion_range[0], len(zs) - 1): min(inclusion_range[1], len(zs))]
    vxs = vxs[min(inclusion_range[0], len(vxs) - 1): min(inclusion_range[1], len(vxs))]
    vys = vys[min(inclusion_range[0], len(vys) - 1): min(inclusion_range[1], len(vys))]
    vzs = vzs[min(inclusion_range[0], len(vzs) - 1): min(inclusion_range[1], len(vzs))]
    
    gal_center = [np.mean(xs), np.mean(ys), np.mean(zs)]
    xs_from_center = [x - gal_center[0] for x in xs]
    ys_from_center = [y - gal_center[1] for y in ys]
    zs_from_center = [z - gal_center[2] for z in zs]
    mean_vel = [np.mean(vxs), np.mean(vys), np.mean(vzs)]

    #We choose to position the galaxy at the same distance from Earth as Fornax.
    if dist_from_sun is None: 
        dist_from_sun = data_archive.getDistanceFromSun(['boylan' + str(halo_number),'MR'])

    viewer_vec_in_gal_coords = [dist_from_sun * np.cos(viewer_phi) * np.sin(viewer_theta), dist_from_sun * np.sin(viewer_phi) * np.sin(viewer_theta), dist_from_sun * np.cos(viewer_theta)]
    gal_vec_in_viewer_coords = [0.0, dist_from_sun, 0.0]

    #We want to rotate the galaxy coordinate frame so that the new y-axis lines up with the viewer axis
    
    # Recalling that the defined phi and theta are the angles of the viewere as seen by the galaxy, we reflect them through the galaxy origin... 
    new_y_phi = (viewer_phi + math.pi) % (2.0 * math.pi)
    new_y_theta = math.pi - viewer_theta
    #...and define the axis (in the current galaxy coords) that we would like for the y-axis to rotate into
    target_gal_new_y_axis = [np.cos(new_y_phi) * np.sin(new_y_theta), np.sin(new_y_phi) * np.sin(new_y_theta), np.cos(new_y_theta)]

    #Now perform two rotations.
    z_rot_angle = (new_y_phi + 3.0 * math.pi / 2.0) % (2.0 * math.pi)
    x_rot_angle = new_y_theta - math.pi / 2.0

    rotate_vector = lambda star_vec: [(star_vec[0] * np.cos(z_rot_angle) + star_vec[1] * np.sin(z_rot_angle)),
                                      (-star_vec[0] * np.sin(z_rot_angle) + star_vec[1] * np.cos(z_rot_angle)) * np.cos(x_rot_angle) - star_vec[2] * np.sin(x_rot_angle),
                                      (-star_vec[0] * np.sin(z_rot_angle) + star_vec[1] * np.cos(z_rot_angle)) * np.sin(x_rot_angle) + star_vec[2] * np.cos(x_rot_angle) ]
    # First, rotate around z so that the new x-axis is orthogonl to the target y-axis
    print ('This should be [0,1,0]: rotate_vector(target_gal_new_y_axis) = ' + str(rotate_vector(target_gal_new_y_axis)))

    rotated_stars = [rotate_vector([xs_from_center[i], ys_from_center[i], zs_from_center[i]]) for i in range(len(xs)) ]
    viewed_stars = [[star[0], star[1] + dist_from_sun, star[2]] for star in rotated_stars]
    mean_viewed_star_positions = np.mean(viewed_stars, axis = 0)
    centered_viewed_stars  = np.array(viewed_stars) - [mean_viewed_star_positions[0],0.0, mean_viewed_star_positions[2]]
    rotated_vels = [rotate_vector([vxs[i], vys[i], vzs[i]]) for i in range(len(xs)) ]

    viewed_dists = [np.sqrt(np.sum(np.array(star) ** 2.0)) for star in viewed_stars]
    viewed_thetas = [(np.arccos(viewed_stars[i][2] / viewed_dists[i])) / deg_to_rad for i in range(len(viewed_stars))]
    viewed_Decs = [(math.pi / 2.0 - viewed_thetas[i] * deg_to_rad) / deg_to_rad for i in range(len(viewed_thetas))]

    viewed_RAs = [np.arccos(viewed_stars[i][0] / (np.sin(viewed_thetas[i] * deg_to_rad) * viewed_dists[i])) / deg_to_rad
                  if abs(viewed_thetas[i]) > 0.0 else 0.0 for i in range(len(viewed_stars))]


    cent_ra=np.mean(viewed_RAs)
    cent_dec=np.mean(viewed_Decs)
    
    corr_ra = np.array(viewed_RAs) - cent_ra
    x1=np.cos(np.array(viewed_Decs) * deg_to_rad) * np.cos(corr_ra * deg_to_rad)
    x2=np.cos(np.array(viewed_Decs) * deg_to_rad) * np.sin(corr_ra * deg_to_rad)
    x3=np.sin(np.array(viewed_Decs) * deg_to_rad)

    x1bar=math.cos(cent_dec * deg_to_rad)*x1+math.sin(cent_dec * deg_to_rad)*x3
    x2bar=x2
    x3bar=-math.sin(cent_dec * deg_to_rad)*x1 + math.cos(cent_dec * deg_to_rad)*x3
    corr_ra = np.angle(x1bar + x2bar*1j, deg = True) #centered RA in degrees
    corr_dec = np.arcsin(x3bar) * (1 / deg_to_rad) #centered Dec in degrees
    radial_dists = (np.sqrt(corr_ra**2 + corr_dec**2) * deg_to_rad * dist_from_sun) #physical distance of star from center, in pc


    #We also want to project these angular positions onto the plane tangent to the center of the galaxy.
    # These are important since we are integrating along a square prism
    theta = -corr_dec + 90.0 #traditional spherical theta coord in degrees
    phi = corr_ra #traditional spherical phi coord in degrees
    proj_x = np.sin(phi * deg_to_rad) / np.cos(phi * deg_to_rad) / deg_to_rad
    proj_y = np.cos(theta * deg_to_rad) / (np.cos(phi * deg_to_rad) * np.sin(theta * deg_to_rad)) / deg_to_rad

    los_vels = [vel[1] for vel in rotated_vels] 
    mean_los_vel = np.mean(los_vels)

    #Note that we don't need to correct the los velocities since they are not spectroscopically measured.  We have the actual, correct velocities! 
    helio_prop_mot = [np.mean([vel[0] for vel in rotated_vels] ), np.mean([vel[2] for vel in rotated_vels] )]
    
    #angMotionConvert = astro_params.getAngularMotionConversionFactor()
    #corr_los_vels = [uncorrected_los_vels[i] -  angMotionConvert * dist_from_sun * (viewed_phis[i] * helio_prop_mot[0] + viewed_thetas[i] * helio_prop_mot[1])
    #                 for i in range(len(uncorrected_los_vels))]
    sigSqr_los_vel = sum([((vel - mean_los_vel) ** 2.0) / len(los_vels) for vel in los_vels])

    if arcmin_limits_R is None: 
        arcmin_limits_R = [-max(np.sqrt(np.array(proj_x) ** 2.0 + np.array(proj_y) ** 2.0)) * 1.1 * 60.0, max(np.sqrt(np.array(proj_x) ** 2.0 + np.array(proj_y) ** 2.0)) * 1.1 * 60.0]
    if arcmin_limits_z is None: 
        arcmin_limits_z = arcmin_limits_R
    all_proj_x = proj_x
    all_proj_y = proj_y
    star_projections_in_lims = [[proj_x[i], proj_y[i]] for i in range(len(viewed_stars))
                                if (abs(proj_x[i]) <= max(np.abs(arcmin_limits_R)) / 60.0 ) and (abs(proj_y[i]) <= max(np.abs(arcmin_limits_z)) / 60.0 )]
    print ('len(star_projections_in_lims) = ' + str(len(star_projections_in_lims))) 
    proj_x = np.array([proj[0] for proj in star_projections_in_lims])
    proj_y = np.array([proj[1] for proj in star_projections_in_lims])

    #return_set = [xs,         ys,          zs,      full_rot_xs, full_rot_ys, full_rot_zs,
    #              viewed_xs, viewed_ys, viewed_zs, 
    #              viewed_RAs, viewed_Decs, corr_ra, corr_dec,    los_vels,
    #              proj_x, proj_y, sigSqr_los_vel, arcmin_limits_R, arcmin_limits_z]
    return_set = [xs,         ys,          zs,
                  [star[0] for star in viewed_stars], [star[1] for star in viewed_stars], [star[2] for star in viewed_stars], 
                  viewed_RAs, viewed_Decs,
                  corr_ra, corr_dec,
                  all_proj_x, all_proj_y, 
                  proj_x, proj_y,
                  los_vels, sigSqr_los_vel,
                  arcmin_limits_R, arcmin_limits_z, dist_from_sun]

    
    return return_set

