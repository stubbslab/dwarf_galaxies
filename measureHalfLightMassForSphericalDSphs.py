from ObservedGalaxyStarData import ObservedGalaxyStarData
from DwarfGalDataArchive import DwarfGalDataArchive 
import math
from AstronomicalParameterArchive import AstronomicalParameterArchive

def measureHalfLightMassForSphericalDSphs(rs, rs_err, M, M_err, c, galaxy, populations, halo_model, pop_selection_method = 'metal_rigid', fixed_half_light_radius = None):
    star_data_by_pop = []
    G =  4.302 * 10 ** -3.0 #Newton's constant in pc (km/s)**2.0 / M_sun
    walker_half_mass_estimator = lambda half_light_radius, sigsqr: 5.0 * half_light_radius * sigsqr / (2.0 * G)
    walker_half_mass_err = lambda half_light_radius, sigsqrE: 5.0 * half_light_radius / 2.0 * G * sigsqrE

    dsph_archive = DwarfGalDataArchive() 
    dist = dsph_archive.getDistanceFromSun([galaxy, 'dummy_pop'])
    astro_params = AstronomicalParameterArchive() 
    deg_to_rad = astro_params.getDegToRad() #math.pi/180.0

    print halo_model.lower() 
    if halo_model.lower() in ['nfw']:
        cut_function = lambda a: math.log(a + 1.0) - a / (1.0 + a)
        cut_function_der = lambda a: 1.0 / (a * (1.0 + a) ** 2.0)
    elif halo_model.lower() in ['cored']:
        cut_function = lambda a: math.log(a + 1.0) - (3.0 * a ** 2.0 + 2.0 * a) / (2.0 * (1.0 + a) ** 2.0)
        cut_function_der = lambda a: 1.0 / (1.0 + a) ** 3.0
    elif halo_model.lower() in ['burkert']:
        cut_function = lambda a: 1.0 / 4.0 * (2.0 * math.log(1.0 + a) + math.log(1.0 + a ** 2.0) - 2.0 * math.atan(a) )
        cut_function_der = lambda a: 1.0 / ((1.0 + a) * (1.0 + a ** 2.0))
    else:
        print 'No recognized halo given.  Returning unit scaling. '
        cut_function = lambda a: 1.0 
    true_half_mass_funct = lambda half_light_radius: M / (4.0 * math.pi) * cut_function(half_light_radius / rs) / cut_function(c)
    true_half_mass_funct_err = lambda half_light_radius: math.sqrt( (M_err * 1.0 / (4.0 * math.pi) * cut_function(half_light_radius / rs) / cut_function(c)) ** 2.0
                                                                    + (rs_err * M / (4.0 * math.pi) * cut_function_der(half_light_radius / rs)/ cut_function(c)
                                                                       * half_light_radius / (rs ** 2.0))  ** 2.0
                                                                  )

    Walker_masses = []
    true_masses = []
    Walker_mass_errs = []
    true_mass_errs = [] 
    for pop in populations:
        print 'For pop ' + pop 
        new_star_data = ObservedGalaxyStarData([galaxy, pop], pop_selection_method)
        angular_half_light_radius = new_star_data.half_light_rad
        print 'angular_half_light_radius = ' + str(angular_half_light_radius)
        if fixed_half_light_radius is None: 
            half_light_radius = angular_half_light_radius * deg_to_rad * dist
        else:
            half_light_radius = fixed_half_light_radius
        print 'half_light_radius = ' + str(half_light_radius) 
        sig_sqr = new_star_data.sigSqr
        sig_sqrE = new_star_data.sigSqrE
        print 'sig_sqr = ' + str(sig_sqr) 
        walker_Mh = walker_half_mass_estimator (half_light_radius, sig_sqr)
        star_data_by_pop = star_data_by_pop + [new_star_data]

        Walker_masses = Walker_masses + [walker_half_mass_estimator(half_light_radius, sig_sqr)]
        Walker_mass_errs = Walker_mass_errs + [walker_half_mass_err (half_light_radius, sig_sqrE)]
        true_masses = true_masses + [true_half_mass_funct(half_light_radius)]
        true_mass_errs = true_mass_errs + [true_half_mass_funct_err(half_light_radius)]
    return [Walker_masses, Walker_mass_errs, true_masses, true_mass_errs]

        
                
