# A hard coded data class that basically includes all of the best parameters describing the galaxies and stellar subpopulations I want to study.
# This file does not read in any data in files, though it does include the name of the files that include data to be read in
# This file is all about the observed properties of galaxies.  Parameters involved in simulations should be accounted for elsewhere.
#For reference, files are defind as:
# McConnachie: 'The Observed Properties of Dwarf Galaxies in and Around the Local Group'
# Walker1: 'Stellar Velocities in the Carina, Fornax, Sculptor, and Sextans dSph Galaxies: Data from the Magellan/MMFS Survey
# Amarisco: A Troublesome Past: Chemodynamics of the Fornax DwarfSpheroidal
# Walker2: A Method for Measuring (Slopes Of) the Mass Profiles
# Lokas1: The mass and velocity anisotropy of the Carina, Fornax, Sculptor and Sextans dwarf spheroidal galaxies
# Piatek1: Proper Motions of Dwarf Spheroidal  Galaxies from Hubble Space Telescope Imaging. V: Final Measurement for Fornax 
# Dinescu1: https://arxiv.org/pdf/astro-ph/0405260.pdf (old Fornax PM) 

import numpy as np
import math
from AstronomicalParameterArchive import AstronomicalParameterArchive
from ComputationalArchive import ComputationalArchive

class DwarfGalDataArchive:

    def inObservation(self,R,z,observation):
        #print 'test'
        if observation['shape'] == 'circle':
            #print 'shape = circle'
            center = observation['center']
            #print 'center = ' +str(center)
            size = observation['size']
            #print 'size = ' + str(size)
            #print np.sqrt((R - center[0] ) ** 2 + (z - center[1]) ** 2)
            return np.sqrt((R - center[0] ) ** 2 + (z - center[1]) ** 2) <= size / 2.0
        elif observation['shape'] == 'square':
            center = observation['center']
            size = observation['size']
            return (center[0] - size / 2.0 <=R) * (center[0] + size / 2.0 >= R) * (center[1] - size / 2.0 <=z) * (center[1] + size / 2.0 >=z)

    def getTotalTime(self,Rmesh,zmesh,log,masking_fields = []):
        obs_time = np.zeros(np.shape(Rmesh))
        for observation in log:
            time_in_field = observation['time']
            if len(masking_fields) == 0 or observation['name'] in masking_fields: 
                obs_time = obs_time +  time_in_field * self.inObservation(Rmesh, zmesh, observation)
        return obs_time

    def getMinObsTime(self,Rmesh,zmesh,log,masking_fields = []):
        default_min_time = 1.0 
        min_time = np.zeros(np.shape(Rmesh)) + default_min_time 
        for observation in log:
            min_time_in_field = observation['min_t']
            if len(masking_fields) == 0 or observation['name'] in masking_fields: 
                min_time = np.minimum( min_time, min_time_in_field * self.inObservation(Rmesh, zmesh, observation) ) ###MUST BE FIXED => SHOULD BE INFINITIES WHEREVER NEW TIME IS NOT 
        return min_time

    def getObsForSufficientTimeMask(self, Rmesh, zmesh, log, masking_fields=[]):
        default_min_time = 1.0
        obs_for_sufficient_time = np.zeros(np.shape(Rmesh))
        for observation in log:
            time_in_field = observation['time']
            min_time_req_in_field = observation['min_t']
            if len(masking_fields) == 0 or observation['name'] in masking_fields:
                obs_for_sufficient_time_in_fields = (time_in_field >= min_time_req_in_field) * self.inObservation(Rmesh, zmesh, observation) 
                obs_for_sufficient_time = np.maximum(obs_for_sufficient_time, obs_for_sufficient_time_in_fields)
        return obs_for_sufficient_time 

    #Should determine probability of not being observed by reading in positions of all CANDIDATES
    # and comparing to positions of all stars with observed velocities stars
    def getProbOfMeasuringVelMask(self, Rmesh, zmesh, log, galaxy_observation, masking_fields = []):
        positions_of_stellar_candidates = getCandidateList(galaxy_observation) #Projected positions of stars on the sky.
        kernal_size = 1.0
        number_of_stars_at_position_funct = lambda xs, ys, star_positions: sum( [np.exp( -1.0 * (((xs - star_coords[0]) ** 2.0 + (ys - star_coords[1])) ** 2.0)/ k ) for star_coords in star_positions] )
        number_of_candidates_mesh = number_of_candidates_at_position_funct (Rmesh, zmesh, positions_of_stellar_candidates)
        number_of_observed_targets_mesh = np.zeros(np.shape( number_of_candidates_mesh )) 
        for observation in log:
            #n_observed_with_vel = observation['n_vel']
            #n_candidates_in_field = observation['n_cands']
            if len(masking_fields) == 0 or observation['name'] in masking_fields:
                field_mask = self.inObservation(Rmesh, zmesh, observation)
                stars_in_observation = observation.stars_with_vels
                number_of_observed_targets_mesh = number_of_observed_targets_mesh + number_of_candidates_at_position_funct (Rmesh, zmesh, [[star[RA], star[Dec]] for star in stars_in_observation])
        prob_of_being_observed = number_of_observed_targets_mesh / number_of_candidates_mesh 
        return prob_of_being_observed

    def getProbOfMissingTrueStarDueToCrowdingMask(self, Rmesh, zmesh, log, stellar_surface_prof, masking_fields = []):
        default_stellar_width = 2
        print ('This function will give you a probability mask based on the chances of missing a star due to it being blocked by surrounding stars. ') 
        print ('The function is not yet completed. ') 
        return 1
            
    #returns a probability mask based on where observations were taken for a given galaxy.
    #R and z must be given in degrees from galaxy center
    def getObservationMask(self,population,R,z,masking_fields = [], mask_types = ['observed'] ):
        #print 'R in degrees ='
        #print R
        #print 'z in degrees = '
        #print z
        Rmesh,zmesh = np.meshgrid(R,z)
        obsLog = self.obsLog[population[0]]
        
        n_obs = self.getNObservedPerFied(Rmesh, zmesh, obsLog, masking_fields = masking_fields)
        #We must return the transpose of the array, since the probability array is also transposed (in order to make the likelihood interpolator function as desired)
        mask = np.zeros(np.shape(Rmesh)) + 1.0
        extra_masks = []

        if 'observed' in mask_types:
            extra_mask = self.getObservedForSufficientTimeMask(Rmesh,zmesh,obsLog, masking_fields = masking_fields)
            #total_time = self.getTotalTime(Rmesh,zmesh,obsLog,masking_fields)
            #min_obs_time = self.getMinObsTime(Rmesh, zmesh, obsLog, masking_fields)
            #extra_mask = (total_time.transpose() >= min_obs_time.transpose())
            extra_masks = extra_masks + [extra_mask]
            
        if 'n_vel_meas' in mask_types:
            #n_obs = self.getNObservedPerFied(Rmesh, zmesh, obsLog, masking_fields)
            #max_pos_n_obs = self.getMaxPosibleNObservedPerFied(Rmesh, zmesh, obsLog, masking_fields)
            extra_mask = getProbOfMeasuringVelMask(Rmesh, zmesh, obsLog, population[0], masking_fields)
            #extra_mask = (n_obs.transpose() >= max_pos_n_obs)
            extra_masks = extra_masks + [extra_mask]

        if 'stellar_crowding' in mask_types:
            print ('I cannot yet correct for stellar crowding.  That is coming soon! ') 
            extra_mask = np.zeros(np.shape(Rmesh)) + 1.0 
        

        for extra_mask in extra_masks:
            mask = mask * extra_mask

        return mask 

    #Returns observation bounds specified unit
    def getObservationBounds(self,population,return_unit = 'degree'):
        obsLog = self.obsLog[population[0]]
        minRA = 0.0
        maxRA = 0.0
        minDec = 0.0
        maxDec = 0.0
        for obs in obsLog:
            center = obs['center']
            size = obs['size']
            shape = obs['shape']
            if shape == 'circle' or shape == 'square':
                if center[0] - size / 2.0 < minRA: minRA = center[0] - size / 2.0
                if center[0] + size / 2.0 > maxRA: maxRA = center[0] + size / 2.0
                if center[1] - size / 2.0 < minDec: minDec = center[1] - size / 2.0
                if center[1] + size / 2.0 > maxDec: maxDec = center[1] + size / 2.0
        if return_unit == 'arcmin':
            return[[minRA*60.0,maxRA*60.0],[minDec*60.0,maxDec*60.0]]
        elif return_unit == 'arcsec':
            return[[minRA*3600.0,maxRA*3600.0],[minDec*3600.0,maxDec*3600.0]]
        elif return_unit == 'degree':
            return[[minRA,maxRA],[minDec,maxDec]]
        else:
            return[[minRA,maxRA],[minDec,maxDec]]

    
    
    def getPopulations(self,galaxy):
        return self.distinct_populations[galaxy]

    def getFile(self,population):
        return self.gal_data_files[population[0]]

    def getCandidateFile(self, galaxy):
        return self.gal_cand_data_files[galaxy]
    
    def getDistanceFromSun(self,population):
        return self.sun_distances[population[0]]
    
    def getDistanceFromSunErr(self,population):
        return self.sun_distance_errs[population[0]] 
    
    def getTotalMass(self,population):
        return self.total_masses[population[0]] 
        
    def getHorizontalBranchVMagnitude(self,population): #V_HB = {"carina" : 20.9, "fornax" : 21.3,"sculptor" : 20.1,"sextans" : 20.25}
        return self.V_horizontal_branch_magnitude[population[0]]
        
    def getSigMgScaleCorr(self,population): #0.079
        return self.sig_mg_scale_corr[population[0]]
        
    def getRaDec(self,population):
        return self.ra_dec_tuples[population[0]]
    
    def getPopParameters(self,population):
        return self.pop_params[population[0]]
    
    def getMetallicityCuts(self,population):
        return self.metallicity_cuts[population[0]]
        
    def getRotationalVelocity(self,population):     
        return self.rotational_velocity[population[0]][population[1]] 
    
    def getVelocityDispersion(self,population):
        return self.velocity_dispersion[population[0]][population[1]]
    
    def getHeliocentricProperMotion(self,population):
        return self.heliocentric_proper_motion[population[0]]

    def getHeliocentricProperMotionError(self,population):
        return self.heliocentric_proper_motion_err[population[0]] 

    def __init__(self, pop_selection_method = ''):
        astro_archive = AstronomicalParameterArchive()
        computational_params_archive = ComputationalArchive()
        deg_to_rad = astro_archive.getDegToRad()
        # number of distinct populations identified in each of the galaxies.  Amarisco for fornax, Walker 1 for others.
        if pop_selection_method in ['none', None, 'None']:
            self.distinct_populations = {"carina":['MP'],"fornax":['MP'],"sculptor":['MP'],"sextans":['MP']}
        else: 
            self.distinct_populations = {"carina":['MP','MR'],"fornax":['MP','IM','MR'],"sculptor":['MP','MR'],"sextans":['MP','MR'],
                                         'boylan1':['all'], 'boylan2':['all'], 'boylan3':['all']}
        #McConnachie: distance of center of galaxy from Sun, in pc
        self.sun_distances = {"carina":105000.0,"fornax":147000.0,"sculptor":86000.0,"sextans":86000.0}
        self.sun_distance_errs = {"carina":0.0,"fornax":0.0,"sculptor":0.0,"sextans":0.0}

        #For the moment, I have pegged some Boylan artificial data to Fornax as reference
        # I then applied some rough rescaling to the distances so that all galaxies would have
        # a similar angular size (~2 deg diameter) on the sky.  
        self.sun_distances['boylan0'] = self.sun_distances['fornax']
        self.sun_distances['boylan1'] = self.sun_distances['fornax'] * 5.0
        self.sun_distances['boylan2'] = self.sun_distances['fornax'] / 5.0
        
        #Lokas1: total galactic mass, in M_sun
        self.total_masses = {"carina":2.3*10**7,"fornax":15.7*10**7,"sculptor":3.1*10**7,"sextans":4.0*10**7}
        #For the moment, I have pegged some Boylan artificial data to Fornax for concreteness of reference
        self.total_masses['boylan0'] = self.total_masses['fornax']
        self.total_masses['boylan1'] = self.total_masses['fornax']
        self.total_masses['boylan2'] = self.total_masses['fornax']

        #Walker 2: typical magnitude of the horizontal branch stars in each galaxy
        self.V_horizontal_branch_magnitude = {"carina" : 20.9, "fornax" : 21.3,"sculptor" : 20.1,"sextans" : 20.25} 
        #Walker 1: correction scale for metallicities 
        self.sig_mg_scale_corr = {"carina" : 0.079, "fornax" : 0.079,"sculptor" : 0.079,"sextans" : 0.079}
        # rotational_velocity about central axis of symmetry , in km/s 
        self.rotational_velocity = {"carina":{"MP": 0.0,"MR": 0.0}, "fornax":{"MP": 0.0,"IM": 0.0, "MR": 0.0},
                                    "sculptor":{"MP": 0.0 ,"MR": 0.0}, "sextans":{"MP": 0.0,"MR": 0.0}}
        
        #For the moment, I have pegged some Boylan artificial data no rotation 
        self.rotational_velocity['boylan0'] = {"all":0.0}
        self.rotational_velocity['boylan1'] = {"all":0.0}
        self.rotational_velocity['boylan2'] = {"all":0.0}

        
        #Amarisco for Fornax, Walker 2 for others: velocity dispersion in km/s
        #I should be able to measure these MYSELF, given the populations that I break the parameters into.
        # But to do that, I must read in data, which is done in another file. 
        self.velocity_dispersion = {"carina":{"MP": 41.70,"MR": 28840.3}, "fornax":{"MP": 14.0**2,"IM": 11.5**2, "MR": 8.5**2},
                                    "sculptor":{"MP": 134.90 ,"MR": 41.70}, "sextans":{"MP": 100.0,"MR": 100.0}}

        self.mask_limit_time = 1.0
        self.mask_limit_n_obs = {"carina" : 256, "fornax" : 256,"sculptor" : 256,"sextans" : 256}
        
        #Hard coded on sky coordinates of galactic centers , in degrees
        #Walker2 for Carina, Fornax, and Sculptor
        #McConnachi for Sextans 
        #Note: another possible set of options available at:
        # http://link.springer.com.ezp-prod1.hul.harvard.edu/article/10.1007/s001590050019
        car_ra= (float(6.0) + float(41.0)/60.0 + float(37.0)/3600.)/24.0*360.0
        car_dec=(float(-50.0) + float(-58.0)/60.0 + float(0.0)/3600.)
        for_ra=(float(2.0) + float(39.0)/60.0 + float(59.0)/3600.)/24.0*360.0
        for_dec=(float(-34.0) + float(-27.0)/60.0 + float(0.0)/3600.)
        scu_ra=(float(1.0) + float(0.0)/60.0 + float(9.0)/3600.)/24.0*360.0
        scu_dec=(float(-33.0) + float(-42.0)/60.0 + float(-30.0)/3600.)
        sex_ra=(float(10.0) + float(13.0)/60.0 + float(3.0)/3600.)/24.0*360.0
        sex_dec=(float(-1.0) + float(-36.0)/60.0 + float(-53.0)/3600.)
        self.ra_dec_tuples = {"carina" : (car_ra,car_dec),"fornax" : (for_ra,for_dec),"sculptor" : (scu_ra,scu_dec),"sextans" : (sex_ra,sex_dec)}

        #Piatek1 for Fornax, not looked up for others
        #GAIA collaboration at: https://arxiv.org/pdf/1804.09381.pdf for all 
        # heliocentric proper motion of dSph center in mas/yr
        # note these are the (mu_alpha,mu_delta) coordinates, which is change is rate of change in right ascension and declination

        self.heliocentric_proper_motion = {"carina":[0.495, 0.143], "fornax":[0.376, -0.413], "sculptor":[0.082, -0.131 ], "sextans":[ -0.496, 0.077]}
        self.heliocentric_proper_motion_err = {"carina":[0.015, 0.014], "fornax":[0.003, 0.003], "sculptor":[0.005, 0.004], "sextans":[0.0025, 0.020]}
        correct_PM_for_spherical_coords = 0

        #print 'self.heliocentric_proper_motion = '
        #print self.heliocentric_proper_motion

        if correct_PM_for_spherical_coords: 
           #correct for fact that change in RA in coordinates is not true angular velocity
            for key in self.heliocentric_proper_motion:
                uncorr_pm = self.heliocentric_proper_motion[key]
                uncorr_pm_err = self.heliocentric_proper_motion_err[key]
                ra_dec = self.ra_dec_tuples[key] 
                self.heliocentric_proper_motion[key] = [uncorr_pm[0] * math.cos(ra_dec[1] * deg_to_rad), uncorr_pm[1]]
                self.heliocentric_proper_motion_err[key] = [uncorr_pm_err[0] * math.cos(ra_dec[1] * deg_to_rad), uncorr_pm_err[1]]
            
        #print 'self.heliocentric_proper_motion = '
        #print self.heliocentric_proper_motion
        
        #Amarisco for Fornax, Walker 2 for others: for each of MP, IM, MR populations, the best average SigMg, StdSigMg, Rh
        #self.pop_params = {"carina": [[0.30,0.09,263.03],[0.044,0.09,218.31],[0.44,0.09,218.31]],"fornax": [[0.27,0.07,888],[0.46,0.058,605],[0.53,0.11,480]],"sculptor":[[0.24,0.01,301.995],[0.24,0.01,301.995],[0.36,0.01,166.1]],"sextans":[[0.24,0.01,301.995],[0.24,0.01,301.995],[0.36,0.01,166.1]]}
        self.pop_params = {"carina": {'MP':{'meanSigMg':0.30,'meanSigMgE':0.09,'Rh':263.03},  'IM':{'meanSigMg':0.44,'meanSigMgE':0.09,'Rh':218.31},'MR':{'meanSigMg':0.44,'meanSigMgE':0.09,'Rh':218.31}},
                           "fornax": {'MP':{'meanSigMg':0.26,'meanSigMgE':0.05,'Rh':935.0,'f':0.27},     'IM':{'meanSigMg':0.453,'meanSigMgE':0.065,'Rh':610.0,'f':0.63},   'MR':{'meanSigMg':0.56,'meanSigMgE':0.10,'Rh':437.0,'f':0.1}},
                           "sculptor":{'MP':{'meanSigMg':0.24,'meanSigMgE':0.01,'Rh':301.995},'IM':{'meanSigMg':0.24,'meanSigMgE':0.01,'Rh':301.995},'MR':{'meanSigMg':0.36,'meanSigMgE':0.01,'Rh':166.1}},
                           "sextans":{'MP':{'meanSigMg':0.24,'meanSigMgE':0.01,'Rh':301.995}, 'IM':{'meanSigMg':0.24,'meanSigMgE':0.01,'Rh':301.995},'MR':{'meanSigMg':0.36,'meanSigMgE':0.01,'Rh':166.1}} }
        #the cuts used to distinguish the various populations for each galaxy.
        #Selected myself: Note that for those galaxies with only 2 populations, the lowest metallicity is just the lowest value of any star in the galaxy
        self.metallicity_cuts={"carina":[0.04,0.57],"fornax":[0.3483,0.5186],"sculptor":[0.3,0.5],"sextans":[-0.9,-0.6]}
        galaxy_data_directory = computational_params_archive.getDSphDataDir() 
        file_name_prefix = "stellarVelocitiesintheCarinaFornaxSculptorandSextansdSphGalaxiesDatafromtheMagellanMMFSSurvey_"
        cand_file_name_suffix = "_rgbcands.res"
        carina_file_name = galaxy_data_directory + file_name_prefix + "CarinaDataFile_justData.txt"  
        fornax_file_name = galaxy_data_directory + file_name_prefix + "FornaxDataFile_justData.txt"  
        sculptor_file_name = galaxy_data_directory + file_name_prefix + "SculptorDataFile_justData.txt"  
        sextans_file_name = galaxy_data_directory + file_name_prefix + "SextansDataFile_justData.txt"
        carina_cand_file_name = galaxy_data_directory + 'car' + cand_file_name_suffix 
        fornax_cand_file_name = galaxy_data_directory + 'for' + cand_file_name_suffix 
        sculptor_cand_file_name = galaxy_data_directory + 'scl' + cand_file_name_suffix 
        sextans_cand_file_name = galaxy_data_directory + 'sex' + cand_file_name_suffix
        boylan_artificial_file_name0 = galaxy_data_directory + 'alexArtificialHalo1.txt'
        boylan_artificial_file_name1 = galaxy_data_directory + 'alexArtificialHalo2.txt'
        boylan_artificial_file_name2 = galaxy_data_directory + 'alexArtificialHalo3.txt'
        
        self.gal_data_files = {"carina":carina_file_name, "fornax":fornax_file_name, "sculptor":sculptor_file_name, "sextans":sextans_file_name,
                               "Boylan0":boylan_artificial_file_name0, "Boylan1":boylan_artificial_file_name1, "Boylan2":boylan_artificial_file_name2}
        self.gal_cand_data_files = {"carina":carina_cand_file_name, "fornax":fornax_cand_file_name, "sculptor":sculptor_cand_file_name, "sextans":sextans_cand_file_name}

        #Should create, for each galaxy, a list of libraries for observation.
        #The format for each observation should be:
        #observation = {'name':THE NAME OF THE OBSERVATION, FOR EACH OF CROSS CHECKING.  eg-Walker38
        #               'time':THE EXPOSURE TIME IN SECONDS,
        #               'shape': THE SHAPE OF THE DETECTION REGION (currently supports squares and circles),
        #               'center': THE CENTER OF THE OBSERVATION, corrected for the galaxy center coordinate.  True for circles and squares
        #               'size': THE SIZE THAT CHARACTERIZES THE }
        #Exotic shapes might have additional parameters.  Not every key necessarily needed for every observation
        #   For exapmle, a rectangular region needs an additional key that a circle and square won't need
        self.obsLog = {"carina"  :[{'names':['WalkerTest0'], 'time':1000.0, 'shape':'circle', 'center': [0.0,0.0], 'size' : 0.2, 'min_t':1.0, 'max_n_obs':256, 'n_targs':256, 'n_candidates':0}
                                  ],
                       "fornax"  :[{'names':['Walker15','15R','15B'] , 'time':23100.0, 'shape':'circle', 'center':[-20.71/60.0, -6.52/60.0 ], 'size':20.0/60.0, 'min_t':1.0, 'n_vel':191, 'n_targs':220, 'n_cands':220},
                                   #{'names':['Walker20','20R','20B']  , 'time':23100.0, 'shape':'circle', 'center' : [-7.5/60.0   , -20.0/60.0]   , 'size':0.04, 'n_cands':0},
                                   {'names':['Walker10','10R','10B'] , 'time':9600.0 , 'shape':'circle', 'center':[-12.72/60.0, 4.65/60.0  ], 'size':20.0/60.0, 'min_t':1.0, 'n_vel':201, 'n_targs':223, 'n_cands':223 },
                                   {'names':['Walker29','29R','29B'] , 'time':9600.0 , 'shape':'circle', 'center':[11.98/60.0 , 9.82/60.0  ], 'size':20.0/60.0, 'min_t':1.0, 'n_vel':177, 'n_targs':216, 'n_cands':216 },
                                   {'names':['Walker34','34R','34B'] , 'time':9600.0 , 'shape':'circle', 'center':[-11.26/60.0, -26.64/60.0], 'size':20.0/60.0, 'min_t':1.0, 'n_vel':147, 'n_targs':163, 'n_cands':163 },
                                   {'names':['Walker592','592R','592B'], 'time':7800.0 , 'shape':'circle', 'center':[1.04/60.0  , -20.64/60.0], 'size':20.0/60.0, 'min_t':1.0, 'n_vel':109, 'n_targs':127, 'n_cands':127 },
                                   {'names':['Walker4','4R','4B'], 'time':52800.0 , 'shape':'circle', 'center':[-2.78/60.0 , 2.39/60.0  ], 'size':20.0/60.0, 'min_t':1.0, 'n_vel':935, 'n_targs':1117, 'n_cands':1117},
                                   {'names':['Walker9','9R','9B'], 'time':8100.0 , 'shape':'circle', 'center':[-11.51/60.0, -11.39/60.0], 'size':20.0/60.0, 'min_t':1.0, 'n_vel':595, 'n_targs':835, 'n_cands':835 },
                                   #{'names':['Walker4_3','15R','15B'], 'time':6000.0 , 'shape':'circle', 'center':[-2.5/60.0 , 2.5/60.0  ], 'size':20.0/60.0, 'min_t':1.0, 'n_vel':218, 'n_targs':224, 'n_cands':893 },
                                   #{'names':['Walker9_2','15R','15B'], 'time':29700.0 , 'shape':'circle', 'center':[-11.0/60.0, -11.5/60.0], 'size':20.0/60.0, 'min_t':1.0, 'n_vel':109, 'n_targs':192, 'n_cands':638 },
                                   {'names':['Walker54','54R','54B'] , 'time':10800.0, 'shape':'circle', 'center':[-29.47/60.0, -23.39/60.0], 'size':20.0/60.0, 'min_t':1.0, 'n_vel':74 , 'n_targs':80 , 'n_cands':80  },
                                   {'names':['Walker81','81R','81B'] , 'time':8100.0 , 'shape':'circle', 'center':[24.38/60.0 , 19.16/60.0 ], 'size':20.0/60.0, 'min_t':1.0, 'n_vel':95 , 'n_targs':101, 'n_cands':101 },
                                   #{'names':['Walker9_3','15R','15B'], 'time':8100.0 , 'shape':'circle', 'center':[-11.0/60.0, -11.5/60.0], 'size':20.0/60.0, 'min_t':1.0, 'n_vel':205, 'n_targs':224, 'n_cands':446 },
                                   {'names':['Walker33','33R','33B'] , 'time':1000.0 , 'shape':'circle', 'center':[17.96/60.0 , -6.07/60.0 ], 'size':20.0/60.0, 'min_t':1.0, 'n_vel':88 , 'n_targs':105, 'n_cands':105 },
                                   {'names':['Walker31','31R','31B'] , 'time':1000.0 , 'shape':'circle', 'center':[-28.24/60.0, -5.19/60.0 ], 'size':20.0/60.0, 'min_t':1.0, 'n_vel':119, 'n_targs':129, 'n_cands':129 },
                                   #{'names':['Walker4_4','15R','15B'], 'time':8100.0 , 'shape':'circle', 'center':[-2.5/60.0 , 2.5/60.0  ], 'size':20.0/60.0, 'min_t':1.0, 'n_vel':216, 'n_targs':224, 'n_cands':669 },
                                   {'names':['Walker70','70R','70B'] , 'time':1000.0 , 'shape':'circle', 'center':[-11.81/60.0, -37.51/60.0], 'size':20.0/60.0, 'min_t':1.0, 'n_vel':60 , 'n_targs':72 , 'n_cands':72  },
                                   #{'names':['Walker4_5','15R','15B'], 'time':9000.0 , 'shape':'circle', 'center':[-2.5/60.0 , 2.5/60.0  ], 'size':20.0/60.0, 'min_t':1.0, 'n_vel':128, 'n_targs':223, 'n_cands':445 },
                                   #{'names':['Walker9_6','15R','15B'], 'time':5400.0 , 'shape':'circle', 'center':[-11.0/60.0, -11.5/60.0], 'size':20.0/60.0, 'min_t':1.0, 'n_vel':96 , 'n_targs':222, 'n_cands':222 },
                                   #{'names':['Walker4_6','15R','15B'], 'time':7200.0 , 'shape':'circle', 'center':[-2.5/60.0 , 2.5/60.0  ], 'size':20.0/60.0, 'min_t':1.0, 'n_vel':165, 'n_targs':222, 'n_cands':222 },
                                   {'names':['Walker47','47R','47B'] , 'time':8100.0 , 'shape':'circle', 'center':[6.37/60.0, 23.59/60.0 ], 'size':20.0/60.0, 'min_t':1.0, 'n_vel':178, 'n_targs':215, 'n_cands':215 },
                                   {'names':['Walker38','38R','38B'] , 'time':8100.0 , 'shape':'circle', 'center':[-9.42/60.0 , 21.97/60.0 ], 'size':20.0/60.0, 'min_t':1.0, 'n_vel':121, 'n_targs':143, 'n_cands':143 },
                                   {'names':['Walker41','41R','41B'] , 'time':7200.0 , 'shape':'circle', 'center':[-26.15/60.0, 15.42/60.0 ], 'size':20.0/60.0, 'min_t':1.0, 'n_vel':60 , 'n_targs':68 , 'n_cands':68  }
                                  ],
                       "sculptor":[{'names':['WalkerTest0'], 'time':1000.0, 'shape':'circle', 'center': [0.0,0.0], 'size':0.2, 'min_t':1.0, 'n_vel':256, 'n_targs':256, 'n_cands':0}
                                  ],
                       "sextans" :[{'names':['WalkerTest0'], 'time':1000.0, 'shape':'circle', 'center': [0.0,0.0], 'size':0.2, 'min_t':1.0, 'n_vel':256, 'n_targs':256, 'n_cands':0}
                                  ]
                      }

                
                                   
