import math
import numpy as np
from DwarfGalDataArchive import DwarfGalDataArchive
from ObservedGalaxyStarData import ObservedGalaxyStarData
from CandidateGalaxyStarData import CandidateGalaxyStarData
from cantrips import getAreaMesh
import scipy.interpolate as interpolate

class GalaxyMask:

    def inObservation(self,R,z,observation):
        #print 'test'
        if observation['shape'] == 'circle':
            center = observation['center']
            size = observation['size']
            return np.sqrt((R - center[0] ) ** 2 + (z - center[1]) ** 2) <= size / 2.0
        elif observation['shape'] == 'square':
            center = observation['center']
            size = observation['size']
            return (center[0] - size / 2.0 <=R) * (center[0] + size / 2.0 >= R) * (center[1] - size / 2.0 <=z) * (center[1] + size / 2.0 >=z)

    #Should determine probability of not being observed by reading in positions of all CANDIDATES
    # and comparing to positions of all stars with observed velocities stars
    def getProbOfMeasuringVelMask(self, xmesh, ymesh, log, galaxy_observation, masking_fields = []):
        kernel_size = 2.0 * 1.0 / 60.0  #Guassian smoothing width, in units of read-in coordinates (degrees by default)
        #print ('[xmesh, ymesh, self.candidate_star_pos] = ' + str([xmesh, ymesh, self.candidate_star_pos]))
        number_of_stars_at_position_funct = lambda xs, ys, star_positions: sum( [np.exp( -1.0 * (( (xs - star_coords[0]) ** 2.0 + (ys - star_coords[1]) ** 2.0))/ (2.0 * kernel_size ** 2.0) ) for star_coords in star_positions] )
        number_of_candidates_mesh = number_of_stars_at_position_funct (xmesh, ymesh, self.candidate_star_pos)
        number_of_observed_targets_mesh = np.zeros(np.shape( number_of_candidates_mesh ))
        for observation in log:
            #n_observed_with_vel = observation['n_vel']
            #n_candidates_in_field = observation['n_cands']
            if len(masking_fields) == 0 or (set(observation['name']).intersection(masking_fields) > 0):
                meas_stars_in_field = self.vel_measured_star_pos_by_observation[observation['names'][0]]
                number_of_observed_targets_mesh = number_of_observed_targets_mesh + number_of_stars_at_position_funct (xmesh, ymesh, meas_stars_in_field)

        prob_of_being_observed = number_of_observed_targets_mesh / number_of_candidates_mesh
        prob_too_high = (prob_of_being_observed > 1)
        prob_of_being_observed = prob_of_being_observed * (1 - prob_too_high) + prob_too_high * 1.0
        valid_prob = prob_of_being_observed > 0.0
        prob_of_being_observed = prob_of_being_observed * valid_prob
        prob_of_being_observed = np.nan_to_num(prob_of_being_observed)
        return prob_of_being_observed


    def getObservedForSufficientTimeMask(self, xmesh, ymesh, log, masking_fields=[]):
        default_min_time = 1.0
        obs_for_sufficient_time = np.zeros(np.shape(xmesh))
        for observation in log:
            time_in_field = observation['time']
            min_time_req_in_field = observation['min_t']
            if len(masking_fields) == 0 or (set(observation['name']).intersection(masking_fields) > 0):
                obs_for_sufficient_time_in_fields = (time_in_field >= min_time_req_in_field) * self.inObservation(xmesh, ymesh, observation)
                obs_for_sufficient_time = np.maximum(obs_for_sufficient_time, obs_for_sufficient_time_in_fields)
        return obs_for_sufficient_time


    #Provide a mask for likelihood of having stars overlap.  Note this is based on a preexisting profile.
    # Masking fields have no impact on this particular calculation.  Computations are done in sqr degrees.
    # Ideally, we would do this with each of the populations in the dSph.  But I have not been able to
    # come up with an efficient way of inverting the algorithm when there are more than 1 unknown variables (
    # (number of stars in each of the populations).
    def getCrowdingMask(self, xmesh, ymesh, n_cand_stars, frac_stars_in_gal = 1.0, star_area_sqrDegrees = (1.0 / (60.0 * 60.0)) ** (2.0), probability_profile = None):

        rs = probability_profile.rs
        dist = probability_profile.dist
        if probability_profile is None:
            print ('Did not get a valid probability profile.  Returning 1. ')
            return np.zeros(np.shape(xmesh)) + 1.0

        n_detected_member_stars = n_cand_stars * frac_stars_in_gal

        n_member_stars = n_detected_member_stars - 1
        print ('n_member_stars starts at ' + str(n_member_stars) )

        area_mesh = getAreaMesh([xmesh, ymesh])

        #print 'area_mesh = ' + str(area_mesh)
        #print '1.0 / (180.0 / math.pi * 60.0) ** 2.0 * probability_profile.onSkyInterpolator((xmesh * math.pi / (180.0) * dist / rs, ymesh * math.pi / (180.0) * dist / rs )) = '
        #print 1.0 / (180.0 / math.pi * 60.0) ** 2.0 * probability_profile.onSkyInterpolator((xmesh * math.pi / (180.0) * dist / rs, ymesh * math.pi / (180.0) * dist / rs ))

        print ('This should be 1: ')
        #print 'np.sum(area_mesh * (math.pi / 180.0) ** 2.0 * probability_profile.onSkyInterpolator((xmesh * math.pi / (180.0) * dist / rs, ymesh * math.pi / (180.0) * dist / rs )) * 1.0 / (180.0 / math.pi * 60.0) ** 2.0 ) = '
        print (np.sum(area_mesh * (math.pi / 180.0) ** 2.0 * probability_profile.onSkyInterpolator((xmesh * math.pi / (180.0) * dist / rs, ymesh * math.pi / (180.0) * dist / rs ))))

        #Probability profile must be evenly spaced so that integral can be done by simple sum
        #print 'n_detected_member_stars = ' + str(n_detected_member_stars)
        max_number_steps = 1000
        arg_accuracy_threshold = 10 ** (-7.0)
        n_steps = 0
        #print 'n_detected_member_stars = ' + str(n_detected_member_stars)
        detection_steps = [100,10,1]
        previous_steps = 0
        for steps in detection_steps:
            predicted_n_detected_member_stars = -1
            n_member_stars = n_member_stars - previous_steps
            while (predicted_n_detected_member_stars < n_detected_member_stars and n_steps <= max_number_steps):
                n_steps = n_steps + 1
                n_member_stars = n_member_stars + steps
                #Note that interpolator takes xmesh, ymesh in units of scale radii,
                # so I must convert xmesh and ymesh into scale radii (degrees -> radians -> scale radii

                arg = (4.0 * star_area_sqrDegrees * (math.pi / 180.0) ** 2.0 * n_member_stars
                       * probability_profile.onSkyInterpolator((xmesh * math.pi / 180.0 * dist / rs, ymesh * math.pi / 180.0 * dist / rs )) )

                arg [arg < arg_accuracy_threshold] = arg_accuracy_threshold

                predicted_n_detected_member_stars = (1.0 / (4.0 * star_area_sqrDegrees) * np.sum(area_mesh * (1.0 - np.exp(-arg))))
                #print 'n_member_stars = ' + str(n_member_stars)
                #print 'n_member_stars = ' + str(n_member_stars) + ' leads to predicted_n_detected_member_stars = ' + str(predicted_n_detected_member_stars)
            previous_steps = steps

        if n_steps > max_number_steps:
            print ('Exceeded number of steps. ' )

        print ('Settling on n_member_stars = ' + str(n_member_stars))
        arg = (4.0 * star_area_sqrDegrees * n_member_stars * 1.0 / (180.0 / math.pi) ** 2.0
                             * probability_profile.onSkyInterpolator((xmesh * math.pi / (180.0) * dist / rs, ymesh * math.pi / (180.0) * dist / rs )))
        print ('np.min(arg) = ' + str(np.min(arg)))
        arg [arg < arg_accuracy_threshold] = arg_accuracy_threshold
        mask = (1.0 - np.exp(-arg)) / arg
        print ('(1.0 - np.exp(-np.min(arg)) ) / (np.min(arg)) = ' + str((1.0 - np.exp(-np.min(arg)) ) / (np.min(arg))) )
        print ('np.min(mask.flatten()) = ' + str(np.min(mask.flatten())) )
        print ('np.max(mask.flatten() ) = ' + str(np.max(mask.flatten())) )
        mask_max_index = np.argmax(mask.flatten())
        print ('mask_max_index = ' + str(mask_max_index) )
        x_at_mask_max_index = arg.flatten()[mask_max_index]
        print ('x that leads to max value is: ' + str(x_at_mask_max_index) )
        print ('And verifying, here is the value directly computed ' + str((1.0 - np.exp(-x_at_mask_max_index)) / x_at_mask_max_index) )

        return mask


    def __init__(self, xmesh, ymesh, galaxy, mask_types = ['observed'], masking_fields = [], vel_meas_stars = None, candidate_stars = None, probability_profile = None, units_of_coord_meshes = 'degrees'):
        #This program does all operations in degrees.  So we need to immediately correct the meshes
        if units_of_coord_meshes in ['arcmin','arcmins']:
            xmesh = xmesh * 1.0 / 60.0
            ymesh = ymesh * 1.0 / 60.0
        if units_of_coord_meshes in ['arcsec','arcsecs']:
            xmesh = xmesh * 1.0 / (60.0) ** 2.0
            ymesh = ymesh * 1.0 / (60.0) ** 2.0
        elif units_of_coord_meshes in ['deg','degrees','degree','degs']:
            xmesh = xmesh
            ymesh = ymesh
        elif units_of_coord_meshes in ['rad', 'radian','rads','radians']:
            xmesh = xmesh * (180.0 / math.pi )
            ymesh = xmesh * (180.0 / math.pi)

        dwarf_arch = DwarfGalDataArchive()
        mask = np.zeros(np.shape(xmesh)) + 1.0
        extra_masks = []
        log = dwarf_arch.obsLog[galaxy]

        if candidate_stars is None:
            candidate_stars = CandidateGalaxyStarData(galaxy)

        if vel_meas_stars is None:
            vel_meas_stars = ObservedGalaxyStarData([galaxy, 'dummy_var'], 'none')

        if 'observed' in mask_types:
           extra_mask = self.getObservedForSufficientTimeMask(xmesh, ymesh, log, masking_fields = masking_fields)
           extra_masks = extra_masks + [extra_mask]

        if 'n_vel_meas' in mask_types:
            self.vel_measured_star_pos_by_observation = {}

            for observation in log:
                self.vel_measured_star_pos_by_observation[observation['names'][0]] = []

            for i in range(len(vel_meas_stars.Field)):
                field = vel_meas_stars.Field[i]
                found_match = 0
                for observation in log:
                    if field in observation['names']:
                        found_match = 1
                        self.vel_measured_star_pos_by_observation[observation['names'][0]] = self.vel_measured_star_pos_by_observation[observation['names'][0]] + [[vel_meas_stars.corrRa[i], vel_meas_stars.corrDec[i]]]
                        break
                if not(found_match):
                    print ('Star with index: ' + str(i) + ' has field: ' + str(field) + ', which cannot be located in the log.  ' )


           #Projected positions of stars on the sky.

            self.candidate_star_pos = [ [candidate_stars.xs_degrees[i], candidate_stars.ys_degrees[i]] for i in range(len(candidate_stars.xs_degrees)) ]

            extra_mask = self.getProbOfMeasuringVelMask(xmesh, ymesh, log, galaxy, masking_fields)
            #extra_mask = (n_obs.transpose() >= max_pos_n_obs)
            extra_masks = extra_masks + [extra_mask]

        if 'crowding' in mask_types:

            n_cand_stars = len(candidate_stars.xs_degrees)
            frac_stars_in_gal = vel_meas_stars.membership_fraction

            extra_mask = self.getCrowdingMask(xmesh, ymesh, n_cand_stars, frac_stars_in_gal = frac_stars_in_gal, probability_profile = probability_profile)

            extra_masks = extra_masks + [extra_mask]

        if 'stellar_crowding' in mask_types:
            print ('I cannot yet correct for stellar crowding.  That is coming soon! ' )
            extra_mask = np.zeros(np.shape(xmesh)) + 1.0


        for extra_mask in extra_masks:
            mask = mask * extra_mask

        self.final_mask = mask
        self.final_mask_interp = interpolate.RegularGridInterpolator((xmesh[0], ymesh[:,0]), np.transpose(self.final_mask), method='linear', bounds_error=False, fill_value=0.0)
