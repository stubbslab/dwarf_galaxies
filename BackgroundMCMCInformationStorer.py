import math

class BackgroundStorer:

    #NOTE: halo_center x and z are now degrees away from light center, rather than physical distance from center in parsecs!!!
    # halo_center y (which is basically never varied) is still distance in parsecs along the line of sight
    def __init__(self):
        self.r_half_light = 791.0
        self.M_star = 10.0 ** 7.39
        typical_sigsqr = 100.0
        Gyr_to_s = 365.2425 * 24 * 3600 * 10.0 ** 9.0
        self.param_ranges = {'M':[1.0 * self.M_star, 600.0 * self.M_star], 'rs':[0.1 * self.r_half_light, 10.0 * self.r_half_light], 'h_x_center':[-0.5,0.5], 'h_y_center':[-0.5,0.5], 'h_z_center':[-0.5,0.5], 'phi':[0.0, math.pi], 'theta':[0.0, 2.0 * math.pi],
                                'Rd':[0.1 * self.r_half_light, 10.0 * self.r_half_light], 'el':[0.1,10.0], 'd_x_center':[-0.5,0.5], 'd_y_center':[-0.5,0.5], 'd_z_center':[-0.5,0.5], 'eps':[0.0, 0.5], 'lam':[0.05, 0.15], 'a':[0.0, math.pi], 'b':[0.0, 2.0 * math.pi],
                                'sigsqr_RR':[1.0 / 4.0 * typical_sigsqr, 4 * typical_sigsqr], 'omega_phi':[-10.0 * 1.0 / Gyr_to_s , 10.0 * 1.0 / Gyr_to_s]}
        self.param_step_sizes = {'M':self.M_star * 6.0, 'rs':self.r_half_light * 0.1, 'h_x_center':0.005, 'h_y_center':0.005, 'h_z_center':0.005,
                                 'Rd':self.r_half_light * 0.1, 'el':0.1, 'd_x_center':0.005, 'd_y_center':0.005, 'd_z_center':0.005, 'eps':0.01, 'lam':0.01,
                                 'sigsqr_RR':4 * typical_sigsqr / 100.0, 'omega_phi':10.0 * 1.0 / Gyr_to_s / 100.0 * 2}
