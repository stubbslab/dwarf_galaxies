import math

class BackgroundStorer:

    #NOTE: halo_center x and z are now degrees away from light center, rather than physical distance from center in parsecs!!!
    # halo_center y (which is basically never varied) is still distance in parsecs along the line of sight
    def __init__(self):
        self.r_half_light = 791.0
        self.M_star = 10.0 ** 7.39
        self.param_ranges = {'M':[1.0 * self.M_star, 600.0 * self.M_star], 'rs':[0.1 * self.r_half_light, 10.0 * self.r_half_light], 'h_x_center':[-0.5,0.5], 'h_y_center':[-0.5,0.5], 'h_z_center':[-0.5,0.5],
                               'sigsqr_rr_0':[25.0, 400.0], 'sigsqr_rr_inf_to_0_rat':[0.1, 10.0], 'r_sigsqr_rr0':[0.1 * self.r_half_light, 10.0 * self.r_half_light], 'alpha_sigsqr_rr':[0.01, 5.0],
                               'gamma_for_beta_inf':[-0.99, 0.99], 'r_beta0':[0.1 * self.r_half_light, 10.0 * self.r_half_light], 'alpha_beta':[0.01, 5.0], }
        self.param_step_sizes = {'M':self.M_star * 6.0, 'rs':self.r_half_light * 0.1, 'h_x_center':0.005, 'h_y_center':0.005, 'h_z_center':0.005,
                                 'sigsqr_rr_0':4.0, 'sigsqr_rr_inf_to_0_rat':self.r_half_light * 0.1, 'r_sigsqr_rr0':0.1, 'alpha_sigsqr_rr':0.06,
                                 'gamma_for_beta_inf':0.02, 'r_beta0':self.r_half_light * 0.1, 'alpha_beta':0.06, }
