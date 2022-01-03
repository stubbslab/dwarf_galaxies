import math 

class BackgroundStorer:

    def __init__(self):
        self.r_half_light = 791.0
        self.M_star = 10.0 ** 7.39 
        self.param_ranges = {'M':[1.0 * self.M_star, 400.0 * self.M_star], 'rs':[0.1 * self.r_half_light, 10.0 * self.r_half_light], 'h_x_center':[-self.r_half_light,self.r_half_light], 'h_z_center':[-self.r_half_light,self.r_half_light], 'phi':[0.0, math.pi], 'theta':[0.0, 2.0 * math.pi], 
                                'Rd':[0.1 * self.r_half_light, 10.0 * self.r_half_light], 'el':[0.1,10.0], 'd_x_center':[-self.r_half_light,self.r_half_light], 'd_z_center':[-self.r_half_light,self.r_half_light], 'eps':[0.0, 0.5], 'lam':[0.05, 0.15], 'a':[0.0, math.pi], 'b':[0.0, 2.0 * math.pi], }
        self.param_step_sizes = {'M':self.M_star * 1.0, 'rs':self.r_half_light * 10.0 * 0.01, 'h_x_center':self.r_half_light/100.0, 'h_z_center':self.r_half_light/100.0,
                                 'Rd':self.r_half_light * 10.0 * 0.01, 'el':10.0 * 0.01, 'd_x_center':self.r_half_light/100.0, 'd_z_center':self.r_half_light/100.0, 'eps':0.01, 'lam':0.01}
