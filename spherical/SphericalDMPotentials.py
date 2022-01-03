import numpy as np
import scipy.interpolate as interpolate

#Returns the potential function, scaled by M_halo / (4 pi r_s ^ 3)
def getPotentialDerivativeFunction(pot_type, c):
    if pot_type in ['nfw','NFW','Nfw']:
        scaling = np.log(1.0 + c) - c / (1.0 + c)
        raw_pot_deriv_funct = lambda xs:  np.where(xs  > 0.0, (np.log(1.0 + xs) - xs / (1.0 + xs)) / (xs ** 2.0), 0.5)

    elif pot_type in ['cored','acored','Cored','Acored', 'CORED', 'ACORED']:
        scaling = np.log(1.0 + c) + (4.0 * c + 3.0) / (2.0 * (1.0 + c) ** 2.0 )
        raw_pot_deriv_funct = lambda xs: np.where(xs  > 0.0, np.divide(1.0, (xs ** 2.0), out=np.zeros_like(xs), where=xs!=0) * (np.log(1.0 + xs) - ((1.5 * xs + 1.0) * xs) / ((xs + 1.0) ** 2.0)), 0.0)

    elif pot_type in ['burkert','burk','Burkert', 'Burk','BURKERT','BURK']:
        scaling = 0.5 * np.log(1.0 + c) + 0.25 * np.log(c ** 2.0 + 1.0) - 0.5 * np.arctan(c)
        #int_mass_funct = lambda xs: 1.0 / xs * (np.log((xs ** 2.0 + 1.0)/ ((xs + 1.0) ** 2.0)) + xs * np.log(xs ** 4.0 / ((xs ** 2.0 + 1.0) * (xs + 1.0) ** 2.0)) - 2.0 * np.arctan(xs) + 2.0 - 2.0 * special.hyp2f1(-0.5, 1.0, 0.5, -xs ** 2.0))
        raw_pot_deriv_funct = lambda xs: np.where(xs > 0.0, 0.5 * (0.5 * np.log(xs ** 2.0 + 1.0) + np.log(xs + 1.0) - np.arctan(xs)) / (xs ** 2.0), 0.0)

    else:
        print ('Gravitationl potential type: "' + pot_type + '" not recognized.  Returning point mass potential. ')
        scaling = np.log(c)
        raw_pot_funct = lambda xs: 1.0 / rs

    potential_funct = lambda rs, c=c: np.where(rs <= c, raw_pot_deriv_funct(rs) / scaling, -np.divide(1.0, rs ** 2.0, out=np.zeros_like(rs), where=rs!=0) )
    return potential_funct

#Returns the potential function, scaled by M_halo / (4 pi r_s ^ 3)
def getPotentialFunction(pot_type, c):
    if pot_type in ['nfw','NFW','Nfw']:
        scaling = np.log(1.0 + c) - c / (1.0 + c)
        raw_pot_funct = lambda xs:  np.where(xs  > 0.0, -(np.log(1.0 + xs)) / xs, -1.0)

    elif pot_type in ['cored','acored','Cored','Acored', 'CORED', 'ACORED']:
        scaling = np.log(1.0 + c) + (4.0 * c + 3.0) / (2.0 * (1.0 + c) ** 2.0 )
        raw_pot_funct = lambda xs: np.where(xs  > 0.0, 1.0 / (2.0 * (1.0 + xs)) - np.log(1.0 + xs) / xs, -0.5)

    elif pot_type in ['burkert','burk','Burkert', 'Burk','BURKERT','BURK']:
        scaling = 0.5 * np.log(1.0 + c) + 0.25 * np.log(c ** 2.0 + 1.0) - 0.5 * np.arctan(c)
        #int_mass_funct = lambda xs: 1.0 / xs * (np.log((xs ** 2.0 + 1.0)/ ((xs + 1.0) ** 2.0)) + xs * np.log(xs ** 4.0 / ((xs ** 2.0 + 1.0) * (xs + 1.0) ** 2.0)) - 2.0 * np.arctan(xs) + 2.0 - 2.0 * special.hyp2f1(-0.5, 1.0, 0.5, -xs ** 2.0))
        raw_pot_funct = lambda xs: np.where(xs > 0.0, 0.5 * np.arctan(xs) - 0.5 * np.log(1.0 + xs) + 0.5 / xs * (np.arctan(xs) - np.log(1.0 + xs) - 0.5 * np.log(1.0 + xs ** 2.0)) + 0.25 * np.log(1.0 + xs ** 2.0), 0.0)
    else:
        print ('Gravitationl potential type: "' + pot_type + '"" not recognized.  Returning point mass potential. ')
        scaling = np.log(c)
        raw_pot_funct = lambda xs: 1.0 / rs

    potential_funct = lambda rs, c=c: np.where(rs <= c, -(raw_pot_funct(c) - raw_pot_funct(rs)) / scaling - 1.0 / c, -np.divide(1.0, rs, out=np.zeros_like(rs), where=rs!=0) )
    return potential_funct
