import numpy as np
import scipy.interpolate as interpolate
import SphericalDMPotentials as sdmp
import scipy.integrate as integrate
import scipy.special as special

    #We define the spherical velocity dispersion profile as: sigsqr_rr(r) = sigsqr_rr_inf - (sigsqr_rr_inf - sigsqr_rr_0) * 1 / (1 + r/r_sig0 ** alpha_rr), alpha_rr > 0.0
    #We define the spherical velocity dispersion profile as: 1/sigsqr_rr(r) = 1/sigsqr_rr_inf - (1/sigsqr_rr_inf - 1/sigsqr_rr_0) * exp(-(r/r0) ** alpha_rr), alpha_rr > 0.0
    #
    # dispersion_rr_params = [sigsqr_rr_0, sigsqr_rr_inf, r0, alpha_rr]
    #We define the velocity beta dispersion profile as: beta(r) = beta_inf - beta_inf * 1 / (1 + r/r_beta0 ** alpha_beta) = beta_inf * (r/r_beta0 ** alpha_beta) / (1 + r/r_beta0 ** alpha_beta)
    # dispersion_beta_params = [gamma_inf, r_beta0, alpha_beta]

def returnSigsqrRRFunct(dispersion_rr_params):
    #sigsqr_rr_0, sigsqr_rr_inf, r0, alpha_rr = dispersion_rr_params
    sigsqr_rr_0, sigsqr_rr_inf, r_sigsqr_rr0, alpha_sigsqr_rr = dispersion_rr_params
    sigsqr_rr_funct = lambda r: (1.0/sigsqr_rr_inf - (1.0/sigsqr_rr_inf - 1.0/sigsqr_rr_0) * np.exp(-(r/r_sigsqr_rr0) ** alpha_sigsqr_rr)) ** -1.0
    return sigsqr_rr_funct

def returnBetaFunct(dispersion_beta_params) :
    gamma_for_beta_inf, r_beta0, alpha_beta = dispersion_beta_params
    beta_inf = 2.0 * gamma_for_beta_inf / (1.0 + gamma_for_beta_inf)
    beta_funct = lambda rs: beta_inf - beta_inf * 1 / (1 + (rs/r_beta0) ** alpha_beta)
    return beta_funct

def returnBetaOverRIntFunct(dispersion_beta_params):
    gamma_for_beta_inf, r_beta0, alpha_beta = dispersion_beta_params
    beta_inf = 2.0 * gamma_for_beta_inf / (1.0 + gamma_for_beta_inf)
    if abs(alpha_beta) > 0.0:
        beta_over_r_int_funct = lambda rs: 2.0 * beta_inf * np.log( (rs / r_beta0) ** alpha_beta + 1.0) / alpha_beta
    else:
        beta_over_r_int_funct = lambda rs: np.log(rs / r_beta0) * beta_inf / 2.0
    return beta_over_r_int_funct


def returnDPotOverSig_rr_Nint_Funct(dispersion_rr_params, halo_type, n_interp_rs, c):
    sigsqr_rr_0, sigsqr_rr_inf, r0, alpha_rr = dispersion_rr_params
    sigsqr_rr_funct = returnSigsqrRRFunct(dispersion_rr_params)
    #over_rsqr_sigsqr_int = lambda r: ((sigsqr_rr_0 - sigsqr_rr_inf) * special.hyp2f1(1, -1/alpha_rr, (alpha_rr - 1) / alpha_rr, - sigsqr_rr_inf * (r/r0) ** alpha_rr / sigsqr_rr_0) - sigsqr_rr_0 ) / (sigsqr_rr_0 * sigsqr_rr_inf * r)
    #over_rsqr_sigsqr_int = lambda r: [integrate.quad(lambda rp: dphi_dr_funct(rp) / sigsqr_rr_funct(rp), 0, r)[0] for r in interp_rs]
    over_rsqr_sigsqr_int = lambda r: np.where(np.abs(r) > 0.0, -(1.0 / sigsqr_rr_inf) / np.array(r) - (1.0 / sigsqr_rr_inf - 1.0 / sigsqr_rr_0) * (special.gammaincc(-1.0 / alpha_rr + 1.0, (c/r0) ** alpha_rr) - special.gammaincc(-1.0 / alpha_rr + 1.0, (np.array(r)/r0) ** alpha_rr)) / (alpha_rr * r0 * special.gamma(-1.0 / alpha_rr + 1.0)), 0.0)
    interp_rs = np.linspace(0.0, c, n_interp_rs)
    dphi_dr_funct = sdmp.getPotentialDerivativeFunction(halo_type, c)
    interp_dphi_drs_over_sigs = [integrate.quad(lambda rp: dphi_dr_funct(rp) / sigsqr_rr_funct(rp), 0, r)[0] for r in interp_rs]
    dphi_dr_over_sigrsqr_interp = interpolate.interp1d(interp_rs, interp_dphi_drs_over_sigs, fill_value = 0.0, bounds_error = False)
    full_interp = lambda rs: np.where(rs < c, dphi_dr_over_sigrsqr_interp (rs), over_rsqr_sigsqr_int(rs) - over_rsqr_sigsqr_int(c) + dphi_dr_over_sigrsqr_interp (c))
    return full_interp
