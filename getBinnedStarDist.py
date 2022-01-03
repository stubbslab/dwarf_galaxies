#Define our function to add up the stars in a set of boxed bins

import numpy as np


def binnedStarDist (ra_range,ra_bin_width,dec_range,dec_bin_width,star_ras,star_decs,weights):
    binned_ras = np.arange(ra_range[0],ra_range[1],ra_bin_width) + ra_bin_width / 2.0
    binned_decs = np.arange(dec_range[0],dec_range[1],dec_bin_width) + dec_bin_width / 2.0
    ra_bins_mesh, dec_bins_mesh = np.meshgrid(binned_ras, binned_decs) 
    star_bins = np.zeros(np.shape(ra_bins_mesh))
    star_ra_bin_num = ((np.array(star_ras) + ra_range[0]) / ra_bin_width).astype(int)
    star_dec_bin_num = ((np.array(star_decs) + dec_range[0]) / dec_bin_width).astype(int)
    
    for i in range(len(star_ra_bin_num)):
        star_bins[star_ra_bin_num[i]][star_dec_bin_num[i]] = star_bins[star_ra_bin_num[i]][star_dec_bin_num[i]] + 1
    return ra_bins_mesh,dec_bins_mesh,star_bins
    