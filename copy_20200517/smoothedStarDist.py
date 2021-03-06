#Return a surface luminosity profile based on an array of stars and some chosen smoothing. 
import numpy as np

def smoothedStarDist (x,y, width, star_ras,star_decs, weights): 
    full_lum_prof = np.zeros(np.shape(x))
    for i in range(len(star_ras)):
        #print i
        full_lum_prof = full_lum_prof + gaussian_lum_prof(x,y,star_ras[i],star_decs[i],width,weights[i])
    return full_lum_prof
