#Return probability of getting a given binned probability distribution from a binned probability array

import numpy as np
from nChooser import nChooser 

def probOfStarDistribution (star_dist, prob_array):
    if np.shape(star_dist) != np.shape(prob_array):
        print 'Star distribution to test and probability array must have same dimensions.'
        print 'Dimensions of given star distribution and probability array are:'
        print  str(np.shape(star_dist)) + ' and ' + str(np.shape(prob_array)) + ', respectively.'
        return 0.0
    N_precision_factors = 0 # number of factors of additional precision needed to tell our result isn't 0
    N_stars = int(np.sum(star_dist))
    print 'Star_dist has a total of ' + str(N_stars) + ' stars.'
    remaining_stars = N_stars
    star_dist_flat = star_dist.flatten().astype(int)
    prob_array_flat = prob_array.flatten()
    total_prob = 1.0
    for i in range (np.size(star_dist_flat)):
        bin_prob = prob_array_flat[i]
        bin_stars = star_dist_flat[i]
        #if bin_stars > 0: print 'bin_stars = ' + str(bin_stars)
        #if bin_stars > 0: print 'remaining_stars = ' + str(remaining_stars)
        #print 'with ' + str(remaining_stars) + ' left, total_prob = ' + str(total_prob)
        total_prob = total_prob * (bin_prob ** bin_stars) * nChooser(remaining_stars,bin_stars)
        if total_prob < 10**(-100):
            total_prob = total_prob * 10**100
            N_precision_factors = N_precision_factors + 1
        remaining_stars = int(remaining_stars - bin_stars)
    return total_prob, N_precision_factors