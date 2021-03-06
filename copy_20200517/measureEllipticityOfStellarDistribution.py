#Measure ellipticity via various methods.

import math
from showDSphProfile import showDSphProfile
import matplotlib
import matplotlib.pyplot as plt
import numpy as np

#the Unweighted Quadrupole moments method
# (see https://www.kaggle.com/c/mdm/details/ellipticity)
        
def measureQuadrupoleMomentsOfStellarDistribution(star_RAs, star_Decs, star_weights = None, mean_be_center = 0, center_shift = [0.0,0.0], inclusion_radius = None,
                                                   show_stars = 0, x_lims = [-0.8,0.8], y_lims = [-0.8,0.8]):
    if not star_weights:
        star_weights = [1.0 for star_RA in star_RAs]

    mean_star_RA = sum( [star_weights[i] * star_RAs[i] for i in range(len(star_RAs))] ) / sum(star_weights)
    mean_star_Dec = sum( [star_weights[i] * star_Decs[i] for i in range(len(star_Decs))] ) / sum(star_weights)
    if mean_be_center:
        center_shift = [mean_star_RA,mean_star_Dec]
    #print 'center_shift = '
    #print center_shift 
    if not inclusion_radius:
        includedDecs = star_Decs
        includedRAs = star_RAs
    else: 
        includedDecs = [star_Decs[i] for i in range(len(star_Decs)) if math.sqrt((star_RAs[i] - center_shift[0]) ** 2 + (star_Decs[i] - center_shift[1]) ** 2) < inclusion_radius ]
        includedRAs = [star_RAs[i] for i in range(len(star_RAs)) if math.sqrt((star_RAs[i] - center_shift[0]) ** 2 + (star_Decs[i] - center_shift[1]) ** 2) < inclusion_radius ]

    
    Q11 = sum( [star_weights[i] * (mean_star_RA - includedRAs[i]) * (mean_star_RA - includedRAs[i]) for i in range(len(includedRAs))] ) / sum(star_weights)
    Q12 = sum( [star_weights[i] * (mean_star_RA - includedRAs[i]) * (mean_star_Dec - includedDecs[i]) for i in range(len(includedRAs))] ) / sum(star_weights)
    Q22 = sum( [star_weights[i] * (mean_star_Dec - includedDecs[i]) * (mean_star_Dec - includedDecs[i]) for i in range(len(includedRAs))] ) / sum(star_weights)

    if Q11 == 0.0 and Q12 == 0.0 and Q22 == 0.0:
        print 'No stars included in inclusion radius.'
        return (0.0,0.0)

    (el1,el2) = ( (Q11-Q22) / (Q11 + Q22 + 2 * (Q11 * Q22 - Q12 ** 2) ** 0.5),
                  (2 * Q12) / (Q11 + Q22 + 2 * (Q11 * Q22 - Q12 ** 2) ** 0.5) ) 
                  

    if show_stars:
        plt.scatter(includedRAs, includedDecs,s=4)
        matplotlib.rcParams['xtick.direction'] = 'out'
        matplotlib.rcParams['ytick.direction'] = 'out'

        plt.ylim = y_lims #ylimits from Walker plot
        plt.xlim = x_lims #xlimits from Walker plot
        plt.xlabel('ra sep')
        plt.ylabel('dec sep')
        plt.title('Fornax Stars Position on Sky')
        plt.axis('scaled')
        plt.axis(x_lims + y_lims) #combined limits from Walker plot 
        plt.show()

    return (el1,el2)

def plotEllipticities(star_RAs, star_Decs, inclusion_radii, star_weights = None, mean_be_center = 0, center_shift = [0.0,0.0],
                                                   show_plot = 0, x_lims = [-0.8,0.8], y_lims = [-0.8,0.8], save_plot = 0,
                                                   measure_method = 'quad',plot_dir = '/Users/sasha/Documents/Harvard/physics/randall/plots/'):
    el1 = []
    el2 = []
    for radii in inclusion_radii: 
        if measure_method == 'quad':
            new_els = measureQuadrupoleMomentsOfStellarDistribution(star_RAs, star_Decs,
                                                                    star_weights = star_weights, mean_be_center = mean_be_center, center_shift = center_shift, inclusion_radius = radii,
                                                                    show_stars = 0, x_lims = x_lims, y_lims = y_lims)
            el1 = el1 + [new_els[0]]
            el2 = el2 + [new_els[1]]
    #print 'el1 = '
    #print el1
    #print 'el2 = '
    #print el2
    #f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex='col', figsize = (15,8))
    #elsqr =[math.sqrt(el1[i]**2 + el2[i]**2) for i in range(len(el1)) ]
    #theta = [math.atan(el2[i] / el1[i]) / 2.0 if abs(el1[i]) > 0.0 else 0.0 for i in range(len(el1)) ]
    #ax1.scatter(inclusion_radii.tolist(), el1)
    #ax1.set_title('el1 (ellipticity in RA)')
    #ax2.scatter(inclusion_radii.tolist(), el2)
    #ax2.set_title('el2 (ellipticity in Dec)')    
    #ax3.scatter(inclusion_radii.tolist(), elsqr)
    #ax3.set_xlabel('Inclusion radius (degrees)' )
    #ax3.set_title('(el1 ^2 + el2 ^ 2)^0.5 (overall ellipticity)') 
    #ax4.scatter(inclusion_radii.tolist(), theta)
    #ax4.set_xlabel('Inclusion radius (degrees)' )
    #ax4.set_title('theta (angle of semimajor axis off of increasing RA)')
    f,ax = plt.subplots(figsize = (9,9))
    ax.plot(el1,el2,'-o')
    ax.set_title('Ellipticity values for various inclusion radii')
    ax.set_xlabel('el1')
    ax.set_ylabel('el2')
    for i in range(1,len(inclusion_radii),10):
        ax.annotate('ri=' + str(inclusion_radii[i]), (el1[i],el2[i]))
    #plt.axis('scaled')
    #plt.axis(x_lims + y_lims)
    ax.set_xbound(lower = -1.0, upper = 1.0 )
    ax.set_ybound(lower = -1.0, upper = 1.0 )
    if save_plot:
        
        center_string = 'centers_'
        if center_shift[0] < 0.0:
            center_string = center_string + 'neg' + '%.4f' % abs(center_shift[0]) + '_'
        else:
            center_string = center_string + 'pos' + '%.4f' % abs(center_shift[0]) + '_'
        if center_shift[1] < 0.0:
            center_string = center_string + 'neg' + '%.4f' % abs(center_shift[1])
        else:
            center_string = center_string + 'pos' + '%.4f' % abs(center_shift[1])   

        plot_name = 'Ellipticities_' + center_string + '.png'
        print 'Saving plot to directory ' + plot_dir + plot_name 
        plt.savefig(plot_dir + plot_name)
    if show_plot:
        plt.show()


    return 0
