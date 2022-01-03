#Define a function to read in dSph data from a file.
#The current format is that given by Walker et al.  It should be updatable to accomodate other files, but that would require a good deal of hard coded change.

import matplotlib
import matplotlib.pyplot as plt
import csv
import math
import numpy as np
import scipy.integrate as integrate
from matplotlib.colors import LogNorm
from logList import logList 
from smoothedStarDist import smoothedStarDist 
from DwarfGalDataArchive import DwarfGalDataArchive
from AstronomicalParameterArchive import AstronomicalParameterArchive
from cantrips import readInFileLineByLine 

def readInGalaxyStarData(population, pop_selection_method='none', plot='no', plot_metallicity='no', save_fig = 0):
    data_archive = DwarfGalDataArchive()
    data_file = data_archive.getFile(population)
    #vhel_index = 7
    corr_vhel_index = 9
    corr_vhel_err_index = 10
    proj_x_index = 19
    proj_y_index = 20 
    #Get distance from sun in parsecs. 
    dist_from_sun = data_archive.getDistanceFromSun(population) #{"carina":105000,"fornax":147000,"sculptor":86000,"sextans":86000} in pc
    dist_from_sun_err = data_archive.getDistanceFromSunErr(population)
    total_mass = data_archive.getTotalMass(population) #{"carina":1.28*10**8,"fornax":1.28*10**8,"sculptor":1.28*10**8,"sextans":1.28*10**8} in M_sun 
    astro_params = AstronomicalParameterArchive() 
    
    #Read the data into various variables 
    Target= []
    Field = []
    Date= []
    RA = []
    DEC = [] 
    Vmag=[]
    VImag =[]
    Vhel=[]
    VhelE =[]
    SigFe=[]
    SigFeE=[]
    SigMg=[]
    SigMgE = []
    Pm=[]
    RFe=[]
    RMg=[]
    RFeE=[]
    RMgE=[]

    #The membership probability threshold to include star in galaxy.  
    inclusion_limit = 0.75

    print ('data_file = ' + str(data_file)) 
    #file = open(data_file, 'rb')
    #Now read in data line by line.  Note that it goes by character in the string.
    lines = readInFileLineByLine(data_file)
    for line in lines: 
        #Tells us whether to include a given line (based on likelihood of membership to dSph)
        #Note that if a given line is not the first instance of a star, it has no pm value.
        #By assigning such lines a pm value of 0, we avoid reading in data on the same star twice. 
        pm = float(line[93:98]) if line[93:98]!='     ' else 0.0 #Population Membership
        #if star inclusion likelihood is above threshold, add to list 
        if pm >=inclusion_limit:
            Target.append(line[0:8])
            Field.append((line[9:12] if line[11] == ',' else line[9:13]).strip()) 
            date = line[14:22]
            Date.append(2450000+float(date))
            #give the right ascension in degrees (RA in archours in file)
            ra = (float(line[23:25])+ float(line[26:28])/60.0 + float(line[29:34])/3600.0)/24.0*360.0
            RA.append(ra)
            if line[35]=='-':
                dec=-1*float(line[36:38])-float(line[39:41])/60.0-float(line[42:46])/3600.0
            else:
                dec=float(line[36:38])+float(line[39:41])/60.0+float(line[42:46])/3600.0
            DEC.append(dec)
            vmag = float(line[47:52]) if line[51]!=' ' else 0.0
            Vmag.append(vmag)
            vimag = float(line[53:58]) if line[57]!=' ' else 0.0
            VImag.append(vimag)
            vhel = float(line[99:104]) if line[99:104]!='     ' else float(line[59:65]) if line[59:65]!='      ' else 0.0
            Vhel.append(vhel)
            vhelE = float(line[105:108]) if line[105:108]!='   ' else float(line[66:70]) if line[66:70]!='    ' else 0.0
            VhelE.append(vhelE)
            sigFe = float(line[71:76]) if line[75]!=' ' else 0.0
            SigFe.append(sigFe)
            sigFeE = float(line[77:81]) if line[80]!=' ' else 0.0
            SigFeE.append(sigFeE)    
            sigMg = float(line[109:114]) if line[109:114]!='     ' else float(line[82:87]) if line[82:87]!='     ' else 0.0
            SigMg.append(sigMg)
            sigMgE = float(line[115:119]) if line[115:119]!='    ' else float(line[88:92]) if line[88:92]!='    ' else 0.0
            SigMgE.append(sigMgE)
    #print 'We read in a total of ' + str(len(SigMg)) + ' stars.  ' 
    #Vega_solar_lums = astro_params.getVegaSolarLums() #40.12 #solar luminosity of Vega
    #solar_lum = 3.825 * 10**23 #solar luminosity in kilowatts
    #Vega_dist_pars = astro_params.getVegaDistPars() #7.68 #distance to Vega in parsecs
    #parsec_to_m = 3.086 * 10**16
    V_intensities = astro_params.convertMagnitudesToIntensities (np.array(Vmag)) #10**(np.array(Vmag)/(-2.5))*(Vega_solar_lums)/(4.0*math.pi*Vega_dist_pars**2)
    #print 'Data from files read in. '

    #To correct metallicities (according to Walker 2009), we need to apply a linear correction.
    # For that, we need the apparent magnitude of the horizontal branch stars (which they generously DO NOT provide for us)
    # SO I'm left to base VHB for:
    # Carina, Fornax, Sculptor from http://adsabs.harvard.edu/abs/2011ApJ...742...20W
    # Fornax on the Battaglia paper (fig 4), and I'm not sure what to do for the others
    # Could read off VHB for Carina from http://arxiv.org/pdf/astro-ph/0301388v1.pdf
    # Couldread off VHB for Sculptor from http://www.aanda.org/articles/aa/pdf/2011/04/aa16398-10.pdf
    # Attempt to read off VHB for Sextans from http://iopscience.iop.org/article/10.1086/379171/pdf
    V_HB = data_archive.getHorizontalBranchVMagnitude(population) #V_HB = {"carina" : 20.9, "fornax" : 21.3,"sculptor" : 20.1,"sextans" : 20.25}
    SigMg_scale_corr = data_archive.getSigMgScaleCorr(population) #0.079
    SigMg_corr = [SigMg[i] + SigMg_scale_corr * (Vmag[i] - V_HB) for i in range(len(SigMg))]

    #We also will need distance from the sun, as measured in Amarisco et al.
    
   

    #We want to convert the RA/DEC coordiantes into more usable coordiantes (ie, just raw angular separation)
    #First, we need the center of the galaxies, in RA and DEC.  These are just hardcoded values pulled from papers
    ra_dec_tuple = data_archive.getRaDec(population) 
    #car_ra= (float(6) + float(41)/60.0 + float(36.7)/3600.)/24.0*360.0
    #car_dec=(float(-50) + float(-57)/60.0 + float(-58)/3600.)
    #for_ra=(float(2) + float(39)/60.0 + float(52)/3600.)/24.0*360.0
    #for_dec=(float(-34) + float(-26)/60.0 + float(-57)/3600.)
    #scu_ra=(float(1) + float(0)/60.0 + float(9.3)/3600.)/24.0*360.0
    #scu_dec=(float(-33) + float(-42)/60.0 + float(-33)/3600.)
    #sex_ra=(float(10) + float(13)/60.0 + float(2.9)/3600.)/24.0*360.0
    #sex_dec=(float(-1) + float(-36)/60.0 + float(-53)/3600.)

    #ra_dec_tuples = {"carina" : (car_ra,car_dec),"fornax" : (for_ra,for_dec),"sculptor" : (scu_ra,scu_dec),"sextans" : (sex_ra,sex_dec)}
    
    
    #Compute spherical coordinates recentered at center of galaxy
    # coordinate transformations pulled from stack exchange at this link: 
    # http://math.stackexchange.com/questions/92301/transforming-from-one-spherical-coordinate-system-to-another
    cent_ra=ra_dec_tuple[0]
    cent_dec=ra_dec_tuple[1]
    corr_ra = np.array(RA) - cent_ra
    deg_to_rad = astro_params.getDegToRad() #math.pi/180.0
    x1=np.cos(np.array(DEC) * deg_to_rad) * np.cos(corr_ra * deg_to_rad)
    x2=np.cos(np.array(DEC) * deg_to_rad) * np.sin(corr_ra * deg_to_rad)
    x3=np.sin(np.array(DEC) * deg_to_rad)

    x1bar=math.cos(cent_dec * deg_to_rad)*x1+math.sin(cent_dec * deg_to_rad)*x3
    x2bar=x2
    x3bar=-math.sin(cent_dec * deg_to_rad)*x1 + math.cos(cent_dec * deg_to_rad)*x3
    corr_ra = np.angle(x1bar + x2bar*1j, deg = True) #centered RA in degrees
    corr_dec = np.arcsin(x3bar) * (1 / deg_to_rad) #centered Dec in degrees
    radial_dists = (np.sqrt(corr_ra**2 + corr_dec**2) * deg_to_rad * dist_from_sun) #physical distance of star from center, in pc
 

    #We also want to project these angular positions onto the plane tangent to the center of the galaxy.
    # These are important since we are integrating along a square prism
    theta = -corr_dec + 90.0 #traditional spherical theta coord in degrees
    phi = corr_ra #traditional spherical phi coord in degrees
    proj_x = np.sin(phi * deg_to_rad) / np.cos(phi * deg_to_rad) / deg_to_rad
    proj_y = np.cos(theta * deg_to_rad) / (np.cos(phi * deg_to_rad) * np.sin(theta * deg_to_rad)) / deg_to_rad
     

    #We also need to correct the heliocentric los velocities to account for the proper motion of the Fornax dSph
    #Gives proper motion in mas/yr; we want motion in km/s 
    helio_prop_mot = data_archive.getHeliocentricProperMotion(population)
    helio_prop_mot_err = data_archive.getHeliocentricProperMotionError(population)
    helio_prop_mot_norm = math.sqrt(helio_prop_mot[0]**2 + helio_prop_mot[1]**2)

    angMotionConvert = astro_params.getAngularMotionConversionFactor()
    corr_Vhel = [Vhel[i] -  angMotionConvert * dist_from_sun * (corr_ra[i] * helio_prop_mot[0] + corr_dec[i] * helio_prop_mot[1]) for i in range(len(Vhel))]
    corr_VhelE = [math.sqrt(
                   VhelE[i]**2
                    + dist_from_sun_err ** 2 * (angMotionConvert * (corr_ra[i] * helio_prop_mot[0] + corr_dec[i] * helio_prop_mot[1]) )**2 
                    + helio_prop_mot_err[0] ** 2 * (angMotionConvert * dist_from_sun * corr_ra[i]) ** 2
                    + helio_prop_mot_err[1] ** 2 * (angMotionConvert * dist_from_sun * corr_dec[i]) ** 2
                    ) for i in range(len(Vhel)) ]
    
    #Need to assign membership probabilities, based on paramaters in this paper: http://arxiv.org/pdf/1206.6691v2.pdf
    pop_params=data_archive.getPopParameters(population) # {"carina": [[0.27,0.07,888],[0.46,0.058,605],[0.53,0.11,480]],"fornax": [[0.27,0.07,888],[0.46,0.058,605],[0.53,0.11,480]],"sculptor":[[0.24,0.022,302],[0.24,0.022,302],[0.36,0.022,166]],"sextans":[[0.36,0.077,166],[0.24,0.022,302],[0.24,0.022,302]]} #for each of MP, IM, MR populations, the best average SigMg, StDSigMg, Rh
    #pop_params[0][1]
    distinct_pops = data_archive.getPopulations(population[0])
    metallicity_cuts=data_archive.getMetallicityCuts(population) #{"carina":[0.04,0.57],"fornax":[0.28,0.5],"sculptor":[0.3,0.5],"sextans":[-0.9,-0.6]}
    
    #Define normal distribution
    #mu is mean of normal, sigma is standard deviation
    normal = lambda x, mu, sigma: 1/math.sqrt(2 * sigma**2 * math.pi) * math.exp(-(x-mu)**2/(2*sigma**2)) 

    #Define the Plummer density profile
    #rh is scale radius of halo
    #3d profile is: (3M)/(4 Pi rh^3)*(1 + r^2/rh^2)^(-5/2)
    #integrate z from -inf to inf, 2d profile is: M/(rh^2 Pi (1+(R/rh)^2)**2)
    #So total luminosity (stand in for total relative probability) in a box of area dA is:
    # M/(rh^2 Pi (1+(R/rh)^2)**2) * dA.
    #So relative probability should just be 1/(rh^2 Pi (1+(R/rh)^2)**2)
    # (and I think we could even do away with the overall constant on the front) 
    #plummerProbability = lambda rh,R: 1/(rh**2 * math.pi ) * 1/((1+(R/rh)**2)**2)

    #According to: https://drive.google.com/drive/folders/0B98bAoTCicyTNVZzb3JMdk9XSHc
    plummerProbability = lambda rh, R: 2*R/(rh**2) * 1/(1 + R**2/(rh**2) )

    #print 'pop_params[i] are: '
    #print pop_params

    #print 'distinct_pops = '
    #print distinct_pops

    #print metallicity cuts for metallicity point method
    
    
    memb_prob_from_metallicity_point = [[normal(sigmaMg, pop_params[pop]['meanSigMg'], pop_params[pop]['meanSigMgE']) for pop in distinct_pops] for sigmaMg in SigMg_corr]
    memb_prob_from_position_point = [[plummerProbability(pop_params[pop]['Rh'], radial_dist) for pop in distinct_pops] for radial_dist in radial_dists.tolist()]
    full_memb_prob_point = np.array(memb_prob_from_position_point)*np.array(memb_prob_from_metallicity_point)

    break_up_populations = 1
    #For now, let's say we always have three populations.
    #We have 4 selection method options:
    # 'metal_rigid' = classify stars by rigid metallicity cuts; is a star below or above specific values
    # 'metal_point' = classify stars by where they most likely belong, based on normal probability distributions of metallicity
    # 'position_point' = classify stars by where they most likely belong, based on normal probability distributions of position
    # 'point' = classify stars by where they most likely belong, based on direct product of normal probability distributions of position and metallicity
    pop_selection_method = pop_selection_method.lower()
    #print 'metallicity_cuts = '
    #print metallicity_cuts
    if pop_selection_method == 'metal_rigid':
        pop_membership=['MP' if metallicity <= metallicity_cuts[0] else 
                        'MR' if metallicity >= metallicity_cuts[1] else
                        'IM'
                        for metallicity in SigMg_corr]
    elif pop_selection_method == 'metal_point':
        pop_membership = ['MP' if probabilities[0] == max(probabilities) else 
                          'MR' if probabilities[2] == max(probabilities) else 
                          'IM'
                          for probabilities in memb_prob_from_metallicity_point]
    elif pop_selection_method == 'position_point':
        pop_membership = ['MP' if probabilities[0] == max(probabilities) else 
                          'MR' if probabilities[2] == max(probabilities) else 
                          'IM'
                          for probabilities in memb_prob_from_position_point]
    elif pop_selection_method == 'point' or pop_selection_method == 'full_point':
        pop_membership =['MP' if probabilities[0] == max(probabilities) else 
                         'MR' if probabilities[2] == max(probabilities) else 
                         'IM'
                         for probabilities in full_memb_prob_point.tolist()]
    else: 
        pop_membership = ['NoSelection' for metallicity in SigMg_corr]
        break_up_populations = 0
    
    #print 'pop_membership = '
    #print pop_membership   
    #pop_membership  = data_archive.getPopMembership(population, pop_selection_method) 
    #print 'len(Target) = ' + str(len(Target))
    #print 'len(Date) = ' + str(len(Date))
    #print 'len(RA) = ' + str(len(RA))
    #print 'len(DEC) = ' + str(len(DEC))
    #print 'len(Vmag) = ' + str(len(Vmag))
    #print 'len(VImag) = ' + str(len(VImag))
    #print 'len(Vhel) = ' + str(len(Vhel))
    #print 'len(VhelE) = ' + str(len(VhelE))
    #print 'len(SigFe) = ' + str(len(SigFe))
    #print 'len(SigFeE) = ' + str(len(SigFeE))
    #print 'len(SigMg) = ' + str(len(SigMg))
    #print 'len(SigMgE) = ' + str(len(SigMgE))
    #print 'len(V_intensities) = ' + str(len(V_intensities))
    #print 'len(SigMg_corr) = ' + str(len(SigMg_corr))
    #print 'len(corr_ra) = ' + str(len(corr_ra))
    #print 'len(corr_dec) = ' + str(len(corr_dec))

    #All of the data that we possibly could return 
    returnSet = [Target,        Field,      Date,    RA,      DEC, Vmag, VImag, Vhel, VhelE, corr_Vhel, corr_VhelE, SigFe, SigFeE, SigMg, SigMgE,
                 V_intensities, SigMg_corr, corr_ra, corr_dec, proj_x, proj_y]
    
    #Here, we prune the to-be-returned data list to select only those populations that we wish to examine.  
    if break_up_populations:
        #print 'breaking up by population '
        returnSetFull = returnSet
        returnSet = []
        #Now for each element in the return set, prune the list to only include those stars that fall in the population of interest. 
        for j in range(len(returnSetFull)):
            #print 'j = ' + str(j)  
            #print 'population[1] = ' + population[1]
            #print 'len(returnSetFull[j]) = ' + str(len(returnSetFull[j])) 
            trimmed_return_elem = np.array( [returnSetFull[j][i] for i in range (len(pop_membership)) if pop_membership[i] == population[1]] )
            #print np.shape(trimmed_return_elem) 
            #print 'trimmed_return_elem = '
            #print trimmed_return_elem
            returnSet = returnSet + [trimmed_return_elem]

    #Estimate half light radius by finding (circular) radius such that half of stars are in circle at nominal galaxy center
    
    star_dists_from_center = np.sqrt(np.array(returnSet[proj_x_index]) ** 2.0 + np.array(returnSet[proj_y_index]) ** 2.0 )
    half_light_radius = np.median(star_dists_from_center) 

    #Now we need to compute the velocity dispersion
    #vhel_index = 7
    trimmed_corr_Vhel = returnSet[corr_vhel_index]
    trimmed_corr_Vhel_err = returnSet[corr_vhel_err_index]
    trimmed_Vhel_err = returnSet[8] 
    #print 'A total of ' + str(len(trimmed_corr_Vhel)) + ' stars in this population. '
    mean_corr_Vhel = sum(trimmed_corr_Vhel) / float(len(trimmed_corr_Vhel))
    mean_corr_Vhel_err = math.sqrt(sum([err ** 2.0 for err in trimmed_corr_Vhel_err]) / len(trimmed_corr_Vhel) )
    sigSqrVhel = sum([(trimmed_corr_Vhel[i] - mean_corr_Vhel) ** 2.0 for i in range(len(trimmed_corr_Vhel)) ] ) / len(trimmed_corr_Vhel)
    sigSqrE = 2.0 / len(trimmed_corr_Vhel) * math.sqrt(sum([ ((trimmed_corr_Vhel[i] - mean_corr_Vhel)) ** 2.0  * (trimmed_corr_Vhel_err[i] ** 2.0) for i in range(len(trimmed_corr_Vhel))]) )

    #Now we need to compute the velocity dispersion
    #vhel_index = 6
    #old_trimmed_corr_Vhel = returnSet[corr_vhel_index - 2]
    #print 'A total of ' + str(len(trimmed_corr_Vhel)) + ' stars in this population. '
    #old_mean_corr_Vhel = sum(old_trimmed_corr_Vhel) / float(len(old_trimmed_corr_Vhel))
    #old_sigSqrVhel = sum([(old_trimmed_corr_Vhel[i] - old_mean_corr_Vhel) ** 2.0 for i in range(len(old_trimmed_corr_Vhel)) ] ) / len(old_trimmed_corr_Vhel) 
    #print 'trimmed_Vhel = '
    #print trimmed_Vhel
    #print 'old_mean_corr_Vhel = '
    #print old_mean_corr_Vhel 
    #print 'old_sigSqrVhel = '
    #print old_sigSqrVhel

    returnSet = returnSet + [half_light_radius, sigSqrVhel, sigSqrE] 
    
        #returnSet = [np.array( [returnSet[j][i] for i in range(len(pop_membership)) if pop_membership[i] == population[1]] ) for j in range(len(returnSet))]
    #print 'len(returnSet) = ' + str(len(returnSet))
    #print 'len(returnSet[0]) = ' + str(len(returnSet[0]))
    #returnSet = returnSet + [dist_from_sun]
    #returnSet = returnSet + [full_mass] 
        
    #corr_ra_one_pop = np.array([corr_ra[i] for i in range(len(pop_membership)) if pop_membership[i] == population[1]])
    #corr_dec_one_pop = np.array([corr_dec[i] for i in range(len(pop_membership)) if pop_membership[i] == pop_of_interest[1]])
     
    
    #[[normal(sigmaMg, pop_params[i][0], pop_params[i][1]) for i in range(3)] for sigmaMg in SigMg_corr]
    #memb_prob_from_position_point = [[plummerProbability(pop_params[i][2], radial_dist) for i in range(3)] for radial_dist in radial_dists.tolist()]
    #full_memb_prob_point = np.array(memb_prob_from_position_point)*np.array(memb_prob_from_metallicity_point)

    #if pop_selection_method == 'metal rigid':
    #    pop_membership=['MP' if metallicity <= metallicity_cuts[0] else 
    #                 'MR' if metallicity >= metallicity_cuts[1] else
    #                 'IM' for metallicity in SigMg_corr]
    #elif pop_selection_method == 'metal point':
    #    pop_membership = ['MP' if probabilities[0] == max(probabilities) else 
    #                 'MR' if probabilities[2] == max(probabilities) else 
    #                 'IM' for probabilities in memb_prob_from_metallicity_point]
    #elif pop_selection_method == 'position point':
    #    pop_membership = ['MP' if probabilities[0] == max(probabilities) else 
    #                 'MR' if probabilities[2] == max(probabilities) else 
    #                 'IM' for probabilities in memb_prob_from_position_point]
    #elif pop_selection_method == 'point' or pop_selection_method == 'full point':
    #    pop_membership =['MP' if probabilities[0] == max(probabilities) else 
    #                 'MR' if probabilities[2] == max(probabilities) else 
    #                 'IM' for probabilities in full_memb_prob_point.tolist()]
    #else pop_membership = ['NoSelection']
    # 
    #pop_member

    #pop_membership_point =['MP' if probabilities[0] == max(probabilities) else 
    #                 'MR' if probabilities[2] == max(probabilities) else 
    #                 'IM' for probabilities in full_memb_prob_point.tolist()]
    #print metallicity_cuts[gal_of_interest]
    #pop_membership_rigid=['MP' if metallicity <= (metallicity_cuts[gal_of_interest])[0] else 
    #                 'MR' if metallicity >= (metallicity_cuts[gal_of_interest])[1] else
    #                 'IM' for metallicity in SigMg_corr]

#################################################################################################################
#Below, we just plot of bunch of parameters, if we're supposed to.  Then we return all of the above lists of values. 
#################################################################################################################

    #Show the metallicity distribution (or don't) 
    if plot_metallicity == 'yes':
        print ('len(SigMg_corr) = ' + str(len(SigMg_corr)) )
        metallicity_binwidth=0.014
        fig=plt.figure(1,figsize=(10,10))
        bins= np.arange(min(SigMg_corr), max(SigMg_corr) + metallicity_binwidth, metallicity_binwidth)
        plt.hist(SigMg_corr,bins=bins,color = 'w')
        #metallicity_pop_curves = []
        for pop in distinct_pops:
            print ('pop = ' + pop )
            color = 'r' if pop == 'MR' else 'g' if pop =='IM' else 'b' if pop =='MP' else  'c'
            new_curve = [normal(bin,pop_params[pop]['meanSigMg'],pop_params[pop]['meanSigMgE']) for bin in bins]
            plt.plot(bins,new_curve) 
        plt.ylabel('number of stars')
        plt.xlabel('Corrected SigMg (Angstroms)')
        font = {'family': 'serif',
                'color':  'red',
                'weight': 'normal',
                'size': 16,
                'rotation': 90}
        #plt.text(0.04,0,'<--- lower metallicity cut',fontdict=font,ha='center',va='bottom' )
        #plt.text(0.57,0,'<--- upper metallicity cut',fontdict=font,ha='center',va='bottom' )
        plt.title(population[0][0].upper() + population[0][1:] + ' Corrected Metallicities')
        pltFileDir="/Users/sasha/Documents/Harvard/physics/randall/plots/"
        plt.show()
        if save_fig:
            fig.savefig(pltFileDir + population[0] + '_metallicity_histogram.pdf')
    if plot == 'yes': 
        #Show the radial distribution of star (or don't) 
        radial_binwidth=50
        fig=plt.figure(1,figsize=(10,10))
        plt.hist(radial_dists,bins=np.arange(min(radial_dists), max(radial_dists) + radial_binwidth, radial_binwidth))
        plt.ylabel('number of stars')
        plt.xlabel('Angular distance from center of galaxy')
        #plt.text(-0.9,0,'<--- lower metallicity cut',fontdict=font,ha='center',va='bottom' )
        #plt.text(-0.6,0,'<--- upper metallicity cut',fontdict=font,ha='center',va='bottom' )
        plt.title(gal_of_interest + ' distance from gal center')
        pltFileDir="/Users/sasha/Documents/Harvard/physics/randall/dwarfGalaxyProject/plots/"
        #Uncomment the following line if you want to save the file
        #plt.savefig(pltFileDir + gal_of_interest + '_metallicity_histogram.png')
        
        plt.show()

        #Show the full galaxy of interest (or don't
        ylims = (-60,60) #ylimits from Walker plot
        xlims=(60,-60) #xlimits from Walker plot
        fig=plt.figure(1,figsize=(10,10))
        plt.scatter(corr_ra * 60.0,corr_dec * 60.0 ,s=2)
        plt.ylim(ylims)
        plt.xlim(xlims)
        plt.xlabel('ra sep')
        plt.ylabel('dec sep')
        plt.title('Stars Positions on Sky of desired Population')
        pltFileDir="/Users/sasha/Documents/Harvard/physics/randall/dwarfGalaxyProject/plots/"
        #Uncomment the following line if you want to save the file
        #plt.savefig(pltFileDir + gal_of_interest + '_all.png')
        
        #Uncomment following line to show this plot
        plt.show()

        #Create the smoothed star distribution 
        ra_mesh,dec_mesh = np.meshgrid(np.arange(-60,60 + 3.0,3.0)[::-1],np.arange(-60,60 + 3.0,3.0))
        #print dec_mesh
        pop_of_interest = 'all' 

        #corr_ra_one_pop = np.array([corr_ra[i] for i in range(len(pop_membership_point)) if pop_membership_point[i] == pop_of_interest])
        #corr_dec_one_pop = np.array([corr_dec[i] for i in range(len(pop_membership_point)) if pop_membership_point[i] == pop_of_interest])
        #print len(corr_dec_one_pop)
        smooth_param = 2.0
        intensities = V_intensities
        smoothed_star_mesh = smoothed_star_dist(ra_mesh,dec_mesh,smooth_param,corr_ra*60.0,corr_dec*60.0,intensities*10**9)

        #Display the smoothed distribution on the sky (or don't) 
        matplotlib.rcParams['xtick.direction'] = 'out'
        matplotlib.rcParams['ytick.direction'] = 'out'
        plt.figure(figsize=(13,11))
        #im=plt.imshow(cross_sect,interpolation='spline16',origin='lower',cmap=cm.cool,extent=(-5.4,5.5,-2.5,2.5)) #not working, for some reason
        CS=plt.contour(ra_mesh,dec_mesh,smoothed_star_mesh,20)
        CB_contour=plt.colorbar(shrink=0.8,extend='both',format='%.2f')
        #CB_gradient=plt.colorbar(im,shrink=0.8, orientation='horizontal')
        l,b,w,h = plt.gca().get_position().bounds
        ll,bb,ww,hh = CB_contour.ax.get_position().bounds
        CB_contour.ax.set_position ([ll,b+0.1*h,ww,h*0.8])
        plt.title('smoothed sky brightness of stars on sky')
        plt.xlabel('R')
        plt.ylabel('z')
        pltFileDir="/Users/sasha/Documents/Harvard/physics/randall/dwarfGalaxyProject/plots/"
        #Uncomment the following line if you want to save the file
        #plt.savefig(pltFileDir + gal_of_interest + '_' + pop_of_interest + '_smoothed_isophotes' + '.png')

        #Uncomment the following line if you want to plot the smoothed star distribution 
        #plt.show()

        #Define the binned distribution of stars 
        ra_mesh,dec_mesh,binned_star_mesh = binned_star_dist ([-60,60],2.0,[-60,60],2.0,corr_ra*60.0,corr_dec*60.0,intensities*10**9)
    
        #Show the binned galaxy of interest (or don't)
        matplotlib.rcParams['xtick.direction'] = 'out'
        matplotlib.rcParams['ytick.direction'] = 'out'
        plt.figure(figsize=(13,11))
        #im=plt.imshow(cross_sect,interpolation='spline16',origin='lower',cmap=cm.cool,extent=(-5.4,5.5,-2.5,2.5)) #not working, for some reason
        CS=plt.contour(ra_mesh,dec_mesh,binned_star_mesh,20)
        CB_contour=plt.colorbar(shrink=0.8,extend='both',format='%.2f')
        #CB_gradient=plt.colorbar(im,shrink=0.8, orientation='horizontal')
        l,b,w,h = plt.gca().get_position().bounds
        ll,bb,ww,hh = CB_contour.ax.get_position().bounds
        CB_contour.ax.set_position ([ll,b+0.1*h,ww,h*0.8])

        plt.title('binned sky brightness')
        plt.xlabel('R')
        plt.ylabel('z')

        pltFileDir="/Users/sasha/Documents/Harvard/physics/randall/dwarfGalaxyProject/plots/"
        #Uncomment the following line if you want to save the file
        #plt.savefig(pltFileDir + gal_of_interest + '_' + pop_of_interest + '_binned_isophotes' + '.png')

        #Uncomment the following line if you want to display the file 
        #plt.show()
    
        #Make HR diagram and show it (or don't) 
        fig=plt.figure(1,figsize=(6,6))
        plt.scatter(VImag,Vmag,s=0.5)
        ylims = (24,17) #ylimits from Walker plot
        xlims=(-0.5,2.5) #xlimits from Walker plot
        plt.ylim(ylims)
        plt.xlim(xlims)
        plt.xlabel('V-I')
        plt.ylabel('V')
        plt.title(gal_of_interest + ' Stars H-R Diagram')
        #Uncomment the following line if you want to display the file 
        plt.show()
    
        #Need to assign membership probabilities, based on paramaters in this paper: http://arxiv.org/pdf/1206.6691v2.pdf
        pop_params={"carina": [[0.27,0.07,888],[0.46,0.058,605],[0.53,0.11,480]],"fornax": [[0.27,0.07,888],[0.46,0.058,605],[0.53,0.11,480]],"sculptor":[[0.24,0.022,302],[0.24,0.022,302],[0.36,0.022,166]],"sextans":[[0.36,0.077,166],[0.24,0.022,302],[0.24,0.022,302]]} #for each of MP, IM, MR populations, the best average SigMg, StDSigMg, Rh
        #pop_params[0][1]
        metallicity_cuts={"carina":[0.04,0.57],"fornax":[0.28,0.5],"sculptor":[0.3,0.5],"sextans":[-0.9,-0.6]}
    
        #Define normal distribution
        normal = lambda x, mu, sigma: 1/math.sqrt(2 * sigma**2 * math.pi) * math.exp(-(x-mu)**2/(2*sigma**2)) 

        #Define the Plummer density profile
        plummerProbability = lambda rh,R: (2*R/(rh**2))/((1+(R/rh)**2)**2)
    
        memb_prob_from_metallicity_point = [[normal(sigmaMg, pop_params[gal_of_interest][i][0], pop_params[gal_of_interest][i][1]) for i in range(3)] for sigmaMg in SigMg_corr]
        memb_prob_from_position_point = [[plummerProbability(pop_params[gal_of_interest][i][2], radial_dist) for i in range(3)] for radial_dist in radial_dists.tolist()]
        full_memb_prob_point = np.array(memb_prob_from_position_point)*np.array(memb_prob_from_metallicity_point)

        pop_membership_point =['MP' if probabilities[0] == max(probabilities) else 
                     'MR' if probabilities[2] == max(probabilities) else 
                     'IM' for probabilities in full_memb_prob_point.tolist()]
        pop_membership_rigid=['MP' if metallicity <= (metallicity_cuts[gal_of_interest])[0] else 
                     'MR' if metallicity >= (metallicity_cuts[gal_of_interest])[1] else
                     'IM' for metallicity in SigMg_corr]
 
        #Now plot again with colorization for different populations (or don't)
        colormap=['b' if flag == 'MR' else 'g' if flag == 'IM' else 'r' for flag in pop_membership_point]
        ylims = (-60,60) #ylimits from Walker plot
        xlims=(60,-60) #xlimits from Walker plot
        fig=plt.figure(1,figsize=(8,8))
        plt.scatter(corr_ra * 60.0,corr_dec * 60.0 ,s=4.,color=colormap)
        plt.ylim(ylims)
        plt.xlim(xlims)
        plt.xlabel('ra sep')
        plt.ylabel('dec sep')
        plt.title('Fornax Stars Position on Sky')
        #Uncomment the following line if you want to display the file 
        #plt.show()

        #Now plot  only one population, specified by the flag specified at beginning of program 
        colormap='b' if pop_of_interest == 'MP' else 'g' if pop_of_interest == 'IM' else 'r'
        ylims = (-60,60) #ylimits from Walker plot
        xlims=(60,-60) #xlimits from Walker plot
        fig=plt.figure(1,figsize=(8,8))
        corr_ra_one_pop = np.array([corr_ra[i] for i in range(len(pop_membership_point)) if pop_membership_point[i] == pop_of_interest])
        corr_dec_one_pop = np.array([corr_dec[i] for i in range(len(pop_membership_point)) if pop_membership_point[i] == pop_of_interest])
        plt.scatter(corr_ra_one_pop * 60.0,corr_dec_one_pop * 60.0 ,s=4.,color=colormap)
        plt.ylim(ylims)
        plt.xlim(xlims)
        plt.xlabel('ra sep')
        plt.ylabel('dec sep')
        plt.title('Fornax Stars Position on Sky')
        pltFileDir="/Users/sasha/Documents/Harvard/physics/randall/dwarfGalaxyProject/plots/"
        #Uncomment the following line if you want to save the file
        #plt.savefig(pltFileDir + gal_of_interest + '_' + pop_of_interest + '_metal_cuts_' + str((metallicity_cuts[gal_of_interest])[0]) + '_' + str((metallicity_cuts[gal_of_interest])[1]) + '.png')
        #plt.savefig(pltFileDir + gal_of_interest + '_' + pop_of_interest + '_full_pop_division' + '.png')
        plt.show()
    
        #Define the binned star mesh for that particular population 
        ra_mesh,dec_mesh = np.meshgrid(np.arange(-60,60 + 0.5,0.5)[::-1],np.arange(-60,60 + 0.5,0.5))
        #print ra_mesh 
        #print dec_mesh
        colormap='b' if pop_of_interest == 'MP' else 'g' if pop_of_interest == 'IM' else 'r'
        corr_ra_one_pop = np.array([corr_ra[i] for i in range(len(pop_membership_point)) if pop_membership_point[i] == pop_of_interest])
        corr_dec_one_pop = np.array([corr_dec[i] for i in range(len(pop_membership_point)) if pop_membership_point[i] == pop_of_interest])
        #print len(corr_dec_one_pop)
        smooth_param = 2.0
        weights = np.array(corr_ra_one_pop * 0.0 + 1.0).tolist()
        smoothed_star_mesh = smoothedStarDist(ra_mesh,dec_mesh,smooth_param,corr_ra_one_pop*60.0,corr_dec_one_pop*60.0,weights)
   
   
        #Plot the smoothed contours for the particular population of interest (or don't)
        matplotlib.rcParams['xtick.direction'] = 'out'
        matplotlib.rcParams['ytick.direction'] = 'out'

        plt.figure(figsize=(13,11))

        #im=plt.imshow(cross_sect,interpolation='spline16',origin='lower',cmap=cm.cool,extent=(-5.4,5.5,-2.5,2.5)) #not working, for some reason

        CS=plt.contour(ra_mesh,dec_mesh,star_mesh,20)
        CB_contour=plt.colorbar(shrink=0.8,extend='both',format='%.2f')

        #CB_gradient=plt.colorbar(im,shrink=0.8, orientation='horizontal')

        l,b,w,h = plt.gca().get_position().bounds
        ll,bb,ww,hh = CB_contour.ax.get_position().bounds

        CB_contour.ax.set_position ([ll,b+0.1*h,ww,h*0.8])

        plt.title('smoothed star position on sky for one population')
        plt.xlabel('R')
        plt.ylabel('z')

        pltFileDir="/Users/sasha/Documents/Harvard/physics/randall/dwarfGalaxyProject/plots/"
        #Uncomment the following line if you want to save the file
        #plt.savefig(pltFileDir + gal_of_interest + '_' + pop_of_interest + '_position_metal_smoothed_isophotes' + '.png')

        #Uncomment the following line if you want to display the figure 
        #plt.show()
    
        #Uncomment the following set of lines if you want to save the smoothed stellar distribution
        full_file_dir = '/Users/sasha/Documents/Harvard/physics/randall/realDataFiles/'
        file_name = 'DG_isophotes_' + gal_of_interest + '_Gsmooth_' + str(smooth_param) + '_' + pop_of_interest +'.txt'
        np.savetxt(full_file_dir + file_name,star_mesh,delimiter=',')#,header=header)
        file_name_R_coords = 'DG_isophotes_RAs_' + gal_of_interest + '_Gsmooth_' + str(smooth_param) +'.txt'
        file_name_z_coords = 'DG_isophotes_Dec_' + gal_of_interest + '_Gsmooth_' + str(smooth_param) +'.txt'
        np.savetxt(full_file_dir + file_name_R_coords,ra_mesh,delimiter=',')
        np.savetxt(full_file_dir + file_name_z_coords,dec_mesh,delimiter=',')
    
        bin_width = 2.0 #in arcmin
        ra_mesh_one_pop,dec_mesh_one_pop,binned_star_mesh = binnedStarDist ([-60,60],2.0,[-60,60],2.0,corr_ra_one_pop*60.0,corr_dec_one_pop*60.0,intensities*10**9)
    
        #Display the binned stellar distribution for the particular population of interest (or don't) 
        matplotlib.rcParams['xtick.direction'] = 'out'
        matplotlib.rcParams['ytick.direction'] = 'out'

        plt.figure(figsize=(13,11))

        #im=plt.imshow(cross_sect,interpolation='spline16',origin='lower',cmap=cm.cool,extent=(-5.4,5.5,-2.5,2.5)) #not working, for some reason

        CS=plt.contour(ra_mesh_one_pop,dec_mesh_one_pop,star_mesh,20)
        CB_contour=plt.colorbar(shrink=0.8,extend='both',format='%.2f')
    
        #CB_gradient=plt.colorbar(im,shrink=0.8, orientation='horizontal')

        l,b,w,h = plt.gca().get_position().bounds
        ll,bb,ww,hh = CB_contour.ax.get_position().bounds

        CB_contour.ax.set_position ([ll,b+0.1*h,ww,h*0.8])

        plt.title('binned star position on sky for one population')
        plt.xlabel('R')
        plt.ylabel('z')

        pltFileDir="/Users/sasha/Documents/Harvard/physics/randall/dwarfGalaxyProject/plots/"
        #Uncomment the following line if you want to save the file
        #plt.savefig(pltFileDir + gal_of_interest + '_' + pop_of_interest + '_smoothed_isophotes' + '.png')

        #Uncomment of the following line if you want to display the figure 
        #plt.show()

####################################################################################################
#End plotting session.
#################################################################################################### 
    
    #print 'Done! Dwarf galaxy data read in.'

    return returnSet 
    #return 0
