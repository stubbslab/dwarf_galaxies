import matplotlib
import matplotlib.pyplot as plt
import csv
import math
import numpy as np
import scipy.integrate as integrate
from matplotlib.colors import LogNorm
from logList import logList 
from smoothedStarDist import smoothedStarDist
from ComputationalArchive import ComputationalArchive 

def readInAllGalaxyData():
    computational_params_archive = ComputationalArchive()
    #These two parameters will be used to define what profile you generate 
    #So they are here at the beginning to define 
    gal_of_interest = "fornax" #carina" "sculptor" "sextans"
    pop_of_interest = 'MR' # can be 'MP', 'IM', 'MR'
    
    #Define files for each of the galaxies to be read in
    fileDir = computational_params_archive.getDSphDataDir() 
    file_name_prefix = "stellarVelocitiesintheCarinaFornaxSculptorandSextansdSphGalaxiesDatafromtheMagellanMMFSSurvey_"
    carina_file_name = file_name_prefix + "CarinaDataFile_justData.txt"  
    fornax_file_name = file_name_prefix + "FornaxDataFile_justData.txt"  
    sculptor_file_name = file_name_prefix + "SculptorDataFile_justData.txt"  
    sextans_file_name = file_name_prefix + "SextansDataFile_justData.txt"
    file_names = {"carina" : carina_file_name,"fornax" : fornax_file_name,"sculptor" : sculptor_file_name,"sextans" : sextans_file_name}
    #print 'Files defined. '
    
    #Define the galaxy that you're interested in examining

    
    #Read the data into various variables 
    Target= []
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
    dist=[]


    file = open(fileDir + file_names[gal_of_interest], 'rb')
    for line in file:
        pm = float(line[93:98]) if line[97]!=' ' else 0
        #print pm #Population Membership
        if pm >=0.8:
            Target.append(line[0:8])
            date = line[14:22]
            Date.append(2450000+float(date))
            ra = (float(line[23:25])+ float(line[26:28])/60.0 + float(line[29:34])/3600.0)/24.0*360.0
            RA.append(ra)
            if line[35]=='-':
                dec=-1*float(line[36:38])-float(line[39:41])/60.0-float(line[42:46])/3600.0
            else:
                dec=float(line[36:38])+float(line[39:41])/60.0+float(line[42:46])/3600.0
            DEC.append(dec)
            vmag = float(line[47:52]) if line[51]!=' ' else 0
            Vmag.append(vmag)
            vimag = float(line[53:58]) if line[57]!=' ' else 0
            VImag.append(vimag)
            vhel = float(line[59:65]) if line[64]!=' ' else 0
            Vhel.append(vhel)
            vhelE = float(line[66:70]) if line[69]!=' ' else 0
            VhelE.append(vhelE)
            sigFe = float(line[71:76]) if line[75]!=' ' else 0
            SigFe.append(sigFe)
            sigFeE = float(line[77:81]) if line[80]!=' ' else 0
            SigFeE.append(sigFeE)    
            sigMg = float(line[82:87]) if line[86]!=' ' else 0
            SigMg.append(sigMg)
            sigMgE = float(line[88:92]) if line[91]!=' ' else 0
            SigMgE.append(sigMgE)
    #print DEC
    #print len(SigMg)
    Vega_solar_lums = 40.12 #solar luminosity of Vega
    #solar_lum = 3.825 * 10**23 #solar luminosity in kilowatts
    Vega_dist_pars = 7.68 #distance to Vega in parsecs
    #parsec_to_m = 3.086 * 10**16
    V_intensities = 10**(np.array(Vmag)/(-2.5))*(Vega_solar_lums)/(4.0*math.pi*Vega_dist_pars**2)
    #print 'Data from files read in. '
    
    #Show the V_intensity distribution (or don't) 
    binwidth=0.1*10**-9
    fig=plt.figure(1,figsize=(10,10))
    plt.hist(V_intensities,bins=np.arange(min(V_intensities), max(V_intensities) + binwidth, binwidth))
    plt.ylabel('number of stars')
    plt.xlabel('V intensities, compared to Vega')
    font = {'family': 'serif',
            'color':  'red',
            'weight': 'normal',
            'size': 16,
            'rotation': 90}
    #plt.text(0.04,0,'<--- lower metallicity cut',fontdict=font,ha='center',va='bottom' )
    #plt.text(0.57,0,'<--- upper metallicity cut',fontdict=font,ha='center',va='bottom' )
    plt.title(gal_of_interest + ' V intensities')
    pltFileDir="/Users/sasha/Documents/Harvard/physics/randall/dwarfGalaxyProject/plots/"
    #Uncomment the following line if you want to save the file
    #plt.savefig(pltFileDir + gal_of_interest + '_metallicity_histogram.png')
    
    #Uncomment the following lines if you want to display the figure 
    #plt.show()
    
     #To correct metallicities (according to Walker 2009), we need to apply a linear correction.
    # For that, we need the apparent magnitude of the horizontal branch stars (which they generously DO NOT provide for us)
    # SO I'm left to base VHB for:
    # Carina, Fornax, Sculptor from http://adsabs.harvard.edu/abs/2011ApJ...742...20W
    # Fornax on the Battaglia paper (fig 4), and I'm not sure what to do for the others
    # Could read off VHB for Carina from http://arxiv.org/pdf/astro-ph/0301388v1.pdf
    # Couldread off VHB for Sculptor from http://www.aanda.org/articles/aa/pdf/2011/04/aa16398-10.pdf
    # Attempt to read off VHB for Sextans from http://iopscience.iop.org/article/10.1086/379171/pdf
    V_HB = {"carina" : 20.9, "fornax" : 21.3,"sculptor" : 20.1,"sextans" : 20.25}
    SigMg_scale_corr = 0.079
    SigMg_corr = [SigMg[i] + SigMg_scale_corr * (Vmag[i] - V_HB[gal_of_interest]) for i in range(len(SigMg))]

    #We also will need distance from the sun, as measured in Amarisco et al.
    dist_from_sun = {"carina":105000,"fornax":147000,"sculptor":86000,"sextans":86000}
    
    #Show the metallicity distribution (or don't) 
    binwidth=0.02
    fig=plt.figure(1,figsize=(10,10))
    plt.hist(SigMg_corr,bins=np.arange(min(SigMg_corr), max(SigMg_corr) + binwidth, binwidth))
    plt.ylabel('number of stars')
    plt.xlabel('Corrected SigMg (Angstroms)')
    font = {'family': 'serif',
            'color':  'red',
            'weight': 'normal',
            'size': 16,
            'rotation': 90}
    #plt.text(0.04,0,'<--- lower metallicity cut',fontdict=font,ha='center',va='bottom' )
    #plt.text(0.57,0,'<--- upper metallicity cut',fontdict=font,ha='center',va='bottom' )
    plt.title(gal_of_interest + ' metallicities')
    pltFileDir="/Users/sasha/Documents/Harvard/physics/randall/dwarfGalaxyProject/plots/"
    #Uncomment the following line if you want to save the file
    #plt.savefig(pltFileDir + gal_of_interest + '_metallicity_histogram.png')
    
    #Uncomment the following lines if you want to display the figure 
    #plt.show()

    #We want to convert the RA/DEC coordiantes into more usable coordiantes (ie, just raw angular separation)
    #First, we need the center of the galaxies, in RA and DEC.  These are just hardcoded values pulled from papers
    car_ra= (float(6) + float(41)/60.0 + float(36.7)/3600.)/24.0*360.0
    car_dec=(float(-50) + float(-57)/60.0 + float(-58)/3600.)
    for_ra=(float(2) + float(39)/60.0 + float(52)/3600.)/24.0*360.0
    for_dec=(float(-34) + float(-26)/60.0 + float(-57)/3600.)
    scu_ra=(float(1) + float(0)/60.0 + float(9.3)/3600.)/24.0*360.0
    scu_dec=(float(-33) + float(-42)/60.0 + float(-33)/3600.)
    sex_ra=(float(10) + float(13)/60.0 + float(2.9)/3600.)/24.0*360.0
    sex_dec=(float(-1) + float(-36)/60.0 + float(-53)/3600.)

    ra_dec_tuples = {"carina" : (car_ra,car_dec),"fornax" : (for_ra,for_dec),"sculptor" : (scu_ra,scu_dec),"sextans" : (sex_ra,sex_dec)}
    
    
    #coordinate transformations pulled from stack exchange at this link: 
    # http://math.stackexchange.com/questions/92301/transforming-from-one-spherical-coordinate-system-to-another
    #print ra_dec_tuples[gal_of_interest][0]
    #print DEC
    cent_ra=ra_dec_tuples[gal_of_interest][0]
    cent_dec=ra_dec_tuples[gal_of_interest][1]
    corr_ra = np.array(RA) - cent_ra
    deg_to_rad = math.pi/180.0
    x1=np.cos(np.array(DEC) * deg_to_rad)*np.cos(corr_ra * deg_to_rad)
    x2=np.cos(np.array(DEC) * deg_to_rad)*np.sin(corr_ra * deg_to_rad)
    x3=np.sin(np.array(DEC) * deg_to_rad)

    x1bar=math.cos(cent_dec * deg_to_rad)*x1+math.sin(cent_dec * deg_to_rad)*x3
    x2bar=x2
    x3bar=-math.sin(cent_dec * deg_to_rad)*x1 + math.cos(cent_dec * deg_to_rad)*x3
    corr_ra = np.angle(x1bar + x2bar*1j, deg = True) 
    #print len(corr_ra)
    corr_dec = np.arcsin(x3bar) * 180.0/math.pi 
    radial_dists = (np.sqrt(corr_ra**2 + corr_dec**2)
                    *math.pi/180.0*dist_from_sun[gal_of_interest]) #physical distance of star from center, in pc

    #Show the radial distribution of star (or don't) 
    binwidth=50
    fig=plt.figure(1,figsize=(10,10))
    plt.hist(radial_dists,bins=np.arange(min(radial_dists), max(radial_dists) + binwidth, binwidth))
    plt.ylabel('number of stars')
    plt.xlabel('Angular distance from center of galaxy')
    #plt.text(-0.9,0,'<--- lower metallicity cut',fontdict=font,ha='center',va='bottom' )
    #plt.text(-0.6,0,'<--- upper metallicity cut',fontdict=font,ha='center',va='bottom' )
    plt.title(gal_of_interest + ' distance from gal center')
    pltFileDir="/Users/sasha/Documents/Harvard/physics/randall/dwarfGalaxyProject/plots/"
    #Uncomment the following line if you want to save the file
    #plt.savefig(pltFileDir + gal_of_interest + '_metallicity_histogram.png')
    
    #Uncomment the following lines if you want to display the figure 
    #plt.show()

    #Show the full galaxy of interest (or don't) 
    ylims = (-60,60) #ylimits from Walker plot
    xlims=(60,-60) #xlimits from Walker plot
    fig=plt.figure(1,figsize=(10,10))
    plt.scatter(corr_ra * 60.0,corr_dec * 60.0 ,s=2)
    plt.ylim(ylims)
    plt.xlim(xlims)
    plt.xlabel('ra sep')
    plt.ylabel('dec sep')
    plt.title('Fornax Stars Position on Sky')
    pltFileDir="/Users/sasha/Documents/Harvard/physics/randall/dwarfGalaxyProject/plots/"
    #Uncomment the following line if you want to save the file
    #plt.savefig(pltFileDir + gal_of_interest + '_all.png')
    
    #Uncomment the following lines if you want to display the figure 
    #plt.show()

    #Create the smoothed star distribution 
    ra_mesh,dec_mesh = np.meshgrid(np.arange(-60,60 + 3.0,3.0)[::-1],np.arange(-60,60 + 3.0,3.0))
    #print ra_mesh
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
    #print metallicity_cuts[gal_of_interest]
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
    #print len(corr_ra_one_pop)
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
    
    
    print 'Done reading in galaxy data! '
    
    


    
    
