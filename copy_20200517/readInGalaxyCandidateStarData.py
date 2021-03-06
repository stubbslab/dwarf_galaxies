#Define a function to read in dSph data from a file.
#The current format is that given by Walker et al.  It should be updatable to accomodate other files, but that would require a good deal of hard coded change.

import matplotlib
import matplotlib.pyplot as plt
import csv
import math
import numpy as np
from DwarfGalDataArchive import DwarfGalDataArchive
from AstronomicalParameterArchive import AstronomicalParameterArchive

def readInGalaxyCandidateStarData(galaxy):
    data_archive = DwarfGalDataArchive()
    data_file = data_archive.getCandidateFile(galaxy)
    #print 'Pulling dSph data from file ' + data_file
    #Get distance from sun in parsecs. 
    dist_from_sun = data_archive.getDistanceFromSun([galaxy, 'dummy_population']) #{"carina":105000,"fornax":147000,"sculptor":86000,"sextans":86000} in pc
    dist_from_sun_err = data_archive.getDistanceFromSunErr([galaxy, 'dummy_population'])
    total_mass = data_archive.getTotalMass([galaxy, 'dummy_population']) #{"carina":1.28*10**8,"fornax":1.28*10**8,"sculptor":1.28*10**8,"sextans":1.28*10**8} in M_sun 
    astro_params = AstronomicalParameterArchive() 
    
    #Read the data into various variables 
    xs_degrees = []
    x_index = 0
    ys_degrees = []
    y_index = 1
    Rs_degrees = []
    R_index = 2

    with open(data_file) as csvfile:
        myreader = csv.reader(csvfile, delimiter = ' ')
        for row in myreader:
            row = [elem for elem in row if elem !='']
            xs_degrees = xs_degrees + [float(row[x_index]) / 60.0]
            ys_degrees = ys_degrees + [float(row[y_index]) / 60.0]
            Rs_degrees = Rs_degrees + [float(row[R_index]) / 60.0]

    return xs_degrees, ys_degrees, Rs_degrees
