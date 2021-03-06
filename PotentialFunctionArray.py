#Define the PotentialFunctionArray object
#This object is used to calculate the potential value of dwarf galaxy profiles, and contains various properties of the simulated values of the gravitational potential.
#For a given gravitational potential, the class reads in the r-like and z-like positions of each point, as well as the value at the simulated point.  R-like must be first column, z-like must be second column, and pot value must be 3rd column
# (it does not matter how the underlying function that gave rise to the potential was defined; only that the file consists of an r-like and z-like coordinate and the value at that point).
#It then stores the range of the calculated r-like and z-like values, and interpolates over the grid of values.
#Note that we call these R-like and z-like axes, since they need not actually correspond to true R and z coordinates.  They just need to be some coordinate along the potential's axis of symmetry, and one perpendicular to the axis of symmetry that are related to tru z and R by some constant scaling (that can be different for each).
#When querying values from the interpolating funciton, they must be entered as (R-like value, z-like value).  

#For large files, this can be time consuming, and should thus be done as infrequently as possible
import csv 
import numpy as np 
from scipy.interpolate import RegularGridInterpolator
import time

class PotentialFunctionArray:
    def __init__(self,potential_file):
        print ('In potential function array, working with potential_file = ' + potential_file )
        Rlikeaxis=[]
        zlikeaxis=[]
        potential_function_values=[]
        Rlike_values=[]
        zlike_values=[]
        potential_function_dict={} 
        with open(potential_file) as csvfile:
            myreader=csv.reader(csvfile,delimiter=',')
            for row in myreader:
                #print 'row = ' + str(row) 
                potential_function_dict[(float(row[0]),float(row[1]))]=float(row[2])
                potential_function_values=potential_function_values + [float(row[2])]
                Rlike_values = Rlike_values + [float(row[0])]
                zlike_values = zlike_values + [float(row[1])]
                if not float(row[0]) in Rlikeaxis:
                    Rlikeaxis=Rlikeaxis+[float(row[0])]
                    
                if not float(row[1]) in zlikeaxis:
                    zlikeaxis=zlikeaxis+[float(row[1])]
        Rlikeaxis.sort()
        zlikeaxis.sort()
        potential_function_array=np.array([Rlike_values,zlike_values,potential_function_values]).transpose()
        #potential_mesh=np.zeros((len(zlikeaxis),len(Rlikeaxis)))
        potential_mesh = np.reshape(potential_function_values,(len(Rlikeaxis),len(zlikeaxis)))
        #for i in range(np.shape(potential_mesh)[0]):
        #    for j in range (np.shape(potential_mesh)[1]):
        #        potential_mesh[i][j]=potential_function_array[i*np.shape(potential_mesh)[1]+j,2]
        self.potential_mesh=potential_mesh
        self.interpolating_function = RegularGridInterpolator((Rlikeaxis,zlikeaxis),potential_mesh,method="linear")
        self.potential_function_array=potential_function_array
        self.Rlikeaxis=np.array(Rlikeaxis)
        self.zlikeaxis=np.array(zlikeaxis)
        self.potential_function_values=potential_function_values
        self.potential_function_dict=potential_function_dict
        self.Rlike_range=[np.min(self.Rlikeaxis),np.max(self.Rlikeaxis)]
        self.zlike_range=[np.min(self.zlikeaxis),np.max(self.zlikeaxis)]
