#Here we define a python object that basically includes all of the information on a dSph observation.
#It is a repository class; the data is actually read in using another function. 

from readInGalaxyCandidateStarData import readInGalaxyCandidateStarData
import math

class CandidateGalaxyStarData:
    def __init__(self, galaxy):
        #Returns the positions of the candidate stars in standard coordinates (arcminutes of separation projected onto plane on sky) 
        
        ( self.xs_degrees,  self.ys_degrees,  self.Rs_degrees
          ) = readInGalaxyCandidateStarData(galaxy)
        
     
