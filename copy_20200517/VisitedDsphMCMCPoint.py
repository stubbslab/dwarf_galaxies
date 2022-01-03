#Defines the VisitedDsphMCMCPoint class.
#This class is used to keep track of each point visited in the parameter space during each MCMC chain.
#It is composed of a few pieces :
#  the parameter storer, which holds all of the parameters
#  the log likelihood of the visited point (based on the normalized surface brightness profile, given the parameters used)
#  the number of times the algorithm stayed at the point (starting at 1, and increasing by 1 every time the algorithm doesn't move from a point)
#Then when, the MCMC series is run, the output consists of an array of instances of this class.  

from DwarfGalaxyParametersStorer import DwarfGalaxyParametersStorer

class VisitedDsphMCMCPoint:

    def __init__(self, parameter_storer,log_likelihood, n_visits):
        self.parameter_storer = parameter_storer
        self.log_likelihood = log_likelihood
        self.n_visits = n_visits
        
