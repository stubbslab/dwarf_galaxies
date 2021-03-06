
class ComputationalArchive:

    def getDSphDataDir(self):
        return self.dSph_data_directory 
    
    def getProjectDir(self):
        return self.project_directory

    def getSpotProbabilityDir(self):
        return self.spot_computed_prob_directory

    def getPlotDir(self):
        return self.plot_directory

    def getPotentialDir(self):
        return self.potential_directory

    def getMCMCOutDir(self):
        return self.MCMC_outputs_directory

    def getRandomizationDir(self):
        return self.randomization_directory 


    def __init__(self):
        self.project_directory = '/Users/sashabrownsberger/Documents/Harvard/physics/randall/'
        self.dSph_data_directory = self.project_directory + 'dSphDataFiles/'
        self.spot_computed_prob_directory = self.project_directory + 'probabilityTables/'
        self.plot_directory = self.project_directory + 'plots/'
        self.potential_directory = self.project_directory + 'potential_tables/'
        self.MCMC_outputs_directory = self.project_directory + 'MCMCOutputs/'
        self.randomization_directory = self.project_directory + 'randomizationResults/'
