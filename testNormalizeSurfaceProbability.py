import numpy as np

def normalizeSurfaceProbability(R_axis, z_axis, surface_prob):

    Rmesh, zmesh = np.meshgrid(R_axis,z_axis)
    normalization_constant = 0.0
    for i in range(np.shape(surface_prob)[0]):
        for j in range(np.shape(surface_prob)[1]):
            prob_at_point = surface_prob[i][j]
            print 'surface_prob[i][j] = ' + str(surface_prob[i][j]) 
            if i == 0:
	        #print 'i = 0'
                #print 'zmesh[i+1][j]-zmesh[i][j] = ' + str(zmesh[i+1][j]-zmesh[i][j])
                prob_at_point = prob_at_point * (zmesh[i+1][j]-zmesh[i][j])
            elif i == np.shape(surface_prob)[0]-1:
	        #print 'i = ' + str(i)
                #print 'zmesh[i][j]-zmesh[i-1][j] = ' + str(zmesh[i][j]-zmesh[i-1][j])
                prob_at_point = prob_at_point * (zmesh[i][j]-zmesh[i-1][j])
            else:
                #print ' i = ' + str(i)
                #print 'zmesh[i+1][j]-zmesh[i-1][j] = ' + str(zmesh[i+1][j]-zmesh[i-1][j])
                prob_at_point = prob_at_point * (zmesh[i+1][j]-zmesh[i-1][j]) / 2.0
            if j == 0:
	        #print ' j = ' + str(j)
                #print 'Rmesh[i][j+1]-Rmesh[i][j] = ' + str(Rmesh[i][j+1]-Rmesh[i][j])
                prob_at_point = prob_at_point * (Rmesh[i][j+1]-Rmesh[i][j])
            elif j == np.shape(surface_prob)[1]-1:
	        #print ' j = ' + str(j)
                #print ' Rmesh[i][j]-Rmesh[i][j-1] = ' + str(Rmesh[i][j]-Rmesh[i][j-1])
                prob_at_point = prob_at_point * (Rmesh[i][j]-Rmesh[i][j-1])
            else:
	        #print ' j = ' + str(j)
                #print 'Rmesh[i][j+1]-Rmesh[i][j-1] = ' + str(Rmesh[i][j+1]-Rmesh[i][j-1])
                prob_at_point = prob_at_point * (Rmesh[i][j+1]-Rmesh[i][j-1]) / 2.0
            print 'prob_at_point = ' + str(prob_at_point)
            print 'area = ' + str(prob_at_point / surface_prob[i][j]) 
            normalization_constant = normalization_constant + prob_at_point

    return  surface_prob / normalization_constant 
