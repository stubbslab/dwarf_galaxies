import math
import numpy as np

def getBinWalls(arcmin_limits,step): 
    if arcmin_limits[0] > 0.0 or arcmin_limits[1] < 0.0:
        return np.arange(arcmin_limits[0] - step / 2.0, arcmin_limits[1] + step ,step)
    else:
        return  np.unique(
                   np.around(
                       np.concatenate(
                          (-1 * np.arange( step / 2.0, abs(arcmin_limits[0]) + step / 2.0 ,step ),
                            np.arange( step / 2.0, arcmin_limits[1] + step / 2.0, step )            )
                                     ),
                            5)
                )

def getBinCenters(arcmin_limits,step):
    if arcmin_limits[0] > 0.0 or arcmin_limits[1] < 0.0:
        return np.arange(acrmin_limits[0] - step / 2.0, arcmin_limits[1] + step, step)
    else:
        return np.unique(
                   np.around(
                       np.concatenate(
                          (-1.0 * np.arange( 0.0, abs(arcmin_limits[0]) + step ,step ),
                            np.arange( 0.0, arcmin_limits[1] + step, step )            )
                                     ),
                            5)
                )
    

def getCombinedBinCenters(arcmin_limits, inner_lim, inner_step, outer_step):
    if inner_step > 0.0:
        inner = getBinCenters([lim * inner_lim for lim in arcmin_limits], inner_step)
        inner_extrema = [np.min(inner),np.max(inner)]
    else:
        inner = []
        in_extrema = [0.0,0.0]
    outer_all = getBinCenters(arcmin_limits, outer_step)
    outer = [outer_elem for outer_elem in outer_all if (outer_elem < inner_extrema[0] or outer_elem > inner_extrema[1])]
    
    return np.unique(np.around(np.concatenate( (inner,outer) ),5))

#Gives you values for R,z and los in measures of the scale radius for calculating the probability density on the sky.
# Desired limits must be given in arcminutes.  
def compRZLos(inner_lim,inner_step,outer_step,arcmin_limits_R,arcmin_limits_z,arcmin_limits_los,rs,dist):
    #print 'inner_lim = ' + str(inner_lim)
    #print 'inner_step = ' + str(inner_step)
    #print 'outer_step = ' + str(outer_step)
    #print 'arcmin_limits_R = '
    #print arcmin_limits_R
    #print 'arcmin_limits_z = '
    #print arcmin_limits_z
    #print 'arcmin_limits_los = '
    #print arcmin_limits_los
    #print 'rs = ' + str(rs)
    #print 'dist = ' + str(dist) 
    R_arcmin = getCombinedBinCenters(arcmin_limits_R, inner_lim, inner_step, outer_step)
    z_arcmin = getCombinedBinCenters(arcmin_limits_z, inner_lim, inner_step, outer_step)
    los_arcmin = getCombinedBinCenters(arcmin_limits_los, inner_lim, inner_step, outer_step)
    #los_bins_outer = np.arange(max(R_arcmin) * -1.0, max(R_arcmin) * 1.0 + outer_step, outer_step, dtype=np.float)
    #los_bins_inner = np.arange(max(R_arcmin) * -1.0 * inner_lim, max(R_arcmin) * 1.0 * inner_lim + iner_step, inner_step, dtype=np.float)
    #los_bins_arcmin = np.unique(np.around(np.concatenate( (los_bins_outer,los_bins_inner) ),5))
    #los_bins_arcmin = np.sort(los_bins_arcmin)
    #returns R,z,los_bins in scale radii 
    #los_bins_arcmin = R_arcmin 
    R=R_arcmin*1/60.0*1/180.0*math.pi*dist/rs
    z=z_arcmin*1/60.0*1/180.0*math.pi*dist/rs
    los=los_arcmin*1/60.0*1/180.0*math.pi*dist/rs
    return R,z,los 
