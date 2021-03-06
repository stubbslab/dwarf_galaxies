import math
import numpy as np
from PotentialArchive import PotentialArchive
from PotentialFunctionArray import PotentialFunctionArray
from FullPotential import FullPotential
from AstronomicalParameterArchive import AstronomicalParameterArchive
from DwarfGalDataArchive import DwarfGalDataArchive
from logList import logList 
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

def showGravitationalPotential(galaxy, rs, el, lam, eps, zeta, c, R_halo, z_halo, R_disk, z_disk, disk_funct = 'None', halo_funct = 'None'):
    e = math.sqrt(1-el**2)
    #f = math.log(1+c) - c/(1+c)
    data_archive = DwarfGalDataArchive()
    astro_archive = AstronomicalParameterArchive()
    potential_archive = PotentialArchive()
    disk_file = potential_archive.getDiskFile(lam)
    halo_file = potential_archive.getHaloFile(el)
    #doesn't matter what population in the galaxy we pick.  We just need to pick one so the archival retrieval functions work properly. 
    population = [galaxy,'MP']
    M = data_archive.getTotalMass(population)  
    if disk_funct == 'None':
        disk_funct = PotentialFunctionArray(disk_file)
    if halo_funct == 'None': 
        halo_funct = PotentialFunctionArray(halo_file)
    gamma = astro_archive.getGamma()
    potentialToPlot = FullPotential(gamma, M, rs, eps, c, lam, zeta, el, disk_interpolating_function = disk_funct, halo_interpolating_function = halo_funct)

    Rmesh_halo, zmesh_halo = np.meshgrid(R_halo,z_halo)
    Rmesh_disk, zmesh_disk = np.meshgrid(R_disk,z_disk) 

    potential_mesh = potentialToPlot.full_value(np.abs(Rmesh_halo),np.abs(zmesh_halo),np.abs(Rmesh_disk),np.abs(zmesh_disk))

    #print potential_mesh
    im=plt.figure(figsize=(9,7))
    print 'potential_mesh = '
    print potential_mesh 
    log_levels=logList(np.min(-1*potential_mesh),np.max(-1*potential_mesh),20)
    print log_levels
    #print Rmesh_halo
    #print zmesh_halo
    CS=plt.contour(Rmesh_halo,zmesh_halo,-1*potential_mesh,levels = log_levels)
    CB_contour=plt.colorbar(shrink=0.8,extend='both',format='%.2e')
    plt.show()
    

    return 0 
    
        
