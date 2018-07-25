import numpy as np
import h5py 
import matplotlib.pyplot as plt
import pandas as pd
import os.path
import matplotlib.dates as mdates
from datetime import datetime
from datetime import date

    
def GenerateDbase(itera,nreaz,perm,th_cond):
#    print(perm)
    filename = "./pflotran/Dbase.h5"
    if itera == 0:
      h5file = h5py.File(filename,'w')
    else:
      h5file = h5py.File(filename,'r+')
      
    variables = []
    variables.append("Permeability")
    variables.append("ThermalConductivity")
    values = []
    values.append(perm[:])
    values.append(th_cond[:])

    for i in range(len(variables)):
        if h5file.get(variables[i]):
          del h5file[variables[i]]    
        h5dset = h5file.create_dataset(variables[i],data=values[i])
    h5file.close() 

def GenerateSimuEnsemble(nobs,obs_coord,z,nreaz,obs_time):
    obs_cell = np.zeros(nobs)
    ntime = len(obs_time)
    for i in range(nobs):
        obs_cell[i] = np.argmin(np.absolute(z-obs_coord[i]))
    obs_cell = obs_cell.astype(int)
    simu_ensemble = np.zeros((nobs*ntime,nreaz))
    for ireaz in range(nreaz):
        obs_temp = np.zeros(nobs*ntime)
        j = 0
        for itime in obs_time:
            h5f = h5py.File("./pflotran/1dthermalR{}.h5".format(ireaz+1),'r')
            grp_time = "Time:"+str(" %12.5E" % itime)+" s"
            dset_temp = "Temperature [C]"
            obs_temp[j*nobs:(j+1)*nobs] = h5f[grp_time][dset_temp][0][0][obs_cell]
            j = j+1
        simu_ensemble[:,ireaz] = obs_temp
        h5f.close()
    
    return simu_ensemble
            

    