# ------------Initialize Python packages to be used------------

import os
import numpy as np
import scipy as sp
import shutil as shutil
from matplotlib.ticker import FormatStrFormatter
import matplotlib
import matplotlib.pyplot as plt
import math
import h5py
import util as util
import subprocess 

#--------Read In User Specified Input------------------

print('\n Reading in User Specified Parameters.')
finput = open("./src/user_specified_parameters.txt", 'r')
fpflotran = open("./dainput/1dthermal.in", 'r')
input_array = finput.readlines()
pflotranin = fpflotran.readlines()
finput.close()

# Run lines of text file to define variables
for line in input_array:
    if "obs_coord" in line:
        new_line = line.replace('= [','= np.array([') + ')' 
        exec(new_line)
    elif "path_to_obs_data" in line:
        new_line = line.replace('path_to_obs_data = ','obs = np.loadtxt(') + ')'
        exec(new_line)
    elif "perm_range" in line:
        new_line = line.replace('= [','= np.array([') +')'
        exec(new_line)
    exec(line)
print(' Done.')

subprocess.call("rm -rf ./pflotran",stdin=None, stdout=None,stderr=None,shell=True)
subprocess.call("mkdir ./pflotran",stdin=None, stdout=None,stderr=None,shell=True)
for ireaz in range(nreaz):
    subprocess.call("mkdir ./pflotran/%d/" % ireaz,stdin=None, stdout=None,stderr=None,shell=True)
    
day_to_sec = 24*3600
nz = int(hz/dz) # number of grid blocks 
z = (np.arange(-hz,0,dz)+np.arange(-hz+dz,dz,dz))/2
init_hy_cond = np.power(np.random.normal(hy_cond_mean,hy_cond_sd,nreaz),10)
init_perm = init_hy_cond*w_vis/w_den/g/day_to_sec
init_perm_sd = init_hy_cond*np.std(np.log10(init_perm))
nobs = len(obs_coord)
ntime = np.shape(obs)[0]-1
obs_time = obs[:,0]
day_to_sec = 24*3600
pflotran_exe = "~/pflotran-cori"

perm = np.zeros((ntime,nreaz))
perm[0,:] = init_perm
simu_time = np.zeros((ntime,2))

for itime in range(0,ntime):
    print("itime is %d" % itime)
    if itime == 0:
        perm_update_index = np.array([0])
    else:
        perm_update_index = np.arange(max(0,itime-mw_length+1),itime)
    print("perm_update_index: ",perm_update_index)
#    print(perm_update_index)
    if itime == 0:
        simu_time[itime,:] = np.array([0,da_interval])
    else:
        simu_time[itime,:] = np.array([perm_update_index[0]*da_interval,(itime+1)*da_interval])
#    print(simu_time[itime-1,:])
    collect_start = 0
    if (itime>0):
        collect_start = simu_time[itime-1,1]
    
    collect_index = np.where((obs_time > collect_start) & (obs_time <= simu_time[itime,1]))
    collect_times = obs_time[collect_index]
    ncollect = len(collect_index)
    simu_ensemble = np.zeros((nreaz,nobs*len(collect_times)))

    for iter in range(0,niter):
       util.MakePflotranInput(pflotranin,simu_time,itime,ncollect,collect_times,perm_update_index,nreaz,da_interval,perm)

       subprocess.call("./src/pflotran.sh {} {} {} ".format(nreaz-1,ncore,pflotran_exe),stdin=None, stdout=None,stderr=None,shell=True)
        
       simu_ensemble = util.GenerateSimuEnsemble(nobs,obs_coord,z,nreaz,collect_times)
       
       obs_sd = obs_sd_ratio*np.delete(np.squeeze(obs[collect_index,:],axis=0),0,1)
       obs_sd = obs_sd*(math.sqrt(alpha))
       obs_ensemble = np.repeat(np.delete(obs[collect_index],0,1),nreaz,0)+np.matmul(np.random.normal(0,1,nreaz*nobs*ncollect).reshape(nreaz,nobs*ncollect),np.diag(np.squeeze(obs_sd)))
       
       #update
       state_vector = np.zeros((nreaz,perm_update_index.shape[0]))
       loc_j = 0
       for jtime in perm_update_index:
           state_vector[:,loc_j] = np.log(perm[jtime,:])
           loc_j = loc_j+1
       
       print("state_vector size: ",state_vector)
       print("simu_ensemble size: ", simu_ensemble)    
       cov_state_simu = np.cov(state_vector,simu_ensemble)
       cov_simu = np.cov(simu_ensemble,simu_ensemble)
       inv_cov_simuAddobs = sp.linalg.inv(cov_simu+np.square(np.diag(np.squeeze(obs_sd))))
       kalman_gain = np.matmul(cov_state_simu,inv_cov_simuAddobs)
       state_vector = state_vector+np.matmul(obs_ensemble-simu_ensemble,np.transpose(kalman_gain))
       
       for jtime in perm_update_index:
           perm[jtime,:] = np.exp(state_vector[:,jtime-perm_update_index[0]]) # no +1
           perm[jtime,:][perm[jtime,:]>perm_range[1]] = perm_range[1]
           perm[jtime,:][perm[jtime,:]<perm_range[0]] = perm_range[0]
       
       # prepare input for next timestep
       perm[itime+1,:] = perm[itime,:]
       
       #disturb perm
       perm_temp = np.log(perm[itime+1,:])
       perm_temp = perm_temp+np.random.normal(nreaz,0,np.sqrt(max(np.square(init_perm_sd)-np.square(np.std(perm_temp)),0)))
       perm[itime+1,:] = np.exp(perm_temp)
       perm[itime+1,:][perm[itime+1,:]>perm_range[1]] = perm_range[1]
       perm[itime+1,:][perm[itime+1,:]<perm_range[0]] = perm_range[0]

np.savetxt("/results/perm.txt",perm)
fpflotran.close()