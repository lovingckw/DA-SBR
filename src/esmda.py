# ------------Initialize Python packages to be used------------

import os
import numpy as np
import scipy.linalg as la
import shutil as shutil
from matplotlib.ticker import FormatStrFormatter
import matplotlib
import matplotlib.pyplot as plt
import math
import h5py
import util as util
import subprocess 
import random
import time
import datetime

#--------Read In User Specified Input------------------
time_start = datetime.datetime.now()
print('\n Reading in User Specified Parameters.')
finput = open("./src/user_specified_parameters.txt", 'r')
fpflotran = open("./dainput/1dthermal.in", 'r')
input_array = finput.readlines()
pflotranin = fpflotran.readlines()
finput.close()

ftest = open("./src/test.txt",'w')

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
    elif "th_cond_range" in line:
        new_line = line.replace('= [','= np.array([') +')'
        exec(new_line)
    exec(line)
print(' Done.')

subprocess.call("rm -rf ./pflotran",stdin=None, stdout=None,stderr=None,shell=True)
subprocess.call("mkdir ./pflotran",stdin=None, stdout=None,stderr=None,shell=True)

    
day_to_sec = 24*3600
nz = int(hz/dz) 
z = (np.arange(-hz,0,dz)+np.arange(-hz+dz,dz,dz))/2
init_logperm = np.random.normal(logperm_mean,logperm_sd,nreaz)
init_perm = np.exp(init_logperm)
init_hy_cond = init_perm*w_den*g*day_to_sec/w_vis

for x in range(nreaz):
    init_hy_cond[x] = random.randint(1,100)
init_perm = init_hy_cond*w_vis/w_den/g/day_to_sec

for iperm,vperm in enumerate(init_perm):
    if vperm < perm_range[0]: init_perm[iperm] = perm_range[0]
    if vperm > perm_range[1]: init_perm[iperm] = perm_range[1]
print(init_perm)
 
ftest.write("init_perm is: {} \n ".format(init_perm))
init_logperm = np.log(init_perm)
init_logperm_sd = np.std(init_logperm)
nobs = len(obs_coord)
ntime = np.shape(obs)[0]-1
obs_time = obs[:,0]
day_to_sec = 24*3600
pflotran_exe = "~/pflotran-cori"

perm = np.zeros((ntime,nreaz))
perm[0,:] = init_perm
th_cond = np.zeros((ntime,nreaz))
th_cond[0,:] = init_th_cond
simu_time = np.zeros((ntime,2))
kalman_gain_output = np.zeros((mw_length,ntime))
#dataset_path = './pflotran/dataset'
FNULL = open(os.devnull,'w')

for itime in range(0,ntime):
    simu_time[itime,:] = np.array([itime*da_interval,(itime+1)*da_interval])
    
    collect_start = 0
    if (itime>0):
        collect_start = simu_time[itime-1,1]
    
    collect_index = np.where((obs_time > collect_start) & (obs_time <= simu_time[itime,1]))
    collect_times = obs_time[collect_index]
    ncollect = len(collect_index)
    simu_ensemble = np.zeros((nobs*len(collect_times),nreaz))

    for itera in range(0,niter):
       util.GenerateDbase(itime,itera,nreaz,perm,th_cond)       
       util.MakePflotranInput(pflotranin,simu_time,itime,ncollect,collect_times,da_interval)
       
       bash_start = datetime.datetime.now()
       subprocess.call("./src/pflotran.sh {} {} {} ".format(nreaz,ncore,pflotran_exe),stdin=None, stdout=FNULL,stderr=None,shell=True)
       bash_end = datetime.datetime.now()
       
       simu_ensemble = util.GenerateSimuEnsemble(nobs,obs_coord,z,nreaz,collect_times)

       obs_sd = obs_sd_ratio*np.delete(np.squeeze(obs[collect_index,:],axis=0),0,1)
       obs_sd = obs_sd*(math.sqrt(alpha))

       obs_ensemble = np.repeat(np.transpose(np.delete(obs[collect_index],0,1)),nreaz,1)+np.diag(obs_sd.flatten('C'))@np.random.normal(0,1,nreaz*nobs*ncollect).reshape(nobs*ncollect,nreaz)       
       
       #update
       state_vector = np.zeros((2,nreaz))
       state_vector[0,:] = np.log(perm[itime,:])
       state_vector[1,:] = th_cond[itime,:]
       
       cov_state_simu = np.zeros((2,nobs))
       cov_state_simu = np.cov(state_vector,simu_ensemble)[0:2,2:]
       cov_simu = np.cov(simu_ensemble)

       if nobs == 1:
         inv_cov_simuAddobs = np.array([1/(cov_simu+np.square(np.diag(obs_sd)))])
       else:
         inv_cov_simuAddobs = la.inv(cov_simu+np.square(np.diag(obs_sd)))    
       
       kalman_gain = cov_state_simu@inv_cov_simuAddobs       
           
       state_vector = state_vector+kalman_gain@(obs_ensemble-simu_ensemble)

       perm[itime,:] = np.exp(state_vector[0,:]) # no +1
       perm[itime,:][perm[itime,:]>perm_range[1]] = perm_range[1]
       perm[itime,:][perm[itime,:]<perm_range[0]] = perm_range[0]
       th_cond[itime,:] = state_vector[1,:]
       th_cond[itime,:][th_cond[itime,:]>th_cond_range[1]] = th_cond_range[1]
       th_cond[itime,:][th_cond[itime,:]<th_cond_range[0]] = th_cond_range[0]  
            
       if itime+1 == ntime: break
 
       # prepare input for next timestep
       perm[itime+1,:] = perm[itime,:]
       th_cond[itime+1,:] = th_cond[itime,:]
    
       #disturb perm
       perm_temp = np.log(perm[itime+1,:])
       perm_temp = perm_temp+np.random.normal(0,np.sqrt(max(np.square(init_logperm_sd)-np.square(np.std(perm_temp)),0)),nreaz)
       perm[itime+1,:] = np.exp(perm_temp)
       perm[itime+1,:][perm[itime+1,:]>perm_range[1]] = perm_range[1]
       perm[itime+1,:][perm[itime+1,:]<perm_range[0]] = perm_range[0]
       
       th_cond[itime+1,:] = th_cond[itime,:]+np.random.normal(0,np.sqrt(max(np.square(th_cond_sd)-np.square(np.std(th_cond[itime+1,:])),0)),nreaz)
       th_cond[itime+1,:][th_cond[itime+1,:]>th_cond_range[1]] = th_cond_range[1]
       th_cond[itime+1,:][th_cond[itime+1,:]<th_cond_range[0]] = th_cond_range[0]  
              
np.savetxt("./figure/perm.txt",perm)
time_end = datetime.datetime.now()
time_cost = time_end-time_start
with open("timecost.txt", mode='w') as file:
    file.write('%s.\n'.format(timecost))
fpflotran.close()
