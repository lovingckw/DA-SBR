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
subprocess.call("cp ./dainput/1dthermal.in ./pflotran/",stdin=None, stdout=None,stderr=None,shell=True)
    
day_to_sec = 24*3600
nz = int(hz/dz) 
z = (np.arange(-hz,0,dz)+np.arange(-hz+dz,dz,dz))/2
init_logperm = np.random.normal(logperm_mean,logperm_sd,nreaz)
init_perm = np.exp(init_logperm)
init_hy_cond = init_perm*w_den*g*day_to_sec/w_vis
init_th_cond = np.random.normal(th_cond_mean,th_cond_sd,nreaz)
 
nobs = len(obs_coord)
ntime = np.shape(obs)[0]
obs_time = obs[:,0]
obs_temp = obs[:,1:]
np.savetxt("./figure/obs_temp.txt",obs_temp)
day_to_sec = 24*3600
pflotran_exe = "~/pflotran-cori"

perm = np.zeros((niter,nreaz))
perm[0,:] = init_perm
th_cond = np.zeros((niter,nreaz))
th_cond[0,:] = init_th_cond
simu_time = np.array([obs_interval,nobs*obs_interval])
FNULL = open(os.devnull,'w')
simu_ensemble = np.array((nobs*ntime,nreaz))

avg_perm = list()
avg_th_cond = list()

for i in range(0,niter):
  print(i)
  util.GenerateDbase(i,nreaz,perm[i,:],th_cond)

  subprocess.call("./src/pflotran.sh {} {} {} ".format(nreaz,ncore,pflotran_exe),stdin=None,stdout=FNULL,stderr=None,shell=True)

  simu_ensemble = util.GenerateSimuEnsemble(nobs,obs_coord,z,nreaz,obs_time)
  np.savetxt("./figure/simu_ensemble{}.txt".format(i),simu_ensemble)
  
  obs_sd = math.sqrt(alpha)*obs_sd_ratio*obs_temp
  obs_ensemble = np.repeat(obs_temp.flatten('C').reshape(ntime*nobs,1),nreaz,1)+np.diag(obs_sd.flatten('C'))@np.random.normal(0,1,nreaz*nobs*ntime).reshape(nobs*ntime,nreaz)

#update
  state_vector = np.zeros((2,nreaz))
  state_vector[0,:] = np.log(perm[i,:])
  state_vector[1,:] = th_cond[i,:]

  cov_state_simu = np.zeros((2,nobs))
  cov_state_simu = np.cov(state_vector,simu_ensemble)[0:2,2:]      

  cov_simu = np.cov(simu_ensemble)
   
  if nobs*ntime == 1:
    inv_cov_simuAddobs = np.array(1/(cov_simu+np.square(np.diag(obs_sd))))
  else:
    inv_cov_simuAddobs = la.inv(cov_simu+np.square(np.diag(obs_sd.flatten('C'))))  
     
  kalman_gain = cov_state_simu@inv_cov_simuAddobs            
  state_vector = state_vector+kalman_gain@(obs_ensemble-simu_ensemble)

  perm[i,:] = np.exp(state_vector[0,:]) # no +1
  th_cond[i,:] = state_vector[1,:] 

  avg_perm.append(np.exp(np.mean(state_vector[0,:])))  
  avg_th_cond.append(np.mean(state_vector[1,:]))
    
  print(avg_perm[i])            
  print(avg_th_cond[i])
  
  if i<niter-1:
    perm[i+1,:] = perm[i,:]
    th_cond[i+1,:] = th_cond[i,:]

np.savetxt("./figure/avg_perm.txt",avg_perm)
np.savetxt("./figure/avg_th_cond.txt",avg_th_cond)
np.savetxt("./figure/perm.txt",perm)
np.savetxt("./figure/th_cond.txt",th_cond)
np.savetxt("./figure/init_perm.txt",init_perm)
np.savetxt("./figure/init_th_cond.txt",init_th_cond) 
        
time_end = datetime.datetime.now()
time_cost = time_end-time_start
with open("timecost.txt", mode='w') as file:
    file.write('%s.\n'.format(time_cost))
fpflotran.close()
