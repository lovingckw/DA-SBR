import numpy as np
import h5py 
import matplotlib.pyplot as plt
import pandas as pd
import os.path
import matplotlib.dates as mdates
from datetime import datetime
from datetime import date

#def GeneratePermDataset(itime,nreaz,perm,nz,dataset_path):
##    for jtime in perm_update_index:
#    name_file = "{}/Perm{}.h5".format(dataset_path,itime+1)
#    if os.path.isfile(name_file):
#        os.remove(name_file) 
#    h5f = h5py.File(name_file,'w') 
#    dset_id = h5f.create_dataset('Cell Ids',(nz,),dtype='i')
#    dset_id[:] = range(1,nz+1) 
#    for kreaz in range(nreaz):
#        name_dset = "Permeability{}".format(kreaz+1)
#        dset_perm = h5f.create_dataset(name_dset,(nz,),dtype='f')
#        dset_perm[:] = perm[itime,kreaz]
#    h5f.close()
#    return            
#            
#def GenerateThermalcondDataset(itime,nreaz,th_cond,nz,dataset_path):
##    for jtime in perm_update_index:
#    name_file = "{}/Th_cond{}.h5".format(dataset_path,itime+1)
#    if os.path.isfile(name_file):
#        os.remove(name_file) 
#    h5f = h5py.File(name_file,'w') 
#    dset_id = h5f.create_dataset('Cell Ids',(nz,),dtype='i')
#    dset_id[:] = range(1,nz+1) 
#    for kreaz in range(nreaz):
#        name_dset = "THERMAL_CONDUCTIVITY_WET{}".format(kreaz+1)
#        dset_perm = h5f.create_dataset(name_dset,(nz,),dtype='f')
#        dset_perm[:] = th_cond[itime,kreaz]
#    h5f.close()
#    return  
    
def GenerateDbase(itime,itera,nreaz,perm,th_cond):
    filename = "./pflotran/Dbase.h5"
    if itime == 0:
      h5file = h5py.File(filename,'w')
    else:
      h5file = h5py.File(filename,'r+')
      
    variables = []
    variables.append("Permeability{}".format(itime))
    variables.append("ThermalConductivity{}".format(itime))
    values = []
    values.append(perm[itime,:])
    values.append(th_cond[itime,:])
    if itera == 0:
      for i in range(len(variables)):
        h5dset = h5file.create_dataset(variables[i],data=values[i])
    else:
      for i in range(len(variables)):
        if h5file.get(variables[i]):
          del h5file[variables[i]]
        h5dset = h5file.create_dataset(variables[i],data=values[i])
    h5file.close() 
            
def MakePflotranInput(pflotranin,simu_time,itime,ncollect,collect_times,da_interval):
    
    # determin whether to use restart
    restart_lindex = [i for i, s in enumerate(pflotranin) if "RESTART" in s][0]
    restart_card = ["  RESTART \n",
                       "    FILENAME 1dthermal-restart.chk \n",
                       "    REALIZATION_DEPENDENT \n",
		       "/ \n"]
    if simu_time[itime,0] == 0:
        pflotranin = pflotranin[0:restart_lindex]+pflotranin[restart_lindex+1:]
    else:
        pflotranin = pflotranin[0:restart_lindex]+restart_card+pflotranin[restart_lindex+1:]
    
    # write the end time of simulation
    finaltime_lindex = [i for i, s in enumerate(pflotranin) if "FINAL_TIME" in s][0]
    pflotranin[finaltime_lindex] = "  FINAL_TIME "+ np.array_str(simu_time[itime,1])+" sec"+"\n"
    
    # generate observation 
    snapshot_lindex =  [i for i, s in enumerate(pflotranin) if "SNAPSHOT_FILE" in s][0]
    obs_card = "    TIMES sec   "
    npl = 3 
    nline = np.floor_divide(ncollect,npl)
    if nline == 0 or nline == 1:
        obs_card = obs_card+" "+np.array_str(collect_times[0])+"\n"
    else:
        for iline in range(1,nline):
            obs_card = obs_card+np.array_str(collect_times[((iline-1)*npl+1):(iline*npl)])+"\\"
        obs_card = obs_card+" "+np.array_str(collect_times[(iline*npl+1):ncollect])+"\n"
    
    # add observation time
    pflotranin = pflotranin[0:snapshot_lindex+1]+[obs_card]+pflotranin[snapshot_lindex+2:]     
    
    permiso_lindex = [i for i, s in enumerate(pflotranin) if "PERM_ISO" in s][0]
    pflotranin[permiso_lindex] = "    PERM_ISO DBASE_VALUE Permeability{} \n".format(itime)
    
    thercondwet_lindex = [i for i, s in enumerate(pflotranin) if "THERMAL_CONDUCTIVITY_WET" in s][0]
    pflotranin[thercondwet_lindex] = "  THERMAL_CONDUCTIVITY_WET DBASE_VALUE ThermalConductivity{} \n".format(itime)
    
#    #prepaer perm_dataset_card
#    perm_dataset_card_lindex = [i for i, s in enumerate(pflotranin) if "perm_dataset" in s][0]+1
#    perm_dataset_card = ["DATASET Perm0\n",
#                    "  FILENAME ./pflotran/dataset/Perm0.h5\n",
#                    "  HDF5_DATASET_NAME Permeability \n",
#                    "  REALIZATION_DEPENDENT\n",
#                    "/ \n",
#                    "\n"]
#    perm_dataset_card_length = len(perm_dataset_card)
#    perm_dataset_card_new = list()
#    
#    #generate new dataset card
##    for jtime in perm_update_index:
#    perm_dataset_card[0] = "DATASET Perm{}\n".format(itime+1)
#    perm_dataset_card[1] = "  FILENAME ./dataset/Perm{}.h5\n".format(itime+1)
#    perm_dataset_card_copy = perm_dataset_card[:]
#    perm_dataset_card_new = perm_dataset_card_new+perm_dataset_card_copy
#        
#    pflotranin = pflotranin[0:perm_dataset_card_lindex]+perm_dataset_card_new+pflotranin[perm_dataset_card_lindex+1:]
#    
#    #prepaer th_cond_dataset_card
#    th_cond_dataset_card_lindex = [i for i, s in enumerate(pflotranin) if "th_cond_dataset" in s][0]+1
#    th_cond_dataset_card = ["DATASET Th_cond0\n",
#                    "  FILENAME ./pflotran/dataset/Th_cond0.h5\n",
#                    "  HDF5_DATASET_NAME THERMAL_CONDUCTIVITY_WET \n",
#                    "  REALIZATION_DEPENDENT\n",
#                    "/ \n",
#                    "\n"]
#    th_cond_dataset_card_length = len(th_cond_dataset_card)
#    th_cond_dataset_card_new = list()
#    
#    #generate new dataset card
##    for jtime in perm_update_index:
#    th_cond_dataset_card[0] = "DATASET Th_cond{}\n".format(itime+1)
#    th_cond_dataset_card[1] = "  FILENAME ./dataset/Th_cond{}.h5\n".format(itime+1)
#    th_cond_dataset_card_copy = th_cond_dataset_card[:]
#    th_cond_dataset_card_new = th_cond_dataset_card_new+th_cond_dataset_card_copy
#        
#    pflotranin = pflotranin[0:th_cond_dataset_card_lindex]+th_cond_dataset_card_new+pflotranin[th_cond_dataset_card_lindex+1:]
        
    #prepare strata_card
    strata_card_lindex = [i for i, s in enumerate(pflotranin) if "STRATA" in s][0]-1
    strata_card = ["STRATA\n",
                   "  REGION all\n",
                   "  MATERIAL Alluvium\n",
                   "  START_TIME 0 sec\n",
                   "  FINAL_TIME 0 sec\n",
                   "END\n",
                   "\n"]
    strata_card_length = len(strata_card)
    strata_card_new = list()
    
    #generate new strata card
#    for jtime in perm_update_index:
#    strata_card[2] = "  MATERIAL Perm"+str(itime+1)+"\n"
#        if jtime == perm_update_index[0]:
#            strata_card[3] = "  START_TIME "+"0 sec"+"\n"
#        else:
#            strata_card[3] = "  START_TIME "+str(da_interval*jtime)+" sec"+"\n"
    strata_card[4] = "  FINAL_TIME "+str(da_interval*(itime+1))+" sec"+"\n"
    strata_card_copy = strata_card[:]
    strata_card_new = strata_card_new+strata_card_copy
        
    pflotranin = pflotranin[0:strata_card_lindex]+strata_card_new+pflotranin[strata_card_lindex+strata_card_length:]
    
#    # prepare perm_card
#    material_card_length = 16
#    material_card_lindex = [i for i, s in enumerate(pflotranin) if "MATERIAL_PROPERTY Alluvium" in s][0]
#    material_card = pflotranin[material_card_lindex:material_card_lindex+material_card_length]
#    perm_value_cardindex = [i for i, s in enumerate(material_card) if "PERM_ISO" in s][0]
#    th_cond_value_cardindex = [i for i, s in enumerate(material_card) if "THERMAL_CONDUCTIVITY_WET" in s][0]
#    material_card_new = list()
#    
#    for ireaz in range(nreaz):
#        del material_card_new[:]
##        for jtime in perm_update_index:
#        material_card[0] = "MATERIAL_PROPERTY Perm"+str(itime+1)+"\n"
#        material_card[1] = "  ID "+str(itime+1)+"\n"
#        material_card[perm_value_cardindex] = "    DATASET Perm{} \n".format(itime+1)
#        material_card[th_cond_value_cardindex] = "  THERMAL_CONDUCTIVITY_WET \n"
#        material_card[th_cond_value_cardindex+1] = "    DATASET Th_cond{} \n".format(itime+1)
#        material_card[th_cond_value_cardindex+2] = "/ \n"
#        material_card_copy = material_card[:] 
#        material_card_new = material_card_new+material_card_copy
#            
#        temp_pflotranin = pflotranin[0:material_card_lindex]+material_card_new+pflotranin[material_card_lindex+material_card_length:]
    new_pflotranin = open("./pflotran/1dthermal.in",'w')
    new_pflotranin.writelines(pflotranin)  
    new_pflotranin.close()
    return

def GenerateSimuEnsemble(nobs,obs_coord,z,nreaz,collect_times):
    obs_cell = np.zeros(nobs)
    for i in range(nobs):
        obs_cell[i] = np.argmin(np.absolute(z-obs_coord[i]))
    obs_cell = obs_cell.astype(int)
    
    simu_ensemble = np.zeros((nobs*len(collect_times),nreaz))
    for ireaz in range(nreaz):
        obs_temp = np.zeros(nobs*len(collect_times))
        jj = 0
        for collect_itime in collect_times:
            h5f = h5py.File("./pflotran/1dthermalR{}.h5".format(ireaz+1),'r')
            grp_time = "Time:"+str(" %12.5E" % collect_itime)+" s"
            dset_temp = "Temperature [C]"
            obs_temp[jj*nobs:(jj+1)*nobs] = h5f[grp_time][dset_temp][0][0][obs_cell]
            jj = jj+1
        simu_ensemble[:,ireaz] = obs_temp
        h5f.close()
    
    return simu_ensemble
            
def PlotFlux(path_to_perm,path_to_true_flux):
    perm = np.loadtxt('../results/perm.txt',dtype=float)
    hy_cond = perm*1e6*9.8*3600*24
    true_flux = pd.read_csv('../dainput/true_flux.csv',sep=',',header=0)
    true_flux['time'] = [datetime.strptime(x,"%m/%d/%y %H:%M") for x in true_flux['time']]
    dt = true_flux['time']
    dt = [mdates.date2num(x) for x in dt]
    ntime = perm.shape[0]
    nreaz = perm.shape[1] 
    
    plt.figure()    
    sp1 = plt.subplot(211)
    line1, = plt.plot(dt[0:ntime],hy_cond[0:ntime,0],'orange',linewidth=0.25)


    
    for i in range(1,nreaz):       
        plt.plot(dt[0:ntime],hy_cond[0:ntime,i],'orange',linewidth=0.25)
    line2, = plt.plot(dt[0:ntime],np.mean(hy_cond[0:ntime,:],axis=1),'k-',linewidth=2)

    plt.ylabel("Hydraulic conductivity (m/d)")
    plt.xlim([date(2017,4,1),date(2017,5,1)])
    plt.ylim([0,60])    
    plt.setp(sp1.get_xticklabels(),visible=False)
    plt.legend((line1,line2),('ES-MDA (realization)','ES-MDA (mean)'))
    plt.legend(frameon=False)

    ax = plt.gca()
    daysFmt = mdates.DateFormatter("%m/%d/%y")
    days = mdates.DayLocator(interval=5)
    ax.xaxis.set_major_locator(days)
    ax.xaxis.set_major_formatter(daysFmt)
    ax.xaxis.set_minor_locator(mdates.DayLocator())
#    ax.set_xlabel("Date")
#    ax.set_ylabel("Hydraulic conductivity (m/d)")



    sp2 = plt.subplot(212, sharex=sp1)
    line1, = plt.plot(dt[0:ntime],hy_cond[0:ntime,0]*(true_flux['hy_grad'].values[0:ntime]),'orange',linewidth=0.25)
    for i in range(1,nreaz):
        plt.plot(dt[0:ntime],hy_cond[0:ntime,i]*(true_flux['hy_grad'].values[0:ntime]),'orange',linewidth=0.25)
    line2, = plt.plot(dt[0:ntime],np.mean(hy_cond[0:ntime,:],axis=1)*(true_flux['hy_grad'].values[0:ntime]),'k-',linewidth=2)       
    line3, = plt.plot(dt[0:ntime],-true_flux['true_flux'].values[0:ntime]*24,'b',linewidth=2)
    
    plt.xlabel("Date")
    plt.ylabel("VHEF (m/d)")
    plt.xlim([date(2017,4,1),date(2017,5,1)])
    plt.ylim([-10,10])
    plt.legend((line1,line2,line3),('ES-MDA (realization)','ES-MDA (mean)','True flux'),frameon=False)
    
    ax = plt.gca()
    ax.xaxis.set_major_locator(days)
    ax.xaxis.set_major_formatter(daysFmt)
    ax.xaxis.set_minor_locator(mdates.DayLocator())
#    ax.set_xlabel("Date")
#    ax.set_ylabel("VHEF (m/d)")

    
    plt.savefig('../figure/flux.png')
    plt.show()
    
def PlotHydraulicConductivity(path_to_permfile):
#    perm = np.loadtxt('../results/perm.txt',dtype=float)
    perm = np.loadtxt(path_to_permfile,dtype=float)
    hy_cond = perm*1e6*9.8*3600*24
    true_flux = pd.read_csv('../dainput/true_flux.csv',sep=',',header=0)
    true_flux['time'] = [datetime.strptime(x,"%m/%d/%y %H:%M") for x in true_flux['time']]
    dt = true_flux['time']
    dt = [mdates.date2num(x) for x in dt]
    ntime = perm.shape[0]
    nreaz = perm.shape[1] 
    
    plt.figure(num=1,dpi=300)    
    line1, = plt.plot(dt[0:ntime],hy_cond[0:ntime,0],'orange',linewidth=0.25)


    
    for i in range(1,nreaz):       
        plt.plot(dt[0:ntime],hy_cond[0:ntime,i],'orange',linewidth=0.25)
    line2, = plt.plot(dt[0:ntime],np.mean(hy_cond[0:ntime,:],axis=1),'k-',linewidth=2)

    plt.xlabel("Date")
    plt.ylabel("Hydraulic conductivity (m/d)")
    plt.xlim([date(2017,4,1),date(2017,5,1)])
    plt.ylim([0,60])    
#    plt.setp(sp1.get_xticklabels(),visible=False)
    lg = plt.legend((line1,line2),('ES-MDA (realization)','ES-MDA (mean)'),frameon=False)
#    lg.get_frame().set_edgecolor('k')

    ax = plt.gca()
    daysFmt = mdates.DateFormatter("%m/%d/%y")
    days = mdates.DayLocator(interval=5)
    ax.xaxis.set_major_locator(days)
    ax.xaxis.set_major_formatter(daysFmt)
    ax.xaxis.set_minor_locator(mdates.DayLocator())   
    plt.savefig('../figure/HydraulicConductivity.png') 
    plt.show()

    
def PlotVerticalFlux(path_to_perm,path_to_true_flux):    
    perm = np.loadtxt(path_to_perm,dtype=float)
    hy_cond = perm*1e6*9.8*3600*24
    true_flux = pd.read_csv('../dainput/true_flux.csv',sep=',',header=0)
    true_flux['time'] = [datetime.strptime(x,"%m/%d/%y %H:%M") for x in true_flux['time']]
    dt = true_flux['time']
    dt = [mdates.date2num(x) for x in dt]
    ntime = perm.shape[0]
    nreaz = perm.shape[1] 
    
    plt.figure(num=2,dpi=300)
    line1, = plt.plot(dt[0:ntime],hy_cond[0:ntime,0]*(true_flux['hy_grad'].values[0:ntime]),'orange',linewidth=0.25)
    for i in range(1,nreaz):
        plt.plot(dt[0:ntime],hy_cond[0:ntime,i]*(true_flux['hy_grad'].values[0:ntime]),'orange',linewidth=0.25)
    line2, = plt.plot(dt[0:ntime],np.mean(hy_cond[0:ntime,:],axis=1)*(true_flux['hy_grad'].values[0:ntime]),'k-',linewidth=2)       
    line3, = plt.plot(dt[0:ntime],-true_flux['true_flux'].values[0:ntime]*24,'b',linewidth=2)
    
    plt.xlabel("Date")
    plt.ylabel("VHEF (m/d)")
    plt.xlim([date(2017,4,1),date(2017,5,1)])
    plt.ylim([-10,10])
    plt.legend((line1,line2,line3),('ES-MDA (realization)','ES-MDA (mean)','True flux'),frameon=False)
    
    ax = plt.gca()
    daysFmt = mdates.DateFormatter("%m/%d/%y")
    days = mdates.DayLocator(interval=5)    
    ax.xaxis.set_major_locator(days)
    ax.xaxis.set_major_formatter(daysFmt)
    ax.xaxis.set_minor_locator(mdates.DayLocator())
#    ax.set_xlabel("Date")
#    ax.set_ylabel("VHEF (m/d)")

    
    plt.savefig('../figure/VerticalFlux.png')
    plt.show()    
    
def PlotVerticalFluxPrior(path_to_priorKfile,path_to_true_flux):
    init_hy_cond = np.loadtxt(path_to_priorKfile,dtype=float)
    true_flux = pd.read_csv('../dainput/true_flux.csv',sep=',',header=0)
    true_flux['time'] = [datetime.strptime(x,"%m/%d/%y %H:%M") for x in true_flux['time']]
    dt = true_flux['time']
    dt = [mdates.date2num(x) for x in dt]
    nreaz = init_hy_cond.shape[0] 
    
    plt.figure(num=3,dpi=300)
#    print("dt is :{} \n".format(len(dt)))
#    print("hy_grad is : {} \n".format(len(true_flux['hy_grad'].values)))
#    print("init_hy_cond is:{} \n".format(len(init_hy_cond)))
    line1, = plt.plot(dt,init_hy_cond[0,0]*(true_flux['hy_grad'].values),'orange',linewidth=0.25)
    for i in range(1,nreaz-1):
        plt.plot(dt,init_hy_cond[0,i]*(true_flux['hy_grad'].values),'orange',linewidth=0.25)
    line2, = plt.plot(dt,np.mean(init_hy_cond[:],axis=1)*(true_flux['hy_grad'].values),'k-',linewidth=2)       
    line3, = plt.plot(dt,-true_flux['true_flux'].values*24,'b',linewidth=2)
    
    plt.xlabel("Date")
    plt.ylabel("VHEF (m/d)")
    plt.xlim([date(2017,4,1),date(2017,5,1)])
    plt.ylim([-10,10])
    plt.legend((line1,line2,line3),('ES-MDA (realization)','ES-MDA (mean)','True flux'),frameon=False)

    
    ax = plt.gca()
    daysFmt = mdates.DateFormatter("%m/%d/%y")
    days = mdates.DayLocator(interval=5)    
    ax.xaxis.set_major_locator(days)
    ax.xaxis.set_major_formatter(daysFmt)
    ax.xaxis.set_minor_locator(mdates.DayLocator())
#    ax.set_xlabel("Date")
#    ax.set_ylabel("VHEF (m/d)")

    
    plt.savefig('../figure/VerticalFluxPrior.png')
    plt.show()
    