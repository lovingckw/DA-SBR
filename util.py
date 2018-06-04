def MakePflotranInput(pflotranin,simu_time,itime,ncollect,collect_times,perm_update_index,nreaz):
    #fpflotran = open("1dthermal.in", 'r')
    #pflotranin = fpflotran.readlines()
    # creat a checkpoint at the end of simulation
    lindex = [i for i, s in enumerate(pflotranin) if "CHECKPOINT" in s][0]
    pflotranin[lindex+1] = "    TIMES sec "+np.array_str(simu_time[itime,1])+"\n"
    
    # determin whether to use restart
    lindex = [i for i, s in enumerate(pflotranin) if "pflotran-restart.chk" in s][0]
    if simu_time[itime,0] == 0:
        pflotranin[lindex] = "# RESTART pflotran-restart.chk"+"\n"
    else:
        pflotranin[lindex] = " RESTART 1dthermal-"+np.array_str(simu_time[itime,0])+"sec.chk"+"\n"
    
    # write the end time of simulation
    lindex = [i for i, s in enumerate(pflotranin) if "FINAL_TIME" in s][0]
    pflotranin[lindex] = "      FINAL_TIME "+ np.array_str(simu_time[itime,1])+" sec"+"\n"
    
    # generate observation 
    snapshot_lindex =  [i for i, s in enumerate(pflotranin) if "SNAPSHOT_FILE" in s][0]-1
    obs_card = "    TIMES sec   "
    npl = 3 
    nline = np.floor_divide(ncollect,npl)
    if nline == 0 or nline == 1:
        obs_card = obs_card+" "+np.array_str(collect_times[0])
    else:
        for iline in range(1,nline):
            obs_card = obs_card+np.array_str(collect_times[((iline-1)*npl+1):(iline*npl)])+"\\"
        obs_card = obs_card+" "+np.array_str(collect_times[(iline*npl+1):ncollect])+"\n"
    
    # add observation time
    pflotranin = pflotranin[0:snapshot_lindex+1]+[obs_card]+pflotranin[snapshot_lindex+2:]     
    
    #prepare strata_card
    strata_card_lindex = [i for i, s in enumerate(pflotranin) if "STRATA" in s][0]-1
    strata_card = ["STRATA\n",
                   "REGION\n",
                   "    MATERIAL Perm 0\n",
                   "    START_TIME 0 sec\n",
                   "    FINAL_TIME 0 sec\n",
                   "END\n",
                   "\n"]
    strata_card_length = len(strata_card)
    strata_card_new = list()
    
    #generate new strata card
    for jtime in perm_update_index:
        strata_card[2] = "  MATERIAL Perm "+str(jtime)+"\n"
        if jtime == perm_update_index[0]:
            strata_card[3] = "  START_TIME "+"0 sec"+"\n"
        else:
            strata_card[3] = "  START_TIME "+str(da_interval*jtime)+" sec"+"\n"
        strata_card[4] = "  FINAL_TIME "+str(da_interval*(jtime+1))+" sec"+"\n"
        strata_card_copy = strata_card[:]
        strata_card_new = strata_card_new+strata_card_copy
        
    pflotranin = pflotranin[0:strata_card_lindex]+strata_card_new+pflotranin[strata_card_lindex+strata_card_length:]
    
    # prepare perm_card
    perm_card_length = 14
    perm_card_lindex = [i for i, s in enumerate(pflotranin) if "MATERIAL_PROPERTY Alluvium" in s][0]
    perm_card = pflotranin[perm_card_lindex:perm_card_lindex+perm_card_length]
    perm_value_cardindex = [i for i, s in enumerate(perm_card) if "PERM_ISO" in s][0]
    perm_card_new = list()
    
    for ireaz in range(nreaz):
        perm_card_new.clear()
        for jtime in perm_update_index:
            perm_card[0] = "MATERIAL_PROPERTY Perm"+str(jtime)+"\n"
            perm_card[1] = "  ID "+str(jtime)+"\n"
            perm_card[perm_value_cardindex] = "  PERM_ISO "+str(perm[jtime,ireaz])+"\n"
            perm_card_copy = perm_card[:] 
            perm_card_new = perm_card_new+perm_card_copy
            
        temp_pflotranin = pflotranin[0:perm_card_lindex]+perm_card_new+pflotranin[perm_card_lindex+perm_card_length:]
        new_pflotranin = open("./pflotran/%d/1dthermal.in" % ireaz,'w')
        new_pflotranin.writelines(temp_pflotranin)  
        new_pflotranin.close()
    return

def GenerateSimuEnsemble(nobs,obs_coord,z,nreaz,collect_times):
    obs_cell = np.zeros(nobs)
    for i in range(nobs):
        obs_cell[i] = np.argmin(np.absolute(z-obs_coord[i]))
    obs_cell = obs_cell.astype(int)
    
    simu_ensemble = np.zeros((nreaz,nobs*len(collect_times)))
    for ireaz in range(nreaz):
        obs_temp = np.zeros(nobs*len(collect_times))
        jj = 0
        for collect_itime in collect_times:
            h5_file = h5py.File("./pflotran/%d/1dthermal.h5" % ireaz,'r')
            key = "Time:"+str(" %12.5E" % collect_itime)+" s/Temperature [C]"
            obs_temp[jj*nobs:(jj+1)*nobs] = h5_file[key][0][0][obs_cell]
            jj = jj+1
        simu_ensemble[ireaz,:] = obs_temp
        h5_file.close()
    
    return simu_ensemble
            
        