#Read Tracer Temp Function
import numpy as np
from datetime import datetime
import numpy.matlib
#from joblib import Parallel, delayed
#import multiprocessing
from temp_parfor import *

def read_tracer_temp(temp_fname_base, plist, time_interp):
    
    num_fields = 12
    
    #Maybe rephrase this if-statement?
    
    if time_interp.size != 0:
        interp_flag=True
    else:
        interp_flag=False
    
    if interp_flag:

        time_start = time_interp[0]
        time_end = [0,0]
    
    # Determine time range of data files and adjust time range to fit data if necessary
    
        for _id in range(0,len(plist)):
        
            p_id = plist[_id]
            
            time_min = time_start
            time_max = time_interp[len(time_interp)-1]
            
            temp_fname = temp_fname_base + str(p_id)
            
            #form_array = np.array(['%f'])
            #temp_fform = np.matlib.repmat(form_array,1,num_fields)
        
            temp_id = open(temp_fname,'r')
            
            
            if temp_id >= 0:
                file_read = np.loadtxt(temp_id)
                ###file_read = cell2mat(file_read) -> this is already an ndarray, not a cell array
                time_init = file_read[0,0]
                
                '''
                offset_val = temp_id.tell()
                offset = -offset_val - 1
                temp_id.seek(offset,2)
                file_read = np.loadtxt(temp_id)
                '''
                reverse_file_read = file_read[::-1]
                time_final = reverse_file_read[0,0]
                
                ii = 0
                jj = (len(time_interp)-1)
                while time_min < time_init and ii < jj:
                    ii = ii + 1
                    time_min = time_interp[ii]
            
                while time_max > time_final and jj > ii:
                    jj = jj - 1
                    time_max = time_interp[jj]
            
                if time_min > time_start:
                    time_start = time_min
                    
                if time_final > time_end[0]:
                    time_end = [time_max, time_final]
                
                temp_id.close()
        
        
       # (_, ii) = np.unique(time_interp, return_index=True)
        ii = np.where(time_interp==time_start)
        ii=int(ii[0])
        jj = np.where(time_interp==time_end[0])
        jj=int(jj[0])
        
        #time_interp = [time_interp(ii:jj-1),time_end];
        #time_interp = [time_interp(ii:jj-1),time_max];
        
        time_interp_slice = time_interp[ii:(len(time_interp))]
        time_interp_slice = time_interp_slice[::-1]
        time_interp = np.append(time_interp_slice,time_init)
        time_interp = time_interp[::-1]
        time_interp = np.unique( time_interp )   
            
        time_     = np.transpose(time_interp)
        
        tic = datetime.today()
        #all code following this tic goes here
        
        (temp,density,ye,flxtot,nu_temp) = temp_parfor(plist,temp_fname_base,time_interp)
        
        toc = datetime.today()
        print toc-tic, 'sec Elapsed'
        
    else:
        p_id = plist[0]
        temp_fname = temp_fname_base + str(p_id)
        temp_id = open(temp_fname,'r')
        
        if temp_id >= 0:
             temp_data = np.loadtxt(temp_id)
             time_from_temp = temp_data[:,0]
             (_, ia) = np.unique(time_from_temp, return_index=True)
             temp_data  = temp_data[ia,:]
             
             time_ = temp_data[:,0]
             temp = temp_data[:,1]
             density = temp_data[:,2]
             ye = temp_data[:,3]
             
             flxtot   = np.zeros((len(time_),4))
             nu_temp  = np.zeros((len(time_),4))
             
             for i in range(0,4):
                 flxtot[:,i] = temp_data[:,4+i]
                 nu_temp[:,i] = temp_data[:,8+i]
                 
             temp_id.close()
    
    eps = np.finfo(float).eps
    
    for k in range(0,len(plist)):
        for j in range(0,len(time_)):
            
            if temp[j,k] < 0.0:
                temp[j,k] = eps
            if density[j,k] < 0.0:
                density[j,k] = eps
            if ye[j,k] < 0.0:
                ye[j,k] = eps
            if ye[j,k] > 1:
                ye[j,k] = 1.0
            for i in range(0,4):
                if nu_temp[j,i,k] < 0.0:
                    nu_temp[j,i,k] = eps
                if np.isnan(nu_temp[j,i,k]) == True:
                    nu_temp[j,i,k] = eps
    
    return (time_, temp, density, ye, flxtot, nu_temp)