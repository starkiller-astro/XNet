#Read Tracer Fate Function
import numpy as np
from datetime import datetime
import numpy.matlib
#from joblib import Parallel, delayed
#import multiprocessing
from fate_parfor import fate_parfor
from math import pi

def read_tracer_fate(fate_fname_base, plist, time_interp):
    
    num_fields = 15
    
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
            
            fate_fname = fate_fname_base + str(p_id)
            form_array = np.array(['%f'])
            fate_fform = np.matlib.repmat(form_array,1,num_fields)
        
            fate_id = open(fate_fname,'r')
            
            
            if fate_id >= 0:
                file_read = np.loadtxt(fate_id)
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
                
                fate_id.close()
        
        
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
        
        (radius,theta,v_rad,v_theta,temp,density,ye,enpy,pe_int,pe_bind,press,lapse,dpe_nuc,dpe_neut) = fate_parfor(plist,fate_fname_base,time_interp)
        
        toc = datetime.today()
        print toc-tic, 'sec Elapsed'
        
    else:
        p_id = plist[0]
        fate_fname = fate_fname_base + str(p_id)
        fate_id = open(fate_fname,'r')
        
        if fate_id >= 0:
             fate_data = np.loadtxt(fate_id)
             time_from_fate = fate_data[:,0]
             (_, ia) = np.unique(time_from_fate, return_index=True)
             fate_data  = fate_data[ia,:]
             
             time_ = fate_data[:,0]
             radius    = fate_data[:,1]
             theta     = fate_data[:,2]
             v_rad     = fate_data[:,3]
             v_theta   = fate_data[:,4]
             temp      = fate_data[:,5]
             density   = fate_data[:,6]
             ye        = fate_data[:,7]
             enpy      = fate_data[:,8]
             pe_int    = fate_data[:,9]
             pe_bind   = fate_data[:,10]
             press     = fate_data[:,11]
             lapse     = fate_data[:,12]
             dpe_nuc   = fate_data[:,13]
             dpe_neut  = fate_data[:,14]

                 
             fate_id.close()
    
    eps = np.finfo(float).eps
    c = 2.99792458E+10
    
    for k in range(0,len(plist)):
        for j in range(0,len(time_)):
            
            if radius[j,k] < 0.0:
                radius[j,k] = eps
            if theta[j,k] > pi:
                theta[j,k] = pi
            if theta[j,k] < 0.0:
                theta[j,k] = 0.0
            if temp[j,k] < 0.0:
                temp[j,k] = eps
            if v_rad[j,k] >= c:
                v_rad[j,k] = eps
            if v_theta[j,k] >= c:
                v_theta[j,k] = eps
            if density[j,k] < 0.0:
                density[j,k] = eps
            if press[j,k] < 0.0:
                press[j,k] = eps
            if ye[j,k] < 0.0:
                ye[j,k] = eps
            if ye[j,k] > 1:
                ye[j,k] = 1.0

    
    return (time_,radius,theta,v_rad,v_theta,temp,density,ye,enpy,pe_int,pe_bind,press,lapse,dpe_nuc,dpe_neut)