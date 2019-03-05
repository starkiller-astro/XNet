import numpy as np
from scipy import interpolate

def fate_parfor(plist,fate_fname_base,time_interp):
    
    radius   = np.zeros((len(time_interp),len(plist)))
    theta    = np.zeros((len(time_interp),len(plist)))
    v_rad    = np.zeros((len(time_interp),len(plist)))
    v_theta  = np.zeros((len(time_interp),len(plist)))
    temp     = np.zeros((len(time_interp),len(plist)))
    density  = np.zeros((len(time_interp),len(plist)))
    ye       = np.zeros((len(time_interp),len(plist)))
    enpy     = np.zeros((len(time_interp),len(plist)))
    pe_int   = np.zeros((len(time_interp),len(plist)))
    pe_bind  = np.zeros((len(time_interp),len(plist)))
    press    = np.zeros((len(time_interp),len(plist)))
    lapse    = np.zeros((len(time_interp),len(plist)))
    dpe_nuc  = np.zeros((len(time_interp),len(plist)))
    dpe_neut = np.zeros((len(time_interp),len(plist)))
    
    for _id in range(0,len(plist)):
        p_id = plist[_id]
        time_range = sorted(time_interp)
        num_fields = 15
        
        fate_fname = fate_fname_base + str(p_id)
        print fate_fname
        
        fate_id = open(fate_fname,'r')
        
        if fate_id >= 0:
             f_data = np.loadtxt(fate_id)
             time_from_fate = f_data[:,0]
             (_, ia) = np.unique(time_from_fate, return_index=True)
             f_data = [f_data[i,:] for i in ia]
             fate_data = np.array(f_data)
             
             if (len(time_from_fate)) > 1: 
                 #p_data = np.zeros((len(time_range),num_fields))
                 
                 interp_func = interpolate.interp1d(fate_data[:,0],fate_data,kind='linear',axis=0,fill_value='extrapolate')
                 p_data = interp_func(time_range)

             else:
                 p_data = np.tile(fate_data[0,:], ((len(time_range)),1))   
                 
                
             radius[:,_id] = p_data[:,1]
             theta[:,_id] = p_data[:,2]
             v_rad[:,_id] = p_data[:,3]
             v_theta[:,_id] = p_data[:,4]
             temp[:,_id] = p_data[:,5]
             density[:,_id] = p_data[:,6]
             ye[:,_id] = p_data[:,7]
             enpy[:,_id] = p_data[:,8]
             pe_int[:,_id] = p_data[:,9]
             pe_bind[:,_id] = p_data[:,10]
             press[:,_id] = p_data[:,11]
             lapse[:,_id] = p_data[:,12]
             dpe_nuc[:,_id] = p_data[:,13]
             dpe_neut[:,_id] = p_data[:,14]
                 
             fate_id.close()
        
        p_data = []
        fate_data = []
        
    return (radius,theta,v_rad,v_theta,temp,density,ye,enpy,pe_int,pe_bind,press,lapse,dpe_nuc,dpe_neut)             