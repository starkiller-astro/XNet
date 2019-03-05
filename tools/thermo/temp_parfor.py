import numpy as np
from scipy import interpolate
'''
from get_time_input import get_time_input
time_start = 0.0
time_bounce  = 2.63176547E-01
time_final   = 1.67357
N = 5000
N_bounce = 10

time_interp = get_time_input(time_start, time_final, time_bounce, N, N_bounce)
path_temp = 'C:\\Users\\taylor\\Documents\\MATLAB\\B12\\B12-WH07\\full_profile\\p_row44\\temp_S12-0'
plist = np.arange(1721,1761)
temp_fname_base = path_temp
'''

def temp_parfor(plist,temp_fname_base,time_interp):
    
    temp     = np.zeros((len(time_interp),len(plist)))
    density  = np.zeros((len(time_interp),len(plist)))
    ye       = np.zeros((len(time_interp),len(plist)))
    flxtot   = np.zeros((len(time_interp),4,len(plist)))
    nu_temp  = np.zeros((len(time_interp),4,len(plist)))
    
    for _id in range(0,len(plist)):
        p_id = plist[_id]
        time_range = sorted(time_interp)
        num_fields = 12
        
        temp_fname = temp_fname_base + str(p_id)
        print temp_fname
        
        temp_id = open(temp_fname,'r')
        
        if temp_id >= 0:
             t_data = np.loadtxt(temp_id)
             time_from_temp = t_data[:,0]
             (_, ia) = np.unique(time_from_temp, return_index=True)
             t_data = [t_data[i,:] for i in ia]
             temp_data = np.array(t_data)
             
             #temp_data = []
             #for i in ia:
                #temp_data  = np.append(temp_data,t_data[i,:])
            ###Change to make temp_data append/replace values.
             
             if (len(time_range)) > 1: 
                 #p_data = np.zeros((len(time_range),num_fields))
                 
                 interp_func = interpolate.interp1d(temp_data[:,0],temp_data,kind='linear',axis=0,fill_value='extrapolate')
                 p_data = interp_func(time_range)

             else:
                 p_data = np.tile(temp_data[0,:], ((len(time_range)),1))   
                
             temp[:,_id] = p_data[:,1]
             density[:,_id] = p_data[:,2]
             ye[:,_id] = p_data[:,3]
             for i in range(0,4):
                 flxtot[:,i,_id] = p_data[:,4+i]
                 nu_temp[:,i,_id] = p_data[:,8+i]
                 
             temp_id.close()
        
        p_data = []
        temp_data = []
        
    return (temp,density,ye,flxtot,nu_temp)             