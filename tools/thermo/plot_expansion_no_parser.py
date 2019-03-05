import os,sys,itertools
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as plticker
import matplotlib.gridspec as gridspec
import matplotlib.axes as plax
from datetime import datetime
from scipy import interpolate
from scipy import signal
import scipy.io as scio
import numpy.matlib
from math import floor
import warnings

with warnings.catch_warnings():
    warnings.filterwarnings("ignore",category=DeprecationWarning)
    import md5, sha

from get_time_input import get_time_input
from read_tracer_temp import *
from read_tracer_fate import *

'''Import all data into variables via parser in MATLAB: Here, we initialize and use functions'''

model_mass = '12'
model_name = 'B'+str(model_mass)+'-WH07'

time_start = 0.0
time_bounce  = 2.63176547E-01
time_final   = 1.67357
N = 5000
N_bounce = 10

time_interp = get_time_input(time_start, time_final, time_bounce, N, N_bounce)

path_temp = 'C:\\Users\\taylor\\Documents\\MATLAB\\B12\\B12-WH07\\full_profile\\p_row44\\temp_S12-0'
path_fate = 'C:\\Users\\taylor\\Documents\\MATLAB\\B12\\B12-WH07\\full_profile\\p_row44\\fate_S12-0'

plist = np.arange(1721,1761)
temp_fname_base = path_temp
fate_fname_base = path_fate

(time_,radius,theta,v_rad,v_theta,temp,density,ye,enpy,pe_int,pe_bind,press,lapse,dpe_nuc,dpe_neut) = read_tracer_fate(fate_fname_base, plist, time_interp)
(time_, temp, density, ye, flxtot, nu_temp) = read_tracer_temp(temp_fname_base, plist, time_interp)

 
time_extrap_min = 0.0025
time_extrap_max = 0.15
adiabatic_tol = 0.05
temp_nse = 8.0
temp_final = 0.5
sgwin_tau = 0.025
sgorder_tau = 6
sgwin_adi = 0.05
sgorder_adi = 2
tau_rho_min = 0.01
tau_rho_max = 1.0
min_jumpdiff = 5
change_min = 0.001
t_extend = 0.0
temp_min = 0.02
output_tdel = 0.0
nu_time_stop = 10.0

#last time temp leaves nse temp
t_start_array = []

#For each particle, we find indices where its temp exceeds the NSE temp, and we
#take the last index and get the final time when it was above that temperature.

for k in range(0,len(plist)):
        temp_per_p = temp[:,k]
        index = np.where(temp_per_p>=temp_nse)[0]
        if type(index) == int:
            if index or index==0:
                i = int(index[len(index)-1])
                t_start = time_[i]
                t_start_array = np.append(t_start_array,t_start)      
            else:
                t_start_array = np.append(t_start_array,0)
                      
        elif type(index) == numpy.ndarray:
            if bool(index.any()):
                i = int(index[len(index)-1])
                t_start = time_[i]
                t_start_array = np.append(t_start_array,t_start)       
            else:
                t_start_array = np.append(t_start_array,0) 


#times from which to start extrapolation if extrap_flag==true:
t_stop_array = [time_[len(time_)-1]] #default

#Initialize more variables:

print_flag = False
write_flag = True
extrap_flag = True
plot_flag = True

prc_min = 25
prc_max = 75

time_extrap_max0 = time_extrap_max;

#Make sure t_stop_array is the final time value.

if t_stop_array[0] == 0.0:
    t_stop_array[0] = time_[len(time_)]

if time_final > 0.0:
    #t_stop_array[t_stop_array > time_final] = time_final
    for x in range(0,len(t_stop_array)):
        if t_stop_array[x] > time_final:
            t_stop_array[x] = time_final


plot_folder = '.\\'  #unnecessary in OG code

if print_flag:
    plot_flag = True
    if not os.path.isdir(plot_folder):
        os.mkdir(plot_folder,0777)

if write_flag:
    profile_folder = os.getcwd()
    if not os.path.isdir(profile_folder):
        os.mkdir(profile_folder)
    th_format = [np.tile('{:15.7e}',(1, 4)), np.tile('{:12.3e}', (1, 8)), '\n' ]
        
    th_fname_base = profile_folder+'\\th_profile\\th-' 

adiabatic_string = "%s%g%s" % (r'${\Delta}(T{\rho}^{-1/3}) <$ ',adiabatic_tol,r'$\times(T{\rho}^{-1/3})_{\mathrm{f}}$' )


cm = plt.get_cmap('gist_rainbow')

#fig = plt.figure()
#ax = fig.add_subplot(111)
#ax.set_color_cycle([cm(1.*i/NUM_COLORS) for i in range(NUM_COLORS)])
#for i in range(NUM_COLORS):
    #ax.plot(np.arange(10)*(i+1))

color_array = []
NUM_COLORS = len(t_stop_array)

#Assign colors to an array for plotting later: 1 color per time value to start 
#extrapolation from.

if len(t_stop_array) > 6:   
    for k in range(NUM_COLORS):
        color_array = np.append(color_array, cm(1.*k/NUM_COLORS))
    color_array = color_array.reshape(NUM_COLORS,4)        
elif len(t_stop_array) < 1:    
    color_array=[]
elif len(t_stop_array)==1:   
    color_array=[('0.','0.','1.','1.')]   
elif len(t_stop_array)==2:
    color_array=[[('0.','0.','1.','1.')], [('1.','0.','0.','1.')]]   
elif len(t_stop_array)==3:    
    color_array=[[('0.','0.','1.','1.')], [('1.','0.','0.','1.')], [('0.','0.5','0.','1.')]]   
elif len(t_stop_array)==4:
    color_array=[[('0.','0.','1.','1.')], [('1.','0.','0.','1.')], [('0.','0.','0.5','1.')], [('0.','0.75','0.75','1.')]]
elif len(t_stop_array)==5:
    color_array=[[('0.','0.','1.','1.')], [('1.','0.','0.','1.')], [('0.','0.','0.5','1.')], [('0.','0.75','0.75','1.')], [('0.75','0.','0.75','1.')]]
elif len(t_stop_array)==6:
    color_array=[[('0.','0.','1.','1.')], [('1.','0.','0.','1.')], [('0.','0.','0.5','1.')], [('0.','0.75','0.75','1.')], [('0.75','0.','0.75','1.')], [('0.75','0.75','0.','1.')]]

a = np.arange(0,1)

#Set t_start_array value equal to the initial for every particle if len==1.
if len(t_start_array) == 1 and len(plist) > 1:
    for i in range(1,len(plist)):
        t_start_array = np.concatenate((t_start_array,a))
        t_start_array[i] = t_start_array[0]        
                                                
i_peak = np.ones(plist.shape)


#temp_peak is an array of the maximum temp values per particle
temp_peak = np.zeros((1,len(plist)))
for _id in range(0,len(plist)):
    temp_peak[0,_id] = np.max(temp[:,_id])

for id_ in range (0,len(plist)):
    
    if len(plist) == 1:
        p_id = 0
    else:
        p_id = id_
    
    if t_start_array[p_id] > 0.0 :
        time_minus_t_start = abs( time_ - t_start_array[p_id] )
        i_nse = np.argmin(time_minus_t_start)
        i_peak[p_id] = i_nse
     
    else:
        i_nse = 0
        if temp_peak[0,p_id] >= temp_nse and temp[len(time_)-1,p_id] < temp_nse:
            i_nse = np.where(temp[:,p_id] > temp_nse)[0]
            if i_nse or i_nse==0:
                i_nse = int(i_nse[len(i_nse)-1])
            
        elif temp[len(temp)-1,p_id] >= temp_nse:
            i_nse = len(time_)-1
            
        elif temp_peak[0,p_id] < temp_nse:
            i_nse = 0
        
        if type(i_nse) == int:
            if i_nse or i_nse==0:    
                t_start_array[p_id] = time_[i_nse]
                i_peak[p_id]       = i_nse
        else:
            if i_nse.any():    
                t_start_array[p_id] = time_[i_nse]
                i_peak[p_id]       = i_nse

tau = np.zeros((len(t_stop_array),len(plist)))
th_out = 0

for _id in range(0,len(plist)):
#for _id in range(3,4):
              
    if len(plist) == 1:
        p_id = 0
    else:
        p_id = _id
        print p_id
    
    #Reset extrap flag to initial value
    extrap  = extrap_flag
    
    t_stop_max = np.max(t_stop_array)
    time_minus_t_stop =  abs( time_ - t_stop_max) 
    istop = np.argmin(time_minus_t_stop)

    if temp_final >= temp[istop,p_id]:
        
        temp_over_temp_final = np.where(temp[:,p_id] >= temp_final)[0]
        temp_over_temp_final = temp_over_temp_final[len(temp_over_temp_final)-1]
        
        if (np.all(temp_over_temp_final==0)):
            istop = temp_over_temp_final[0]       

        t_stop   = time_[istop]
        extrap   = False
        
    elif temp[istop,p_id] >= temp_nse:
        t_stop   = time_[istop]
        extrap   = False
        
    else:
        t_stop   = time_[istop]
        
    
    #This is the time which is written to the profile to begin post-processing
    t_start = t_start_array[p_id]
    
    if istop <= i_peak[p_id]:
        peak = 1
    else:
        peak = i_peak[p_id]
    
    #This is the time from which we may extrapolate
    t_peak  = time_[peak]

    if extrap or plot_flag:
        rho = density[:,p_id]
        max_peak_times = max(t_stop - time_extrap_max, t_peak)
        time_difference = abs(time_-max_peak_times)
        iwin_min = np.argmin(time_difference)
        span_tau = floor( sgwin_tau * (istop - iwin_min + 1) / ( time_[istop] - time_[iwin_min] ) )
        span_tau = int(span_tau)
        if (span_tau%2) == 0:
            span_tau = span_tau + 1
           
        span_tau = max( span_tau, sgorder_tau+1 )
        iwin_min = max( iwin_min-span_tau, 1 )
        iwin_max = min( istop+span_tau, len(time_) )
        rho[iwin_min:iwin_max] = signal.savgol_filter(density[iwin_min:iwin_max, p_id], span_tau, sgorder_tau)
        #rho(iwin_min:iwin_max) = smooth( time(iwin_min:iwin_max), density(iwin_min:iwin_max,p_id), span_tau, 'sgolay', sgorder_tau );
        
        rhodot = (np.gradient(rho))/(np.gradient(time_))
        
        tau_rho = (-rho) / rhodot

        adiabatic_raw = 0.34 * (temp[:,p_id]**3) / density[:,p_id]
        adiabatic = adiabatic_raw
        iwin_min = np.argmin(time_difference)
        span_adi = floor( sgwin_adi * (istop - iwin_min + 1) / ( time_[istop] - time_[iwin_min] ) )
        span_adi = int(span_adi)
        if (span_adi%2) == 0:
            span_adi = span_adi + 1

        
        span_adi = max( span_adi, sgorder_adi+1 )
        iwin_min = max( iwin_min-span_adi, 1 )
        iwin_max = min( istop+span_adi, len(time_) )
        adiabatic[iwin_min:iwin_max] = signal.savgol_filter(adiabatic_raw[iwin_min:iwin_max], span_adi, sgorder_adi )

#         rhodot = gradient( density(:,p_id), time );
#         rhodot = [ 0; diff(density(:,p_id)) ./ diff(time) ];
#         rhodot = deriv_3pt( density(:,p_id), time );
#         [~,span_tau] = min( abs( time - t_stop + sgwin_tau ) );
#         span_tau = max( istop - span_tau + 1, sgorder_tau + 1 )
#         tic;
#         rho = smooth( time, density(:,p_id), span_tau, 'sgolay', sgorder_tau );
#         toc;
#         rhodot = gradient( rho, time );
#         rhodot = deriv_5pt( rho, time );
#         
#         tau_rho = -density(:,p_id) ./ rhodot;
#         tau_rho = -rho ./ rhodot;
#         
#         temp_smooth = smooth( time, temp(:,p_id), sgwin_adi, 'sgolay', sgorder_adi );
# 
#         [~,span_adi] = min( abs( time - t_stop + sgwin_adi ) );
#         span_adi = max( istop - span_adi + 1, sgorder_adi + 1 )
#         
#         adiabatic_raw = 0.34 * temp(:,p_id).^3 ./ density(:,p_id);
#         adiabatic_raw = temp_smooth .* rho.^(-1/3);
# 
#         tic;
#         adiabatic   = smooth( time, adiabatic_raw, span_adi, 'sgolay', sgorder_adi );
#        toc;
#         adiabatic   = smooth( adiabatic_raw, 'rlowess', sgwin_adi );
#         adiabatic   = adiabatic_raw;

    
    if plot_flag:
#         tau_rho_avg = smooth( tau_rho, 'rlowess' );
#         tau_rho_avg = smooth( tau_rho, sgwin_tau, 'sgolay', sgorder_tau);
        tau_rho_avg = tau_rho
     
        
        h_ratio = [3,1]
        total_height = h_ratio[0] + h_ratio[1]
        gs = gridspec.GridSpec(total_height, 1)

        fig_h = plt.figure()

        axis_h1 = plt.subplot(gs[:h_ratio[0], :])
        axis_h1 = plt.gca()
        plax.Axes.set_yscale(axis_h1,'linear')
        ylim1 = plax.Axes.set_ylim(axis_h1,0,temp_nse)
        plt.yticks(np.arange(0,9))
        plax.Axes.set_xscale(axis_h1,'linear')
        xlim1 = plax.Axes.set_xlim(axis_h1,t_peak-time_bounce, t_stop-time_bounce)

        axis_h1.grid(True)
        #plt.minorticks_on()
        #plt.tick_params(axis='y', which='minor', direction='out')


        axis_h2 = plt.subplot(gs[h_ratio[0],:],sharex=axis_h1)
        plt.subplots_adjust(wspace=0, hspace=0)
        axis_h2 = plt.gca()
        plax.Axes.set_yscale(axis_h2,'linear') 
        axis_h2 = plt.gca()     
        plax.Axes.set_ylim(axis_h2,-2.9999,2.9999)
        plt.yticks([-2.0,2.0])
        plax.Axes.set_xscale(axis_h2,'linear')
        plax.Axes.set_xlim(axis_h2,t_peak-time_bounce, t_stop-time_bounce)

        axis_h2.grid(True)
        #plt.minorticks_on()
        #plt.tick_params(axis='x', which='minor', direction='out')     
        
        if time_bounce != 0.0:
            axis_h2.set_xlabel('Time after bounce [s]' )
        else:
            axis_h2.set_xlabel('Elapsed time [s]' )
        
        axis_h1.set_ylabel('Temperature [GK]')
        axis_h2.set_ylabel('Expansion Timescale [s]')
        
        #Plot the non-extrapolated temperature profile
        handle1, = axis_h1.plot((time_[peak:istop] - time_bounce),temp[peak:istop,p_id],'k-', linewidth = 1.5, label='Temperature')
        
        ##############
        raw_data_dict = {}
        data_path = 'C:\\Users\\taylor\\Documents\\MATLAB\\B12\\B12-WH07\\raw_data'
        raw_data_dict = scio.loadmat(data_path, appendmat=True, variable_names=('temp','time'))
        temp_from_raw = raw_data_dict['temp']
        time_from_raw = raw_data_dict['time']
        '''
        axis_h1.plot((time_from_raw[peak:istop]-time_bounce),temp_from_raw[peak:istop,p_id],'r--',linewidth=1.5)
        #^We will leave this up for now, but this reads and plots MATLAB data for comparison.
        '''##############

        #Plot the expansion timescale: change back to 'k--'
        handle2, = axis_h2.plot((time_[(peak+5):istop]-time_bounce),tau_rho_avg[(peak+5):istop],'k--',label=r'$\tau_{\mathrm{exp}}$')
        
        '''#########
        #This reads tau_rho_avg per particle from the MATLAB script and plots for comparison.
        
        tau_rho_avg_dict = {}
        data_path = 'C:\\Users\\taylor\\Documents\\MATLAB\\B12\\B12-WH07\\rho_values'+str(p_id+1)
        tau_rho_avg_dict = scio.loadmat(data_path, appendmat=True, variable_names=('tau_rho_avg'))
        tau_rho_avg_tracer1 = tau_rho_avg_dict['tau_rho_avg']
        
        axis_h2.plot((time_from_raw[(peak+5):istop]-time_bounce),tau_rho_avg_tracer1[(peak+5):istop],'r--')
        '''###########

        axis_h3 = axis_h1.twinx()

        ylim3=plax.Axes.set_ylim(axis_h3,min(adiabatic[peak:istop]),max(adiabatic[peak:istop]))

        #Plot the isentropic curve. Change back to 'k-.' when no longer comparing!
        handle3, = axis_h3.plot(time_[peak:istop]-time_bounce, adiabatic[peak:istop],'k-.',linewidth=0.5, label=r'$0.34 T_{9}^{3} / \rho_{5}$')
        
        num_ticks = len(axis_h1.get_yticks())
        ytick = np.linspace( ylim3[0], ylim3[1], num_ticks - 2)
        plax.Axes.set_yticks(axis_h3,ytick)
        plax.Axes.set_ylim(axis_h3, ylim3[0] - (ytick[1]-ytick[0]), ylim3[1] + (ytick[len(ytick)-1]-ytick[len(ytick)-2]) )
        axis_h3.set_ylabel(r'$0.34 T_{9}^{3} / \rho_{5}$')
        
        #plax.Axes.set_yticks(axis_h3,np.linspace(axis_h3.get_yticks()[0],axis_h3.get_yticks()[-1],len(axis_h1.get_yticks())))
        
        '''###########
        #This reads adiabatic data from MATLAB script per particle and plots for comparison.
        
        adiabatic_dict = {}
        data_path = 'C:\\Users\\taylor\\Documents\\MATLAB\\B12\\B12-WH07\\adiabatic_values'+str(p_id+1)
        adiabatic_dict = scio.loadmat(data_path, appendmat=True, variable_names=('adiabatic'))
        adiabatic_tracer1 = adiabatic_dict['adiabatic']
        
        axis_h3.plot(time_from_raw[peak:istop]-time_bounce, adiabatic_tracer1[peak:istop],'r-.',linewidth=0.5)
        '''#############
                        
                
        fig_h2 = plt.figure()
        axis_h2_1 = plt.subplot(111)
        axis_h2_1 = plt.gca()
        plax.Axes.set_yscale(axis_h2_1,'log')
        #ylim1 = plax.Axes.set_ylim(axis_h2_1,0,temp_nse)
        #plt.yticks(np.arange(0,9))
        plax.Axes.set_xscale(axis_h2_1,'linear') ############
        #xlim1 = plax.Axes.set_xlim(axis_h1,t_peak-time_bounce, t_stop-time_bounce)
         
        handle_flx1,= axis_h2_1.plot((time_[peak:istop] - time_bounce),flxtot[peak:istop,0,p_id],'k-', linewidth = 1.5, label='Electron Neutrino Flux')   
        handle_flx2,= axis_h2_1.plot((time_[peak:istop] - time_bounce),flxtot[peak:istop,1,p_id],'r-', linewidth = 1.5, label='Electron Anti-Neutrino Flux')        
        #handle_nutemp1,= axis_h2_1.plot((time_[peak:istop] - time_bounce),nu_temp[peak:istop,0,p_id],'b-', linewidth = 1.5, label='Column 1 Neutrino Temp')        
        #handle_nutemp2,= axis_h2_1.plot((time_[peak:istop] - time_bounce),nu_temp[peak:istop,1,p_id],'g-', linewidth = 1.5, label='Column 2 Neutrino Temp')        
        axis_h2_1.set_xlabel('Time (s)')
        axis_h2_1.set_ylabel('Neutrino Flux (neutrinos/s)')
        
        
        '''        
        ##############
        neutrino_dict = {}
        data_path = 'C:\\Users\\taylor\\Documents\\MATLAB\\B12\\B12-WH07\\nu_values'
        neutrino_dict = scio.loadmat(data_path, appendmat=True, variable_names=('flxtot','nu_temp'))
        flxtot_MAT = neutrino_dict['flxtot']
        nu_temp_MAT = neutrino_dict['nu_temp']
        
        axis_h2_1.plot((time_from_raw[peak:istop] - time_bounce),flxtot_MAT[peak:istop,0,p_id],'c--', linewidth = 1.5, label='Column 1 Neutrino Flux')   
        axis_h2_1.plot((time_from_raw[peak:istop] - time_bounce),flxtot_MAT[peak:istop,1,p_id],'y--', linewidth = 1.5, label='Column 2 Neutrino Flux')        
        axis_h2_1.plot((time_from_raw[peak:istop] - time_bounce),nu_temp_MAT[peak:istop,0,p_id],'k--', linewidth = 1.5, label='Column 1 Neutrino Temp')        
        axis_h2_1.plot((time_from_raw[peak:istop] - time_bounce),nu_temp_MAT[peak:istop,1,p_id],'m--', linewidth = 1.5, label='Column 2 Neutrino Temp')
        ##############  
        '''      
        
    #The for loop lets us plot a line of a different color for each possible time to extrapolate from.
    for n in range(0,len(t_stop_array)):
        istop = np.argmin( abs( time_ - t_stop_array[n] ) )
#       istop = find( time > t_stop_array(n), 1, 'first' )
        extrap = extrap_flag
        extrap_window = []
               
        if temp_final >= temp[istop,p_id]:
            find_val = np.where(temp[:,p_id] > temp_final)[0]
            if find_val.any():
                istop = find_val[len(find_val)-1]
            
            tau[n,p_id] = 0.0
            t_stop   = time_[istop]
            t_exp    = t_stop
            t_extrap = t_stop
            extrap   = False
            
        elif temp_nse <= temp[istop,p_id]:
            tau[n,p_id] = 0.0
            t_stop   = time_[istop]
            t_exp    = t_stop
            t_extrap = t_stop
            extrap   = False
        else:
            t_stop   = time_[istop]
            
        
        if extrap:
## Determine the time window to be used to estimate the expansion timescale

            # Set the maximum size of the window
            time_extrap_max = time_extrap_max0
            if t_stop <= (t_peak + time_extrap_max):
                time_extrap_max = t_stop - t_peak
                warnings.warn('Time range since last peak is less than time_extrap_max for t_stop='+str(t_stop))
            elif t_stop <= t_peak + time_extrap_min:
                time_extrap_max = 0.0
                warnings.warn('Time range since last peak is less than time_extrap_min for t_stop='+str(t_stop))
            

            # Set the time index for the end of the time window
            iwin_stop = istop
            iwin_stop = int(iwin_stop)
            
            # Determine the time index for the maximum window size
            iwin_max = np.argmin( abs( time_[peak:istop] - t_stop + time_extrap_max ) )
            iwin_max = iwin_max + peak
            iwin_max = int(iwin_max)

            # Determine the time index for the minimum window size
            iwin_min = np.argmin( abs( time_[iwin_max:istop] - t_stop + time_extrap_min ) )
            iwin_min = iwin_min + iwin_max
            iwin_min = int(iwin_min)

            # Determine how far back from the end of the simulation the constant isentropic/adibatic criteria holds
            adiabatic_istop = adiabatic[istop]
            
            adi_value = abs((adiabatic[iwin_max:istop]-adiabatic_istop)/(adiabatic_istop))

            iwin_adi = np.where(adi_value >= adiabatic_tol)[0]

            if type(iwin_adi) == numpy.ndarray:
                if iwin_adi.any():
                    iwin_adi = iwin_adi[len(iwin_adi)-1]
                    
            iwin_adi = iwin_adi + iwin_max
            
            # Set the beginning of the window from isentropic/adiabatic criteria if between iwin_max and iwin_min
            if not iwin_adi:
                iwin_start = iwin_max
            elif iwin_adi <= iwin_min:
                iwin_start = iwin_adi
            else:
                iwin_start = istop
                warnings.warn('Isentropic conditions not met')
            
            iwin_start = int(iwin_start)
            
            # First consider time ranges where constant isentropic/adiabatic criteria is true
            adi_window = np.arange(iwin_start,iwin_stop+1)

            # Restrict the window to time ranges where expansion timescale is within limits
            tau_window = []
            for i in range(0,len(adi_window)):
                if  tau_rho[adi_window[i]] > tau_rho_min and tau_rho[adi_window[i]] < tau_rho_max:
                    tau_window = np.append(tau_window, adi_window[i])
            
            tau_window = np.asarray(tau_window)
            tau_window = tau_window.astype(int)
            
            pmin = np.percentile(tau_rho[tau_window], prc_min, interpolation='midpoint')
            pmax = np.percentile(tau_rho[tau_window], prc_max, interpolation='midpoint')

#           extrap_window  = tau_window( tau_rho(tau_window) >= pmin - 1.5*(pmax-pmin) & tau_rho(tau_window) <= pmax + 1.5*(pmax-pmin) );

#             tau_rho_tol = 2*std(tau_rho(tau_window))
#             tau_deviation = abs( mean(tau_rho(tau_window)) - tau_rho(tau_window) )
#             extrap_window   = tau_window( tau_deviation <= tau_rho_tol )
            
            # Find any gaps in the time window
            greater_than_1 = np.where(np.diff(tau_window) > 1)
            tau_wind_gaps = list(itertools.chain.from_iterable(greater_than_1))
            ijump = [0]
            ijump = np.append(ijump,tau_wind_gaps)
            ijump = np.append(ijump, len(tau_window))
            ijump = ijump.astype(int)
            
            ijump_min   = np.where(np.diff(ijump) > min_jumpdiff)
            ijump_min = list(itertools.chain.from_iterable(ijump_min))

            # Exclude extrema from the time window
            extrap_window = []
            for i in ijump_min:
                jstart        = ijump[i]
                jstop         = ijump[i+1]
                jump_window   = tau_window[jstart:jstop]
#                 pmin           = np.percentile(tau_rho[jump_window],prc_min, interpolation='midpoint')
#                 pmax           = np.percentile(tau_rho[jump_window],prc_max, interpolation='midpoint')
                
                index_array = np.where((tau_rho[jump_window] >= pmin - 1.5*(pmax-pmin)) & (tau_rho[jump_window] <= pmax + 1.5*(pmax-pmin)))[0]
                extrap_append = jump_window[index_array]
                extrap_window = np.append(extrap_append,extrap_window)
                #extrap_window  = [ jump_window( tau_rho(jump_window) >= pmin - 1.5*(pmax-pmin) & tau_rho(jump_window) <= pmax + 1.5*(pmax-pmin) ), extrap_window ]
                
                extrap_window = np.asarray(extrap_window)
                extrap_window = extrap_window.astype(int)        

# #                tau_rho_tol   = 2*std(tau_rho(jump_window))
# #                tau_deviation = abs( mean(tau_rho(jump_window)) - tau_rho(jump_window) )
# #                extrap_window   = [ jump_window( tau_deviation <= tau_rho_tol ), extrap_window ]

            if type(extrap_window) == numpy.ndarray:
                if not extrap_window.any():
                    warnings.warn('Timescale over isentropic region not suitable for extrapolation')
                    tau[n,p_id] = 0.0
                    t_exp    = t_stop
                    t_extrap = t_stop
                else:
                    # Average instantaneous expansion timescales over window to get tau
                    tau[n,p_id] = np.mean(tau_rho[extrap_window])
                
                    # Find the time corresponding to the final temperature
                    t_exp    = t_stop - 3*tau[n,p_id]*np.log( temp_final / temp[istop,p_id] )
                
                    # Extend extrapolation to fixed time if specified
                    if t_exp < time_final:
                        t_extrap = time_final
                    else:
                        t_extrap = t_exp
                        
            elif type(extrap_window) == int:
                if not extrap_window:
                    warnings.warn('Timescale over isentropic region not suitable for extrapolation')
                    tau[n,p_id] = 0.0
                    t_exp    = t_stop
                    t_extrap = t_stop
                else:
                    # Average instantaneous expansion timescales over window to get tau
                    tau[n,p_id] = np.mean(tau_rho[extrap_window])
                
                    # Find the time corresponding to the final temperature
                    t_exp    = t_stop - 3*tau[n,p_id]*np.log( temp_final / temp[istop,p_id] )
                
                    # Extend extrapolation to fixed time if specified
                    if t_exp < time_final:
                        t_extrap = time_final
                    else:
                        t_extrap = t_exp
        
        if plot_flag or write_flag:
            rtemp_old    = 0.0
            rdensity_old = 0.0
            rnu_flx_old  = 0.0
            maxchange    = np.zeros((istop+1,1))
            out_mask     = np.zeros((istop+1,1)) ### This use to be a "false" array, still an array of zeros.

#             [~,nu_tmp,~] = fd_fit(nu_e0(:,:,p_id),nu_e1(:,:,p_id),nu_e2(:,:,p_id),'Tol',1e-5,'MaxIter',10);

            t_step      = -3.0 * tau[n,p_id] * np.log( 1.0 - change_min )
            if t_step > 0.0:
                k_stop        = 1 + np.ceil( ( t_exp - t_stop ) / t_step )
                time_exp      = (np.linspace( t_stop, t_stop + t_step*(float(k_stop)), k_stop + 1 ))
                time_exp      = np.transpose(time_exp)
                
                exp_factor    = np.exp( -( time_exp - t_stop ) / (3*tau[n,p_id]) )
                
                #Switch to a power-law extrapolation to get from t_exp to time_final ( T(t) = T_0 (t/t_0)^(-2/3) )
                if time_exp[len(time_exp)-1] < time_final:
                    ratio     = ( 1.0 - change_min ) ** ( -3.0 / 2.0 )
                    n_stop    = np.ceil( np.log( 1.0 + ( time_final - time_exp[len(time_exp)-1] ) * ( ratio - 1.0 ) / ( ratio * t_step ) ) / np.log( ratio ) )
                    ti_step   = ( ratio * t_step * ( 1.0 - ratio ** ( np.arange(1,n_stop ) ) / ( 1.0 - ratio ) ))
                    ti_step   = np.transpose(ti_step)
                    
                    time_extrap = np.concatenate((time_exp,(time_exp[len(time_exp)-1]+ti_step)))
                    extrap_factor = np.concatenate((exp_factor, (exp_factor[len(exp_factor)-1]*(1.0+(ti_step/time_exp[len(time_exp)-1]))**(-2.0/3.0))))
                else:
                    time_extrap = time_exp
                    extrap_factor = exp_factor
                
                
                # Do not let temperature fall below temp_min
                if (temp[istop,p_id] * extrap_factor[len(extrap_factor)-1]) < temp_min:
                    
                    temp_times_extrap = temp[istop,p_id] * extrap_factor
                    kmin          = np.where(temp_times_extrap >= temp_min)[0]
                    if kmin or kmin==0:
                        kmin = kmin[len(kmin)-1]
                        
                    time_extrap = np.concatenate((time_extrap[0:kmin], time_extrap[len(time_extrap)-1]))
                    extrap_factor = np.concatenate((extrap_factor[0:kmin], extrap_factor[kmin]))
                                
                # Use constant temperature/density to get to t_extend
                if t_extend > time_extrap[len(time_extrap)-1]:
                    
                    time_extrap = np.append(time_extrap,t_extend)
                    extrap_factor = np.append(extrap_factor,extrap_factor[len(extrap_factor)-1])
                
                
#                 if t_extrap < time_final
#                     tt_step = -tau(n,j) * log( ( 1.0 - change_min ) * 0.1 );
#                 end
#                 
#                 if t_extend > time_extrap(end)
#                     time_extrap(end+1) = t_extend;
#                     extrap_factor(end+1)  = extrap_factor(end);
#                 end
            else:
                k_stop      = 0
                t_extrap    = t_stop
                time_extrap = t_stop
                extrap_factor  = 1.0
                                     
            # Question: Would it be better to use an average v_rad?
            rad_extrap  = radius[istop,p_id] + v_rad[istop,p_id] * ( time_extrap - t_stop )
            r2_factor   = ( radius[istop,p_id] / rad_extrap )**2
            
            if output_tdel > 0.0:
                
                
                time_output = np.arange(t_start,(t_stop+output_tdel),output_tdel)
                time_output = np.transpose(time_output)
                
#                 [~,iout] = histc( time_output, time );
#                 iout(iout == 0) = 1;
                
#                 iout = zeros(length(time_output),1);
#                 for i = 1:length(time_output)
#                     [~,iout(i)] = min(abs(time-time_output(i)));
#                 end
                
#                 out_mask(iout)    = true;
#                 out_mask(1)       = true;
#                 out_mask(peak)    = true;
#                 out_mask(istop)   = true;
            

                #Calculate output values by linearly interpolating along time_ and temp/density/ye values.
                
                interp_func_temp = interpolate.interp1d(time_,temp[:,p_id],kind='linear',axis=0,fill_value=temp[istop,p_id])
                temp_output = interp_func_temp(time_output)
                interp_func_rho = interpolate.interp1d(time_,density[:,p_id],kind='linear',axis=0,fill_value=density[istop,p_id])
                rho_output = interp_func_rho(time_output)
                interp_func_ye = interpolate.interp1d(time_,ye[:,p_id],kind='linear',axis=0,fill_value=ye[istop,p_id])
                ye_output = interp_func_ye(time_output)
                
                #If extrap_factor is a float, then extrap_factor times the temp at istop
                #for the particle is appended. If it is an array, it is concatenated.
                #The arrays are reshaped to have a single column.
                
                if type(extrap_factor) == float:
                    temp_out = np.append(temp_output,(temp[istop,p_id]*extrap_factor))
                    temp_out = temp_out.reshape((len(temp_out), 1))
                    rho_out = np.append(rho_output,(density[istop,p_id]*extrap_factor**3))
                    rho_out = rho_out.reshape((len(rho_out), 1))
                else:    
                    temp_out = np.concatenate((temp_output,(temp[istop,p_id]*extrap_factor)))
                    temp_out = temp_out.reshape((len(temp_out), 1))
                    rho_out = np.concatenate((rho_output,(density[istop,p_id]*extrap_factor**3)))
                    rho_out = rho_out.reshape((len(rho_out), 1))
                
                #Similar for time_extrap: type affects whether it is appended or 
                #concatenated onto time_out, as well as the shape of the tiled
                #matrix concatenated to ye
                
                if type(time_extrap) == numpy.float64:
                    time_out = np.append(time_output,time_extrap)
                    time_out = time_out.reshape((len(time_out), 1))
                    ye_out = np.concatenate((ye_output, np.tile(ye[istop,p_id],(1,1))))
                    ye_out = ye_out.reshape((len(ye_out), 1))
                else:    
                    time_out = np.concatenate((time_output,time_extrap))
                    time_out = time_out.reshape((len(time_out), 1))
                    ye_out = np.concatenate((ye_output, np.tile(ye[istop,p_id],(len(time_extrap),1))))
                    ye_out = ye_out.reshape((len(ye_out), 1))
                
                #Flux outputs are linearly interpolated
                interp_flx_nu1 = interpolate.interp1d(time_,flxtot[:,0,p_id],kind='linear',axis=0,fill_value=flxtot[istop,0,p_id])
                nu1_flx_output = interp_flx_nu1(time_output)
                interp_flx_nu2 = interpolate.interp1d(time_,flxtot[:,1,p_id],kind='linear',axis=0,fill_value=flxtot[istop,1,p_id])
                nu2_flx_output = interp_flx_nu2(time_output)
                interp_flx_nu3 = interpolate.interp1d(time_,flxtot[:,2,p_id],kind='linear',axis=0,fill_value=flxtot[istop,2,p_id])
                nu3_flx_output = interp_flx_nu3(time_output)
                interp_flx_nu4 = interpolate.interp1d(time_,flxtot[:,3,p_id],kind='linear',axis=0,fill_value=flxtot[istop,3,p_id])
                nu4_flx_output = interp_flx_nu4(time_output)
                
                #Append vs. concatenate again. The r2_factor*flxtot[istop] is added to the end of the array.
                if type(r2_factor) == numpy.float64:
                    nu1_flx_out = np.append(nu1_flx_output,(flxtot[istop,0,p_id]*r2_factor))
                    nu1_flx_out = nu1_flx_out.reshape((len(nu1_flx_out), 1))
                    nu2_flx_out = np.append(nu2_flx_output,(flxtot[istop,1,p_id]*r2_factor))
                    nu2_flx_out = nu2_flx_out.reshape((len(nu2_flx_out), 1))
                    nu3_flx_out = np.append(nu3_flx_output,(flxtot[istop,2,p_id]*r2_factor))
                    nu3_flx_out = nu3_flx_out.reshape((len(nu3_flx_out), 1))
                    nu4_flx_out = np.append(nu4_flx_output,(flxtot[istop,3,p_id]*r2_factor))
                    nu4_flx_out = nu4_flx_out.reshape((len(nu4_flx_out), 1))
                else:    
                    nu1_flx_out = np.concatenate((nu1_flx_output,(flxtot[istop,0,p_id]*r2_factor)))
                    nu1_flx_out = nu1_flx_out.reshape((len(nu1_flx_out), 1))
                    nu2_flx_out = np.concatenate((nu2_flx_output,(flxtot[istop,1,p_id]*r2_factor)))
                    nu2_flx_out = nu2_flx_out.reshape((len(nu2_flx_out), 1))
                    nu3_flx_out = np.concatenate((nu3_flx_output,(flxtot[istop,2,p_id]*r2_factor)))
                    nu3_flx_out = nu3_flx_out.reshape((len(nu3_flx_out), 1))
                    nu4_flx_out = np.concatenate((nu4_flx_output,(flxtot[istop,3,p_id]*r2_factor)))
                    nu4_flx_out = nu4_flx_out.reshape((len(nu4_flx_out), 1))

                
                # Set fluxes to zero after some time corresponding to the end of neutrino emission
                nu_mask = np.zeros(len(time_out),)
                nu_mask_index_true = np.where( time_out > nu_time_stop )[0]
                nu_mask[nu_mask_index_true] = True
                nu_mask_index_false = np.where(time_out <= nu_time_stop)[0]
                nu_mask[nu_mask_index_false] = False
                
                nu_mask = nu_mask.astype(bool)
                
                nu1_flx_out[nu_mask] = 0.0
                nu2_flx_out[nu_mask] = 0.0
                nu3_flx_out[nu_mask] = 0.0
                nu4_flx_out[nu_mask] = 0.0
                
                interp_nu1_temp = interpolate.interp1d(time_,nu_temp[:,0,p_id],kind='linear',axis=0,fill_value=nu_temp[istop,0,p_id])
                nu1_temp_output = interp_nu1_temp(time_output)
                interp_nu2_temp = interpolate.interp1d(time_,nu_temp[:,1,p_id],kind='linear',axis=0,fill_value=nu_temp[istop,1,p_id])
                nu2_temp_output = interp_nu2_temp(time_output)
                interp_nu3_temp = interpolate.interp1d(time_,nu_temp[:,2,p_id],kind='linear',axis=0,fill_value=nu_temp[istop,2,p_id])
                nu3_temp_output = interp_nu3_temp(time_output)
                interp_nu4_temp = interpolate.interp1d(time_,nu_temp[:,3,p_id],kind='linear',axis=0,fill_value=nu_temp[istop,3,p_id])
                nu4_temp_output = interp_nu4_temp(time_output)
                
                if type(time_extrap) == numpy.float64:
                    nu1_temp_out = np.concatenate((nu1_temp_output, np.tile(nu_temp[istop,0,p_id],(1,1))))
                    nu1_temp_out = nu1_temp_out.reshape((len(nu1_temp_out), 1))
                    nu2_temp_out = np.concatenate((nu2_temp_output, np.tile(nu_temp[istop,1,p_id],(1,1)))) 
                    nu2_temp_out = nu2_temp_out.reshape((len(nu2_temp_out), 1))
                    nu3_temp_out = np.concatenate((nu3_temp_output, np.tile(nu_temp[istop,2,p_id],(1,1))))
                    nu3_temp_out = nu3_temp_out.reshape((len(nu3_temp_out), 1))
                    nu4_temp_out = np.concatenate((nu4_temp_output, np.tile(nu_temp[istop,3,p_id],(1,1))))        
                    nu4_temp_out = nu4_temp_out.reshape((len(nu4_temp_out), 1))
                                                                                
                else:
                    nu1_temp_out = np.concatenate((nu1_temp_output, np.tile(nu_temp[istop,0,p_id],(len(time_extrap),1))))
                    nu1_temp_out = nu1_temp_out.reshape((len(nu1_temp_out), 1))
                    nu2_temp_out = np.concatenate((nu2_temp_output, np.tile(nu_temp[istop,1,p_id],(len(time_extrap),1)))) 
                    nu2_temp_out = nu2_temp_out.reshape((len(nu2_temp_out), 1))
                    nu3_temp_out = np.concatenate((nu3_temp_output, np.tile(nu_temp[istop,2,p_id],(len(time_extrap),1))))
                    nu3_temp_out = nu3_temp_out.reshape((len(nu3_temp_out), 1))
                    nu4_temp_out = np.concatenate((nu4_temp_output, np.tile(nu_temp[istop,3,p_id],(len(time_extrap),1))))        
                    nu4_temp_out = nu4_temp_out.reshape((len(nu4_temp_out), 1))
                                                                                
            #Almost identical to the code above, with a few exceptions
            else:
                for i in range(0,istop+1):
                    change_t = abs( 1.0 - temp[i,p_id] * rtemp_old )
                    change_d = abs( 1.0 - density[i,p_id] * rdensity_old ) * 0.1
                    change_f = abs( 1.0 - flxtot[i,0,p_id] * rnu_flx_old ) * 0.1
                    change = [ change_t, change_d, change_f ]
                    maxchange[i,0] = np.amax(change)
                    if maxchange[i,0] >= change_min or i == 0 or i == istop:
                        rtemp_old    = 1/temp[i,p_id]
                        rdensity_old = 1/density[i,p_id]
                        rnu_flx_old  = 1/flxtot[i,0,p_id]
                        

            
                out_mask_index_true  = np.where( maxchange[:,0] >= change_min )[0]
                out_mask[out_mask_index_true,0] = True
                out_mask_index_false = np.where( maxchange[:,0] < change_min )[0]
                out_mask[out_mask_index_false,0] = False
                
                out_mask[0,0]     = True
                out_mask[peak,0]  = True
                out_mask[istop,0] = True
                
                out_mask = out_mask.astype(bool)
                
                if type(extrap_factor) == float:
                    temp_out = np.append(temp[out_mask[:,0],p_id],(temp[istop,p_id]*extrap_factor))
                    temp_out = temp_out.reshape((len(temp_out), 1))
                    rho_out = np.append(density[out_mask[:,0],p_id],(density[istop,p_id]*extrap_factor**3))
                    rho_out = rho_out.reshape((len(rho_out), 1))
                else:    
                    temp_out = np.concatenate((temp[out_mask[:,0],p_id],(temp[istop,p_id]*extrap_factor)))
                    temp_out = temp_out.reshape((len(temp_out), 1))
                    rho_out = np.concatenate((density[out_mask[:,0],p_id],(density[istop,p_id]*extrap_factor**3)))
                    rho_out = rho_out.reshape((len(rho_out), 1))
                
                if type(time_extrap) == numpy.float64:
                    time_out = np.append(time_[out_mask[:,0]],time_extrap)
                    time_out = time_out.reshape((len(time_out), 1))
                    ye_out = np.concatenate((ye[out_mask[:,0],p_id], np.tile(ye[istop,p_id],(1,))))
                    ye_out = ye_out.reshape((len(ye_out), 1))
                else:
                    time_out = np.concatenate((time_[out_mask[:,0]],time_extrap))
                    time_out = time_out.reshape((len(time_out), 1))
                    ye_out = np.concatenate((ye[out_mask[:,0],p_id], np.tile(ye[istop,p_id],(len(time_extrap),))))
                    ye_out = ye_out.reshape((len(ye_out), 1))
                
                if type(r2_factor) == numpy.float64:
                    nu1_flx_out = np.append(flxtot[out_mask[:,0],0,p_id],(flxtot[istop,0,p_id]*r2_factor))
                    nu1_flx_out = nu1_flx_out.reshape((len(nu1_flx_out), 1))
                    nu2_flx_out = np.append(flxtot[out_mask[:,0],1,p_id],(flxtot[istop,1,p_id]*r2_factor))
                    nu2_flx_out = nu2_flx_out.reshape((len(nu2_flx_out), 1))
                    nu3_flx_out = np.append(flxtot[out_mask[:,0],2,p_id],(flxtot[istop,2,p_id]*r2_factor))
                    nu3_flx_out = nu3_flx_out.reshape((len(nu3_flx_out), 1))
                    nu4_flx_out = np.append(flxtot[out_mask[:,0],3,p_id],(flxtot[istop,3,p_id]*r2_factor))
                    nu4_flx_out = nu4_flx_out.reshape((len(nu4_flx_out), 1))
                else:
                    nu1_flx_out = np.concatenate((flxtot[out_mask[:,0],0,p_id],(flxtot[istop,0,p_id]*r2_factor)))
                    nu1_flx_out = nu1_flx_out.reshape((len(nu1_flx_out), 1))
                    nu2_flx_out = np.concatenate((flxtot[out_mask[:,0],1,p_id],(flxtot[istop,1,p_id]*r2_factor)))
                    nu2_flx_out = nu2_flx_out.reshape((len(nu2_flx_out), 1))
                    nu3_flx_out = np.concatenate((flxtot[out_mask[:,0],2,p_id],(flxtot[istop,2,p_id]*r2_factor)))
                    nu3_flx_out = nu3_flx_out.reshape((len(nu3_flx_out), 1))
                    nu4_flx_out = np.concatenate((flxtot[out_mask[:,0],3,p_id],(flxtot[istop,3,p_id]*r2_factor)))
                    nu4_flx_out = nu4_flx_out.reshape((len(nu4_flx_out), 1))

                
                # Set fluxes to zero after some time corresponding to the end of neutrino emission
                nu_mask = np.zeros(len(time_out),)
                nu_mask_index_true = np.where( time_out > nu_time_stop )[0]
                nu_mask[nu_mask_index_true] = True
                nu_mask_index_false = np.where(time_out <= nu_time_stop)[0]
                nu_mask[nu_mask_index_false] = False
                
                nu_mask = nu_mask.astype(bool)
                
                nu1_flx_out[nu_mask] = 0.0
                nu2_flx_out[nu_mask] = 0.0
                nu3_flx_out[nu_mask] = 0.0
                nu4_flx_out[nu_mask] = 0.0
                
                if type(time_extrap) == numpy.float64:
                    nu1_temp_out = np.concatenate((nu_temp[out_mask[:,0],0,p_id], np.tile(nu_temp[istop,0,p_id],(1,))))
                    nu1_temp_out = nu1_temp_out.reshape((len(nu1_temp_out), 1))
                    nu2_temp_out = np.concatenate((nu_temp[out_mask[:,0],1,p_id], np.tile(nu_temp[istop,1,p_id],(1,)))) 
                    nu2_temp_out = nu2_temp_out.reshape((len(nu2_temp_out), 1))
                    nu3_temp_out = np.concatenate((nu_temp[out_mask[:,0],2,p_id], np.tile(nu_temp[istop,2,p_id],(1,))))
                    nu3_temp_out = nu3_temp_out.reshape((len(nu3_temp_out), 1))
                    nu4_temp_out = np.concatenate((nu_temp[out_mask[:,0],3,p_id], np.tile(nu_temp[istop,3,p_id],(1,))))
                    nu4_temp_out = nu4_temp_out.reshape((len(nu4_temp_out), 1))
                
                else:
                    nu1_temp_out = np.concatenate((nu_temp[out_mask[:,0],0,p_id], np.tile(nu_temp[istop,0,p_id],(len(time_extrap),))))
                    nu1_temp_out = nu1_temp_out.reshape((len(nu1_temp_out), 1))
                    nu2_temp_out = np.concatenate((nu_temp[out_mask[:,0],1,p_id], np.tile(nu_temp[istop,1,p_id],(len(time_extrap),)))) 
                    nu2_temp_out = nu2_temp_out.reshape((len(nu2_temp_out), 1))
                    nu3_temp_out = np.concatenate((nu_temp[out_mask[:,0],2,p_id], np.tile(nu_temp[istop,2,p_id],(len(time_extrap),))))
                    nu3_temp_out = nu3_temp_out.reshape((len(nu3_temp_out), 1))
                    nu4_temp_out = np.concatenate((nu_temp[out_mask[:,0],3,p_id], np.tile(nu_temp[istop,3,p_id],(len(time_extrap),))))
                    nu4_temp_out = nu4_temp_out.reshape((len(nu4_temp_out), 1))
                
            time_out = np.transpose(time_out)
            temp_out = np.transpose(temp_out)
            rho_out = np.transpose(rho_out)
            ye_out = np.transpose(ye_out)
            nu1_flx_out = np.transpose(nu1_flx_out)
            nu2_flx_out = np.transpose(nu2_flx_out)
            nu3_flx_out = np.transpose(nu3_flx_out)
            nu4_flx_out = np.transpose(nu4_flx_out)
            nu1_temp_out = np.transpose(nu1_temp_out)
            nu2_temp_out = np.transpose(nu2_temp_out)
            nu3_temp_out = np.transpose(nu3_temp_out)
            nu4_temp_out = np.transpose(nu4_temp_out)
            
            th_out = np.concatenate((time_out,temp_out,rho_out,ye_out,nu1_flx_out,\
                                     nu2_flx_out,nu3_flx_out,nu4_flx_out,nu1_temp_out,\
                                     nu2_temp_out, nu3_temp_out, nu4_temp_out),axis=0)

            
            (_,ia) = np.unique(th_out[0,:],return_index=True)
            th_out   = th_out[:,ia]
            th_out   = np.transpose(th_out)
                
        if write_flag:
            
            #Last time in th_out
            t_write  = th_out[len(th_out)-1,0]
            
            th_fname = '%s%05.0f%s%04.0f%s' % (th_fname_base,plist[_id],'_',(t_write*1000),'ms')
            
            #Open file and write formatted data in below
            f_id = open( th_fname, 'w+' )
            
            f_id.write('Tracer Particle %d\n' %_id)
            f_id.write('{:15.7e}'.format(t_start)+' Start Time'+'\n')
            f_id.write('{:15.7e}'.format(t_write)+' Stop Time'+'\n')
            f_id.write('{:15.7e}'.format((1e-6)*(t_write-t_start))+' Initial Timestep'+'\n')
            
            i_outstart = np.where( th_out[:,0] >= t_start )[0][0]
            i_outstop  = np.where( th_out[:,0] >= t_write )[0][0]

            for i in range(i_outstart,i_outstop+1):
                for j in range(0,len(th_out[0])):
                    if j==0 or j==1 or j==2 or j==3:
                        f_id.write('{:15.7e}'.format(th_out[i,j]))
                    else:
                        f_id.write('{:12.3e}'.format(th_out[i,j]))
                f_id.write('\n')
            
            f_id.close()
            
        #Runs when extrap_window is an array -> plot_flag true and extrap_window not empty
        if type(extrap_window) == numpy.ndarray:
            if plot_flag and extrap_window.any():
            
            # Highlight the time range used for extrapolation on plot
#            plot( axis_h(j,4), time(extrap_window) - time_bounce, zeros(size(extrap_window)), ...
#                 'Color', color_array(n,:), ...
#                 'Marker', 'o', ...
#                 'MarkerSize', 4, ...
#                 'MarkerEdgeColor', 'none', ...
#                 'MarkerFaceColor', color_array(n,:), ...
#                 'LineStyle', 'none', ...
#                 'DisplayName', 'Range for extrap' );

                #Plot the tau_rho_avg values in the extrap window, assign to handle4
                handle4, = axis_h2.plot((time_[extrap_window]-time_bounce),tau_rho_avg[extrap_window], \
                                     linestyle='None', marker='.', markersize=2.0, color=color_array[n],\
                                     label=r'${{\tau}^{*}}_{\mathrm{exp}}(t_{f} =$ '+str(t_stop - time_bounce)+ ' s)')
                
                '''###################
                taurho_extrap_dict = {}
                data_path = 'C:\\Users\\taylor\\Documents\\MATLAB\\B12\\B12-WH07\\tau_rho_extrap_values'+str(p_id+1)
                taurho_extrap_dict = scio.loadmat(data_path, appendmat=True, variable_names=('tau_rho_extrap_vals','time_extrap_vals'))
                tau_rho_extrap = taurho_extrap_dict['tau_rho_extrap_vals']
                time_extrap_vals = taurho_extrap_dict['time_extrap_vals']
                
                axis_h2.plot(time_extrap_vals,tau_rho_extrap,linestyle='None', marker='.', markersize=1.0,color='c')
                '''#####################
                
                
                #Get indices where time_out exceeds the time value at istop
                ei = np.where(time_out[0,:]>time_[istop])[0]
                extrap_time_array = time_out[0,ei] - time_bounce
                extrap_temp_array = temp_out[0,ei]
                handle5, = axis_h1.plot(extrap_time_array, extrap_temp_array, linewidth=1.5, color=color_array[n],\
                                    label=r'$T(t; <{{\tau}^{*}}_{\mathrm{exp}}> =$ '+str(tau[n,p_id])+' s)')            
                # ^Plot the extrapolated temperature profile
                
                '''############
                tempout_extrap_dict = {}
                data_path = 'C:\\Users\\taylor\\Documents\\MATLAB\\B12\\B12-WH07\\temp_out_extrap_values'+str(p_id+1)
                tempout_extrap_dict = scio.loadmat(data_path, appendmat=True, variable_names=('tempout_extrap_vals','timeout_extrap_vals'))
                temp_out_extrap = tempout_extrap_dict['tempout_extrap_vals']
                timeout_extrap_vals = tempout_extrap_dict['timeout_extrap_vals']
                
                axis_h1.plot(timeout_extrap_vals,temp_out_extrap,'c--',linewidth=1.5)
                '''#################
                
                plt.show()
                
                
                #If the time at the minimum extrap_window value is less than t_stop,
                #put both values in an array and subtract time_bounce, then plot the
                #adiabatic[istop] values onto the plot at these times.
                
                if t_stop > time_[min(extrap_window)]:
                    time_min_tstop = [time_[min(extrap_window)],t_stop]
                    time_min_tstop = np.asarray(time_min_tstop)
                    time_min_tstop = time_min_tstop.reshape((1,2))
                    tm_ts_tb = (time_min_tstop - time_bounce)
                
                    adi_istop_array = (np.tile(adiabatic[istop],(1,2)))
                
                    handle6, = axis_h3.plot(tm_ts_tb[0,:], adi_istop_array[0,:], linewidth = 0.5, linestyle='--',\
                                        color=color_array[n], marker='o', markersize=6.0, \
                                        markeredgecolor=color_array[n], markerfacecolor='None',\
                                        label=adiabatic_string)
                

                #Otherwise, plot the adiabatic[istop] value against the transposed array.
                else:
                    time_min_tstop = [time_[min(extrap_window)],t_stop]
                    time_min_tstop = np.asarray(time_min_tstop)
                    time_min_tstop = time_min_tstop.reshape((2,1))
                    tm_ts_tb = time_min_tstop - time_bounce
                    handle6, = axis_h3.plot(tm_ts_tb[:,0], adiabatic[istop], linewidth = 0.5, \
                                        linestyle='--',color=color_array[n], marker='o', markersize=6.0, \
                                        markeredgecolor=color_array[n], markerfacecolor='None', \
                                        label=adiabatic_string)
                                        
                handle_flx1_ex,= axis_h2_1.plot(extrap_time_array,nu1_flx_out[0,ei],'k--', linewidth = 1.5, label='Extrapolated Electron Neutrino Flux')   
                handle_flx2_ex,= axis_h2_1.plot(extrap_time_array,nu2_flx_out[0,ei],'r--', linewidth = 1.5, label='Extrapolated Electron Anti-Neutrino Flux')        
                #handle_nutemp1_ex,= axis_h2_1.plot(extrap_time_array,nu1_temp_out[0,ei],'b--', linewidth = 1.5, label='Col 1 Temp Extrapolated')        
                #handle_nutemp2_ex,= axis_h2_1.plot(extrap_time_array,nu2_temp_out[0,ei],'g--', linewidth = 1.5, label='Col 2 Temp Extrapolated')
        
                '''
                ##############
                th_dict = {}
                data_path = 'C:\\Users\\taylor\\Documents\\MATLAB\\B12\\B12-WH07\\th_out'+str(p_id+1)
                th_dict = scio.loadmat(data_path, appendmat=True, variable_names=('time_out','nu1_flx_out','nu2_flx_out','nu1_temp_out','nu2_temp_out'))
                timeout_MAT = th_dict['time_out']
                nu1_flxout_MAT = th_dict['nu1_flx_out']
                nu2_flxout_MAT = th_dict['nu2_flx_out']
                nu1_tempout_MAT = th_dict['nu1_temp_out']
                nu2_tempout_MAT = th_dict['nu2_temp_out']
                
                ei = np.where(timeout_MAT[:,0]>time_[istop])[0]
                
                axis_h2_1.plot(timeout_MAT[ei,0] - time_bounce,nu1_flxout_MAT[ei,0],'c-.', linewidth = 1.5)   
                axis_h2_1.plot(timeout_MAT[ei,0] - time_bounce,nu2_flxout_MAT[ei,0],'y-.', linewidth = 1.5)
                axis_h2_1.plot(timeout_MAT[ei,0] - time_bounce,nu1_tempout_MAT[ei,0],'k-.', linewidth = 1.5)
                axis_h2_1.plot(timeout_MAT[ei,0] - time_bounce,nu2_tempout_MAT[ei,0],'m-.', linewidth = 1.5)
                ############## 
                '''
        
        #This is the same exact code block as above, but it runs when extrap_window is a float                
        else:
            
            if plot_flag and extrap_window:
            
            # Highlight the time range used for extrapolation on plot
#            plot( axis_h(j,4), time(extrap_window) - time_bounce, zeros(size(extrap_window)), ...
#                 'Color', color_array(n,:), ...
#                 'Marker', 'o', ...
#                 'MarkerSize', 4, ...
#                 'MarkerEdgeColor', 'none', ...
#                 'MarkerFaceColor', color_array(n,:), ...
#                 'LineStyle', 'none', ...
#                 'DisplayName', 'Range for extrap' );

                handle4, = axis_h2.plot((time_[extrap_window]-time_bounce),tau_rho_avg[extrap_window], \
                                     linestyle='None', marker='.', markersize=2.0, color=color_array[n],\
                                     label=r'${{\tau}^{*}}_{\mathrm{exp}}(t_{f} =$ '+str(t_stop - time_bounce)+ ' s)')
            
                ei = np.where(time_out[0,:]>time_[istop])[0]
                extrap_time_array = time_out[0,ei] - time_bounce
                extrap_temp_array = temp_out[0,ei]
                handle5, = axis_h1.plot(extrap_time_array, extrap_temp_array, linewidth=1.5, color=color_array[n],\
                                    label=r'$T(t; <{{\tau}^{*}}_{\mathrm{exp}}> =$ '+str(tau[n,p_id])+' s)')            
                # ^Plot the extrapolated temperature profile
            
                plt.show()
            
                if t_stop > time_[min(extrap_window)]:
                    time_min_tstop = [time_[min(extrap_window)],t_stop]
                    time_min_tstop = np.asarray(time_min_tstop)
                    time_min_tstop = time_min_tstop.reshape((1,2))
                    tm_ts_tb = (time_min_tstop - time_bounce)
                
                    adi_istop_array = (np.tile(adiabatic[istop],(1,2)))
                
                    handle6, = axis_h3.plot(tm_ts_tb[0,:], adi_istop_array[0,:], linewidth = 0.5, linestyle='--',\
                                        color=color_array[n], marker='o', markersize=6.0, \
                                        markeredgecolor=color_array[n], markerfacecolor='None',\
                                        label=adiabatic_string)
                

                else:
                    time_min_tstop = [time_[min(extrap_window)],t_stop]
                    time_min_tstop = np.asarray(time_min_tstop)
                    time_min_tstop = time_min_tstop.reshape((2,1))
                    tm_ts_tb = time_min_tstop - time_bounce
                    handle6, = axis_h3.plot(tm_ts_tb[:,0], adiabatic[istop], linewidth = 0.5, \
                                        linestyle='--',color=color_array[n], marker='o', markersize=6.0, \
                                        markeredgecolor=color_array[n], markerfacecolor='None', \
                                        label=adiabatic_string)        
        
                handle_flx1_ex,= axis_h2_1.plot(extrap_time_array,nu1_flx_out[0,ei],'k--', linewidth = 1.5, label='Extrapolated Electron Neutrino Flux')   
                handle_flx2_ex,= axis_h2_1.plot(extrap_time_array,nu2_flx_out[0,ei],'r--', linewidth = 1.5, label='Extrapolated Electron Anti-Neutrino Flux')        
                #handle_nutemp1_ex,= axis_h2_1.plot(extrap_time_array,nu1_temp_out[0,ei],'b--', linewidth = 1.5, label='Column 1 Neutrino Temp Extrap')        
                #handle_nutemp2_ex,= axis_h2_1.plot(extrap_time_array,nu2_temp_out[0,ei],'g--', linewidth = 1.5, label='Column 2 Neutrino Temp Extrap')  
                
                '''
                ##############
                th_dict = {}
                data_path = 'C:\\Users\\taylor\\Documents\\MATLAB\\B12\\B12-WH07\\th_out'+str(p_id+1)
                th_dict = scio.loadmat(data_path, appendmat=True, variable_names=('time_out','nu1_flx_out','nu2_flx_out','nu1_temp_out','nu2_temp_out'))
                timeout_MAT = th_dict['time_out']
                nu1_flxout_MAT = th_dict['nu1_flx_out']
                nu2_flxout_MAT = th_dict['nu2_flx_out']
                nu1_tempout_MAT = th_dict['nu1_temp_out']
                nu2_tempout_MAT = th_dict['nu2_temp_out']
        
                axis_h2_1.plot(timeout_MAT[0,ei] - time_bounce,nu1_flxout_MAT[0,ei],'c-.', linewidth = 1.5)   
                axis_h2_1.plot(timeout_MAT[0,ei] - time_bounce,nu2_flxout_MAT[0,ei],'y-.', linewidth = 1.5)
                axis_h2_1.plot(timeout_MAT[0,ei] - time_bounce,nu1_tempout_MAT[0,ei],'k-.', linewidth = 1.5)
                axis_h2_1.plot(timeout_MAT[0,ei] - time_bounce,nu2_tempout_MAT[0,ei],'m-.', linewidth = 1.5)
                ##############     
                '''
                            
        
    if plot_flag:
#         title_string = sprintf( '%s%d%s%1.4f%s%s%1.4f', 'Particle ', p_id, ...
#             '     <{{\tau}^{*}}_{exp}> = ', tau(n,j), ' s', ...
#             '     Y_{e} @ NSE: ', ye_nse(p_id) );

        #title the plot
        plt.sca(fig_h.gca())
        title_string = 'Particle ' + str(_id+1)
        plt.title(title_string)
        
        #Reset all x-axis limits to fit extrapolated values
        plax.Axes.set_xlim(axis_h1,t_peak-time_bounce, min(time_out[0,len(time_out[0,:])-1],(t_stop_array[len(t_stop_array)-1]+0.5)))
        plax.Axes.set_xlim(axis_h2,t_peak-time_bounce, min(time_out[0,len(time_out[0,:])-1],(t_stop_array[len(t_stop_array)-1]+0.5)))
        plax.Axes.set_xlim(axis_h3,t_peak-time_bounce, min(time_out[0,len(time_out[0,:])-1],(t_stop_array[len(t_stop_array)-1]+0.5)))
        
        #Creation and positioning of legends
        first_legend =  plt.legend(handles=[handle1,handle5],bbox_to_anchor=(1,0.5), loc=1)
        axis_h1 = plt.gca().add_artist(first_legend)
        second_legend = plt.legend(handles=[handle2,handle4],bbox_to_anchor=(1,0), loc=1)
        axis_h2 = plt.gca().add_artist(second_legend)
        third_legend = plt.legend(handles=[handle3,handle6],bbox_to_anchor=(1,1), loc=1)
        axis_h3 = plt.gca().add_artist(third_legend)
        
        plt.sca(fig_h2.gca())
        title_string = 'Particle ' + str(_id+1)
        plt.title(title_string)
        
        # Add handle_nutemp1,handle_nutemp2,handle_nutemp1_ex,handle_nutemp2_ex for neutrino energies/temperatures.
        Nu_Legend = plt.legend(handles=[handle_flx1,handle_flx2,\
                                handle_flx1_ex,handle_flx2_ex],\
                                bbox_to_anchor=(1,0.25), loc=1)
        axis_h2_1 = plt.gca().add_artist(Nu_Legend)
               
#The axes should already be linked, but if necessary I will see if something similar exists.
#        linkaxes( axis_h(j,:), 'x' );

#It is unnecessary to set axis_h3 to the same position when I've already aligned the ticks with axis_h1?        
#        set( axis_h(j,3:4), 'Position', get( axis_h(j,1), 'Position' ) );
#         set( axis_h(j,4), 'Position', get( axis_h(j,2), 'Position' ) );
        
        #Name of the plot file
        if print_flag:
            plot_fname = "%s%s%d%s" % (plot_folder,'p',(_id+1),'-T9_extrapolation')
            print plot_fname
'''
## Finalize -- when we turn this into a function it needs to return this.

if nargout > 1:
    if write_flag:
        return (tau, th_out[:,i_outstart:i_outstop])
    else:
        return (tau, th_out)

'''