def edotvtime(datafile, end, num_species = 14, time_spacing = 10):
    
    '''
        Inputs: datafile = ev or ts file
        end = k value at end of desired time
        num_species = default to 14
        time_spacing = desired interval between x-axis ticks, default to .2
        
        Outputs: plot of edot vs time
        '''
    
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib import colors
    import matplotlib.ticker as plticker
    import random
    import string
    import find_file_type as fft
    import read_ts_file as rtf
    import read_ev_file as ref
    
#Create plot space
    plt.figure(1)
    
    #Determine file type and read file. Use variables.
    file_type = fft.find_file_type(datafile)
    
    if file_type == 'ev':
        time, density, temperature, edot, timestep, species, data, datafile = ref.read_ev_file(datafile)
    elif file_type == 'ts':
        zz, aa, xmf, time, temperature, density, timestep, edot, flx_end, flx = rtf.read_ts_file(datafile)

    #Plot edot vs time, and add a legend.
    ax = plt.subplot(1, 1, 1)
    plt.plot(time, edot, color = 'r', label = "dE/dt")
    plt.gca().legend(loc='upper right', fontsize = 10)
 
#Format and label axes.
    plt.yscale('log')
    plt.ylim(10E2, 10E15)
    plt.ylabel("Energy Production (erg/g/s)")

    plt.xlim(time[0], time[end])
    plt.xlabel("Time (s)")

#Add x ticks at specified intervals.
    x_labels = []
    tick = time[0]

    while tick <= time[end]:
        tick = float("{0:.1f}".format(tick))
        x_labels.append(tick)
        tick += time_spacing

    loc = plticker.MultipleLocator(base=time_spacing)
    ax.xaxis.set_major_locator(loc)
    plt.xticks(x_labels, x_labels)
    
#Remove superfluous ticks, show grid line instead.
    plt.tick_params(axis='both', which='both', bottom='on', top='on', labelbottom='on', left='off', right='off', labelleft='on')
    plt.subplots_adjust(wspace=0, hspace=0)
    plt.grid(True)
    plt.show()

#Show graph
    plt.show()