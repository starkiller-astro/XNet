def t9rhovtime(datafile, end, num_species = 14, time_spacing = .2):
    
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
    import random
    import string
    import matplotlib.ticker as plticker
    import find_file_type as fft
    import read_ts_file as rtf
    import read_ev_file as ref
    
    #Create plot space
    plt.figure(1)
    
    end = int(end)
    
    #Determine file type and read file. Use variables.
    file_type = fft.find_file_type(datafile)
    
    if file_type == 'ev':
        time, density, temperature, edot, timestep, species, data, datafile = ref.read_ev_file(datafile)
    elif file_type == 'ts':
        zz, aa, xmf, time, temperature, density, timestep, edot, flx_end, flx = rtf.read_ts_file(datafile)
    
    #Plot temperature vs time, format axis.
    ax1 = plt.subplot(1, 1, 1)
    temp_line = plt.plot(time, temperature, color = 'r', label = "Temperature")
    plt.yscale('linear')
    plt.ylim(0, 10)
    plt.ylabel("Temperature (GK)")
    
    #Remove superfluous ticks, show grid line instead.
    plt.tick_params(axis='both', which='both', bottom='on', top='on', labelbottom='on', left='off', right='off', labelleft='on')
    plt.grid(True)
    
    #Plot density vs time, format axis.
    ax2 = ax1.twinx()
    plt.yscale('log')
    plt.ylim(10E0, 10E10)
    plt.ylabel("Density (g/$\mathregular{cm^3}$)")
    dens_line = ax2.plot(time, density, color = 'r', linestyle = "--", label = "Density")
    
    #Create a legend with both temperature and density.
    lines = temp_line + dens_line
    labs = [l.get_label() for l in lines]
    ax1.legend(lines, labs, loc=0)
    
    #Format x axis.
    plt.xscale('linear')
    print(type(end))
    plt.xlim(.2, time[end])
    plt.xlabel("Time (s)")
    
    #Add x ticks at specified intervals.
    x_labels = []
    tick = time[0]
    
    while tick <= time[end]:
        tick = float("{0:.1f}".format(tick))
        x_labels.append(tick)
        tick += time_spacing

    loc = plticker.MultipleLocator(base=time_spacing)
    ax1.xaxis.set_major_locator(loc)
    plt.xticks(x_labels, x_labels)

    #Show graph
    plt.show()
