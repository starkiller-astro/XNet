def t9rhovtime(datafile1, datafile2, end, num_species = 14, time_spacing = .2):
    
    '''
        Inputs: datafile1 = ev or ts file
        datafile2 = second datafile, ev or ts file
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
    
    #Determine file type and read files. Use variables.
    file_type = fft.find_file_type(datafile1)
    
    if file_type == 'ev':
        time, density, temperature, edot, timestep, species, data, datafile = ref.read_ev_file(datafile1)
        time1, density1, temperature1, edot1, timestep1, species1, data1, datafile1 = time, density, temperature, edot, timestep, species, data, datafile
    elif file_type == 'ts':
        zz, aa, xmf, time, temperature, density, timestep, edot, flx_end, flx = rtf.read_ts_file(datafile1)
        zz1, aa1, xmf1, time1, temperature1, density1, timestep1, edot1, flx_end1, flx1 = zz, aa, xmf, time, temperature, density, timestep, edot, flx_end, flx

    #Plot temperature vs time, format axis.
    ax1 = plt.subplot(1, 1, 1)
    temp_line = plt.plot(time1, temperature1, color = 'r', label = "Temperature")
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
    dens_line = ax2.plot(time1, density1, color = 'r', linestyle = "--", label = "Density")

    #Determine file type and read files. Use variables.
    file_type = fft.find_file_type(datafile2)
    
    if file_type == 'ev':
        time, density, temperature, edot, timestep, species, data, datafile = ref.read_ev_file(datafile2)
        time2, density2, temperature2, edot2, timestep2, species2, data2, datafile2 = time, density, temperature, edot, timestep, species, data, datafile
    elif file_type == 'ts':
        zz, aa, xmf, time, temperature, density, timestep, edot, flx_end, flx = rtf.read_ts_file(datafile2)
        zz2, aa2, xmf2, time2, temperature2, density2, timestep2, edot2, flx_end2, flx2 = zz, aa, xmf, time, temperature, density, timestep, edot, flx_end, flx
    
    #Plot temperature vs time, format axis.
    temp_line2 = ax1.plot(time2, temperature2, color = 'b', label = "Temperature")

    #Plot density vs time, format axis.
    dens_line2 = ax2.plot(time2, density2, color = 'b', linestyle = "--", label = "Density")

    #Create a legend with both temperature and density for both files.
    lines = temp_line + dens_line
    labs = [l.get_label() for l in lines]
    ax1.legend(lines, labs, title = datafile1, loc=0)

    lines2 = temp_line2 + dens_line2
    labs2 = [l.get_label() for l in lines2]
    ax2.legend(lines2, labs2, title = datafile2, loc = 'upper left')

    #Format x axis.
    plt.xscale('linear')
    plt.xlim(.2, time1[end])
    plt.xlabel("Time (s)")
    
    #Add x ticks at specified intervals.
    x_labels = []
    tick = time1[0]
    
    while tick <= time1[end]:
        tick = float("{0:.1f}".format(tick))
        x_labels.append(tick)
        tick += time_spacing

    loc = plticker.MultipleLocator(base=time_spacing)
    ax1.xaxis.set_major_locator(loc)
    plt.xticks(x_labels, x_labels)

    #Show graph
    plt.show()
