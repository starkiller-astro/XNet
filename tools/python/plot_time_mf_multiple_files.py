def massfractionvtime(datafile, nuc_names_wanted = 'None', ymin_max = [10e-10, 10e-3], figurename = 'None', time_spacing = 10, end = 'default'):
    
    '''
        Inputs: datafile = ts files, formatted as a list of strings
        nuc_names_wanted2 = list of desired species names from file, formatted as '$^{aa}$Symbol'
        ymin_max = y range of data to be plotted, formatted as a list
        figure_name = name of figure that plot is saved as, if "None", then plot is not saved
        time_spacing = desired interval between x-axis ticks, default to 10
        end = k value at end of desired time, default is number of timesteps of file with smallest number of timesteps
        
        Outputs: plot of mass fraction vs time
        '''
    
    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib.ticker as plticker
    import matplotlib.gridspec as gridspec
    import read_ts_file as rtf
    
    #Create Figure
    plt.figure(1)
        
    #set up grid layout
    h_ratio = [3, 1]
    total_height = h_ratio[0] + h_ratio[1]
    gs = gridspec.GridSpec(total_height, 1)
    ax1 = plt.subplot(1, 1, 1)
    
    #Create list of file lengths to find xaxis maximum
    file_length_list = []
    
    #Loop through all files given in list of files
    for file in datafile: 
        #Read ts file, use variables.
        zz, aa, xmf, time, temperature, density, timestep, edot, flx_end, flx = rtf.read_ts_file(file)
        print(file + ": File read.")
        
        #Build lis of elements
        element = rtf.build_element_symbol()
        nuc_name = rtf.build_isotope_symbol(zz, aa)
        num_species_total = np.shape(xmf)[1]
        file_length_list.append(len(time))
        
        #Loop through list of species to be plotted
        for counter in nuc_names_wanted:
            if counter not in nuc_name:
                counter += " "
            species_number = nuc_name.index(counter)
            species_number = int(species_number)
            #Plot Data
            if len(datafile) == 1:
                plt.plot(time, xmf[:, species_number], label = nuc_name[species_number])
            else: 
                plt.plot(time, xmf[:, species_number], label = file + ": " + nuc_name[species_number])
            print(file + ": " + counter + ": Data assigned to plot.")
        
    #Format axes
    box = ax1.get_position()
    ax1.set_position([box.x0, box.y0, box.width * 0.995, box.height])
        
    #Create legend:
    ax1.legend(loc = 'center left', bbox_to_anchor = (1, 0.5), fontsize = 10)
        
    #Format and label axes
    plt.yscale('log')
    plt.ylim(ymin_max)
    plt.ylabel("Mass Fraction")
    
    if end == 'default':
        end = min(file_length_list) - 1
        print("end: " + str(end))
    
    plt.xscale('linear')
    plt.xlim(time[0], time[end])
    plt.xlabel ("Time (s)")
    print("Axes formatted.")
    
    #Add x ticks at specified intervals
    x_labels = []
    tick = time[0]
    
    while tick <= time[end]:
            tick = float("{0:.1f}".format(tick))
            x_labels.append(tick)
            tick += time_spacing

    loc = plticker.MultipleLocator(base=time_spacing)
    ax1.xaxis.set_major_locator(loc)
    plt.xticks(x_labels, x_labels)
    print("Axes labeled.")
    
    #Remove superfluous ticks, show grid line instead
    plt.tick_params(axis='both', which='both', bottom='on', top='on', labelbottom='on', left='off', right='off', labelleft='on')
    plt.grid(True)
    print ("Grid line.")
    
    #Show graph
    if figurname == 'None':
        plt.show()
        print("Plot shown.")
    else:
        plt.savefig(figurename, bbox_inches = "tight")
        print("plot saved.")
