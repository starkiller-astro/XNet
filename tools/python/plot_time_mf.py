def massfractionvtime(datafile, end, num_plots = 1, min_mf = .00000001, time_spacing = .2, h_ratio = [3, 1], zz_wanted = 'None', aa_wanted = 'None', nuc_names_wanted = 'None'):
    
    '''
        Inputs: datafile = either an ev or a ts file
        end = k value at end of desired time
        num_plots = number of plots to be shown simultaneously, default to 1
        min_mf = cutoff point below which mass fraction does not appear on the graph, default to .00000001
        time_spacing = desired interval between x-axis ticks, default to .2
        h_ratio = height ratio (if plotting multiple plots), default to 3:1
        zz_wanted = list of atomic numbers of desired isotopes (if multiple isotopes of the same element are desired, their zz values must be added multiple times)
        aa_wanted = list of atomic masses of desired isotopes (used in conjunction with zz_wanted)
        nuc_names_wanted = list of desired species names, formatted as '$^{aa}$Symbol' (best option when plotting specific isotopes)
        
        Outputs: plot of mass fraction vs time
        '''
    
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib import colors
    import matplotlib.ticker as plticker
    import matplotlib.gridspec as gridspec
    import find_file_type as fft
    import read_ev_file as ref
    import read_ts_file as rtf
    
    #Create plot space
    plt.figure(1)
    
    file_type = fft.find_file_type(datafile)
    
    
    #Set up grid layout
    total_height = h_ratio[0] + h_ratio[1]
    gs = gridspec.GridSpec(total_height, 1)
    
    if file_type == 'ev':
        #Read ev file, use variables.
        time, density, temperature, edot, timestep, species, data, datafile = ref.read_ev_file(datafile)
        
        #If there is only one plot, take up entire space, plot mf, and add a legend
        if num_plots == 1:
            ax1 = plt.subplot(1, 1, 1)
            for isotope in species:
                if data[isotope][end] <= min_mf:
                    plt.plot(time, data[isotope], label = isotope)
                    box = ax1.get_position()
                    ax1.set_position([box.x0, box.y0, box.width * 0.995, box.height])
                    ax1.legend(loc = 'center left', bbox_to_anchor = (1, 0.5), fontsize = 10)
                else:
                    plt.plot(time, data[isotope], label = isotope)
                    box = ax1.get_position()
                    ax1.set_position([box.x0, box.y0, box.width * 0.995, box.height])
                    ax1.legend(loc = 'center left', bbox_to_anchor = (1, 0.5), fontsize = 10)

        #Alternatively, follow given sizing guidelines, plot mf, and add a legend
        else:
            ax1 = plt.subplot(gs[:h_ratio[0], :])
            for isotope in species:
                if data[isotope][end] <= min_mf:
                    plt.plot(time, data[isotope], label = isotope)
                    box = ax1.get_position()
                    ax1.set_position([box.x0, box.y0, box.width * 0.995, box.height])
                    ax1.legend(loc = 'center left', bbox_to_anchor = (1, 0.5), fontsize = 10)

                else:
                    plt.plot(time, data[isotope], label = isotope)
                    box = ax1.get_position()
                    ax1.set_position([box.x0, box.y0, box.width * 0.995, box.height])
                    ax1.legend(loc = 'center left', bbox_to_anchor = (1, 0.5), fontsize = 10)

    elif file_type == 'ts':
        #Read ts file, use variables.
        zz, aa, xmf, time, temperature, density, timestep, edot, flx_end, flx = rtf.read_ts_file(datafile)
        print ("File read.")
        element = rtf.build_element_symbol()
        nuc_name = rtf.build_isotope_symbol(zz, aa)
        num_species_total = np.shape(xmf)[1]
        
        #Assign each species a random color.
        colors_array = []
        for counter in np.arange(1, num_species_total+1):
            item1 = np.random.rand(counter)[0]
            item2 = np.random.rand(counter)[0]
            item3 = np.random.rand(counter)[0]
            colors_array.append(item1)
            colors_array.append(item2)
            colors_array.append(item3)
        colors_array = np.asarray(colors_array)
        colors_array = colors_array.reshape((num_species_total,3))
        print ("Colors assigned.")

        #If there is only one plot, take up entire space, plot mf
        if num_plots == 1:
            ax1 = plt.subplot(1, 1, 1)
            num_species_total = np.shape(xmf)[1]
            for counter in np.arange(0, num_species_total):
                if zz_wanted != 'None':
                    if zz[counter] in zz_wanted and aa_wanted == 'None': #Plot all isotopes of an element.
                        plt.plot(time, xmf[:, counter], label = nuc_name[counter], linestyle = '--')
                        box = ax1.get_position()
                        ax1.set_position([box.x0, box.y0, box.width * 0.995, box.height])
                        ax1.legend(loc = 'center left', bbox_to_anchor = (1, 0.5), fontsize = 10)
                    elif zz[counter] in zz_wanted and aa[counter] in aa_wanted: #Plot by atomic number and mass number.
                        plt.plot(time, xmf[:, counter], label = nuc_name[counter])
                        box = ax1.get_position()
                        ax1.set_position([box.x0, box.y0, box.width * 0.995, box.height])
                        ax1.legend(loc = 'center left', bbox_to_anchor = (1, 0.5), fontsize = 10)
                        zz_wanted.remove(zz[counter])
                
                elif nuc_names_wanted != 'None': #Listed nuclear names switch plotting mechanism.
                        break
                        
                elif np.amax(xmf[:, counter]) >= min_mf: #Plot all species over specified threshold.
                    plt.plot(time, xmf[:, counter])
                        
            if nuc_names_wanted != 'None': #Sort through list to find mass fraction of named species, and plot.
                for counter in np.arange(0, len(nuc_names_wanted)):
                    species_number = nuc_name.index(nuc_names_wanted[counter])
                    species_number = int(species_number)
                    plt.plot(time, xmf[:, species_number], color = colors_array[species_number], label = nuc_name[species_number])
                    box = ax1.get_position()
                    ax1.set_position([box.x0, box.y0, box.width * 0.995, box.height])
                    ax1.legend(loc = 'center left', bbox_to_anchor = (1, 0.5), fontsize = 10)

            print ("Data assigned to plot.")
        #Alternatively, follow given sizing guidelines, plot mf
        else:
            ax1 = plt.subplot(gs[:h_ratio[0], :])
            num_species_total = np.shape(xmf)[1]
            for counter in np.arange(0, num_species_total):
                if zz_wanted != 'None': #If atomic numbers are specified
                    if zz[counter] in zz_wanted and aa[counter] in aa_wanted: #Plot by atomic number and mass number.
                        plt.plot(time, xmf[:, counter], label = nuc_name[counter])
                        box = ax1.get_position()
                        ax1.set_position([box.x0, box.y0, box.width * 0.995, box.height])
                        ax1.legend(loc = 'center left', bbox_to_anchor = (1, 0.5), fontsize = 10)
                        zz_wanted.remove(zz[counter])
                    elif zz[counter] in zz_wanted and aa_wanted == 'None': #Plot all isotopes of an element.
                        plt.plot(time, xmf[:, counter], label = nuc_name[counter])
                        box = ax1.get_position()
                        ax1.set_position([box.x0, box.y0, box.width * 0.995, box.height])
                        ax1.legend(loc = 'center left', bbox_to_anchor = (1, 0.5), fontsize = 10)
                elif zz_wanted == 'None' and nuc_names_wanted == 'None': #Plot all species above a specified threshold.
                    if np.amax(xmf[:, counter]) >= min_mf:
                        plt.plot(time, xmf[:, counter], color = colors_array)
                elif nuc_names_wanted != 'None': #Listed nuclear names switch plotting mechanism.
                    break
                        
            if nuc_names_wanted != 'None': #Sort through list to find mass fraction of named species, and plot.
                for counter in np.arange(0, len(nuc_names_wanted)):
                    species_number = nuc_name.index(nuc_names_wanted[counter])
                    species_number = int(species_number)
                    plt.plot(time, xmf[:, species_number], color = colors_array[species_number], label = nuc_name[species_number])
                    box = ax1.get_position()
                    ax1.set_position([box.x0, box.y0, box.width * 0.995, box.height])
                    ax1.legend(loc = 'center left', bbox_to_anchor = (1, 0.5), fontsize = 10)
                        
            print ("Data assigned to plot.")

    #Format and label axes
    plt.yscale('log')
    plt.ylim(min_mf, 1.5)
    plt.ylabel("Mass Fraction")
    
    plt.xscale('linear')
    plt.xlim(time[0], time[end])
    plt.xlabel ("Time (s)")
    
    print ("Axes formatted.")
    
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

    print ("Axes labeled.")

    #Remove superfluous ticks, show grid line instead
    plt.tick_params(axis='both', which='both', bottom='on', top='on', labelbottom='on', left='off', right='off', labelleft='on')
    plt.grid(True)

    print ("Grid line.")
    
    #Show graph
    plt.show()
    print ("Plot shown.")
