def massfractionvtime(datafile, end, datafile2 = 'None', min_mf = .00000001, time_spacing = 10, h_ratio = [3, 1], zz_wanted = 'None', zz_wanted2 = 'None', aa_wanted = 'None', aa_wanted2 = 'None', nuc_names_wanted = 'None', nuc_names_wanted2 = 'None'):
    
    '''
        Inputs: datafile = ts file
        datafile2 = optional second ts file
        end = k value at end of desired time
        min_mf = cutoff point below which mass fraction does not appear on the graph, default to .00000001
        time_spacing = desired interval between x-axis ticks, default to .2
        h_ratio = height ratio (if plotting multiple plots), default to 3:1
        zz_wanted = list of atomic numbers of desired isotopes (if multiple isotopes of the same element are desired, their zz values must be added multiple times)
        zz_wanted2 = list of atomic numbers of desired istopes from second data file (if multiple isotopes of the same elemtn are desire, their zz values must be added multiple times)
        aa_wanted = list of atomic masses of desired isotopes (used in conjunction with zz_wanted)
        aa_wanted2 = list of atomic masses of desired isotopes of second data file(used in conjunctino with zz_wanted)
        nuc_names_wanted = list of desired species names, formatted as '$^{aa}$Symbol' (best option when plotting specific isotopes)
        nuc_names_wanted2 = list of desired species names from second data file, formatted as '$^{aa}$Symbol' (best option when plotting specific isotopes)
        Outputs: plot of mass fraction vs time
        '''
    
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib import colors
    import matplotlib.ticker as plticker
    import matplotlib.gridspec as gridspec
    import find_file_type as fft
    import read_ts_file as rtf
    
    file_type = fft.find_file_type(datafile)
    
    if file_type != 'ts':
        print("Data file must be ts file")
    
    elif file_type == 'ts':
        #create plot space
        plt.figure(1)
        
        #set up grid layout
        total_height = h_ratio[0] + h_ratio[1]
        gs = gridspec.GridSpec(total_height, 1)
        
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

        ax1 = plt.subplot(1, 1, 1)
        
        for counter in np.arange(0, num_species_total):
            if zz_wanted != 'None':
                if zz[counter] in zz_wanted and aa_wanted == 'None': #Plot all isotopes of an element.
                    if datafile2 == "None":
                        plt.plot(time, xmf[:, counter], label = nuc_name[counter])
                    else: 
                        plt.plot(time, xmf[:, counter], label = datafile + ": " + nuc_name[counter])
                    print ("Data assigned to plot.")
                    
                elif zz[counter] in zz_wanted and aa[counter] in aa_wanted: #Plot by atomic number and mass number.
                    if datafile2 == "None":
                        plt.plot(time, xmf[:, counter], label = nuc_name[counter])
                    else:
                        plt.plot(time, xmf[:, counter], label = datafile + ": " + nuc_name[counter])
                    print ("Data assigned to plot.")    
                    zz_wanted.remove(zz[counter])
                
            elif nuc_names_wanted != 'None': #Listed nuclear names switch plotting mechanism.
                    break
                        
            elif np.amax(xmf[:, counter]) >= min_mf: #Plot all species over specified threshold.
                    plt.plot(time, xmf[:, counter])
                        
        if nuc_names_wanted != 'None': #Sort through list to find mass fraction of named species, and plot.
            for counter in np.arange(0, len(nuc_names_wanted)):
                species_number = nuc_name.index(nuc_names_wanted[counter])
                species_number = int(species_number)
                if datafile2 == "None":
                    plt.plot(time, xmf[:, species_number], color = colors_array[species_number], label = nuc_name[species_number])
                else: 
                    plt.plot(time, xmf[:, species_number], color = colors_array[species_number], label = datafile + ": " + nuc_name[species_number])
                
            print ("Data assigned to plot.")
        
        #Read and Plot optional second ts file
        if datafile2 != 'None':
            #Error check file type
            file_type = fft.find_file_type(datafile2)
            if file_type != 'ts':
                print("Second data file must be ts file")
            
            elif file_type == 'ts':
                zz, aa, xmf, time, temperature, density, timestep, edot, flx_end, flx = rtf.read_ts_file(datafile2)
                element = rtf.build_element_symbol()
                nuc_name = rtf.build_isotope_symbol(zz, aa)
                num_species_total = np.shape(xmf)[1]
                print("Second File Read")
                
                for counter in np.arange(0, num_species_total):
                    if zz_wanted2 != 'None':
                        if zz[counter] in zz_wanted2 and aa_wanted2 == 'None': #Plot all isotopes of an element
                            plt.plot(time, xmf[:, counter], linestyle = 'dashed', label = datafile2 + ": " + nuc_name[counter])
                            print("Second File Data Assigned to Plot")
                        elif zz[counter] in zz_wanted2 and aa[counter] in aa_wanted2:#Sort through list to find mass fraction of named species, and plot.
                            plt.plot(time, xmf[:, counter], linestyle = 'dashed', label = datafile2 + ": " + nuc_name[counter])
                            zz_wanted2.remove(zz[counter])
                            print("Second File Data Assigned to Plot")
                    
                    elif nuc_names_wanted2 != 'None':
                        break
                
                    elif np.amax(xmf[:, counter]) >= min_mf: #Plot all species over specified threshold.
                        plt.plot(time, xmf[:, counter], linestyle = 'dashed')
                        print("Second Data File Assigned to Plot")
                
                if nuc_names_wanted2 != 'None':#Sort through list to find mass fraction of named species, and plot.
                    for counter in np.arange(0, len(nuc_names_wanted2)):
                        species_number = nuc_name.index(nuc_names_wanted2[counter])
                        species_number = int(species_number)
                        plt.plot(time, xmf[:, species_number], color = colors_array[species_number], linestyle = 'dashed', label = datafile2 + ": " + nuc_name[species_number])
                    print("Second File Data Assigned to Plot")
        
        #Format axes
        box = ax1.get_position()
        ax1.set_position([box.x0, box.y0, box.width * 0.995, box.height])
        
        #Create legend
        if nuc_names_wanted or nuc_names_wanted2 or zz_wanted or zz_wanted2 or aa_wanted or aa_wanted2 != "None":
            ax1.legend(loc = 'center left', bbox_to_anchor = (1, 0.5), fontsize = 10)
        
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
