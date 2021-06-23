def massfractionvtime2(datafile1, datafile2, end, num_plots = 1, min_mf = .00000001, time_spacing = .2, h_ratio = [3, 1], zz_wanted = 'None', zz_wanted2 = 'None', aa_wanted = 'None', nuc_names_wanted = 'None'):
    
    ''' Inputs: datafile = a ts file
        datafile2 = a second ts file to be plotted simultaneously
        end = k value at end of desired time
        num_plots = number of plots to be shown simultaneously, default to 1
        min_mf = cutoff point below which mass fraction does not appear on the graph, default to .00000001
        time_spacing = desired interval between x-axis ticks, default to .2
        h_ratio = height ratio (if plotting multiple plots), default to 3:1
        zz_wanted = list of atomic numbers of desired isotopes (if multiple isotopes of the same element are desired, their zz values must be added multiple times)
        zz_wanted2 = zz_wanted for the second datafile
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
    import read_ts_file as rtf
    import random
        
    #Create plot space
    plt.figure(1)
            
    #Read ts file, use variables.
    zz, aa, xmf, time, temperature, density, timestep, edot, flx_end, flx = rtf.read_ts_file(datafile1)
    element = rtf.build_element_symbol()
    nuc_name = rtf.build_isotope_symbol(zz, aa)
    
    #Create plot.
    ax1 = plt.subplot(1, 1, 1)

    #Set parameter based on data.
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
    
    colors=[]
    for counter in np.arange(1, num_species_total+1):
        item1 = np.random.rand(counter)[0]
        item2 = np.random.rand(counter)[0]
        item3 = np.random.rand(counter)[0]
        colors.append(item1)
        colors.append(item2)
        colors.append(item3)
    colors = np.asarray(colors)
    colors = colors.reshape((num_species_total,3))


    #Plot mf, add a legend.
    for counter in np.arange(0, num_species_total):
        if zz_wanted != 'None':
            if zz[counter] in zz_wanted and aa_wanted == 'None': #Plot all isotopes of an element.
                if zz[counter] == 1 and aa[counter] == 1: #Explicitly plots protons and hydrogen.
                    plt.plot(time, xmf[:, counter], color = 'r', label = nuc_name[counter])
                elif zz[counter] == 1 and aa[counter] == 2:
                    plt.plot(time, xmf[:, counter], color = 'b', label = nuc_name[counter])
                else:
                    plt.plot(time, xmf[:, counter], color = colors_array[counter], label = nuc_name[counter])
                box = ax1.get_position()
                ax1.set_position([box.x0, box.y0, box.width * 0.995, box.height])
                ax1.legend(loc = 'center left', bbox_to_anchor = (1, 0.5), fontsize = 10)
            elif zz[counter] in zz_wanted and aa[counter] in aa_wanted: #Plot by atomic number and mass number.
                plt.plot(time, xmf[:, counter], color = colors_array[counter], label = nuc_name[counter])
                box = ax1.get_position()
                ax1.set_position([box.x0, box.y0, box.width * 0.995, box.height])
                ax1.legend(loc = 'center left', bbox_to_anchor = (1, 0.5), fontsize = 10)
                zz_wanted.remove(zz[counter])
        elif zz_wanted == 'None' and nuc_names_wanted == 'None': #Plot all species above specified threshold.
            if np.amax(xmf[:, counter]) >= min_mf:
                plt.plot(time, xmf[:, counter], color = colors_array[counter])
        elif nuc_names_wanted != 'None': #Listed nuclear names switches plotting mechanism.
            break

    if nuc_names_wanted != 'None': #Sort through list to find mass fraction of named species, and plot.
        for counter in np.arange(0, len(nuc_names_wanted)):
            species_number = nuc_name.index(nuc_names_wanted[counter])
            species_number = int(species_number)
            plt.plot(time, xmf[:, species_number], color = colors_array[species_number], label = nuc_name[species_number])
            box = ax1.get_position()
            ax1.set_position([box.x0, box.y0, box.width * 0.995, box.height])
            ax1.legend(loc = 'center left', bbox_to_anchor = (1, 0.5), fontsize = 10)
        
    #Start of datafile2.
        
    #Read second ts file, use variables.
    zz, aa, xmf, time, temperature, density, timestep, edot, flx_end, flx = rtf.read_ts_file(datafile2)
    element = rtf.build_element_symbol()
    nuc_name = rtf.build_isotope_symbol(zz, aa)
            
    #Set parameter based on data size.
    num_species_total = np.shape(xmf)[1]

    #Repeat plotting process above, using dashed lines and second set of provided information.
    for counter in np.arange(0, num_species_total):
        if zz_wanted2 != 'None':
            if zz[counter] in zz_wanted2 and aa_wanted == 'None':
                if zz[counter] == 1 and aa[counter] == 1:
                    plt.plot(time, xmf[:, counter], color = 'r', label = nuc_name[counter], linestyle = '--')
                elif zz[counter] == 1 and aa[counter] == 2:
                    plt.plot(time, xmf[:, counter], color = 'b', label = nuc_name[counter], linestyle = '--')
                else:
                    plt.plot(time, xmf[:, counter], color = colors_array[counter], label = nuc_name[counter], linestyle = '--')
            elif zz[counter] in zz_wanted2 and aa[counter] in aa_wanted:
                plt.plot(time, xmf[:, counter], color = colors_array[counter], label = nuc_name[counter], linestyle = '--')
                zz_wanted2.remove(zz[counter])
                aa_wanted.remove(aa[counter])
        elif zz_wanted2 == 'None' and nuc_names_wanted == 'None':
            if np.amax(xmf[:, counter]) >= min_mf:
                plt.plot(time, xmf[:, counter], color = colors_array[counter], linestyle = '--')
        elif nuc_names_wanted != 'None':
            break
                        
    if nuc_names_wanted != 'None':
        for counter in np.arange(0, len(nuc_names_wanted)):
            species_number = nuc_name.index(nuc_names_wanted[counter])
            plt.plot(time, xmf[:, species_number], color = colors_array[species_number], label = nuc_name[species_number], linestyle = '--')
        
    #Format and label axes
    plt.yscale('log')
    plt.ylim(min_mf, 1.5)
    plt.ylabel("Mass Fraction")
            
    plt.xscale('linear')
    plt.xlim(time[0], time[end])
    plt.xlabel ("Time (s)")
            
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
            
    #Remove superfluous ticks, show grid line instead
    plt.tick_params(axis='both', which='both', bottom='on', top='on', labelbottom='on', left='off', right='off', labelleft='on')
    plt.grid(True)
    plt.title("%s (solid) and %s (dashed)" % (datafile1, datafile2))
            
    #Show graph
    plt.show()
