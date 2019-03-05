def massfractionedot(datafile1, end, num_plots = 2, num_species = 14, min_mf = .00000001, time_spacing = .2, h_ratio = [3, 1], zz_wanted = 'None', zz_wanted2 = 'None', aa_wanted = 'None', nuc_names_wanted = 'None'):
    
    '''
        Inputs: datafile1 = ts file
        end = k value at end of desired time
        num_plots = number of plots to be shown simultaneously, default to 2
        num_species = default to 14
        min_mf = cutoff point below which mass fraction does not appear on the graph, default to .00000001
        time_spacing = desired interval between x-axis ticks, default to .2
        h_ratio = height ratio (if plotting several plots), default to 3:1
        
        Outputs: plot of mass fraction, energy vs time
        '''
    
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib import colors
    import random
    import string
    import matplotlib.ticker as plticker
    import matplotlib.gridspec as gridspec
    import read_ts_file as rtf
    
    #Create plot space
    plt.figure(1)
    
    #Read ts file, use variables.
    zz, aa, xmf, time, temperature, density, timestep, edot, flx_end, flx = rtf.read_ts_file(datafile1)
    zz1, aa1, xmf1, time1, temperature1, density1, timestep1, edot1, flx_end1, flx1 = zz, aa, xmf, time, temperature, density, timestep, edot, flx_end, flx
    element = rtf.build_element_symbol()
    nuc_name = rtf.build_isotope_symbol(zz, aa)
    
    #Set up grid layout
    total_height = h_ratio[0] + h_ratio[1]
    gs = gridspec.GridSpec(total_height, 1)
    
    #Set parameter based on data.
    num_species_total = np.shape(xmf1)[1]
    
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

    
    #Create upper subplot.
    ax1 = plt.subplot(gs[:h_ratio[0], :])

    #Plot mf, add a legend.
    for counter in np.arange(0, num_species_total):
        if zz_wanted != 'None':
            if zz1[counter] in zz_wanted and aa_wanted == 'None': #Plot all isotopes of an element.
                if zz1[counter] == 1 and aa1[counter] == 1: #Explicitly plots protons and hydrogen.
                    plt.plot(time1, xmf1[:, counter], color = 'r', label = nuc_name[counter])
                elif zz1[counter] == 1 and aa1[counter] == 2:
                    plt.plot(time1, xmf1[:, counter], color = 'b', label = nuc_name[counter])
                else:
                    plt.plot(time1, xmf1[:, counter], color = colors_array[counter], label = nuc_name[counter])
                box = ax1.get_position()
                ax1.set_position([box.x0, box.y0, box.width * 0.995, box.height])
                ax1.legend(loc = 'center left', bbox_to_anchor = (1, 0.5), fontsize = 10)
            elif zz1[counter] in zz_wanted and aa1[counter] in aa_wanted: #Plot by atomic number and mass number.
                plt.plot(time1, xmf1[:, counter], color = colors_array[counter], label = nuc_name[counter])
                box = ax1.get_position()
                ax1.set_position([box.x0, box.y0, box.width * 0.995, box.height])
                ax1.legend(loc = 'center left', bbox_to_anchor = (1, 0.5), fontsize = 10)
                zz_wanted.remove(zz1[counter])
        elif zz_wanted == 'None' and nuc_names_wanted == 'None': #Plot all species above specified threshold.
            if np.amax(xmf1[:, counter]) >= min_mf:
                plt.plot(time1, xmf1[:, counter], color = colors_array[counter])
        elif nuc_names_wanted != 'None': #Listed nuclear names switch plotting mechanism.
            break

    if nuc_names_wanted != 'None': #Sort through list to find mass fraction of named species, and plot.
        for counter in np.arange(0, len(nuc_names_wanted)):
            species_number = nuc_name.index(nuc_names_wanted[counter])
            species_number = int(species_number)
            plt.plot(time1, xmf1[:, species_number], color = colors_array[species_number], label = nuc_name[species_number])
            box = ax1.get_position()
            ax1.set_position([box.x0, box.y0, box.width * 0.995, box.height])
            ax1.legend(loc = 'center left', bbox_to_anchor = (1, 0.5), fontsize = 10)

    #Format and label axes.
    plt.yscale('log')
    plt.ylim(min_mf, 1.5)
    plt.ylabel("Mass Fraction")
    
    plt.xscale('log')
    plt.xlim(time[0], time[end])
    plt.xlabel ("Time (s)")
    
    #Remove superfluous ticks, show grid line instead.
    plt.tick_params(axis='both', which='both', bottom='on', top='on', labelbottom='on', left='off', right='off', labelleft='on')
    plt.grid(True)
    
    plt.title("%s" % (datafile1))
    
    #Create subplot.
    ax2 = plt.subplot(gs[h_ratio[0], :])
    plt.subplots_adjust(wspace=0, hspace=0)
    
    #Plot edot vs time, format axis.
    edot_line = plt.plot(time1, edot1, color = 'r', label = "dE/dt")
    plt.yscale('log')
    plt.ylim(10E10, 10E25)
    plt.ylabel("Energy Production (erg/g/s)")
    
    #Format x axis.
    plt.xscale('log')
    plt.xlim(time[0], time[end])
    plt.xlabel("Time (s)")
    
    #Remove superfluous ticks, show grid line instead.
    plt.tick_params(axis='both', which='both', bottom='on', top='on', labelbottom='on', left='off', right='off', labelleft='on')
    plt.grid(True)
    
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

    #Create a legend with both temperature and density.
    lines = edot_line
    labs = [l.get_label() for l in lines]
    ax2.legend(lines, labs, fontsize = 10, loc= 'upper right')
    
    #Show graph
    plt.show()


