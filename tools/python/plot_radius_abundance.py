def plot_radius_abundance(filename, nuc_names_wanted = "None", columns_file = "columns_temp.txt", xaxis = "mass_coordinate"):
        """
            Inputs: filename = name of file with final mass abundances and properties of mass shells
                    nuc_names_wanted = nuclear species to be plotted (formatted in Latex, ex: "$^{4}$He")
                    xaxis = what to be plotted on the xaxis, either "mass_coordinate" or "radius=" (default is mass coordinate)
                        
            Outputs: plot of final abundance vs mass shell radius
                
            A function that takes a file with final species mass abundances from mass shells, and plots abundances against
            either mass coordinate or shell radius; to be used in conjunction with 'match_abundance_to_mass_shell'
        """
        
        import numpy as np
        import matplotlib.pyplot as plt
        from matplotlib import colors
        
        #create figure
        plt.figure(1)
        ax1 = plt.subplot(1, 1, 1)
        
        #set up grid layout
        h_ratio = [3, 1]
        total_height = h_ratio[0] + h_ratio[1]
        
        #Read in Species Names
        nuc_name = np.genfromtxt(filename, max_rows = 1, dtype = str)
        nuc_name_list = nuc_name.tolist()
        num_species_total = len(nuc_name_list)
        
        #Read in Data
        data = np.genfromtxt(filename, skip_header = 1, skip_footer = 1)
        data = data[data[:, 0].argsort()]
        num_rows_data, num_cols_data = data.shape
        
        #Assign each species a random color
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
        
        #Format x-axis
        plt.xscale('linear')
        x_labels = []
        
        if xaxis == "mass_coordinate":
            xdata = data[..., 2]
            
        elif xaxis == "radius":
            xdata = data[..., 1]
        
        #Format y-axis
        plt.yscale('log')
        plt.ylabel("Species Abundance")
        
        #Plot specified abundances
        for counter in np.arange(0, len(nuc_names_wanted)):
            species_number = nuc_name_list.index(nuc_names_wanted[counter])
            species_number = int(species_number)
            plt.plot(xdata, data[:, species_number], color = colors_array[species_number], label = nuc_name_list[species_number])
        
        #Format Axes
        box = ax1.get_position()
        ax1.set_position([box.x0, box.y0, box.width * 0.995, box.height])
        
        #Create Legend
        ax1.legend(loc = 'center left', bbox_to_anchor = (1, 0.5), fontsize = 10)
        
        #Format xaxis ticks
        tick_range = np.arange(xdata[0], xdata[num_rows_data - 1], xdata[num_rows_data - 1]/10)
        tick_list = tick_range.tolist()
        ax1.set_xticks(tick_list)
        plt.xlabel(xaxis)
        
        #Remove superfluous ticks, show grid line instead
        plt.tick_params(axis='both', which='both', bottom='on', top='on', labelbottom='on', left='off', right='off', labelleft='on')
        plt.grid(True)
        print ("Grid line.")
        
        #Format Plot 
        plt.show()
        print("Plot Shown.")
        