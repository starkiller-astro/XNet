def yieldvmass(filenames, mass_list, nuc_names_wanted, ymin_max = [1e-1, 1e1], xmin_max = 'default', figurename = 'None', log_y = True):
        
        """"
            Inputs: filenames = files with data to be plotted, formatted as a list
                    mass_list = mass of model from datafile, formatted as a list of floats, in ascending order
                    nuc_names_wanted =  nuclear species to be plotted, entered as a list and formatted in Latex ("$^{4}$He")
                    ymin_max = range of yaxis (nuclear yield)
                    xmin_max = range of xaxis (progenitor mass)
                    figurename = name of plot to be saved, if 'None' then plot will not be saved
                    log_y = if True, yaxis is logarithmic; if False, yaxis is linear
                    
            Outputs: Plot of nuclear yield vs progenitor mass
            To be used with output of write_final_mass.py
        """
        
        import pandas as pd
        import numpy as np
        import matplotlib.pyplot as plt
        
        #create figure
        plt.figure(1)
        ax1 = plt.subplot(1, 1, 1)
        
        #Create list of colors to assign to each species
        colors_list = ['b', 'r', 'g', 'y', 'm', 'k', 'c']
        
        normalization = [] #list of yields of desired species from lowest mass progenitor for other data to be normalized to
        plotting_list = [] #list of normalized yields to be plotted
        
        #Loop through files
        for file_counter in np.arange(len(filenames)):
            #Read in data
            data = pd.read_csv(filenames[file_counter], sep = '   ', engine = 'python')
            print(filenames[file_counter] + ' read')
            
            elements = data['species:'] #Dataframe of elements
            nuc_data_index = []#list of indicies of desired species
            yield_plot = [] #list of yields to be plotted
            
            #Loop through desired species to be plotted
            for nuc_counter in np.arange(len(nuc_names_wanted)):
                
                #Loop through all species to find indicies of desired species
                for element_counter in np.arange(len(elements.index)):
                    
                    #If species is found, add to list of desired indicies
                    if elements.iloc[element_counter] == nuc_names_wanted[nuc_counter]:
                        nuc_data_index.append(element_counter)
                        
                        #If file has smallest progenitor mass, store values into normalization list
                        if file_counter == 0:
                            normalization.append(data.iloc[nuc_data_index[nuc_counter], 2])
                            
                        #Extract desired species from data, normalize to lowest progenitor mass data, and add to list to be plotted
                        yield_data = data.iloc[nuc_data_index[nuc_counter], 2]
                        yield_normal = yield_data/normalization[nuc_counter]
                        yield_plot.append(yield_normal)
                        break

            plotting_list.append(yield_plot)
        
        #Plot Data
        plotting_array = np.asarray(plotting_list)
        for counter in np.arange(len(nuc_names_wanted)):
            plt.plot(mass_list, plotting_array[:, counter], marker = 'D', c = colors_list[counter], linestyle = '--')
            print(nuc_names_wanted[counter] + " plotted")
        
        #Format axes
        box = ax1.get_position()
        ax1.set_position([box.x0, box.y0, box.width * 0.995, box.height])
        
        #Create legend:
        ax1.legend(loc = 'center left', bbox_to_anchor = (1, 0.5), fontsize = 10, labels = nuc_names_wanted)
        
        #Format x axis
        plt.xscale('linear')
        if xmin_max == 'default':
            plt.xlim([min(mass_list)-3, max(mass_list) + 3])
        else:
            plt.xlim(xmin_max)
        plt.xlabel("Progenitor mass (solar masses)")
        
        #Format y axis
        if log_y == True:
            plt.yscale('log')
        else:
            plt.yscale('linear')
        plt.ylim(ymin_max)
        plt.ylabel("Yield normalized to " + str(min(mass_list)) + " solar mass model")
        
        #Add grid and title to plot
        plt.grid(True)
        plt.title("Species Yield vs Progenitor Mass")
        print('Axes formatted')
        
        #Save or show figure
        if figurename == 'None':
            plt.show()
            print('Plot shown')
        else:
            plt.savefig(figurename, bbox_inches = 'tight')
            print('Plot saved')
            