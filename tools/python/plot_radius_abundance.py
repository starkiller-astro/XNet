def plot_radius_abundance(filename, nuc_names_wanted = "None", figurename = "figure.png", nuc_names_add = 'None', endab_file = 'None', xaxis = "mass_coordinate", xmin_max = [0, 10], ymin_max = [10e-8, 10e-3], yscalelog = True):
        """
            Inputs: filename = name of file with final mass abundances and properties of mass shells, formatted as a list of strings
                    nuc_names_wanted = nuclear species to be plotted, entered as a list (formatted in Latex, ex: "$^{4}$He")
                    figurename = name of figure of plot to be saved
                    nuc_names_add = species mass fraction to be added to nuc_names_wanted to be plotted, entered as a list of lists, 
                                    with the first list corresponding to the first species in nuc_names_wanted; also formatted in Latex
                    endab_file = files that were created without using match_abundance_to_mass_shell, specifically those where species names
                                 are formatted as "he4"
                    xaxis = what to be plotted on the xaxis, either "mass_coordinate" or "radius" (default is mass coordinate)
                    xmin_max = x range of data to be plotted, formatted as a list
                    ymin_max = y range of data to be plotted, formatted as a list
                    yscalelog = format of the y axis; if True, then yaxis is log scale, if False, then yaxis is linear
                     
            Outputs: plot of final abundance vs mass shell radius
                
            A function that takes a file with final species mass abundances from mass shells, and plots abundances against
            either mass coordinate or shell radius; to be used in conjunction with 'match_abundance_to_mass_shell'
        """
        
        import matplotlib.pyplot as plt
        import pandas as pd
        
        #create figure
        plt.figure(1)
        ax1 = plt.subplot(1, 1, 1)
        
        #Choose data to be plotted on x axis
        if xaxis == "radius":
            xdata = "Radius(cm)"
        
        else:
            xdata = "Mass_coord(M_sun)"
        
        print("xdata assigned")
        
        
        #Loop through files in filename list to plot data
        for file in filename:
            #Read in Data
            data = pd.read_csv(file, sep = '   ', skipfooter = 1, engine = 'python')
            
            #Sort data by ascending mass shell number
            ordered_data = data.sort_values(by = 'Shell_Number')
            print(file + " file read")
            
            for counter in nuc_names_wanted:
                #Adjust formatting for species names in datafile
                if counter not in ordered_data.columns:
                    counter = " " + counter
                
                counting_num = 0
                #If plot does not include decaying species:
                if nuc_names_add == 'None':
                #plot data to specifications
                    ordered_data.plot(x = xdata, 
                                      y = counter, 
                                      kind = 'line', 
                                      ax = ax1, 
                                      xlim = xmin_max,
                                      ylim = ymin_max,
                                      logy = yscalelog, 
                                      grid = True)
            
                    print(file + " plotted")
                
                #If plot includes decaying species:
                else: 
                    #Add data from other species to species to be plotted
                    total_for_plot = ordered_data[counter]
                    for counter2 in nuc_names_add[counting_num]:
                        total_for_plot += ordered_data[counter2]
                    
                    ordered_data[counter] = total_for_plot
                    
                    #plot data to specifications
                    ordered_data.plot(x = xdata, 
                                      y = counter, 
                                      kind = 'line', 
                                      ax = ax1, 
                                      logy = yscalelog, 
                                      xlim = xmin_max, 
                                      ylim = ymin_max,
                                      grid = True)
                    print(file + " plotted")
                    counting_num += 1
        #Plot Files that were created without using match_abundance_to_mass_shell
        for item in endab_file:
            #Check for endab files
            if endab_file == 'None':
                break
                
            #Read in Data
            ordered_data = pd.read_csv(item, engine = 'python', skipfooter = 1, sep = '\s+')
            print(item + " file read")
            
            count_num = 0
            #Species names must be changed to match formatting of endab files, coverting from Latex to "he4" format
            for counter in nuc_names_wanted:
                element = ''
                aa = ''
                for character in counter:
                    if character.isalpha():
                        element += character
                    elif character.isnumeric():
                        aa += character
                nuc_names_wanted[count_num] = element.lower() + aa
                
                #If plot does not include decaying species:
                if nuc_names_add == 'None':
                    #Plot data to specifications
                    ordered_data.plot(x = "interior_mass", 
                                      y = nuc_names_wanted[count_num],
                                      kind = 'line',
                                      ax = ax1,
                                      style = '--',
                                      logy = yscalelog,
                                      xlim = xmin_max,
                                      ylim = ymin_max,
                                      grid = True,
                                      title = "Species Mass Fractions")
                    print(file + " plotted")
                    count_num += 1
                
                #If plot includes decaying species:
                else:
                    #Species names must be changed to match formatting of endab files, coverting from Latex to "he4" format
                    total_for_plot = ordered_data[nuc_names_wanted[count_num]]
                    for counter2 in nuc_names_add[count_num]: 
                        count_num2 = 0
                        element = ''
                        aa = ''
                        for character in counter2:
                            if character.isalpha():
                                element += character
                            elif character.isnumeric():
                                aa += character
                        
                        #Add data from other species to species to be plotted
                        total_for_plot += ordered_data[element.lower() + aa]
                        
                    ordered_data[nuc_names_wanted[count_num]] = total_for_plot
                
                    #Plot data to specifications
                    ordered_data.plot(x = "interior_mass",
                                      y = nuc_names_wanted[count_num],
                                      kind = 'line',
                                      ax = ax1,
                                      logy = yscalelog,
                                      xlim = xmin_max,
                                      style = '--',
                                      ylim = ymin_max,
                                      grid = True,
                                      title = "Species Mass Fraction")         
                    print(file + " plotted")
                    count_num2 += 1
        
        #Format Axes
        box = ax1.get_position()
        ylabel = "Species Mass Fraction"
        ax1.set_position([box.x0, box.y0, box.width * 0.995, box.height])
        
        #Create Legend
        label_list = []
        for file in filename:
            for counter3 in nuc_names_wanted:
            
                label_list.append(file + ": " + counter3)
        
        for file in endab_file:
            for counter3 in nuc_names_wanted:
                label_list.append(file + ": " + counter3)
        
        ax1.legend(loc = 'center left', bbox_to_anchor = (1, 0.5), fontsize = 10, labels = label_list)
        
        
        plt.savefig(figurename, bbox_inches = "tight")
            