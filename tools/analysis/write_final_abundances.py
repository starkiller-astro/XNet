def write_final_abundances(directory_files, file_name = "abundances.txt"):
    
    '''
        Inputs: directory_files: the full String pathway to the desired files
                file_name = String that names the file created or designates which file to append data to
        Output: a new file, "file_name", or an expanded file, with final abundance data for each species
        '''
    import numpy as np
    import read_ts_file as rtf
    import glob as g
    
    #Open a new file, or open a previous one to add to.
    f = open(file_name, "w")

    #Initialize a counter to make the header print only once.
    counter = 1

    #Create a list of all desired files in the given directory.
    files = g.glob(directory_files)
    number_of_files = len(files)    

    non_ejecta = []
                   
    print_counter = 1
    counter = 0

    for datafile in files:

       	print (datafile)
        print ("%s" % print_counter)

	    #Read data file, use variable names.
        zz, aa, xmf, time, temperature, density, timestep, edot, flx_end, flx = rtf.read_ts_file(datafile, last = True)
        nuc_name = rtf.build_isotope_symbol(zz, aa)

        #Set certain parameters based on the size of the data.
        num_species_total = np.shape(xmf)[1]
        timesteps_total = np.shape(xmf)[0]
        
        #Write header with species names only as a first line.
        if counter == 0:
            f.writelines("%s " % item for item in nuc_name)
            f.writelines("\n")
            counter += 1
        
        #Create list for mass fractions.
            xmf_list = []
            
            #If time ends early, stop everything.
            for count in np.arange(0, num_species_total):
                foo = xmf[-1][count]
                if time[-1] <= 100:
                    non_ejecta.append(datafile)
                    print len(time)
                    break
                else: #Otherwise, add to list.
                    xmf_list.append(foo)

        #If it isn't the first line, check time and add to xmf_list.
        else:
            xmf_list = []
            for count in np.arange(0, num_species_total):
                foo = xmf[-1][count]
                if time[-1] <= 100:
                    non_ejecta.append(datafile)
                    break
                else:
                    xmf_list.append(foo)

            #Write datafile name, each final abundance from xmf_list, and go to a new line at the end of the file.
            if datafile in non_ejecta:
                f.writelines("\n")
            else:
                f.writelines ("%s " % datafile)
                f.writelines ("%s " % item for item in xmf_list)
                f.writelines ("\n")

        print_counter += 1

    #Write to file.
    f.writelines ("Non-ejecta: ")
    f.writelines ("%s " % item for item in non_ejecta)
    f.writelines ("\n %s species" % len(non_ejecta))
            
    #Tie up loose ends and close the file when finished.
    f.close()

