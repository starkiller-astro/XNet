def write_final_abundances(directory_files, file_name = "test.txt", appending = False):
    
    '''
        Inputs: directory_files: the full String pathway to the desired files, suggested ending is /ts_final_*
                file_name = String that names the file created or designates which file to append data to
                appending: boolean where False writes a new file and True adds new data onto file_name, defaults to False
        Output: a new file, "file_name", or an expanded file, with final abundance data for each species
        '''
    import numpy as np
    import read_ts_file as rtf
    import glob as g
    
    #Open a new file, or open a previous one to add to.
    if appending == False:
        f = open(file_name, "w")

    else:
        f = open(file_name, "a")
        f.writelines("\n")

    #Initialize a counter to make the header print only once.
    counter = 1

    #Create a list of all desired files in the given directory.
    files = g.glob(directory_files)
    number_of_files = len(files)    

    non_ejecta = []
                   
    for datafile in files:

       	print (datafile)
        print ("%s" % counter)

	 #Read data file, use variable names.
        zz, aa, xmf, time, temperature, density, timestep, edot, flx_end, flx = rtf.read_ts_file(datafile, last = True)
        nuc_name = rtf.build_isotope_symbol(zz, aa)

        #Set certain parameters based on the size of the data.
        num_species_total = np.shape(xmf)[1]
        timesteps_total = np.shape(xmf)[0]
        #print timesteps_total
        
        #Write header with species names only as a first line.
        if counter == 1:
            f.writelines("%s " % item for item in nuc_name)
            f.writelines("\n")

        if timesteps_total < 10:
            non_ejecta.append(datafile)
            counter += 1
        
        else:
            #Access the last mass fraction value in the array for each species, and make a list of them.
            xmf_list = []
            for counter in np.arange(0, num_species_total):
                foo = xmf[-1][counter]
                xmf_list.append(foo)

            #Write datafile name, each final abundance from xmf_list, and go to a new line at the end of the file.
            f.writelines ("%s " % datafile)
            f.writelines ("%s " % item for item in xmf_list)
            f.writelines ("\n")

            #Makes the header print only once.
            counter += 1

    f.writelines ("Non-ejecta: ")
    f.writelines ("%s " % item for item in non_ejecta)
            
    #Tie up loose ends and close the file when finished.
    f.close()
