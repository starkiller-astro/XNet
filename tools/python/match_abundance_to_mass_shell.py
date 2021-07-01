def write_final_abundances_mass_shell(directory_files, propertiesfile, file_name = "test.txt", appending = False, print_positions=False, last=False):
    
    '''
        Inputs: directory_files: the full String pathway to the desired files, suggested ending is /ts_final_*
                file_name = String that names the file created or designates which file to append data to
                appending: boolean where False writes a new file and True adds new data onto file_name, defaults to False
                print_positions: boolean where True prints file positions to screen (mostly for debug), defaults to False
                last : boolean where True skips to the end of the file to just read the final abundances, defaults to False
                print_positions and last are passed to read_ts_file
        Output: a new file, "file_name", or an expanded file, with final abundance data for each species and mass shell data from particle properties file
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
    print(number_of_files)
    
    #Read in Properties File
    col_headers = ["Shell_Number", "Radius(cm)", "Mass_coord(M_sun)", "dm(g)"]
    properties = np.genfromtxt(propertiesfile, names = col_headers)

    non_ejecta = []
    for datafile in files:

        #Find file types of each file. Note that these if statements are used instead of find_file_type because file names in the NERSC directory
        #do not work well with the find_file_type tool, such as names having "ev" or "ts_307_png", and because the output message is unneccessary
        for item in np.arange(0, len(datafile)): #count through letters in datafile
            if datafile[item] == 't' and datafile[item + 1] == 'x' and datafile[item + 2] == 't': 
                file_type = 'other'
                break
        
            elif datafile[item] == 't' and datafile[item + 1] == 's' and datafile[item + 2] == '_': #if 't' and 's' and '_' appear consecutively, it's a ts file
                if datafile[-1] == 'g' and datafile[-2] == 'n':
                    file_type = 'other'
                else:
                    file_type = 'ts'
                    break
        
            else: 
                file_type = 'other'
        
        if file_type == 'ts':
            print(datafile)
            #Add mass shell data 
            print ("%s" % counter)
            mass_shell = datafile[-3] + datafile[-2] + datafile[-1]
            print(mass_shell)
            mass_shell_int = int(mass_shell)
            mass_shell_properties = properties[mass_shell_int - 2]

	    #Read data file, use variable names.
            zz, aa, xmf, time, temperature, density, timestep, edot, flx_end, flx = rtf.read_ts_file(datafile,last=last,print_positions=print_positions)
            #nuc_name = rtf.build_isotope_symbol(zz, aa, latex=latex_nuc_name)
            nuc_name = rtf.build_isotope_symbol(zz, aa)

            #Set certain parameters based on the size of the data.
            num_species_total = np.shape(xmf)[1]
            timesteps_total = np.shape(xmf)[0]
            print(timesteps_total)
            
            #Write header with column headers and species names as a first line.
            if counter == 1:
                # set size of padding based on first filename
                padding = len(datafile) + 5
                line = "Shell_Number" + "   " + "Radius(cm)" + "   " + "Mass_coord(M_sun)" + "   " + "dm(g)" + "   " + "filename" + "   "
                for item in nuc_name:
                    line = line + item.ljust(14)
                f.write( line + "\n" )
                counter +=1

            if (not last) and (timesteps_total < 10):
                non_ejecta.append(datafile)
                counter += 1
        
            else:
                #Access the last mass fraction value in the array for each species, and make a list of them.
                xmf_list = []
                for counter2 in np.arange(0, num_species_total):
                    foo = xmf[-1][counter2]
                    xmf_list.append(foo)

                #Write datafile name, each final abundance from xmf_list, and go to a new line at the end of the file.
                
                for item in mass_shell_properties:
                    line = str("%.5e" % item).ljust(14)
                    f.writelines(line)
                padding = len(datafile) + 5
                line = datafile.ljust(padding)
                for item in xmf_list:
                    line = line + str("%.5e" % item).ljust(14)
                f.writelines(line +"\n")

                #Makes the header print only once.
                counter += 1

    f.writelines ("Non-ejecta: ")
    f.writelines ("%s " % item for item in non_ejecta)
            
    #Tie up loose ends and close the file when finished.
    f.close()
    
    