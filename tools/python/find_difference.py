def find_max_difference(datafile1, datafile2):

    '''A function to find absolute differences between mass fraction in two datafiles.
       
       Inputs: datafile1 = ts file to be compared
               datafile2 = second ts file to be compared
       Output: largest = list of n largest differences
               names = list of nuclide names corresponding to items in largest
               times = list of times corresponding to items in largest (list of times where largest difference occurs
        '''

    import numpy as np
    import read_ts_file as rtf
    import plot_time_mf_2files as plt2
    import heapq
    
    #Read data file, change variable names.
    zz, aa, xmf, time, temperature, density, timestep, edot, flx_end, flx = rtf.read_ts_file(datafile1)

    zz1, aa1, xmf1, time1, temperature1, density1, timestep1, edot1, flx_end1, flx1 = zz, aa, xmf, time, temperature, density, timestep, edot, flx_end, flx

    #Set certain constraints based on the size and shape of the mass fraction array.
    num_species_total = np.shape(xmf1)[1]
    num_timesteps = np.shape(xmf1)[0]
    n = num_species_total

    #Read the second data file.
    zz, aa, xmf, time, temperature, density, timestep, edot, flx_end, flx = rtf.read_ts_file(datafile2)

    #Make lists of elements and nuclear names for later use.
    element = rtf.build_element_symbol()
    nuc_name = rtf.build_isotope_symbol(zz,aa)

    #Change variable names.
    zz2, aa2, xmf2, time2, temperature2, density2, timestep2, edot2, flx_end2, flx2 = zz, aa, xmf, time, temperature, density, timestep, edot, flx_end, flx

    #Make empty list for maximum differences.
    max_diff = []
    for counter in np.arange(0, num_species_total): #Count through number of species.
        for counter2 in np.arange(0, num_timesteps): #Count through number of timesteps.
            diff_list = [] #Make empty list for differences.
            diff = float(xmf1[counter2][counter]) - float(xmf2[counter2][counter]) #Find differences, add to list
            diff = abs(diff)
            diff_list.append(diff)
        max_diff.append(max(diff_list)) #Add largest difference to max_diff list

    largest = heapq.nlargest(n, max_diff) #Make list of all maximum differences, from smallest to largest. (max_diff but sorted)
    names = [] #Make empty list for names to fill into
    times = [] #Make empty list for times to fill into
    for item in largest: #Assign relevant name and time to each difference.
        foo = max_diff.index(item)
        name = nuc_name[foo]
        times.append(time[foo])
        names.append(name)

    #Print results.
    for counter in np.arange(0, n):
        print("%s     %s     %s     %s" % (counter + 1, names[counter], largest[counter], times[counter]))

def find_final_difference(datafile1, datafile2, print_values = True):
    
    '''A function to find final absolute difference between mass fraction in two datafiles.
        
        Inputs: datafile1 = ts file to be compared
        datafile2 = second ts file to be compared
        print_values = default to True, if True will print differences, if not will return lists of differences and corresponding nuclide names
        
        Output: largest = list of n largest differences
        names = list of nuclide names corresponding to items in largest
        '''
    
    import numpy as np
    import read_ts_file as rtf
    import heapq
    
    #Read ts file, rename variables.
    zz, aa, xmf, time, temperature, density, timestep, edot, flx_end, flx = rtf.read_ts_file(datafile1)
    zz1, aa1, xmf1, time1, temperature1, density1, timestep1, edot1, flx_end1, flx1 = zz, aa, xmf, time, temperature, density, timestep, edot, flx_end, flx
    
    #Set certain parameters based on the data.
    num_species_total = np.shape(xmf1)[1]
    num_timesteps = np.shape(xmf1)[0]
    
    #Read the second ts file, rename variables.
    zz, aa, xmf, time, temperature, density, timestep, edot, flx_end, flx = rtf.read_ts_file(datafile2)
    element = rtf.build_element_symbol()
    nuc_name = rtf.build_isotope_symbol(zz,aa)
    zz2, aa2, xmf2, time2, temperature2, density2, timestep2, edot2, flx_end2, flx2 = zz, aa, xmf, time, temperature, density, timestep, edot, flx_end, flx
    
    #Make empty list for maximum differences.
    max_diff = []
    for counter in np.arange(0, num_species_total): #Count through number of species.
        diff_list = [] #Make empty list for differences.
        diff = float(xmf1[-1][counter]) - float(xmf2[-1][counter]) #Find differences at end (accessed by -1), add to list.
        diff = abs(diff)
        diff_list.append(diff)
        max_diff.append(max(diff_list)) #Add largest difference to max_diff list.
    
    largest = heapq.nlargest(n, max_diff) #list of final absolute differences, from largest to smallest
    names = [] #Make empty list for names to fill in to.
    for item in largest: #Assign relevant name to each difference.
        foo = max_diff.index(item)
        name = nuc_name[foo]
        names.append(name)

    #Either print or return largest and names.
    if print_values == True:
        for counter in np.arange(0, n):
            print("%s     %s     %s" % (counter + 1, names[counter], largest[counter]))
    else:
        return largest, names


def find_point_difference(datafile1, datafile2, t):
    
    '''A function to find differences between mass fraction in two datafiles at a specific timestep.
        
        Inputs: datafile1 = ts file to be compared
        datafile2 = second ts file to be compared
        t = timestep desired
        
        Output: largest = list of n largest differences
        names = list of nuclide names corresponding to items in largest
    '''
    
    import numpy as np
    import read_ts_file as rtf
    import heapq
    
    #Read datafile, rename variables.
    zz, aa, xmf, time, temperature, density, timestep, edot, flx_end, flx = rtf.read_ts_file(datafile1)
    zz1, aa1, xmf1, time1, temperature1, density1, timestep1, edot1, flx_end1, flx1 = zz, aa, xmf, time, temperature, density, timestep, edot, flx_end, flx
    
    #Set certain parameters based on the data.
    num_species_total = np.shape(xmf1)[1]
    num_timesteps = np.shape(xmf1)[0]
    
    #Read second ts file, rename and use variables.
    zz, aa, xmf, time, temperature, density, timestep, edot, flx_end, flx = rtf.read_ts_file(datafile2)
    element = rtf.build_element_symbol()
    nuc_name = rtf.build_isotope_symbol(zz,aa)
    zz2, aa2, xmf2, time2, temperature2, density2, timestep2, edot2, flx_end2, flx2 = zz, aa, xmf, time, temperature, density, timestep, edot, flx_end, flx

    diff_list = [] #Make empty list for differences.
    for counter in np.arange(0, num_species_total): #Count through number of species.
        diff = float(xmf1[t][counter]) - float(xmf2[t][counter]) #Find each difference at specified timestep, add to list.
        diff = abs(diff)
        diff_list.append(diff)

    largest = heapq.nlargest(n, diff_list) #Rearrange diff_list into largest to smallest order.
    names = [] #Make empty list for names to fill in to.
    for item in largest: #Assign relevant name to each difference.
        foo = diff_list.index(item)
        name = nuc_name[foo]
        names.append(name)
    return largest, names #Return lists of differences and names.

    #Print results.
    for counter in np.arange(0, n):
        print("%s     %s     %s" % (counter + 1, names[counter], largest[counter]))

def per_diff(datafile1, datafile2):

    '''A function to find final percent difference between two data sets of mass fraction.
    Inputs: datafile1 = ts file to be compared
            datafile2 = second ts file to be compared
    Output: smallest = list of differences from smallest to largest
            names = list of nuclide names corresponding to items in largest
    '''

    import numpy as np
    import read_ts_file as rtf
    import heapq

    #Take in data from the first file, give it variable names.
    zz, aa, xmf, time, temperature, density, timestep, edot, flx_end, flx = rtf.read_ts_file(datafile1)

    zz1, aa1, xmf1, time1, temperature1, density1, timestep1, edot1, flx_end1, flx1 = zz, aa, xmf, time, temperature, density, timestep, edot, flx_end, flx

    #Set certain constraints based on the input data. (How many species, how many timesteps, etc)
    num_species_total = np.shape(xmf1)[1]
    num_timesteps = np.shape(xmf1)[0]
    n = num_species_total
    
    #Take in data from the second file.
    zz, aa, xmf, time, temperature, density, timestep, edot, flx_end, flx = rtf.read_ts_file(datafile2)
    
    #Make lists of elements and nuclear names for later use.
    element = rtf.build_element_symbol()
    nuc_name = rtf.build_isotope_symbol(zz,aa)

    #Change variable names.
    zz2, aa2, xmf2, time2, temperature2, density2, timestep2, edot2, flx_end2, flx2 = zz, aa, xmf, time, temperature, density, timestep, edot, flx_end, flx

    #Set certain constraints based on the new set of input data.
    num_species_total2 = np.shape(xmf2)[1]
    num_timesteps2 = np.shape(xmf2)[0]

    #Create an empty list for the maximum percentage differences.
    max_per_diff = []
    
    for counter in np.arange(0, num_species_total): #count through each column (species)
        avg = (float((xmf1[-1][counter]) + float(xmf2[-1][counter]))/2) #take the average of each item, -1 as an index accesses the last entry, so it's end time
        if avg == 0.0: #if average is zero, make it something calculable
            avg = 10E-20
        diff_list = [] #make empty list for differences
        diff = float(xmf1[-1][counter]) - float(xmf2[-1][counter]) #calculate differences and add to list
        diff = abs(diff)
        diff_list.append(diff)
        per_diff_list = [] #create list of percent differences
        for item in diff_list: #calculate percent differences, add to list
            per_diff = item / avg
            per_diff_list.append(per_diff)
            max_per_diff.append(max(per_diff_list)) #find largest percentage differences and add to list

    smallest = heapq.nsmallest(n, max_per_diff) #make list of smallest percentage differences (remember that n = num_species_total, so this is a list of all perdiffs in ascending order)
    names = [] #make empty list of names
    for item in smallest: #assign relevant name to each percentage difference
        foo = max_per_diff.index(item)
        name = nuc_name[foo]
        names.append(name)

    for counter in np.arange(0, n):
        print("%s     %s     %s" % (counter + 1, names[counter], smallest[counter]))

def compare_final(datafile1, datafile2):

    '''A function to compare one datafile relative to another
    Inputs: datafile1 = ts file to be compared (reported errors are relative to this file)
            datafile2 = second ts file to be compared
    Output: rerr_sort = list of relative errors from smallest to largest
            names = list of nuclide names corresponding to items in rerr_sort
    '''

    import numpy as np
    import read_ts_file as rtf
    import heapq

    #Take in data from the first file, give it variable names.
    zz1, aa1, xmf1, time1, temperature1, density1, timestep1, edot1, flx_end1, flx1 = rtf.read_ts_file(datafile1)
    en1 = np.multiply(edot1,timestep1)
    enuc1 = np.cumsum(en1)

    #Set certain constraints based on the input data. (How many species, how many timesteps, etc)
    num_species_total = np.shape(xmf1)[1]
    num_timesteps = np.shape(xmf1)[0]
    
    #Take in data from the second file.
    zz2, aa2, xmf2, time2, temperature2, density2, timestep2, edot2, flx_end2, flx2 = rtf.read_ts_file(datafile2)
    en2 = np.multiply(edot2,timestep2)
    enuc2 = np.cumsum(en2)
    
    #Make lists of elements and nuclear names for later use.
    element = rtf.build_element_symbol()
    nuc_name = rtf.build_isotope_symbol(zz2,aa2)

    #Set certain constraints based on the new set of input data.
    num_species_total2 = np.shape(xmf2)[1]

    #Create lists for the differences
    aerr = np.abs(np.subtract(xmf2[-1][:],xmf1[-1][:]))
    rerr = np.divide(aerr,np.abs(xmf1[-1][:]))

    isort = np.argsort(rerr)
    aerr_sort = aerr[isort]
    rerr_sort = rerr[isort]
    xmf1_sort = xmf1[-1][isort]
    xmf2_sort = xmf2[-1][isort]
    name_sort = [nuc_name[i] for i in isort]

    print("%5s %5s\t%10s\t%16s\t%16s\t%16s\t%16s" % ('i','isort','name','X1','X2','|dX|','|dX| / |X1|'))

    fmt = "%5s %5s\t%10s\t%16.8e\t%16.8e\t%16.8e\t%16.8e"
    for i in np.arange(0, num_species_total):
        print(fmt % (i+1, isort[i]+1, name_sort[i], xmf1_sort[i], xmf2_sort[i], aerr_sort[i], rerr_sort[i]))

    print("")

    fmt = "%5s %5s\t%10s\t%16s\t%16s\t%16.8e\t%16.8e"
    aerr_norm = np.linalg.norm(aerr,ord=2)
    rerr_norm = np.divide(aerr_norm,np.linalg.norm(xmf1[-1][:],ord=2))
    print(fmt % ('', '', '2-norm', '', '', aerr_norm, rerr_norm))

    fmt = "%5s %5s\t%10s\t%16.8e\t%16.8e\t%16.8e\t%16.8e"
    aerr = np.abs(np.subtract(temperature2[-1], temperature1[-1]))
    rerr = np.divide(aerr, np.abs(temperature1[-1]))
    print(fmt % ('', '', 'T', temperature1[-1], temperature2[-1], aerr, rerr))

    aerr = np.abs(np.subtract(enuc2[-1], enuc1[-1]))
    rerr = np.divide(aerr, np.abs(enuc1[-1]))
    print(fmt % ('', '', 'E_nuc', enuc1[-1], enuc2[-1], aerr, rerr))
