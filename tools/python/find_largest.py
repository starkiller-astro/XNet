def absolute_largest(datafile1, datafile2, print_values = True):
    '''
        Inputs: datafile1 = ts file
                datafile2 = second ts file
                print_values = default to True, if True will print differences, if not will return lists of differences and corresponding nuclide names
                
        Outputs: names= list of nuclide names corresponding to largest species in datafile1
                 largest = final abundances of species in datafile1 in descending order
                 names2 = list of nuclide names corresponding to largest species in datafile2
                 largest2 = final abundances of species in datafile2 in descending order
        '''
    
    import numpy as np
    import read_ts_file as rtf
    import heapq

    #Read datafile1, use variable names.
    zz, aa, xmf, time, temperature, density, timestep, edot, flx_end, flx = rtf.read_ts_file(datafile1)
    nuc_name = rtf.build_isotope_symbol(zz, aa)

    #Set certain parameters based on the shape of the data.
    num_species_total = np.shape(xmf)[1]
    n = num_species_total

    xmf_list = [] #Make empty list of abundances.
    for counter in np.arange(0, num_species_total): #Count through all species.
        foo = xmf[-1][counter] #Find final abundances (shown by -1)
        xmf_list.append(foo) #Add final abundance to list.

    largest = heapq.nlargest(n, xmf_list) #Rearrange abundance list into descending order.
    names = [] #Create empty list for names to fill into.
    for item in largest: #Add relevant names to correct position in list.
        bar = xmf_list.index(item)
        name = nuc_name[bar]
        names.append(name)

    #Read second data file, use variable names.
    zz, aa, xmf, time, temperature, density, timestep, edot, flx_end, flx = rtf.read_ts_file(datafile2)
    nuc_name = rtf.build_isotope_symbol(zz, aa)

    #Set parameter according to data.
    num_species_total = np.shape(xmf)[1]

    xmf_list2 = [] #Make second empty list of abundances.
    for counter in np.arange(0, num_species_total): #Count through all species, find final abundances (shown by -1), add to list.
        foo = xmf[-1][counter]
        xmf_list2.append(foo)

    largest2 = heapq.nlargest(n, xmf_list2) #Rearrange xmf_list2 into descending order
    names2 = [] #Create empty list for names to fill in to.
    for item in largest2: #Add relevant names to correct position in list.
        bar = xmf_list2.index(item)
        name = nuc_name[bar]
        names2.append(name)

    #Either print or return results.
    if print_values == True:
        for counter in np.arange(0, n):
            print("%s      %s      %s" % (counter + 1, names[counter], largest[counter], names2[counter], largest2[counter]))

    else:
        return names, largest, names2, largest2
