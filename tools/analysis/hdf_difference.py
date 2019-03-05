def hdf_difference(file1, file2, column = False, column_number = 0, number_species = 0):

    '''
        Inputs: file1 = an hdf file of abundances
                file2 = a second hdf file of abundances
                column = option to analyze a specific column, defaults to False
                column_number = requires column to be True, specifies which column to analyze, defaults to 0
                number_species = number of species that constitute percentage (x out of total), defaults to 0
                
        Outputs: hdf file of the abundance differences between the two given files
        
        WHEN RETURNING PARTICLE IDS WITH THE LARGEST DIFFERENCES, YOU MUST ADD ONE TO THE ID TO GET THE CORRECT ID, SINCE THERE IS AN OFFSET.
        '''
    
    import h5py
    import numpy as np
    import heapq

    #Create hdf file, structure and name it.
    r = h5py.File('ts_abundances_difference.h5', 'w')
    group = r.create_group("differences")
    array_name = '%s and %s' % (file1, file2)

    #Use the given hdf files.
    with h5py.File(file1, 'r') as f1:
        with h5py.File(file2, 'r') as f2:
            array1 = f1['abundance_arrays/local_abundances'][()]
            array2 = f2['abundance_arrays/local_abundances'][()]

            #Subtract the two arrays.
            result = np.subtract(array1, array2)
            
            #Restore the first column to a numbered label of each file.
            for counter in np.arange(0, 4720):
                result[counter][0] = counter + 1

    #Create the dataset and fill it in with the subtracted array. maxshape = (None, None) makes the size mutable.
    dset = group.create_dataset(array_name, data = result, maxshape = (None, None))

    #Create list for largest particle names.
    largest_files = list()
    
    #Format desired column so that it can be accessed.
    column_chosen = result[:, column_number]

    #If selecting a column, organize that column from largest to smallest.
    if column == True:
        for counter in np.arange(0, len(column_chosen)):
            largest = heapq.nlargest(len(column_chosen), column_chosen)

        #Add each largest difference's corresponding name to the largest files list (same as position in list, initially)
        for item in largest:
            position = np.where(column_chosen == item)
            largest_files.append(position)

        #Format and print largest file names (will be ugly, but legible enough I suppose) YOU MUST ADD ONE TO EACH FILE NAME IN ORDER TO GET THE CORRECT FILE, SINCE PYTHON STARTS COUNTING AT 0
        for item in largest_files:
            item = item[0]
        print largest_files
        
        #Take the specified number of species and find what fraction of the total difference it constitutes.
        if number_species != 0:
            sum = 0
            partial_sum = 0
            for counter in np.arange(0, 4720): #Calculate total sum.
                sum += abs(largest[counter])
            for counter in np.arange(0, number_species - 1): #Calculate partial sum.
                partial_sum += abs(largest[counter])
    #Print results.
            print "Total sum: %s" % sum
            print "Partial sum: %s" % partial_sum
            print "Percentage of total difference in species: %s" % (partial_sum/sum)

    #Close the file.
    r.close()
