def find_file_type(datafile):
    
    '''A function to determine the type of file used. Assumes that file type is described by two letters in file name.
        
        Inputs: datafile = ev or ts file
        Output: file_type = either ev or ts
    '''

    import numpy as np
    
    for counter in np.arange(0, len(datafile)): #count through letters in datafile
        if datafile[counter] == 'e' and datafile[counter + 1] == 'v': #if 'e' and 'v' appear consecutively, it's an ev file
            file_type = 'ev'
            break
        elif datafile[counter] == 't' and datafile[counter + 1] == 's': #if 't' and 's' appear consecutively, it's a ts file
            file_type = 'ts'
            break
        else: #if neither of these conditions are met, print a warning message and assume it's a ts file
            print("In the future, please include 'ev' or 'ts' somewhere in the file name. We're just going to assume that you've given us a ts file now.")
            file_type = 'ts'

    return file_type
