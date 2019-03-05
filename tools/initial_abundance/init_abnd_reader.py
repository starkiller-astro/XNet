def init_abnd_reader(init_file):
    
    import numpy as np
    import re as re
    import imp
    #from check_mass_number import check_mass_number
    #from isfloat import isfloat
    
    cmn = imp.load_source('check_mass_number.py', 'C:\\Users\\taylor\\Canopy\\scripts\\check_mass_number.py')
    isf = imp.load_source('isfloat.py', 'C:\\Users\\taylor\\Canopy\\scripts\\isfloat.py')
    
    with open(init_file, 'r') as f:
        line_1 = f.readline()
        print line_1
    line_list = line_1.split()
    version = int(line_list[1])
    print version
    
    #Reads WH07 file format.
    if version == 10104:
        
        #Open the large file of abundances
        f_pathname = init_file
        f_opener = open(f_pathname,'r')
        #Delimiter is whitespace by default.

        #Use genfromtxt to read from file (skipping header and wind row). Missing_values
        #is an option this way, but it does not work.
        abnd_list = np.genfromtxt(f_opener,dtype=None,skip_header=1,skip_footer=1)
        f_opener.close()

        
    #Reads different WH15 file format.
    elif version == 10201:
        
        f_pathname = init_file
        f_opener = open(f_pathname,'r')
        #Change Delimiter?
        abnd_list = np.genfromtxt(f_opener,dtype=None,skip_header=1)
        f_opener.close()   
                
    isotopes_for_mass = abnd_list[0,13:len(abnd_list[0,:])]
    molar_masses = np.zeros((len(isotopes_for_mass),))
    molar_masses[0] = 1
    for ii in range(1,len(isotopes_for_mass)):
        inputString = isotopes_for_mass[ii]
        mass_bool = cmn.check_mass_number(inputString)
        if mass_bool == True:
            mass_val = int(re.search(r'\d+', inputString).group())
            print mass_val
            molar_masses[ii] = mass_val
            
    if version == 10104:
               
        #Replace necessary values
        for i in range(0,len(abnd_list[:,0])):
            for j in range(0,len(abnd_list[0,:])):
                if abnd_list[i,j] == '---':
                    abnd_list[i,j] = 1.0e-50
                if abnd_list[i,j] == '0.0000000000000000E+00':
                    abnd_list[i,j] = 1.0e-50
                if abnd_list[i,j] == 'nt1':
                    abnd_list[i,j] = 'n'
                if abnd_list[i,j] == 'h1':
                    abnd_list[i,j] = 'p'
                if abnd_list[i,j] == 'h2':
                    abnd_list[i,j] = 'd'
                if abnd_list[i,j] == 'h3':
                    abnd_list[i,j] = 't'   

    elif version == 10201:     
                
        #Replace necessary values
        for i in range(0,len(abnd_list[:,0])):
            for j in range(0,len(abnd_list[0,:])):
                if abnd_list[i,j] == '0':
                    abnd_list[i,j] = 1.0e-50
                if abnd_list[i,j] == '0.00000000000000000e+00':                
                    abnd_list[i,j] = 1.0e-50
                if abnd_list[i,j] == 'H1':
                    abnd_list[i,j] = 'p'
                if abnd_list[i,j] == 'H2':
                    abnd_list[i,j] = 'd'
                if abnd_list[i,j] == 'H3':
                    abnd_list[i,j] = 't' 
                    
                    
    #Turn into a list for easier manipulation
    test_list = abnd_list.tolist()

    #If the strings are floats (using isfloat fxn), turn them into floats.
    #If they are not, leave them be.
    for i in range(0,len(test_list)): #len(abnd_list[:,0])
        for j in range(0,len(test_list[0])): #len(abnd_list[0,:])
            float_bool = isf.isfloat(test_list[i][j])
            if float_bool == True:
                floated = float(test_list[i][j])
                test_list[i][j] = floated

    if version == 10104:
        ### Remove he4 and put behind he3
        he4_column = []
        for row in test_list:
            he4_column.append(row[15])
            del row[15]

        i=0
        for row in test_list:
            row.insert(18, he4_column[i])
            i=i+1

        #Get wind values separately (different # of columns, so genfromtxt won't read it)
        with open(f_pathname, 'rb') as fh:

            wind = None
            for line in (line for line in fh if line.rstrip('\n')):
                wind = line

        wind = wind.split()
        for i in range(0,len(wind)):
            if (wind[i].isalpha() == False):
                if (":" not in wind[i]):
                    wind[i] = float(wind[i])

        print len(wind)
        
    return version,test_list,molar_masses