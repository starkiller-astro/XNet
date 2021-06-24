"""
    Module to read ev files
    """
def read_ev_file(datafile):
    """
        Inputs > datafile, name of input file (ev file)
        Outputs< time, elapsed time as a function of time
        temperature, temperature (in GK) as a function of time
        density, mass density (in g/cm3) as a function of time
        timestep, timestep as a function of time
        edot, energy generation rate as a function of time
        species, a list of each species
        data, a shortcut to np.genfromtxt that reads the given datafile
        datafile, the same ev file from input
        """
    import numpy as np
    
    #Read the data file
    col_headers = np.genfromtxt(datafile, max_rows = 1, dtype = str)
    
    if len(col_headers) != 22:
        col_headers[20] = 'I'
        col_headers = np.append(col_headers, 't')
        col_headers_list = col_headers.tolist()
        data = np.genfromtxt(datafile, skip_header = 1, dtype = float, names = col_headers_list)
    
    else:
        data = np.genfromtxt(datafile, names=True)
    
    #Sort columns of data by name. Store in a list.
    columns = []
    for name in data.dtype.names:
        columns.append(name)
    
    #Save columns as variables.
    time = data['Time']
    density = data['Density']
    temperature = data['TGK']
    edot = data['dEdt']
    timestep = data['Timestep']

    #Make a list of species mass fraction data (xmf).
    species = []
    alphabet = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't', 'v', 'w', 'x', 'y', 'z']
    
    for item in columns:
        if item != 'Time' and item != 'TGK' and item != 'Density' and item != 'dEdt' and item != 'Timestep':
            species.append(item)

    for item in species:
        if item in alphabet:
            species.remove(item)
    
    return time, density, temperature, edot, timestep, species, data, datafile

