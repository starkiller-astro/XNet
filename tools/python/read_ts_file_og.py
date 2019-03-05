"""
Module to read XNet ts files 
"""
def read_ts_file(file_name, debug=False):
    """
    Inputs > file_name, name of input file 
             debug, boolean that turns debug mode on or off, defaults to False
    Outputs< zz, proton number of every species in network 
             aa, atomic mass number of every species in the network
             xmf, mass fraction as a function of time, for every species
             time, elapsed time as a function of time
             temperature, temperature (in GK) as a function of time
             density, mass density (in g/cm3) as a function of time
             timestep, timestep as a function of time
             edot, energy generation rate as a function of time
             flx_end, initial and final species for each reaction flux in the network
             flx, reaction fluxes as a function of time
    """
    print "If plotting fails, try uncommenting all lines marked MODIFICATION FOR AUSTIN'S FILES."

    import numpy as np
    
    # Open file_name

    with open(file_name, 'rb') as file_id:
    
    # Read Run Descriptions
    
        read_desc=np.empty(80,dtype='str')
    
        record_length1 =np.fromfile(file_id,dtype='uint32',count=1)
        read_desc      =np.fromfile(file_id,dtype='uint8',count=80).tostring()
        desc = read_desc
        read_desc      =np.fromfile(file_id,dtype='uint8',count=80).tostring()
        desc=desc+read_desc
        read_desc      =np.fromfile(file_id,dtype='uint8',count=80).tostring()
        desc=desc+read_desc
        read_desc      =np.fromfile(file_id,dtype='uint8',count=80).tostring()
        data_desc= read_desc
        print desc, data_desc
        record_length2 =np.fromfile(file_id,dtype='uint32',count=1)
    
    #Read Run Settings
    
        record_length1 =np.fromfile(file_id,dtype='uint32',count=1)
        kstmx          =np.fromfile(file_id,dtype='int32',count=1)[0]
        kitmx          =np.fromfile(file_id,dtype='int32',count=1)
        iweak          =np.fromfile(file_id,dtype='int32',count=1)
        iscrn          =np.fromfile(file_id,dtype='int32',count=1)
        iconvc         =np.fromfile(file_id,dtype='int32',count=1)
        changemx       =np.fromfile(file_id,dtype='float64',count=1)
        tolm           =np.fromfile(file_id,dtype='float64',count=1)
        tolc           =np.fromfile(file_id,dtype='float64',count=1)
        yacc           =np.fromfile(file_id,dtype='float64',count=1)
        ymin           =np.fromfile(file_id,dtype='float64',count=1)
        tdel_mm        =np.fromfile(file_id,dtype='float64',count=1)
        record_length2 =np.fromfile(file_id,dtype='uint32',count=1)
    
    # Read Abundance Info
    
        record_length1 =np.fromfile(file_id,dtype='uint32',count=1)
        read_desc      =np.fromfile(file_id,dtype='uint8',count=80).tostring()
        abund_file     =read_desc
        read_desc      =np.fromfile(file_id,dtype='uint8',count=80).tostring()
        abund_desc     =read_desc
        record_length2 =np.fromfile(file_id,dtype='uint32',count=1)
    
    # Read Thermodynamic Info
    
        record_length1 =np.fromfile(file_id,dtype='uint32',count=1)
        read_desc      =np.fromfile(file_id,dtype='uint8',count=80).tostring()
        thermo_file    =read_desc
        read_desc      =np.fromfile(file_id,dtype='uint8',count=80).tostring()
        thermo_desc    =read_desc
        record_length2 =np.fromfile(file_id,dtype='uint32',count=1)
    
    # Read Nuclear Info
    
        record_length1 =np.fromfile(file_id,dtype='uint32',count=1)
        ny             =np.fromfile(file_id,dtype='int32',count=1)[0]
        zz             =np.fromfile(file_id,dtype='float64',count=ny)
        aa             =np.fromfile(file_id,dtype='float64',count=ny)
        record_length2 =np.fromfile(file_id,dtype='uint32',count=1)
    
    # Read Flux Info
    
        record_length1   =np.fromfile(file_id,dtype='uint32',count=1)
        nflx             =np.fromfile(file_id,dtype='int32',count=1)[0]
        if nflx>0:
            flx_end      =np.zeros((nflx,2),dtype='int32')
            flx_end[:,0] =np.fromfile(file_id,dtype='int32',count=nflx)
            flx_end[:,1] =np.fromfile(file_id,dtype='int32',count=nflx)
        record_length2   =np.fromfile(file_id,dtype='uint32',count=1)
    
    # Size data arrays
    
        time        = np.zeros(kstmx,dtype='float64')
        temperature = np.zeros(kstmx,dtype='float64')
        density     = np.zeros(kstmx,dtype='float64')
        timestep    = np.zeros(kstmx,dtype='float64')
        edot        = np.zeros(kstmx,dtype='float64')
        xmf         = np.zeros((kstmx,ny),dtype='float64')
        if nflx>0 :
            flx     = np.zeros((kstmx,nflx))
    
    # Read data from each timestep
        
        for k in range(kstmx-1):
            record_length1 =np.fromfile(file_id,dtype='uint32',count=1)
    
    # If end of file, exit
            if record_length1.size == 0:
                break
    
    # Otherwise read data
    
            kstep          =np.fromfile(file_id,dtype='int32',count=1)
            time[k]        =np.fromfile(file_id,dtype='float64',count=1)
            temperature[k] =np.fromfile(file_id,dtype='float64',count=1)
            density[k]     =np.fromfile(file_id,dtype='float64',count=1)
            timestep[k]    =np.fromfile(file_id,dtype='float64',count=1)
            edot[k]        =np.fromfile(file_id,dtype='float64',count=1)
            xmf[k,:]       =np.fromfile(file_id,dtype='float64',count=ny)
                #if nflx>0 : #MODIFICATION FOR AUSTIN'S FILES
                #flx[k,:]   =np.fromfile(file_id,dtype='float64',count=nflx) #MODIFICATION FOR AUSTIN'S FILES
            record_length2 =np.fromfile(file_id,dtype='uint32',count=1)
            if debug:
                print k, kstep, time[k]
    
    # Calculate # of timesteps
    
    nstep = k+1
    
    # Remove empty timesteps
    
    time        = time[0:nstep-1]
    temperature = temperature[0:nstep-1]
    density     =density[0:nstep-1]
    timestep    =timestep[0:nstep-1]
    edot        =edot[0:nstep-1]
    xmf         =xmf[0:nstep-1,:]
    if nflx>0 :
        flx     =flx[0:nstep-1,:]

    # Convert abundance to mass fraction
    
#xmf = xmf * aa #MODIFICATION FOR AUSTIN'S FILES

    return zz, aa, xmf, time, temperature, density, timestep, edot, flx_end, flx
    
def build_element_symbol():
    """
    Builds symbols for Chemical Elements
    Input < None
    Output> element
    """
    element = ('n ','H ','He','Li','Be','B ','C ','N ','O ','F ','Ne',
               'Na','Mg','Al','Si','P ','S ','Cl','Ar','K ','Ca','Sc','Ti',
               'V ','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se',
               'Br','Kr','Rb','Sr',' Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd',
               'Ag','Cd','In','Sn','Sb','Te',' I','Xe','Cs','Ba','La','Ce',
               'Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb',
               'Lu','Hf','Ta','W ','Re','Os','Ir','Pt','Au','Hg','Tl','Pb',
               'Bi','Po','At','Rn','Fr','Ra','Ac','Th','Pa','U ','Np','Pu',
               'Am','Cm','Bk','Cf','Es','Fm','Md','No','Lr','Rf','Db','Sg',
               'Bh','Hs','Mt','Ds','Rg','Cn')

    return element
    
def build_isotope_symbol(zz,aa):
    """
    Builds isotope symbols, with superscripts
    Input < zz, proton number for each isotope in the network
            aa, atomic mass number of each isotope in the network
    Ouputs> nuc_name, isotopic symbol
    """
    import numpy as np
    element = build_element_symbol()
    zint=zz.astype(np.intp)
    aint=aa.astype(np.intp)
    n_iso=np.size(zint)
    nuc_name=[None] * n_iso
    for k in range(n_iso):
        nuc_name[k] = "$^{%u}$" % aint[k]
        nuc_name[k] = nuc_name[k] + element[zint[k]]

    return nuc_name
