def pinput_radii_reader(pinput_file):
    
    import numpy as np
    import imp

    #Read file of radii for each row of particles
    radius_pathname = pinput_file
    fradius = open(radius_pathname,'r')
    rad_list = np.loadtxt(fradius,skiprows=1,usecols=(0,))
    
    with open(radius_pathname, 'r') as f_pnum:
        line_1 = f_pnum.readline()
        print line_1
    line_list = line_1.split()
    n_part_x = int(line_list[1])
    n_part_y = int(line_list[3])
    n_part_z = int(line_list[5])
    
    return rad_list,n_part_x,n_part_y,n_part_z