def calc_abnd(pinput_file,init_file):

    import numpy as np
    from scipy import interpolate
    import re as re
    import matplotlib.pyplot as plt
    import warnings
    import imp
    #from log_interp1d import log_interp1d
    #from pinput_radii_reader import pinput_radii_reader
    #from init_abnd_reader import init_abnd_reader

    prr = imp.load_source('pinput_radii_reader.py', 'C:\\Users\\taylor\\Canopy\\scripts\\pinput_radii_reader.py')
    iar = imp.load_source('init_abnd_reader.py', 'C:\\Users\\taylor\\Canopy\\scripts\\init_abnd_reader.py')

    rad_list,n_part_x,n_part_y,n_part_z = prr.pinput_radii_reader(pinput_file)
    version,test_list,molar_masses = iar.init_abnd_reader(init_file)
    
            
    #ROTATE TEST_LIST TO GET ACCESSIBLE ROW PER ISOTOPE.
    rotated = map(list, zip(*test_list))
    
    #Take rows for all abundances. Store 0th element in list and delete from main
    abundance_rows=rotated[13:] 
    isotope_list = []
    for row in abundance_rows:
        isotope_list.append(row[0])
        del row[0]
    
    #Different WH07 format.
    if version == 10104:
        
        #Customize (no clear way to read from file)
        model_mass = 12
        model_letter = 'B'
        model_name = model_letter+str(model_mass)+'-WH07'
        
        #Take radius row and delete 0th element 'cell_outer_radius'
        radius_vals=rotated[2]
        del radius_vals[0]
        
        mass_coords = rotated[1]
        del mass_coords[0]
        
        cell_Y_e = rotated[11]
        del cell_Y_e[0]
    
    elif version == 10201:
        
        model_mass = 9.6
        model_letter = 'Z'
        model_name = model_letter+str(model_mass)+'-WH15'
            
        radius_vals=rotated[3]
        del radius_vals[0]       

        mass_coords = rotated[2]
        del mass_coords[0]

        cell_Y_e = rotated[12]
        del cell_Y_e[0]

        

    #Same for both:
    mass_list = np.arange(0,1)
    for ii in range(0,len(rad_list)):
        for jj in range(1,len(radius_vals)):
            if rad_list[ii] > radius_vals[jj-1] and rad_list[ii] < radius_vals[jj]:
                mass_approx = mass_coords[jj-1]
                mass_list = np.append(mass_list,mass_approx)
                    
    mass_array = mass_list[1:101]


    for k in range(0,len(rad_list)):

        row_radius = rad_list[k]
    
        '''        
        log_interp_func = log_interp1d(radius_vals,abundance_rows,kind='linear')
        abundances = log_interp_func(row_radius)
        '''
        interp_func = interpolate.interp1d(radius_vals,abundance_rows,kind='linear',bounds_error=False,fill_value="extrapolate")
        abundances = interp_func(row_radius)
    
        isotope_list = np.asarray(isotope_list)
        abundances = np.asarray(abundances)
    
        indices_list = np.where(abundances > 1.0e-30)[0]
    
        for q in range(0,n_part_y):
        
            particle_name = (k*n_part_y)+q
        
            filestring = '%s%05.0f%s' % ("part",particle_name,"abundances.txt")
            f= open(filestring,"w+")
        
            ### Write a comment with f.write for the beginning of the file.
            comment_string = "# Model "+model_name+" Row "+str(k+1)+" Particle "+str(particle_name)+", Approximate Mass Coordinate = "+str(mass_array[k])+"\n"
            f.write(comment_string)
            ye_string = "Ye   "+str(cell_Y_e[k])+"\n"
            f.write(ye_string)
        
                    
            X_tot = sum(abundances)
            if X_tot < 0.99:
                warnings.warn('sum(A*Y) < 0.99 at Row '+str(k+1))
        
            for index in indices_list:
            
                Y_values = (abundances[index]/molar_masses[index])/X_tot
            
                f.writelines(["%-5s %24.15e\n" % (isotope_list[index],Y_values)])
            
            f.close()      