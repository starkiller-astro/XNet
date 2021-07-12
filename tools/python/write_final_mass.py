def write_final_mass(final_mf_file, final_mass_file, nuc_names_wanted = "None"):
        """
            Inputs: final_mf_file = text file with final mass fractions of species with mass shell data
                    final_mass_file = name of file to be written
                    nuc_names_wanted = desired nuclear species to write final masses for
                    
            Outputs: A text file with the final mass abundances of ever species in the model
            
            A Function that takes a file with final mass fractions of species and mass shell data, and 
            outputs a file with the final masses of all species, to be used with match_abundance_to_mass_shell
        """
        
        
        import pandas as pd
        
        #Read in Data from text file
        data = pd.read_csv(final_mf_file, sep = '   ', skipfooter = 1, engine = 'python')
        
        #Create a new file to write final species masses in
        file = open(final_mass_file, 'w')
        
        #If all final species masses are wanted, loop through every column of data:
        if nuc_names_wanted == "None":
            for column in data.columns[5:]: 
                #Multiply mass fractions by shell mass, then sum the shells together
                data["tempcol"] = data[column]*data["dm(g)"]
                mass_sum = data["tempcol"].sum()
                
                #Write final mass to file, along with species name
                file.writelines(column + ":    " + str(mass_sum) + '\n')
        
        #Loop through only desired nuclear species
        else: 
            for item in nuc_names_wanted:
                #Adjust formatting for species names in datafile
                if nuc_names_wanted not in ordered_data.columns:
                    nuc_names_wanted = " " + nuc_names_wanted
                
                #Multiply mass fractions by shell mass, then sum the shells together
                data["tempcol"] = data[item]*data["dm(g)"]
                mass_sum = data["tempcol"].sum()
                
                #Write final mass to file, along with species name
                file.writelines(item + ":   " + str(mass_sum) + '\n')
        
        #Close File
        file.close()
        