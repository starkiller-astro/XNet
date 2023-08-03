'''
This code is used to open the control file and input the training data files into the control file
'''
def append_text_to_file(file_path, text_to_append):
    with open(file_path, 'a') as file:
        file.write(text_to_append + '\n')

for counter in range (1,12):
    index = f"{counter:02d}"
    for i in range (1,27):
        density = 4.0 + 0.2 * (i - 1)
        #print(f"Run = {i} temperature = {temperature:.1f}")
        
        


        for j in range (1,41):
                file_path = '/Users/alancangas/XNet/test/control' #this is for opening the control file to write to
                temperature = 0.1 * j
                text_to_add = (f"Data_alpha/ab_co_trn_{index}")
                append_text_to_file(file_path, text_to_add)
                text_to_add = (f"trng_data_4_AE/training_data_dens{density:.1f}_temp{temperature:.1f}.txt")
                append_text_to_file(file_path, text_to_add)


