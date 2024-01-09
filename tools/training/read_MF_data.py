'''
 This program reads the data from the neon data files and prints the last entry of the temperature and mass fraction
 It also plots the temperature vs mass fraction graph
 Note: The code assumes that the data files are present in a specific directory and follow the naming convention 
 "ev_co_alpha_control_data_{index}_{file_number:02d}{file_extension}". The elmnt_index value must be valid and 
 correspond to the desired element's mass fraction column in the data files.
'''
import os
import numpy as np
import matplotlib.pyplot as plt

data_directory = ''
elmnt_index = int(input('Choose the index for the element you want to analyze\n 9 = Neon\n 12 = Aluminum\n 13 = Argon\n 14 = Calcium\n 17 = iron\n'))
index = input("Choose increment 01=0.05 , 02=0.01 , 03=0.5 :\n")
file_prefix = f"/Users/alancangas/XNet/test/training_AE_summ2023/new_Tests/ev_co_alpha_control_data_{index}_"
file_extension = ''
num_files = int(input("How many files do you want to read? "))

def read_last_entry(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()
        if lines:
            columns = lines[-1].split()
            if len(columns) >= 9:
                return float(columns[2]), float(columns[elmnt_index])
    return None, None

print("Temperature\tNeon Mass Fraction")
print("--------------------------------")
for file_number in range(1, num_files + 1):
    file_name = f"{file_prefix}{file_number:02d}{file_extension}"
    file_path = os.path.join(data_directory, file_name)

    column2, column9 = read_last_entry(file_path)

    if column2 is not None and column9 is not None:
        print(f"{column2}\t\t{column9}")
print("--------------------------------")

print("Temperature vs Mass Fraction graph")

for file_number in range(1, num_files + 1):
    file_name = f"{file_prefix}{file_number:02d}{file_extension}"
    file_path = os.path.join(data_directory, file_name)

    column2, column9 = read_last_entry(file_path)

    if column2 is not None and column9 is not None:
        plt.scatter(column2, column9)
plt.grid(True)
plt.xscale('linear')
plt.yscale('log')
plt.xlabel('Temperature')
plt.ylabel('Mass Fraction')
plt.title('Mass Fraction vs Temperature')
plt.show()

