'''
This script is used to generate training data for the autoencoder. It takes a file with a single line of data and 
replaces the temperature and density values with the desired values. It then saves the new file with the desired name.
'''
import os

original_file_path = '/Users/alancangas/XNet/test/Test_Problems/target.txt'
output_directory = '/Users/alancangas/XNet/test/trng_data_4_AE/'
line_numbers = [5, 6]

increment_index1 = 0.1
start_exp = 4.0
end_exp = 9.0

for increment_index2_for_run in range(1, 27):
    exponent = start_exp + (increment_index2_for_run - 1) * 0.2
    value2_base = 10 ** exponent
    print(f"Run {increment_index2_for_run}: exponent = {exponent:.1f}, value2_base = {value2_base:.1f}")


    for increment_index1_for_run in range(1, 42):
        value1_base = (increment_index1_for_run - 1) * increment_index1
        
        for i in range(40):
            temperature_value = value1_base + (i * increment_index1)

            if temperature_value > 4.0 or temperature_value == 0.0:
                break

            value2 = value2_base * (10 ** (i * increment_index1))

            updated_file_path = f"{output_directory}training_data_dens{exponent:.1f}_temp{temperature_value:.1f}.txt"
            with open(original_file_path, 'r') as original_file:
                lines = original_file.readlines()

            updated_lines = lines.copy()

            for line_number, line in enumerate(updated_lines, 1):
                if line_number in line_numbers:
                    values = line.split()
                    if len(values) >= 4:
                        try:
                            values[1] = "{:.6E}".format(temperature_value)
                            values[2] = "{:.6E}".format(value2)
                            updated_lines[line_number - 1] = ' '.join(values) + "\n"
                        except ValueError:
                            print(f"Line {line_number} does not contain valid numbers at index 1 and/or index 2.")

            with open(updated_file_path, 'w') as updated_file:
                updated_file.writelines(updated_lines)

            print(f"Updated data saved to {updated_file_path}.")
            
