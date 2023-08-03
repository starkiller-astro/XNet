'''
This script renames the files in the directory to the format "train_co_alpha_{ab_label}_{density_label}_{temp_label}.txt"
where ab_label is the alpha-beta ratio, density_label is the density, and temp_label is the temperature.
Note: it assumes that the files are named "train_co_alpha_{file_num}.txt" where file_num is a number from 00001 to 11440.
'''
import os

def extract_numeric_part(filename):
    numeric_part = filename.split('_')[-1]  # Get the last part after the last underscore
    numeric_part = numeric_part.replace('.txt', '')  # Remove the '.txt' extension
    return int(numeric_part)

def rename_files():
    base_path = "/Users/alancangas/XNet/test/trng_ts_4_AE"  # Replace with the actual path to your files

    # Loop through each file number (1 to 11440)
    for file_num in range(1, 11441):
        # Get the numeric part of the original filename
        numeric_part = extract_numeric_part(f"train_co_alpha_{str(file_num).zfill(5)}")

        # Calculate ab, density, and temp from the numeric part
        ab = (numeric_part - 1) // (26 * 40) + 1
        density = ((numeric_part - 1) // 40) % 26 + 1
        temp = (numeric_part - 1) % 40 + 1

        # Format ab, density, and temp labels
        ab_label = f"ab{str(ab).zfill(2)}"
        density_label = f"dens{str(density).zfill(2)}"
        temp_label = f"temp{str(temp).zfill(2)}"

        # Construct the original and new file names
        old_file_name = f"train_co_alpha_{str(file_num).zfill(5)}"
        new_file_name = f"train_co_alpha_{ab_label}_{density_label}_{temp_label}"

        # Construct the full original and new file paths
        old_file_path = os.path.join(base_path, f"{old_file_name}")
        new_file_path = os.path.join(base_path, f"{new_file_name}")

        # Rename the file
        os.rename(old_file_path, new_file_path)

rename_files()
