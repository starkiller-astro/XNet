# Summer 2023 Internship Project

This repository contains a collection of Python scripts and data used during my summer internship project in 2023. The scripts are designed to perform various tasks related to calculating molar fractions, generating training data for an autoencoder, and analyzing temperature vs mass fraction data.

## Requirements

- Python 3.x
- NumPy
- Matplotlib

## Usage Instructions

1. **Calculate Molar Fractions**

   The script `co_ratio_script.py` calculates the molar fraction of C12 and O16 in the C12 + O16 reaction. It generates CO ratios for the autoencoder. Run the script, and it will print the C12 and O16 molar fractions for different values of 'i' ranging from 0 to 10.

2. **Input Training Data**

   The script `place_file_names.py` opens the control file and inputs the training data files into the control file. Modify the `file_path` and `output_directory` variables with the appropriate paths. Then, run the script, and it will append the training data filenames to the control file.

3. **Generate Training Data**

   The script `file_data_gener.py` generates training data for the autoencoder. It takes a file with a single line of data and replaces the temperature and density values with desired values. The output files are saved in the specified output directory. Update the `original_file_path`, `output_directory`, and `line_numbers` variables as needed. Run the script, and it will generate multiple training data files with different temperature and density values.

4. **Rename Files**

   The script `rename_files.py` renames the files in the directory based on the format "train_co_alpha_{ab_label}_{density_label}_{temp_label}.txt". It assumes that the files are named "train_co_alpha_{file_num}.txt" where file_num is a number from 00001 to 11440. Run the script, and it will rename the files in the specified directory.

5. **Analyze Neon Data**

   The script `read_MF_data.py` reads the final molar mass fractions from the training files and prints the last entry of the temperature and mass fraction. It also plots the temperature vs mass fraction graph. Modify the `data_directory`, `elmnt_index`, `index`, `file_prefix`, `file_extension`, and `num_files` variables as needed. Then, run the script to view the temperature vs mass fraction graph and print the last entries.

6. **Data Directory**

   The "new_Tests" directory stores the data used by the script (read_MF_data.py) to create the temperature vs mass fraction graph. Ensure that this directory contains the necessary data files with the specific naming convention "ev_co_alpha_control_data_{index}_{file_number:02d}{file_extension}" for proper execution of the script.

## Note

Before running any of the scripts, ensure that you have the required Python libraries (NumPy and Matplotlib) installed. Modify the file paths and variables in each script as per your system setup and specific use case. Additionally, make sure that the file names are the same or feel free to change them at your liking.

## Acknowledgments

I would like to extend my thanks and appreciation to my mentor Dr. Harris, for their invaluable guidance and support throughout this project.


