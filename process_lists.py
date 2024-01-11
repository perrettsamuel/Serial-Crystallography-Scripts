import os
import re

def process_files(output_file='partialtor_pulse.lst'):
    # Create or clear the output file
    with open(output_file, 'w') as f:
        pass

    # Pattern to match the files and extract the number
    file_pattern = re.compile(r'split_list_pulse_(\d+)\.lst')

    # Loop through files in the current directory
    for filename in os.listdir('.'):
        match = file_pattern.match(filename)
        if match:
            number = match.group(1)  # Extract the number from the filename

            # Read the current file
            with open(filename, 'r') as file:
                lines = file.readlines()

            # Modify and append the lines to the output file
            with open(output_file, 'a') as output:
                for line in lines:
                    modified_line = line.strip() + f' P{number}\n'
                    output.write(modified_line)

# Call the function
process_files()

