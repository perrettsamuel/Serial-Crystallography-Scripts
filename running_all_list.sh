#!/bin/bash

# Directory containing the .lst files
lst_dir="./"

# Directory to save the output .stream files
output_stream_dir="stream_output"

# Create the output directory if they don't exist
mkdir -p "$output_stream_dir"

# Path to the original .toml file
toml_file="/gpfs/exfel/u/scratch/FXE/202202/p002808/sperrett/Lysozyme_October_Final/indexing_v4_listing/xwiz_conf.toml"

# Get the list of all .lst files in the directory
lst_files=($(ls ${lst_dir}split_list_pulse_*.lst))

# Get the total number of .lst files
total_lst_files="${#lst_files[@]}"

# Function to display a progress bar
display_progress() {
    local progress="${1}"
    local total="${2}"
    local filled_length=$(($progress*40/$total))
    local bar_fill=$(printf "%0.s#" $(seq 1 $filled_length))
    local bar_empty=$(printf "%0.s-" $(seq 1 $(($total-$filled_length))))
    printf "\rProgress: [${bar_fill}${bar_empty}] ${progress}/${total}"
}

# Iterate through each .lst file
for lst_file in "${lst_files[@]}"; do
    # Extract the list name from the .lst file name
    list_name=$(basename $lst_file .lst)
    
    # Update the .toml file to use the current .lst file and the extracted list name
    sed -i "s|frames_list_file = \".*\"|frames_list_file = \"${lst_file}\"|" "$toml_file"
    sed -i "s/list_prefix = .*/list_prefix = \"${list_name}\"/" "$toml_file"

    # Run the xwiz-workflow -a command with the modified .toml file
    xwiz-workflow -a 
    
    # Move the resulting .stream file to the stream output directory
    mv "${list_name}.stream" "${output_stream_dir}/${list_name}.stream"
    
    # Display the progress bar
    display_progress $((++count)) $total_lst_files
done

# Print a newline for clarity after progress bar completion
echo ""

