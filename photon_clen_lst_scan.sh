#!/bin/bash

gaussian_sample() {
  local mean=$1
  local std_dev=$2
  local num_samples=$3

  awk -v mean=$mean -v std_dev=$std_dev -v num_samples=$num_samples \
    'BEGIN {
       srand()
       for (i=1; i<=num_samples; i++) {
         u1 = rand(); u2 = rand();
         z0 = sqrt(-2 * log(u1)) * cos(2 * 3.14159265359 * u2);
         print mean + z0 * std_dev;
       }
     }'
}

# Directory containing .lst files
lst_dir="."

# Parameters for Gaussian sampling
mean_clen=0.234
std_dev_clen=0.001
samples_clen=15

mean_photon_energy=9214
std_dev_photon_energy=15
samples_photon_energy=15

# Generate values
clen_values=($(gaussian_sample $mean_clen $std_dev_clen $samples_clen))
photon_energy_values=($(gaussian_sample $mean_photon_energy $std_dev_photon_energy $samples_photon_energy))

# Paths to the .toml and .geom files (unchanged from your script)
toml_file="/gpfs/exfel/u/scratch/FXE/202202/p002808/sperrett/Lysozyme_October_Final/indexing_v7_3dscan_pulse0/xwiz_conf.toml"
geom_file="/gpfs/exfel/u/scratch/FXE/202202/p002808/sperrett/Lysozyme_October_Final/indexing_v7_3dscan_pulse0/output_v17_clen_refines_s.geom"

# Rest of your script (unchanged)
for lst_file in ${lst_dir}/*.lst; do
    x=$(basename $lst_file .lst | cut -d'_' -f4)
    
    for clen in "${clen_values[@]}"; do
        for photon_energy in "${photon_energy_values[@]}"; do
            sed -i "s|frames_list_file = \".*\"|frames_list_file = \"$lst_file\"|" "$toml_file"
            sed -i "s|^list_prefix = .*|list_prefix = \"split_list_pulse_${x}_clen_${clen}_photon_energy_${photon_energy}\"|" "$toml_file"
            sed -i "s|^clen = .*|clen = ${clen}|" "$geom_file"
            sed -i "s|^photon_energy = .*|photon_energy = ${photon_energy}|" "$geom_file"
            echo "Processing ${lst_file} with clen=${clen} and photon_energy=${photon_energy}"
            xwiz-workflow -a 
        done
    done
done
