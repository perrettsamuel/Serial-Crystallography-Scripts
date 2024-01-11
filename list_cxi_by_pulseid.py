"""
Generates cxi file list by pulse ID, that can be separately fed to crystfel if indexing separate pulses needs original cxi file 
python list_cxi_by_pulseid.py p2808_r048.cxi  --start-pulseId 8 --end-pulseId 1488 --factor 96 
""""


import os
import h5py
import argparse
import numpy as np

def get_pulse_id_from_cxi(cxi_filepath):
    with h5py.File(cxi_filepath, 'r') as cxi_file:
        return np.array(cxi_file['entry_1/pulseId'])

def print_progress_bar(position, total, bar_length=50):
    progress = position / total
    arrow = '-' * int(round(progress * bar_length) - 1)
    spaces = ' ' * (bar_length - len(arrow))
    print(f'\rProgress: [{arrow + spaces}] {progress * 100:.2f}%', end='')

def main():
    parser = argparse.ArgumentParser(description='List CXI files based on Pulse ID')
    parser.add_argument('cxi_file', help='Path to the example CXI file')
    parser.add_argument('--start-pulseId', type=int, help='The starting pulse ID')
    parser.add_argument('--end-pulseId', type=int, help='The ending pulse ID')
    parser.add_argument('--factor', type=int, help='The increment factor for pulse ID')

    args = parser.parse_args()
    
    pulse_ids_in_file = get_pulse_id_from_cxi(args.cxi_file)
    
    matching_frames_dict = {}
    
    for i, pulse_id_in_file in enumerate(pulse_ids_in_file):
        if args.start_pulseId <= pulse_id_in_file <= args.end_pulseId and (pulse_id_in_file - args.start_pulseId) % args.factor == 0:
            if pulse_id_in_file not in matching_frames_dict:
                matching_frames_dict[pulse_id_in_file] = []
            matching_frames_dict[pulse_id_in_file].append(i)

        print_progress_bar(i, len(pulse_ids_in_file))
    
    output_directory = 'output_lists'
    os.makedirs(output_directory, exist_ok=True)
    
    for index, frames in enumerate(matching_frames_dict.values()):
        output_filename = os.path.join(output_directory, f"split_list_pulse_{index}.lst")
        
        with open(output_filename, 'w') as output_file:
            for frame in frames:
                output_file.write(f"{args.cxi_file} //{frame}\n")

    print("\nDone.")

if __name__ == "__main__":
    main()

