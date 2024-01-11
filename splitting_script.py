"""
Run with python splitting_script.py --input mosflm.stream  --start-pulseId 8 --end-pulseId 1488 --factor 96

Give input stream file, typically output from crystfel, then pulse ID start and end of sequence along with multiplying factor, this can be found in the raw or h5 files for a run if needed 
"""



import argparse
import os
import sys

def print_progress_bar(position, total, bar_length=50):
    progress = position / total
    arrow = '-' * int(round(progress * bar_length) - 1)
    spaces = ' ' * (bar_length - len(arrow))
    sys.stdout.write(f'\rProgress: [{arrow + spaces}] {progress * 100:.2f}%')
    sys.stdout.flush()

def main():
    parser = argparse.ArgumentParser(description='Extract specific pulseID chunks from stream file')
    parser.add_argument('--input', help='Path to the input stream file')
    parser.add_argument('--start-pulseId', type=int, help='The starting pulse ID')
    parser.add_argument('--end-pulseId', type=int, help='The ending pulse ID')
    parser.add_argument('--factor', type=int, help='The increment factor for pulse ID')

    args = parser.parse_args()

    output_dir = 'output_streams'
    os.makedirs(output_dir, exist_ok=True)

    pulse_ids = list(range(args.start_pulseId, args.end_pulseId + 1, args.factor))

    with open(args.input, 'r') as stream_file:
        total_size = os.path.getsize(args.input)  # Get the total size of the file for the progress bar

        header = stream_file.readline()
        while "----- Begin chunk -----" not in header:
            header += stream_file.readline()

        content = stream_file.readlines()

        for index, pulse_id in enumerate(pulse_ids):
            print(f"Processing Pulse ID: {pulse_id}")
            with open(os.path.join(output_dir, f'filtered_pulse_{index}.stream'), 'w') as output_file:
                output_file.write(header)

                for line_index, line in enumerate(content):
                    if f"header/int//entry_1/pulseId = {pulse_id}" in line:
                        chunk_start = line_index
                        while "----- Begin chunk -----" not in content[chunk_start]:
                            chunk_start -= 1

                        chunk_end = line_index
                        while "----- End chunk -----" not in content[chunk_end]:
                            chunk_end += 1

                        output_file.writelines(content[chunk_start:chunk_end + 1])

            position = index + 1
            print_progress_bar(position, len(pulse_ids))

    print("\nDone.")

if __name__ == "__main__":
    main()

