import os
import re
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import warnings
from scipy.optimize import curve_fit
import random
from matplotlib.ticker import MaxNLocator
def gaussian(x, mean, amplitude, standard_deviation):
    if standard_deviation <= 0:  # Prevent non-positive standard deviation
        return np.inf
    return amplitude * np.exp(- ((x - mean)**2 / (2 * standard_deviation ** 2)))

# Function to extract cell parameters from a line
def extract_cell_parameters(line):
    pattern = r'Cell parameters (\d+\.\d+) (\d+\.\d+) (\d+\.\d+) nm, (\d+\.\d+) (\d+\.\d+) (\d+\.\d+) deg'
    match = re.search(pattern, line)
    if match:
        a, b, c, alpha, beta, gamma = map(float, match.groups())
        return a * 10, b * 10, c * 10, alpha, beta, gamma
    return None

def find_and_extract_cell_parameters_updated(stream_file):
    with open(stream_file, 'r') as file:
        for line in file:
            if 'Cell parameters' in line:
                cell_params = extract_cell_parameters(line)
                if cell_params:
                    yield cell_params  # Using yield to make this function a generator


# New function to plot histogram with Gaussian fit
def plot_histogram_with_fit(ax, values, fit_params, bins, param_label):
    ax.hist(values, bins=bins, color='skyblue', edgecolor='black', alpha=0.7, label=param_label)
    x = np.linspace(min(values), max(values), 1000)
    ax.plot(x, gaussian(x, *fit_params), color='red', label='Gaussian Fit')
    ax.legend()

# Function to extract unit cell parameters from .cell file
def extract_unit_cell_parameters(cell_file):
    params = {}
    with open(cell_file, 'r') as file:
        for line in file:
            if '=' in line and 'A' in line:
                key, value = line.split('=')
                params[key.strip()] = float(value.split()[0].strip())
            elif '=' in line and 'deg' in line:
                key, value = line.split('=')
                params[key.strip()] = float(value.split()[0].strip())
    return params

# Initialization
directory = os.getcwd()
num_cases = []
cell_file_params = None
# Initialize lists to store the standard deviations from Gaussian fitting
std_devs = [[], [], [], [], [], []]
all_pulse_params = []  # Store parameters for each pulse separately

# Read .cell file parameters
for filename in os.listdir(directory):
    if filename.endswith('.cell'):
        cell_file = os.path.join(directory, filename)
        cell_file_params = extract_unit_cell_parameters(cell_file)
        break

if cell_file_params:
    print("Unit cell parameters extracted from .cell file:", cell_file_params)
else:
    warnings.warn("No .cell file found in the directory.")

# Process .stream files
a_avgs, b_avgs, c_avgs, alpha_avgs, beta_avgs, gamma_avgs = ([] for _ in range(6))
avg_pct_diffs = []
a_errs, b_errs, c_errs, alpha_errs, beta_errs, gamma_errs = ([] for _ in range(6))

# Collect all stream files and sort them numerically
stream_files = sorted([filename for filename in os.listdir(directory) if filename.startswith('filtered_pulse_') and filename.endswith('.stream')],
                      key=lambda x: int(re.search(r'filtered_pulse_(\d+)_corrected\.stream', x).group(1)))

for filename in stream_files:
    stream_file = os.path.join(directory, filename)
    pulse_number = int(re.search(r'filtered_pulse_(\d+)_corrected\.stream', filename).group(1))  # Extract pulse number
    pulse_ID = None  # Initialize pulse_ID to None

    with open(stream_file, 'r') as file:
        for line in file:
            if 'header/int//entry_1/pulseId' in line:
                pulse_ID = int(line.split()[-1])  # Extract pulse_ID
                break  # Break the loop after extracting the pulse_ID

    pulse_params = list(find_and_extract_cell_parameters_updated(stream_file))  # Collect parameters for this pulse
    if pulse_params:
        all_pulse_params.append(pulse_params)  # Add this pulse's parameters to the list
        print(f"Found {len(pulse_params)} cases in pulse {pulse_number}, corresponding to stream file {filename}, Pulse ID: {pulse_ID}")  # Print correspondence

        if pulse_params:  
            a_values, b_values, c_values, alpha_values, beta_values, gamma_values = zip(*pulse_params)

            for values, avgs, std_dev_list in zip(
                [a_values, b_values, c_values, alpha_values, beta_values, gamma_values],
                [a_avgs, b_avgs, c_avgs, alpha_avgs, beta_avgs, gamma_avgs],
                std_devs
            ):
                num_bins = max(1, int((max(values) - min(values)) / 0.01) + 1)  # Ensure at least 1 bin
                bins = np.linspace(min(values), max(values), num_bins)

                y, _ = np.histogram(values, bins=bins)
                x = (bins[:-1] + bins[1:]) / 2

                bounds = ([0, 0, 0], [np.inf, np.inf, np.inf])
                try:
                    popt, _ = curve_fit(gaussian, x, y, p0=[np.mean(values), max(y), np.std(values)], bounds=bounds)
                    avgs.append(popt[0])
                    std_dev_list.append(popt[2])
                except Exception as e:
                    print(f"Failed to fit Gaussian for values {values}. Error: {e}")
                    avgs.append(np.mean(values))
                    std_dev_list.append(np.std(values))

            num_cases.append(len(pulse_params))  # Corrected to pulse_params

            pct_diffs = [
                abs((a_avgs[-1] - cell_file_params['a']) / cell_file_params['a']) * 100,
                abs((b_avgs[-1] - cell_file_params['b']) / cell_file_params['b']) * 100,
                abs((c_avgs[-1] - cell_file_params['c']) / cell_file_params['c']) * 100,
                abs((alpha_avgs[-1] - cell_file_params['al']) / cell_file_params['al']) * 100,
                abs((beta_avgs[-1] - cell_file_params['be']) / cell_file_params['be']) * 100,
                abs((gamma_avgs[-1] - cell_file_params['ga']) / cell_file_params['ga']) * 100
            ]
            avg_pct_diffs.append(np.mean(pct_diffs))




# Calculate average standard deviations for each pulse over the six parameters
avg_std_devs = [np.mean(std_dev_list) for std_dev_list in zip(*std_devs)]

# Plotting the six separate cell parameters with the .cell values as a red dashed line
param_avgs = [a_avgs, b_avgs, c_avgs, alpha_avgs, beta_avgs, gamma_avgs]
param_labels = ['a', 'b', 'c', 'α', 'β', 'γ']
 # Corrected the labels
param_units = ['Å', 'Å', 'Å', 'deg', 'deg', 'deg']
colors = sns.color_palette("husl", 8)  # Using seaborn for better colors

fig, axes = plt.subplots(2, 3, figsize=(15, 10), facecolor='w', edgecolor='k')
fig.subplots_adjust(hspace=0.4, wspace=0.4)
plt.suptitle('Average Cell Parameters vs. Pulse Number', fontsize=20)

for i, ax in enumerate(axes.flat):
    x_values = range(1, len(param_avgs[i]) + 1)  # Shift x-values by 1
    ax.errorbar(x_values, param_avgs[i], yerr=std_devs[i], fmt='o-', 
                color=colors[i], ecolor='lightgrey', elinewidth=1, capsize=3, label=param_labels[i])
    ax.set_ylim(cell_file_params[list(cell_file_params.keys())[i]] - 0.2, cell_file_params[list(cell_file_params.keys())[i]] + 0.2)
    #ax.axhline(y=cell_file_params[list(cell_file_params.keys())[i]], color='r', linestyle='--')
    ax.set_xlabel('Pulse Number', fontsize=14)
    ax.set_ylabel(f'{param_labels[i]} ({param_units[i]})', fontsize=14)
   
    ax.set_xlim(0, len(param_avgs[i]) + 1)
    ax.legend(fontsize=14)

    # Add colored background every 4th pulse number
    for pulse_number in range(4, len(param_avgs[i]) + 1, 4):
        ax.axvspan(pulse_number - 3, pulse_number, facecolor='lightblue', alpha=0.3)



normalized_num_cases = [(x/15935)*100 for x in num_cases]

# Plotting total hits per pulse
fig, ax1 = plt.subplots(figsize=(12, 6), facecolor='w', edgecolor='k')
ax1.plot(range(1, len(normalized_num_cases) + 1), normalized_num_cases, 'o-', color=colors[7], label=f'Index Crystals (Total: {int(np.sum(num_cases))})', markersize=8)
ax1.set_xlabel('Pulse Number', fontsize=16)
ax1.set_ylabel('Total Indexed Rate of All Frames (%)', fontsize=16)
# Setting x-axis ticks to be every integer value
ax1.set_xticks(range(len(normalized_num_cases)))
ax1.set_xlim(0.5, len(normalized_num_cases) + 0.5)  # Adjusting the x-axis limits to start from 0
ax1.set_ylim(0, 90)
ax1.legend(fontsize=14)

# Add colored background every 4th pulse number
for pulse_number in range(4, len(normalized_num_cases) + 1, 4):
    ax1.axvspan(pulse_number - 3, pulse_number, facecolor='lightblue', alpha=0.3)



# Normalizing num_cases.
normalized_num_cases = [(x/15935)*100 for x in num_cases]

# Aspect ratio calculation for the new plot
aspect_ratio = 861.7 / 269.3
# Calculate the figure height to match the aspect ratio with a reasonable width
fig_width = 10  # You can adjust this width as needed.
fig_height = fig_width / aspect_ratio

# Creating a separate figure for the pulses with the desired aspect ratio.
fig2, ax2 = plt.subplots(figsize=(fig_width, fig_height), facecolor='w', edgecolor='k')

# Plot for the first 4 pulses
ax2.plot(range(1, 5), normalized_num_cases[:4], 'o-', color='darkred', label='1$^{st}$ Droplet', markersize=8)

# Plot for the pulses 5-9, renumbered to 1-5, using squares as markers
renumbered_pulses = range(1, 5)  # Renumbered hit numbers
ax2.plot(renumbered_pulses, normalized_num_cases[4:8], 's-', color='darkblue', label='2$^{nd}$ Droplet', markersize=8)

# Set the axis labels and ticks
ax2.set_xlabel('Hit Number in Droplet', fontsize=16)
ax2.set_ylabel('Total Indexed Rate (%)', fontsize=16)
ax2.set_ylim(0)  # Setting the lower y-axis limit to 0

# Setting integer ticks on the x-axis
ax2.xaxis.set_major_locator(MaxNLocator(integer=True))

# Add a legend
ax2.legend(fontsize=14)

# Ensure everything fits well without cutting off labels
plt.tight_layout()

# Display the plot
plt.show()

