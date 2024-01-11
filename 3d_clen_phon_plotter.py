import os
import numpy as np
from scipy.stats import skew
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import re
from scipy.interpolate import griddata
from matplotlib import colors
from numpy.polynomial.polynomial import polygrid2d
from scipy.interpolate import bisplrep, bisplev, Rbf
import itertools
import scipy.spatial
from sklearn.preprocessing import StandardScaler


def gaussian_params(data):
    """Extract Gaussian parameters from data."""
    mu = np.mean(data)
    sigma = np.std(data)
    skewness = skew(data)
    gaussianity = skewness * sigma  # Calculate Gaussianity
    return mu, sigma, skewness, gaussianity


def extract_cell_parameters(line):
    """Extract cell parameters from a line."""
    pattern = r'Cell parameters (\d+\.\d+) (\d+\.\d+) (\d+\.\d+) nm, (\d+\.\d+) (\d+\.\d+) (\d+\.\d+) deg'
    match = re.search(pattern, line)
    if match:
        a, b, c, alpha, beta, gamma = map(float, match.groups())
        return a * 10, b * 10, c * 10, alpha, beta, gamma
    return None


def find_and_extract_cell_parameters_updated(stream_file):
    """Find and extract cell parameters from a stream file."""
    with open(stream_file, 'r') as file:
        for line in file:
            if 'Cell parameters' in line:
                cell_params = extract_cell_parameters(line)
                if cell_params:
                    yield cell_params


def extract_number_of_crystals(stream_file):
    """Count the number of occurrences of '--- Begin crystal' in the stream file."""
    with open(stream_file, 'r') as file:
        return sum(1 for line in file if '--- Begin crystal' in line)


def generate_smooth_surface(ax, data, xlim, ylim, zlim, cmap="cividis"):
    """Generate a smooth surface plot."""
    # Create a grid of linearly spaced points
    grid_x, grid_y = np.mgrid[xlim[0]:xlim[1]:100j, ylim[0]:ylim[1]:100j]
    grid_z = griddata((data[:, 2], data[:, 1]), data[:, 3], (grid_x, grid_y), method='linear')
    ax.plot_surface(grid_x, grid_y, grid_z, cmap=cmap, vmin=0.65, vmax=1)


def process_stream_file(filename):
    """Process a stream file, extracting cell parameters and counting crystals."""
    cell_params_list = []
    crystal_count = 0
    
    with open(filename, 'r') as file:
        for line in file:
            if 'Cell parameters' in line:
                cell_params = extract_cell_parameters(line)
                if cell_params:
                    cell_params_list.append(cell_params)
            elif '--- Begin crystal' in line:
                crystal_count += 1
                
    return cell_params_list, crystal_count


colors = np.loadtxt('batlow.txt')
batlow_cmap = plt.cm.colors.ListedColormap(colors)

def plot_smooth_surface(data, title, fig, ax, vmin=0.1, vmax=1, method='rbf'):
    # Load the color palette from the text file
    colors = np.loadtxt('batlow.txt')
    batlow_cmap = plt.cm.colors.ListedColormap(colors)

    # 1. Scale the data
    scaler = StandardScaler()
    scaled_data = scaler.fit_transform(data[:, 1:3])  # We only scale CLEN and Photon Energy
    X_scaled, Y_scaled = scaled_data[:, 0], scaled_data[:, 1]
    
    Z = data[:, 3]
    
    # Plot the raw data as small transparent grey dots
    ax.scatter(data[:, 1], data[:, 2], Z, color='grey', alpha=0.2, s=1)
    
    # Define interpolation grid using the range of the scaled data
    x_min, x_max = X_scaled.min(), X_scaled.max()
    y_min, y_max = Y_scaled.min(), Y_scaled.max()
    grid_x_scaled, grid_y_scaled = np.mgrid[x_min:x_max:1000j, y_min:y_max:1000j]
    
    if method == 'rbf':
        rbf = Rbf(X_scaled, Y_scaled, Z, function='multiquadric', smooth=0.2)
        grid_z = rbf(grid_x_scaled, grid_y_scaled)
    else:
        grid_z = griddata((X_scaled, Y_scaled), Z, (grid_x_scaled, grid_y_scaled), method='cubic')
    
    # Handle NaN values in grid_z and cap values
    grid_z[np.isnan(grid_z)] = np.nanmin(grid_z)
    grid_z[grid_z > 1] = 1
    grid_z[grid_z < 0] = 0
    # Convert grid_x and grid_y back to original scale
    grid_x, grid_y = scaler.inverse_transform(np.vstack([grid_x_scaled.ravel(), grid_y_scaled.ravel()]).T).T
    grid_x = grid_x.reshape(grid_x_scaled.shape)
    grid_y = grid_y.reshape(grid_y_scaled.shape)
    
    surf = ax.plot_surface(grid_x, grid_y, grid_z, cmap=batlow_cmap, vmin=vmin, vmax=vmax, alpha=1)  # Using batlow color palette and setting transparency
    weights = grid_z[grid_z > 0.98]
    centroid_x = np.sum(grid_x[grid_z > 0.98] * weights) / np.sum(weights)
    centroid_y = np.sum(grid_y[grid_z > 0.98] * weights) / np.sum(weights)
    centroid_z = np.sum(grid_z[grid_z > 0.98] * weights) / np.sum(weights)
    
    # Mark the weighted centroid with a red cross
    ax.scatter(0.23572, 9295, 1.001, color='red',edgecolor='black', linewidth=1.5, marker='o', s=100, label='Weighted Centroid', zorder=10000000)
    print(f"Weighted Centroid: Photon Energy = {centroid_x}, CLEN = {centroid_y}, Value = {centroid_z}")



    # Add colorbar
    fig.colorbar(surf, ax=ax, label='Indexing Figure of Merit', shrink=0.7)
    ax.set_zlim(0,1)
    ax.set_xlabel('Detector Length (m)')
    ax.set_ylabel('Photon Energy (ev)')
    ax.set_zlabel('Indexing Figure of Merit')
    ax.legend()  # Display the legend to differentiate max points


# Initialization
dir_path = os.getcwd()
gaussianity_factors = []
max_crystals_per_pulse = {}
max_FoM_per_pulse = {}


# Data collection loop
for filename in os.listdir(dir_path):
    if filename.startswith("split_list_pulse") and filename.endswith(".stream"):
        filepath = os.path.join(dir_path, filename)
        all_params, num_crystals = process_stream_file(filepath)

        split_name = filename.split("_")
        if len(split_name) < 9:
            print(f"Filename {filename} doesn't match the expected format.")
            continue

        pulse_id = float(split_name[3])
        clen = float(split_name[5])
        photon_energy = float(split_name[8].replace(".stream", ""))
        


        if pulse_id not in max_crystals_per_pulse or num_crystals > max_crystals_per_pulse[pulse_id]:
            max_crystals_per_pulse[pulse_id] = num_crystals

        if all_params:
            params = [gaussian_params([item[i] for item in all_params]) for i in range(6)]
            averages = np.mean(params, axis=0)
            gaussianity = averages[3]
            gaussianity_factors.append((pulse_id, photon_energy, clen, gaussianity))

indexing_FoM_data = []
for (pulse_id, photon_energy, clen, gaussianity) in gaussianity_factors:
    max_crystals = max_crystals_per_pulse[pulse_id]
    normalized_crystals = extract_number_of_crystals(os.path.join(dir_path,f"split_list_pulse_{int(pulse_id)}_clen_{clen}_photon_energy_{int(photon_energy) if photon_energy.is_integer() else photon_energy}.stream")) 
    indexing_FoM = normalized_crystals * (1 - abs(gaussianity))
    indexing_FoM_data.append((pulse_id, clen, photon_energy, indexing_FoM))
    # Track the max FoM for each pulse ID
    if pulse_id not in max_FoM_per_pulse or indexing_FoM > max_FoM_per_pulse[pulse_id]:
        max_FoM_per_pulse[pulse_id] = indexing_FoM

# Normalize the Indexing FoM based on max FoM for each pulse ID
for i, data in enumerate(indexing_FoM_data):
    pulse_id, clen, photon_energy, indexing_FoM = data
    indexing_FoM_data[i] = (pulse_id, clen, photon_energy, indexing_FoM/max_FoM_per_pulse[pulse_id])

# Convert to numpy array
indexing_FoM_data_array = np.array(indexing_FoM_data)
unique_pulse_ids = np.unique(indexing_FoM_data_array[:, 0])

# Specify the pulse ID you want to plot
target_pulse_id = 0.0

# Extract the data for the desired Pulse ID
filtered_data = indexing_FoM_data_array[indexing_FoM_data_array[:, 0] == target_pulse_id]

# Get the index of the highest Figure of Merit in the filtered data
highest_idx = np.argmax(filtered_data[:, 3])

# Extract min and max values for energy and CLEN
min_energy = np.min(filtered_data[:, 2])
max_energy = np.max(filtered_data[:, 2])
min_clen = np.min(filtered_data[:, 1])
max_clen = np.max(filtered_data[:, 1])


# Zoomed data extraction
zoomed_data = filtered_data[(filtered_data[:, 2] >= 9250) & (filtered_data[:, 2] <= 9350) & (filtered_data[:, 1] >= 0.233) & (filtered_data[:, 1] <= 0.239)]

# Get the index of the highest Figure of Merit in the zoomed data
highest_idx_zoomed = np.argmax(zoomed_data[:, 3])


# Full range plots
fig1 = plt.figure(figsize=(10, 10))
ax1 = fig1.add_subplot(111, projection='3d')
plot_smooth_surface(filtered_data, f'Pulse ID {target_pulse_id} Full Range', fig1, ax1)

# Figure 1: Original scatter plot
fig2, ax2 = plt.subplots(figsize=(10, 10), subplot_kw={'projection': '3d'})
ax2.scatter(filtered_data[:, 2], filtered_data[:, 1], filtered_data[:, 3], c=filtered_data[:, 3], cmap=batlow_cmap, s=20, alpha=0.8)
ax2.scatter(9295, 0.23572, filtered_data[highest_idx, 3], c='red', s=100, marker='o')
ax2.set_title(f'Pulse ID {target_pulse_id} - Scatter Plot')
ax2.set_xlabel('Photon Energy')
ax2.set_ylabel('CLEN')
ax2.set_zlabel('Indexing Figure of Merit')
ax2.set_zlim(0, 1)

# Zoomed data extraction
zoomed_data = filtered_data[(filtered_data[:, 2] >= 9250) & (filtered_data[:, 2] <= 9350) & (filtered_data[:, 1] >= 0.233) & (filtered_data[:, 1] <= 0.239)]
highest_idx_zoomed = np.argmax(zoomed_data[:, 3])

# Zoomed-in range smooth surface plot
fig3 = plt.figure(figsize=(10, 10))
ax3 = fig3.add_subplot(111, projection='3d')
plot_smooth_surface(zoomed_data, f'Pulse ID {target_pulse_id} Zoomed Range', fig3, ax3)

# Figure 3: Zoomed-in scatter plot
fig4, ax4 = plt.subplots(figsize=(10, 10), subplot_kw={'projection': '3d'})
ax4.scatter(zoomed_data[:, 2], zoomed_data[:, 1], zoomed_data[:, 3], c=zoomed_data[:, 3], cmap=batlow_cmap, s=5)
ax4.scatter(zoomed_data[highest_idx_zoomed, 2], zoomed_data[highest_idx_zoomed, 1], zoomed_data[highest_idx_zoomed, 3], c='red', s=100, marker='*')
ax4.set_title(f'Zoomed View: Pulse ID {target_pulse_id} - Scatter Plot')
ax4.set_xlabel('Photon Energy')
ax4.set_ylabel('CLEN')
ax4.set_zlabel('Indexing Figure of Merit')
ax4.set_xlim(9250, 9350)
ax4.set_ylim(0.233, 0.239)
ax4.set_zlim(0.7, 1)

plt.tight_layout()
plt.show()



