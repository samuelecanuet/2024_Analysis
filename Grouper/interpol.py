import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import make_interp_spline

# Read data from the file
filename = "../../../../../../../mnt/hgfs/shared-2/2024_DATA/R-MATRIX/18N_points.txt"
data = np.loadtxt(filename, delimiter='\t')

# Separate the data into x (energy) and y (counts)
x = data[:, 0]
y = data[:, 1]


# Sort the data by x
sorted_indices = np.argsort(x)
x_sorted = x[sorted_indices]
y_sorted = y[sorted_indices]


# Identify unique elements and their indices
unique, indices = np.unique(x_sorted, return_index=True)

# Sort indices to maintain the original order
sorted_indices = np.sort(indices)

# Select the elements from x and y based on the sorted indices
x_unique = x_sorted[sorted_indices]
y_unique = y_sorted[sorted_indices]

x_sorted = x_unique
y_sorted = y_unique


# Create spline interpolation
spline = make_interp_spline(x_sorted, y_sorted)  # k=3 for cubic spline

# Define the range for interpolation
x_new = np.linspace(0, 5000, 10000)
y_new = spline(x_new)

# replace negative values with 0
y_new[y_new < 0] = 0

# Plot the original data points
plt.figure(figsize=(12, 6))
plt.plot(x_sorted, y_sorted, 'o', label='Data Points')

# Plot the spline interpolation
plt.plot(x_new, y_new, '-', label='Spline Interpolation', color='red')

# Customize the plot
plt.title('Data Points with Spline Interpolation')
plt.xlabel('Energy')
plt.ylabel('Counts')
plt.legend()
plt.yscale('log')
plt.grid(True)
# Show the plot
plt.show()


# convert this intrpolation into a histogram
# Create a histogram
plt.figure(figsize=(12, 6))
plt.hist(x_new, bins=10000, weights=y_new, histtype='step', color='blue', label='Histogram')
plt.yscale('log')
plt.show()

# #write this histogram in file 
# np.savetxt('histogram.txt', np.column_stack((x_new, y_new)), delimiter='\t', header='Energy\tCounts', fmt='%1.4e')

# #save hist in root file histogram
# import ROOT
# from ROOT import TFile, TH1F
# import array

# # Create a new ROOT file
# root_file = TFile('histogram.root', 'RECREATE')

# # Create a histogram
# histogram = TH1F('histogram', 'Energy vs Counts', 10000, 0, 5000)

# # Fill the histogram
# for i in range(len(x_new)):
#     histogram.Fill(x_new[i], y_new[i])

# # Write the histogram to the ROOT file
# histogram.Write()

# # Close the ROOT file
# root_file.Close()