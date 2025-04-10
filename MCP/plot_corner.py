import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import time

from ROOT import TFile

# Enable interactive mode
plt.ion()

# Define colors
colors = list(mcolors.TABLEAU_COLORS)
del colors[0]  # Remove white color

while True:
    plt.clf()  # Clear previous plot
    f = TFile("/mnt/hgfs/shared-2/2025_DATA/MCP_DATA/03_TEST/005_MCP_1p9kV_BeamScan.fast/005_MCP_1p9kV_BeamScan_0001.root", "READ")
    # DisplayTH2D(f.Get("h_Image"))  # Display TH2D
    
    with open("guess_corner_2025.txt", "r") as f:
        counter =0
        for line in f:
         
            try:
                line = line.strip().split(" ")
                # Ensure enough values exist in the line
                if len(line) < 8:
                    print(f"Skipping invalid line: {line}")
                    continue

                # Convert string values to float
                points = [float(x) for x in line]

                # Select color
                color = colors[counter % len(colors)]
                counter += 1

                # Plot scatter points
                for i in range(1, 9, 2):
                    plt.scatter(points[i], points[i+1], color=color, s=5)
                    if i == 1:
                        plt.text(points[i], points[i+1], f"{int(points[0])}", color=color)
            
            except ValueError as e:
                print(f"Error parsing line: {line}, {e}")

    plt.pause(1)  # Pause to update the plot
