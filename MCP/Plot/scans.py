from ROOT import *
from root2mpl import *
import matplotlib.pyplot as plt
import os

values = {}

color_list = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf"]

with open("../Config_Files/MCP_config.txt") as f:
    lines = f.readlines()
    for line in lines:
        if (line[0] == "#"):
            break
        line = line.split(" ")
        values[line[0]] = [float(line[1]), float(line[2])]

    


fig, ax = plt.subplots(figsize=(10, 8))
i=0
for run, value in values.items():
    print(f"Run: {run}, Value: {value}")
    filename = ""
    # get filename of a file containing run
    for file in os.listdir("/mnt/hgfs/shared-2/2025_DATA/MCP_DATA/GROUPED/"):
        if run in file:
            filename = file
            break

    if filename == "":
        print(f"No file found for run {run}")
        continue
    print(f"Filename: {filename}")

    # open the file
    file_path = f"/mnt/hgfs/shared-2/2025_DATA/MCP_DATA/GROUPED/{filename}"
    f = TFile(file_path, "READ")

    if (f.Get("H_Channel_5").GetEntries() < 1e3):
        continue
        
    if (value[1] != 2.):
        continue

    DisplayTH1D(f.Get("H_Channel_5"), ax=ax,
                title=f"Run {run} - Channel 5",
                xlabel="X-axis Label",
                ylabel="Y-axis Label", normalized=True, label=f"B = {value[1]} T, HV = {value[0]} V, (run {run})",
                xlim = (0, 40e3), color=color_list[i], rebin=10)
    i+=1
plt.legend()
plt.show()