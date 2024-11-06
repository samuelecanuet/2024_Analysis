from root2mpl import *
from ROOT import *
import numpy as np
dir = "../../../../../../../mnt/hgfs/shared-2/2024_DATA/DETECTOR_DATA/CALIBRATED/"


file = TFile(dir + "Calibrated.root", "READ")

sizex = 14
sizey = 4
ylim = 50e3

fig1, ax1 = plt.subplots(figsize = (sizex, sizey))
DisplayCanvas(file.Get("D5.1_Calibration"), ax=ax1, ylog = None, xlim = (3200, 3450), ylim = (1, ylim), lw=0.5)
plt.savefig("D5.1.png", dpi=300)


fig, ax  = plt.subplots(figsize = (sizex, sizey))
DisplayCanvas(file.Get("32Ar_SUM"), ax=ax, ylog = True, xlim = (0, 6500), ylim = (1, ylim), lw = 0.5, legend = None, title="", color="black")
plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.2)
plt.savefig("test.png", dpi=300) 