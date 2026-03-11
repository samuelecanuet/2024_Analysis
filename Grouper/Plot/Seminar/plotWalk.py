from ROOT import *
import numpy as np
import matplotlib.pyplot as plt
from root2mpl import *
import os
import matplotlib as mpl



file= TFile("/mnt/hgfs/shared-2/2025_DATA/DETECTOR_DATA/GROUPED/run_098_data_32Ar_grouped_full.root", "READ")


hNOcorrected2D = file.Get("SiPMLow/SiPMLow_Time/BetaLo2/SiPM_ChannelTime_Cutted_BetaLo2")
hNOcorrectedx = hNOcorrected2D.ProjectionX("hNOcorrectedx")
hcorrected2D = file.Get("SiPMLow/SiPMLow_Time/BetaLo2/SiPM_ChannelTime_Walk_Corrected_BetaLo2")
hcorrectedx = hcorrected2D.ProjectionX("hcorrectedx")

fig, ax = plt.subplots(figsize=(10,6))
h1 = DisplayTH1D(hNOcorrectedx, ax=ax, label="Before walk correction", color="red", title="", ylog=True)
h2 = DisplayTH1D(hcorrectedx, ax=ax, label="After walk correction", color="blue", title="", ylog=True)
ax.set_ylim(1, 1e5)
plt.show()
