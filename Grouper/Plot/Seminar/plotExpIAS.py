from root2mpl import *
from ROOT import *
import numpy as np

import matplotlib as mpl

YEAR = 2025


file  =TFile(f"/mnt/hgfs/shared-2/{YEAR}_DATA/DETECTOR_DATA/ANALYSED/32Ar_{YEAR}_analysed_all.root", "READ")

Hs = None
Hc = None
Hnc = None
for det in range(5, 9):
    for strip in range(1, 6):
        if det == 5 and strip == 1:
            Hc = file.Get(f"RandomCorrection/D{det}.{strip}/Write/H_Coinc_CorrectedD{det}.{strip}_3")
            Hnc = file.Get(f"RandomCorrection/D{det}.{strip}/Write/H_NoCoinc_CorrectedD{det}.{strip}_3")

        else:
            Hc.Add(file.Get(f"RandomCorrection/D{det}.{strip}/Write/H_Coinc_CorrectedD{det}.{strip}_3"))
            Hnc.Add(file.Get(f"RandomCorrection/D{det}.{strip}/Write/H_NoCoinc_CorrectedD{det}.{strip}_3"))

Hs = Hc.Clone("H_Sum_Corrected_3")
Hs.Add(Hnc)

fig, ax = plt.subplots(figsize=(8,6))
h1 = DisplayTH1D(Hs, ax=ax, label="Sum", color="black", title="")
h2 = DisplayTH1D(Hc, ax=ax, label="Coincidence", color="red", title="")
h3 = DisplayTH1D(Hnc, ax=ax, label="No Coincidence", color="blue", title="")
ax.set_xlim(3320, 3380)
ax.set_ylim(0, 31000)

## blank
h1.set_visible(False)
h2.set_visible(False)
h3.set_visible(False)
plt.savefig(f"Down_{YEAR}_blank.png")
# coinc
h2.set_visible(True)
plt.savefig(f"Down_{YEAR}_coinc.png")
# nocoinc and singles
h3.set_visible(True)
h1.set_visible(True)
plt.savefig(f"Down_{YEAR}.png")
# plt.show()