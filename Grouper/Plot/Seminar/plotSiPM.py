from root2mpl import *
from ROOT import *
import numpy as np

import matplotlib as mpl
# mpl.rc('font', family='serif', serif='Linguistics Pro')
# mpl.rc('text', usetex=TRue)
mpl.rc('mathtext', fontset='custom', bf='Linguistics Pro:bold')


## full
file = TFile("/mnt/hgfs/shared-2/2025_DATA/DETECTOR_DATA/GROUPED/merged_2025_Grouped.root", "READ")

file = TFile("/mnt/hgfs/shared-2/2025_DATA/DETECTOR_DATA/CALIBRATED/SiPM_Calibrated_2025_all.root", "READ")


rebin = 25
peaklist = [1, 7, 14, 27]
H = {}
for peak in peaklist:
    c = file.Get("32Ar/SiPM_105/cExp_Sim_Low32Ar_SiPM5")
    for pad in range(c.GetListOfPrimitives().GetEntries()):
        prim = c.GetListOfPrimitives().At(pad)
        for h in prim.GetListOfPrimitives():
            if h.GetName() == f"H_SiPMLow_Calibrated_32Ar_Peak_{peak}_SiPM5":
                H[peak] = h

## rebin
for peak in peaklist:
    H[peak].Rebin(rebin)

### H[0] sum of all
H[0] = H[14].Clone(f"H_SiPMLow_Calibrated_32Ar_Peak_sum_SiPM5")
deletepeak=[]
for peak in peaklist:
    if peak != 14:
        if H[peak].GetNbinsX() == H[0].GetNbinsX():
            H[0].Add(H[peak], 1.0)
        else:
            deletepeak.append(peak)
    
peaklist.insert(0, 0)

blues = plt.get_cmap("Blues")
C = ["black"]
C.extend([blues(i) for i in np.linspace(0.5, 1.0, len(peaklist)-1)])
h= {}
fig, ax = plt.subplots(figsize=(12,6))
for peak in peaklist:
    label = rf"$p_{{{peak}}}$"
    if peak == 0:
        label = r"All p"
    if (peak not in deletepeak):
        h[peak] = DisplayTH1D(H[peak], ax=ax, label=label, lw=1.5, color=C[peaklist.index(peak)], title="", ylog=True)

ax.set_ylim(1, 1e5)
ax.set_xlim(0, 8000)
ax.set_xticklabels([int(i) for i in ax.get_xticks()], fontsize=14)
ax.set_yticks([])
ax.set_ylabel("Normalized counts", fontsize=14)
lg1 = ax.legend(title=r"$\textbf{Coincidence with}$", fontsize=14)
plt.setp(lg1.get_title(),fontsize=14,fontweight='bold')	# legend title bold
plt.gca().add_artist(lg1)



for peak, hh in h.items():
    if (peak != 0):
        hh.set_visible(False)
plt.savefig("SIPM_2025_single.png",bbox_inches='tight',pad_inches=0, dpi=300)
for peak, hh in h.items():
    hh.set_visible(True)
plt.savefig("SIPM_2025.png",bbox_inches='tight',pad_inches=0, dpi=300)
# plt.show()