from root2mpl import *
from ROOT import *
import numpy as np
dir = "../../../../../../../mnt/hgfs/shared-2/2024_DATA/DETECTOR_DATA/ANALYSED/"

plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=-0.2)


# file = TFile(dir+"32Ar_analysed.root")
# plt.subplots_adjust(hspace=0.05)

# fig, ax = plt.subplots(figsize=(10, 5))
# ylabel = "Counts/2ns"
# xlabel = r"$\Delta \mathrm{t\,\,[ns]}$"
# DisplayTH1D(file.Get("FakeCorrection/D1.1/Write/H_SiPM_Time_D1.1_9"), ax, color='black', xlabel=xlabel, ylabel=ylabel, title="Time coincidence "+r"$\beta p$", ylog=True, xlim = (-300, 200), lw=0.5, ylim=(1, 1e5))
# fig.savefig("beta_p_time.png", dpi=300)


# file = TFile(dir+"/32Ar_analysed.root")
# DisplayCanvas(file.Get("FakeCorrection/D5.1/Corrected_D5.1_3"), xlim = (3280, 3380), lw = 0.8, rebin=2)
# plt.show()


file = TFile(dir+"../CALIBRATED/SiPM_Calibrated.root")
c = file.Get("32Ar/Peak_14/cExp_Sim1_32Ar_Peak_14_SiPM7")
# loop on T1D on the canvas
h = c.GetPrimitive("H_SiPM_Calibrated_32Ar_Peak_14_SiPM7")
DisplayTH1D(h, xlim = (0, 7000), lw = 0.8, rebin=2)

plt.show()