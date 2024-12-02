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


file = TFile(dir+"../CALIBRATED/Calibrated_Down.root")
# DisplayCanvas(file.Get("32Ar_SUM"), ylog=True, xlim = (0, 6500), lw = 0.8, rebin=5)
# plt.savefig("32Ar_calib.png", dpi=300)


DisplayCanvas(file.Get("D5.1/D5.1_Calibration"))
plt.show()