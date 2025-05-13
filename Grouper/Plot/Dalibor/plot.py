from root2mpl import *
from ROOT import *
import numpy as np
dir = "../../../../../../../mnt/hgfs/shared-2/"

plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=-0.2)


file = TFile(dir+"merged_runs.rootanalysed.root")
plt.subplots_adjust(hspace=0.05)

fig, ax = plt.subplots(figsize=(10, 5))
ylabel = "Counts/2ns"
xlabel = r"$\Delta \mathrm{t\,\,[ns]}$"
c = file.Get("SiPM_Time_Coinc")
#loop on primitives
H=0
for i in range(c.GetListOfPrimitives().GetSize()):
    h = c.GetListOfPrimitives().At(i)
    #print(h.GetName())
    if h.GetName() == "H_SiPM_Time_Coinc_SUM_9":
        H= h
DisplayTH1D(H, ax, color='black', xlabel=xlabel, ylabel=ylabel, title="Time coincidence "+r"$\beta p$", ylog=True, xlim = (-300, 200), lw=0.5, ylim=(1, 1e6))
fig.savefig("beta_p_time.png", dpi=300)


# file = TFile(dir+"../CALIBRATED/Calibrated_Down.root")
# # DisplayCanvas(file.Get("32Ar_SUM"), ylog=True, xlim = (0, 6500), lw = 0.8, rebin=5)
# # plt.savefig("32Ar_calib.png", dpi=300)


# DisplayCanvas(file.Get("D5.1/D5.1_Calibration"))
# plt.show()