from root2mpl import *
import ROOT
import numpy as np
dir = "../../../../../../../mnt/hgfs/shared-2/2024_DATA/DETECTOR_DATA/ANALYSED/"

plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=-0.2)


fig, ax = plt.subplots(figsize=(16, 5))

file = ROOT.TFile("/mnt/hgfs/shared-2/2025_DATA/DETECTOR_DATA/CALIBRATED/Calibrated_2025.root", "READ")
c = file.Get("32Ar_SUM_Down")
for primitive in c.GetListOfPrimitives():
    if primitive.InheritsFrom("TH1D") and primitive.GetName() == "32Ar_Exp_All_Down":
        integral = primitive.Integral(33000, 33800)
        h1 = primitive.Clone("h1")
        

file2 = ROOT.TFile("/mnt/hgfs/shared-2/2024_DATA/DETECTOR_DATA/CALIBRATED/Calibrated_2024.root", "READ")
c = file2.Get("32Ar_SUM_Down")
for primitive in c.GetListOfPrimitives():
    if primitive.InheritsFrom("TH1D") and primitive.GetName() == "32Ar_Exp_All_Down":
        integral1 = primitive.Integral(33000, 33800)
        h2 = primitive.Clone("h2")
        DisplayTH1D(primitive, ax, rebin=10, color='black', scaler=integral1/integral, label = r"2024 ($\times${:.1f})".format(integral/integral1), title="")

DisplayTH1D(h1, ax, rebin=10, color='#851CCA', xlabel="Energy [keV]", ylabel="Counts / keV", ylog=True, label="2025", title="")

ax.set_xlabel("Energy [keV]", fontsize=20)
ax.set_ylabel("Counts / 0.5keV", fontsize=20)
ax.legend()
ax.set_xlim(650, 6500)
ax.set_ylim(1, 200000)
# fig.savefig("32Ar_2024_2025_comparison.png", dpi=300)
# fig.show()




fig, ax = plt.subplots(figsize=(8, 4))
plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.25)
file = ROOT.TFile("/mnt/hgfs/shared-2/2025_DATA/DETECTOR_DATA/ANALYSED/32Ar_analysed_2025.root", "READ")
h=0
for det in range(1, 9):
    for strip in range(1, 6):
        if h == 0:
            h = file.Get(f"RandomCorrection/D{det}.{strip}/Write/H_SiPM_Time_D{det}.{strip}_3")
            print(h.GetEntries())
        else:
            h.Add(file.Get(f"RandomCorrection/D{det}.{strip}/Write/H_SiPM_Time_D{det}.{strip}_3"), 1)
intensity = h.Integral(h.GetXaxis().FindBin(-200), h.GetXaxis().FindBin(200))

file1 = ROOT.TFile("/mnt/hgfs/shared-2/2024_DATA/DETECTOR_DATA/ANALYSED/32Ar_analysed_2024.root", "READ")
h1=0
for det in range(1, 9):
    for strip in range(1, 6):
        if h1 == 0:
            h1 = file1.Get(f"RandomCorrection/D{det}.{strip}/Write/H_SiPM_Time_D{det}.{strip}_3")
        else:
            h1.Add(file1.Get(f"RandomCorrection/D{det}.{strip}/Write/H_SiPM_Time_D{det}.{strip}_3"), 1)
intensity1 = h1.Integral(h1.GetXaxis().FindBin(-200), h1.GetXaxis().FindBin(200))
DisplayTH1D(h1, ax, color='black', xlabel="SiPM Time [ns]", ylabel="Counts / 0.5ns", ylog=True, title="", label =r"2024 ($\times${:.1f})".format(integral/integral1), scaler=integral1/integral, lw=1.5)
DisplayTH1D(h, ax, color='#851CCA', xlabel="SiPM Time [ns]", ylabel="Counts / 0.5ns", ylog=True, title="", label ="2025", lw=1.5)


ax.legend()
ax.set_xlabel(r"$\Delta t _{\beta p}$ [ns]", fontsize=20)
ax.set_ylabel("Counts / 2ns", fontsize=20)
ax.set_xlim(-200, 200)
ax.set_ylim(1, 1100000)
fig.savefig("betap_2024_2025_comparison.png", dpi=300)
plt.show()