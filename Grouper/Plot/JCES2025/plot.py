from root2mpl import *
from ROOT import *
import numpy as np
dir = "/mnt/hgfs/shared-2/2025_DATA/DETECTOR_DATA/CALIBRATED/"

plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=-0.2)
from matplotlib.ticker import LogFormatter
from matplotlib.ticker import FuncFormatter
def scientific_notation(x, pos):
    
    if x > 0:
        if x / 10**np.floor(np.log10(x)) == 1:
            return r'$10^{{{:d}}}$'.format(int(np.floor(np.log10(x))))
        return r'${:.0f}.10^{{{:d}}}$'.format(x / 10**np.floor(np.log10(x)), int(np.floor(np.log10(x))))
    else:
        return r'$0$'

formatter = FuncFormatter(scientific_notation)

file = TFile(dir + "Calibrated_2025test.root", "READ")


ylim = 50e3

# fig, ax = plt.subplots(figsize = (sizex, sizey))
# DisplayCanvas(file.Get("D4.1"), ax=ax, ylog = None, xlim = (3200, 3450), ylim = (1, ylim), lw = 0.5)
# plt.savefig("D4.1.png", dpi=300)

# fig1, ax1 = plt.subplots(figsize = (sizex, sizey))
# DisplayCanvas(file.Get("D5.1"), ax=ax1, ylog = None, xlim = (3200, 3450), ylim = (1, ylim))
# plt.savefig("D5.1.png", dpi=300)

# fig2, ax2 = plt.subplots(figsize = (sizex, sizey))
# DisplayCanvas(file.Get("D4.1_Calibration1"), ax=ax2, ylog = None, xlim = (3200, 3450), ylim = (1, ylim), title = "D4.1")
# plt.savefig("D4.1_c.png", dpi=300)

# fig3, ax3 = plt.subplots(figsize = (sizex, sizey))
# DisplayCanvas(file.Get("D5.1_Calibration1"), ax=ax3, ylog = None, xlim = (3200, 3450), ylim = (1, ylim), title="D5.1")
# plt.savefig("D5.1_c.png", dpi=300)

# fig4, ax4 = plt.subplots(figsize = (sizex, sizey))
# DisplayCanvas(file.Get("D4.1"), ax=ax4, ylog = None, xlim = (0, 6500), ylim = (1, ylim), lw = 0.5)
# plt.savefig("D4.1_0_6500.png", dpi=300)

# fig5, ax5 = plt.subplots(figsize = (sizex, sizey))
# DisplayCanvas(file.Get("D5.1"), ax=ax5, ylog = None, xlim = (0, 6500), ylim = (1, ylim))
# plt.savefig("D5.1_0_6500.png", dpi=300)

# fig6, ax6 = plt.subplots(figsize = (sizex, sizey))
# DisplayCanvas(file.Get("D4.1_Calibration1"), ax=ax6, ylog = None, xlim = (0, 6500), ylim = (1, ylim), title = "D4.1")
# plt.savefig("D4.1_c_0_6500.png", dpi=300)

# fig7, ax7 = plt.subplots(figsize = (sizex, sizey))
# DisplayCanvas(file.Get("D5.1_Calibration1"), ax=ax7, ylog = None, xlim = (0, 6500), ylim = (1, ylim), title="D5.1")
# plt.savefig("D5.1_c_0_6500.png", dpi=300)


c = file.Get("32Ar_SUM_Down")
# loop on primitive
for i in range(c.GetListOfPrimitives().GetEntries()):
    obj = c.GetListOfPrimitives().At(i)
    if obj.GetName() == "32Ar_Exp_All_Down":
        h = obj

fig0, ax0 = plt.subplots(figsize = (14, 4))
h.Rebin(10)
DisplayTH1D(h, ax=ax0, ylog = True, xlim = (0, 6500), ylim = (1, 3*ylim), lw = 0.5, legend = None, color="black", title="")
plt.savefig("32Ar_SUM_Down_log.png", dpi=300)


sizex = 4
sizey = 4
file = TFile(dir + "../ANALYSED/32Ar_analysed.root", "READ")
fig, ax  = plt.subplots(figsize = (sizex, sizey))
DisplayCanvas(file.Get("RandomCorrection/D5.1/Corrected_D5.1_6"), ax=ax, xlim = (3310, 3390), lw = 0.5, legend = None, title="",rebin =2, ylabel="Counts / 0.2 keV")
plt.savefig("D5.1_32Ar_peak.png", dpi=300)

# same for D1.1
fig1 , ax1  = plt.subplots(figsize = (sizex, sizey))
DisplayCanvas(file.Get("RandomCorrection/D1.1/Corrected_D1.1_6"), ax=ax1, xlim = (3310, 3390), lw = 0.5, legend = None, title="",rebin =2, ylabel="Counts / 0.2 keV")
plt.savefig("D1.1_32Ar_peak.png", dpi=300)
# plt.show()




# fig1, ax1 = plt.subplots(figsize = (sizex, sizey))
# DisplayCanvas(file.Get("D5.1"), ax=ax1, ylog = True, xlim = (0, 6500), ylim = (1, ylim), legend = None)
# # add an inset with a zoom on the peak
# axins1 = ax1.inset_axes([0.62, 0.6, 0.35, 0.35])
# axins1.yaxis.set_major_formatter(formatter)
# DisplayCanvas(file.Get("D5.1"), ax=axins1, ylog = None, xlim = (3300, 3400), ylim = (1, 12e3), legend = None, title=None, ylabel = "", xlabel = "")
# ax1.indicate_inset_zoom(axins1)
# plt.savefig("D5.1_log.svg")