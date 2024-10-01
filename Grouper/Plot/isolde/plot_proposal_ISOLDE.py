from root2mpl import *
from ROOT import *
import numpy as np
dir = "../../../../../../../mnt/hgfs/shared-2/2024_DATA/DETECTOR_DATA/CALIBRATED/"

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

file = TFile(dir + "Calibrated.root", "READ")

sizex = 8
sizey = 5
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



fig, ax  = plt.subplots(figsize = (sizex, sizey))
DisplayCanvas(file.Get("D4.1"), ax=ax, ylog = True, xlim = (0, 6500), ylim = (1, ylim), lw = 0.5, legend = None)
# add an inset with a zoom on the peak
axins = ax.inset_axes([0.62, 0.6, 0.35, 0.35])
axins.yaxis.set_major_formatter(formatter)
DisplayCanvas(file.Get("D4.1"), ax=axins, ylog = None, xlim = (3300, 3400), ylim = (1, 12e3), lw = 0.5, legend = None, title=None, ylabel = "", xlabel = "")
ax.indicate_inset_zoom(axins)
plt.savefig("D4.1_log.svg")

fig1, ax1 = plt.subplots(figsize = (sizex, sizey))
DisplayCanvas(file.Get("D5.1"), ax=ax1, ylog = True, xlim = (0, 6500), ylim = (1, ylim), legend = None)
# add an inset with a zoom on the peak
axins1 = ax1.inset_axes([0.62, 0.6, 0.35, 0.35])
axins1.yaxis.set_major_formatter(formatter)
DisplayCanvas(file.Get("D5.1"), ax=axins1, ylog = None, xlim = (3300, 3400), ylim = (1, 12e3), legend = None, title=None, ylabel = "", xlabel = "")
ax1.indicate_inset_zoom(axins1)
plt.savefig("D5.1_log.svg")