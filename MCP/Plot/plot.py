from root2mpl import *
from ROOT import *
from mpl_toolkits.mplot3d import axes3d
from matplotlib.colors import LinearSegmentedColormap
from mpl_toolkits.axes_grid1 import make_axes_locatable, axes_size
import matplotlib.gridspec as gridspec


#fucntion to diplay TF2 in matpltllib
def DisplayTF2(f, ax, levels=[0.5]):
    x = np.linspace(f.GetXmin(), f.GetXmax(), 1000)
    y = np.linspace(f.GetYmin(), f.GetYmax(), 1000)
    X, Y = np.meshgrid(x, y)
    Z = np.array([[f.Eval(xi, yi) for xi in x] for yi in y])

    zmax = max(Z.flatten())
    
    im = ax.contour(X, Y, Z, colors="red", levels=[i*zmax for i in levels], linewidths=1.5, linestyles="dashed")

color=["#08348C", "#89BAD9", "#021526", "#376B8B"]
##################### Coincidence Time #####################

# file = TFile("/mnt/hgfs/shared-2/2025_DATA/MCP_DATA/GROUPED/run_001_StableBeamScan_grouped.root", "READ")
# fig, ax = plt.subplots(figsize=(9, 3), constrained_layout=True)

# c = file.Get("TIME_CORNERS_FULL")
# counter=0
# ax.set_ylim(1, 3e6)
# for h in c.GetListOfPrimitives():
#     if h.ClassName() == "TH1D":
#         DisplayTH1D(h, ax, color=color[counter], xlabel=r"$\Delta$T (ns)", ylabel="Counts / 4 ns", label = f"SA{counter+1}", ylog=True, title="", labelsize=15, rebin=2)
#         counter+=1

# plt.xlim(-300, 1001)
# plt.xlabel(r"$\Delta$T (ns)", fontsize=20)
# plt.ylabel("Counts / 4 ns", fontsize=20)
# plt.legend(loc="upper right", ncol=2)
# plt.savefig("Corner_Coincidence_Time.pdf")
# #Adding window coinc
# ax.axvline(x=-100, color="red", linestyle="--", linewidth=1., label="Coincidence window")
# ax.axvline(x=-190, color="red", linestyle="--", linewidth=1.)
# handles, labels = ax.get_legend_handles_labels()
# handles[2], handles[4] = handles[4], handles[2]
# handles[3], handles[4] = handles[4], handles[3]
# labels[2], labels[4] = labels[4], labels[2]
# labels[3], labels[4] = labels[4], labels[3]
# plt.legend(handles, labels, loc="upper right", ncol=2, columnspacing=-2.9)
# plt.savefig("Corner_Coincidence_Time_with_window_2.pdf")

##################### MCP PICTURES #####################
#personnalized cmap Blues + (only for z=-1 color is red)
cmap = plt.cm.Blues
# cmap.set_under('red')

# #RAW
file = TFile("/mnt/hgfs/shared-2/2025_DATA/MCP_DATA/GROUPED/run_001_StableBeamScan_grouped.root", "READ")
# fig, ax = plt.subplots(figsize = (10, 8), constrained_layout=True)
# h2 = file.Get("H_RAW_2D")
# h2_uncleaned = file.Get("H_RAW_2D_uncleaned")

# maxi=10
    
# # DisplayTH2D(h2, ax, color="Blues", vmax=maxi, xlabel=r"$X$ (a.u)", ylabel=r"$Y$ (a.u)", title="", zlog=False, visible=False)


# plt.savefig("2D_Scan_RAW_wo_fit.pdf")

# xl = []
# yl = []
# for binx in range(1, h2.GetNbinsX()+1):
#     for biny in range(1, h2.GetNbinsY()+1):
#         bin = h2.GetBin(binx, biny)
#         if (h2.GetBinContent(bin) == 0 and h2_uncleaned.GetBinContent(bin) != 0):
#             x = h2.GetXaxis().GetBinCenter(binx)
#             y = h2.GetYaxis().GetBinCenter(biny)
#             xl.append(x)
#             yl.append(y)

#             h2.SetBinContent(bin, -1)

# DisplayTH2D(h2, ax, color=cmap, vmax=maxi, xlabel=r"$X_d$ (a.u)", ylabel=r"$Y_d$ (a.u)", title="", zlog=False)
# ax.set_xlabel(r"$X_d$ (a.u)", fontsize=25)
# ax.set_ylabel(r"$Y_d$ (a.u)", fontsize=25)
# ax.set_xlim(-0.6, 0.6)
# ax.set_ylim(-0.6, 0.6)
# print("bin size = ", h2.GetXaxis().GetBinWidth(1), "x", h2.GetYaxis().GetBinWidth(1))
# cbar = plt.colorbar(ax.images[0], ax=ax, ticks=np.arange(0, maxi+1, 2), label="Counts per bin")
# cbar.ax.set_yticklabels([str(int(i)) if i != cbar.get_ticks()[-1] else r'$>$' + str(int(i)) for i in cbar.get_ticks()])
# Area in points^2
# plt.savefig("Reviewed/2D_Scan_RAW_wo_fit_w_bkg.pdf")


file1 = TFile("../Calibration_2025 copy.root", "READ")

# # LOG
# c = file1.Get("c")
# h2 = None
# for h in c.GetListOfPrimitives():
#     if h.ClassName() == "TF2":
#         f = h

# fig, ax = plt.subplots(figsize = (10, 8), constrained_layout=True)
# h2 = file.Get("h_Image")
# maxi=14
# DisplayTH2D(h2, ax, color="Blues", vmax=maxi, xlabel=r"$X$ (a.u)", ylabel=r"$Y$ (a.u)", title="")
# ax.set_xlim(-0.4, 0.4)
# ax.set_ylim(-0.4, 0.4)
# ax.set_xlabel(r"$X$ (a.u)", fontsize=25)
# ax.set_ylabel(r"$Y$ (a.u)", fontsize=25)

# cbar = plt.colorbar(ax.images[0], ax=ax, ticks=np.arange(0, maxi+2, 2))
# cbar.ax.set_yticklabels([str(int(i)) if i != cbar.get_ticks()[-1] else r'$>$' + str(int(i)) for i in cbar.get_ticks()])
# cbar.set_label("Counts per bin", fontsize=25)
# plt.savefig("2D_Scan_wo_fit.pdf")

# DisplayTF2(f, ax)
# plt.savefig("Reviewed/2D_Scan_w_fit.pdf")

### CALIBRATED

# file = TFile("../Calibration_2025.root", "READ")

# fig, ax = plt.subplots(figsize = (10, 8), constrained_layout=True)
# c = file.Get("c_reconstruction_interpolated")
# h2 = None
# for h in c.GetListOfPrimitives():
#     if h.ClassName() == "TH2D":
#         h2 = h
# h2.Rebin2D(2,2)
# maxi = 16
# DisplayTH2D(h2, ax, color="Blues", vmax=maxi, xlabel=r"$x$ (mm)", ylabel=r"$y$ (mm)", title="")
# ax.set_xlabel(r"$x$ (mm)", fontsize=25)
# ax.set_ylabel(r"$y$ (mm)", fontsize=25)
# ax.set_xlim(-8, 8)
# ax.set_ylim(-8, 8)
# ax.set_yticks(ax.get_xticks())
# cbar = plt.colorbar(ax.images[0], ax=ax, ticks=np.arange(0, maxi+1, 2))
# cbar.ax.set_yticklabels([str(int(i)) if i != cbar.get_ticks()[-1] else r'$>$' + str(int(i)) for i in cbar.get_ticks()])
# cbar.set_label("Counts per bin", fontsize=25)
# plt.savefig("Reviewed/2D_Scan_Calibrated_wo_fit.pdf")

# ### CALIBRATED SECOND STEP
# ## 0T
# file = TFile("/mnt/hgfs/shared-2/2025_DATA/MCP_DATA/CALIBRATED/run_001_StableBeamScan_secondstep_calibrated.root", "READ")

# fig, ax = plt.subplots(figsize = (10, 8), constrained_layout=True)
# h2 = file.Get("H_reconstruction_interpolated_second")
# h2.Rebin2D(2,2)
# maxi = 16
# DisplayTH2D(h2, ax, color="Blues", vmax=maxi, xlabel=r"$x$ (mm)", ylabel=r"$y$ (mm)", title="")
# ax.set_xlabel(r"$x$ (mm)", fontsize=25)
# ax.set_ylabel(r"$y$ (mm)", fontsize=25)
# ax.set_xlim(-4, 4)
# ax.set_ylim(-4, 4)
# ax.set_yticks(ax.get_xticks())
# cbar = plt.colorbar(ax.images[0], ax=ax, ticks=np.arange(0, maxi+1, 2))
# cbar.ax.set_yticklabels([str(int(i)) if i != cbar.get_ticks()[-1] else r'$>$' + str(int(i)) for i in cbar.get_ticks()])
# cbar.set_label("Counts per bin", fontsize=25)
# # plt.savefig("2D_Scan_Calibrated_Second_wo_fit.pdf")

# # adding fillbetween on the plot 
# color = "grey"
# alpha = 0.5
# # for i in range(-3, 3):
# #     ax.fill_between([-4, 4], 2*(i+1) - 0.7, 2*(i) + 0.7, color=color, alpha=alpha, edgecolor=None, linewidth=0)
# #     for j in range(-3, 3):
# #         ax.fill_between([2*(i+1) - 0.7, 2*(i) + 0.7], -0.7 + 2*j, 0.7+2*j, color=color, alpha=alpha, edgecolor=None, linewidth=0)
# # plt.savefig("2D_Scan_Calibrated_Second_w_grid.pdf")

# # adding squares on the plot
# color = "red"
# alpha = 1.0
# length = 1.4
# for i in range(-3, 3):
#     for j in range(-3, 3):
#         square = plt.Rectangle((2*i - length/2, 2*j - length/2), length, length, fill=False, color=color, alpha=alpha, linewidth=2, linestyle="--")
#         ax.add_patch(square)
# plt.savefig("Reviewed/2D_Scan_Calibrated_Second_w_gridlines.pdf")

### 1T
# file = TFile("../MSecondStepCalibrated_006.root", "READ")

# fig, ax = plt.subplots(figsize = (10, 8), constrained_layout=True)
# h2 = file.Get("H_reconstruction_interpolated_second_measurement")
# h2.Rebin2D(2,2)
# maxi = 20
# DisplayTH2D(h2, ax, color="Blues", vmax=maxi, xlabel=r"$x$ (mm)", ylabel=r"$y$ (mm)", title="", log=True)
# ax.set_xlim(-4, 4)
# ax.set_ylim(-4, 4)
# cbar = plt.colorbar(ax.images[0], ax=ax, ticks=np.arange(0, maxi+1, 2))
# cbar.ax.set_yticklabels([str(int(i)) if i != cbar.get_ticks()[-1] else r'$>$' + str(int(i)) for i in cbar.get_ticks()])
# plt.savefig("2D_Scan_Calibrated_Second_wo_fit_1T.pdf")

# # adding fillbetween on the plot 
# color = "grey"
# alpha = 0.5
# # for i in range(-3, 3):
# #     ax.fill_between([-4, 4], 2*(i+1) - 0.7, 2*(i) + 0.7, color=color, alpha=alpha, edgecolor=None, linewidth=0)
# #     for j in range(-3, 3):
# #         ax.fill_between([2*(i+1) - 0.7, 2*(i) + 0.7], -0.7 + 2*j, 0.7+2*j, color=color, alpha=alpha, edgecolor=None, linewidth=0)
# plt.savefig("2D_Scan_Calibrated_Second_w_grid_1T.pdf")

# # adding squares on the plot
# color = "red"
# alpha = 1.0
# length = 1.4
# for i in range(-3, 3):
#     for j in range(-3, 3):
#         square = plt.Rectangle((2*i - length/2, 2*j - length/2), length, length, fill=False, color=color, alpha=alpha, linewidth=2, linestyle="--")
#         ax.add_patch(square)
# plt.savefig("2D_Scan_Calibrated_Second_w_gridlines_1T.pdf")

# ### RESIDUALS

# fig, axx = plt.subplots(figsize=(9, 4), constrained_layout=True, nrows=1, ncols=3, sharex=True, sharey=True)

# for i, l in enumerate(["A", "B", "C"]):
#     ax = axx[i]
#     x = []
#     y = []
#     file = TFile("../Calibration_2025_"+l+".root", "READ")
#     c = file.Get("c_residuals")
#     g = None
#     for h in c.GetListOfPrimitives():
#         if h.ClassName() == "TGraphErrors":
#             g = h

#     for j in range(g.GetN()):
#         x.append(g.GetX()[j])
#         y.append(g.GetY()[j]*g.GetX()[j])

#     mean = np.std([y[j] for j in range(len(x)) if x[j] < 4 and (x[j] < 3.3 or x[j] > 3.6)])
#     ax.scatter(x, y, label=f"Cell {mean}", color=color[i], s=15)
#     ax.fill_between(x, -mean, mean, color=color[i], alpha=0.2)
#     ax.legend()
# plt.show()


# """
# """
# #RESULT ################ 2025
#     #### 2D

# cmap = plt.get_cmap('Blues')

# fig, ax = plt.subplots(figsize = (10, 8), constrained_layout=True)
# file = TFile("run_079_MCP_32Ar_Heinz14kV_grouped.root_fitres_calibrated.root", "READ")
# c = file.Get("MeassurementFitted_2D_View")
# for h in c.GetListOfPrimitives():
#     if h.ClassName() == "TF2":
#         f = h

# c = file.Get("MeassurementFitted_2D_View")
# h2 = None
# for h in c.GetListOfPrimitives():
#     if h.ClassName() == "TH2D":
#         h2 = h

# maxi = 20
# rebin = 1
# lim = 4
# DisplayTH2D(h2, ax, color=cmap, xlabel=r"$x$ (mm)", ylabel=r"$y$ (mm)", title="", zlog=False, vmax = maxi, rebinx= rebin, rebiny=rebin)
# ax.set_xlim(-lim, lim)
# ax.set_ylim(-lim, lim)
# ax.set_yticks(ax.get_xticks())
# ax.set_yticklabels([int(i) for i in ax.get_yticks()])
# ax.set_xlabel(r"$x$ (mm)", fontsize=25)
# ax.set_ylabel(r"$y$ (mm)", fontsize=25)

# cbar = plt.colorbar(ax.images[0], ax=ax, ticks=np.arange(0, maxi+1, maxi/10))
# cbar.ax.set_yticklabels([str(int(i)) if i != cbar.get_ticks()[-1] else r'$>$' + str(int(i)) for i in cbar.get_ticks()])
# cbar.set_label("Counts per bin", fontsize=25)<
# # plt.savefig("2D_Beam_Calibrated_wo_fit.pdf")
# DisplayTF2(f, ax, [0.15, 0.25, 0.5, 0.75, 0.9])
# ax.set_xlim(-lim, lim)
# ax.set_ylim(-lim, lim)
# ax.set_xlabel(r"$x$ (mm)", fontsize=30)
# ax.set_ylabel(r"$y$ (mm)", fontsize=30)
# plt.savefig("Reviewed/2D_Beam_Calibrated_with_fit_2.pdf")
# plt.show()

# #RESULT ################ 2024
#     #### 2D

# cmap = plt.get_cmap('Blues')

# fig, ax = plt.subplots(figsize = (10, 8), constrained_layout=True)
# file = TFile("../Calibration_2024_Saved.root", "READ")
# c = file.Get("MeassurementFitted_2D_View")
# for h in c.GetListOfPrimitives():
#     if h.ClassName() == "TF2":
#         f = h

# c = file.Get("Measurement_2D_View")
# h2 = None
# for h in c.GetListOfPrimitives():
#     if h.ClassName() == "TH2D":
#         h2 = h

# maxi = 20
# DisplayTH2D(h2, ax, color=cmap, xlabel=r"$x$ (mm)", ylabel=r"$y$ (mm)", title="", zlog=False, vmax = maxi)
# ax.set_xlim(-2, 2)
# ax.set_ylim(-2, 2)

# cbar = plt.colorbar(ax.images[0], ax=ax, ticks=np.arange(0, maxi+1, maxi/10))
# cbar.ax.set_yticklabels([str(int(i)) if i != cbar.get_ticks()[-1] else r'$>$' + str(int(i)) for i in cbar.get_ticks()])
# plt.savefig("2D_Beam_Calibrated_wo_fit_2024.pdf")
# DisplayTF2(f, ax, [0.1, 0.25, 0.5, 0.75, 0.9])
# ax.set_xlim(-2, 2)
# ax.set_ylim(-2, 2)
# plt.savefig("2D_Beam_Calibrated_with_fit_2024_2.pdf")
# """

# """
# ###################### RESOLUTION #####################
# fig, ax = plt.subplots(figsize = (9, 3), constrained_layout=False, nrows=1, ncols=3, sharex=True, sharey=True)
# plt.subplots_adjust(wspace=0.15, bottom=0.2, left = 0.05)#, right=1.05)
# file = TFile("../Calibration_2025.root", "READ")
# c = file.Get("c_Resolution_XY")
# hxy, hx, hy = None, None, None
# cmap = plt.get_cmap('Blues')

# for p in c.GetListOfPrimitives():
#     for h in p.GetListOfPrimitives():
            
#         if h.GetName() == "H_Resolution_X_interpolated":
#             hx = h
#         if h.GetName() == "H_Resolution_Y_interpolated":
#             hy = h
#         if h.GetName() == "H_Resolution_XY_interpolated":
#             hxy = h

# xlim = 7
# max = 0.6
# rebin = 100
# DisplayTH2D(hx, ax[0], color=cmap, xlabel=r"$x$ (mm)", ylabel=r"$y$ (mm)", xlim=(-xlim, xlim), ylim=(-xlim, xlim), vmax=max, rebinx = rebin, rebiny = rebin)
# ax[0].set_title(r"$\delta{x}$")
# DisplayTH2D(hy, ax[1], color=cmap, xlabel=r"$x$ (mm)", ylabel="", xlim=(-xlim, xlim), ylim=(-xlim, xlim), vmax=max, rebinx = rebin, rebiny = rebin)
# ax[1].set_title(r"$\delta{y}$")
# DisplayTH2D(hxy, ax[2], color=cmap, xlabel=r"$x$ (mm)", ylabel="", xlim=(-xlim, xlim), ylim=(-xlim, xlim), vmax=max, rebinx = rebin, rebiny = rebin)
# ax[2].set_title(r"$\delta{xy}$")
# ax[0].set_xticks(np.arange(-xlim, xlim + xlim/4, 2*xlim/4))
# ax[0].set_yticks(np.arange(-xlim, xlim + xlim/4, 2*xlim/4))
# ax[1].set_xticks(np.arange(-xlim, xlim + xlim/4, 2*xlim/4))
# ax[1].set_yticks(np.arange(-xlim, xlim + xlim/4, 2*xlim/4))
# ax[2].set_xticks(np.arange(-xlim, xlim + xlim/4, 2*xlim/4))
# ax[2].set_yticks(np.arange(-xlim, xlim + xlim/4, 2*xlim/4))

# ## add subplot for colorbar
# divider = make_axes_locatable(ax[2])
# cax = divider.append_axes("right", size="5%", pad=0.05)
# cbar = plt.colorbar(ax[2].images[0], cax=cax, ticks=np.arange(0, max+max/6, max/6))
# cbar.ax.set_yticklabels([str(int(i)) if i != cbar.get_ticks()[-1] else r'$>$' + str(int(i)) for i in cbar.get_ticks()])
# #ticks labels of the color bar
# labels = [str(round(i, 1)) for i in np.arange(0, max+max/6, max/6)]
# labels[-1] = r'$>$ 0.6'
# cax.set_yticklabels(labels)
# cax.set_ylabel(r"$\sigma$ (mm)", fontsize=15)

# plt.savefig("Resolution_2025.pdf")


# only resolution x and y
# fig, axx = plt.subplots(figsize = (9, 4), constrained_layout=True, nrows=1, ncols=2, sharex=True, sharey=True)
# # plt.subplots_adjust(wspace=0.15, bottom=0.25, left = 0.1, right=0.9)
# file = TFile("/mnt/hgfs/shared-2/2025_DATA/MCP_DATA/CALIBRATED/run_001_StableBeamScan_secondstep_calibrated.root", "READ")
# c = file.Get("c_Resolution_XY_proj")
# hxy, hx, hy = None, None, None
# cmap = plt.get_cmap('Blues')

# for p in c.GetListOfPrimitives():
#     for h in p.GetListOfPrimitives():
            
#         if h.GetName() == "H_Resolution_X_interpolated_proj":
#             hx = h
#         if h.GetName() == "H_Resolution_Y_interpolated_proj":
#             hy = h
#         if h.GetName() == "H_Resolution_XY_interpolated_proj":
#             hxy = h

# xlim = 7
# max = 0.5
# rebin = 1
# DisplayTH2D(hx, axx[0], color=cmap, xlabel=r"$x$ (mm)", ylabel=r"$y$ (mm)", xlim=(-xlim, xlim), ylim=(-xlim, xlim), vmax=max, rebinx = rebin, rebiny = rebin)
# axx[0].set_title(r"$\delta_{x}$", fontsize = 24)
# axx[0].set_xlabel(r"$x$ (mm)", fontsize=25)
# axx[0].set_ylabel(r"$y$ (mm)", fontsize=25)
# DisplayTH2D(hy, axx[1], color=cmap, xlabel=r"$x$ (mm)", ylabel="", xlim=(-xlim, xlim), ylim=(-xlim, xlim), vmax=max, rebinx = rebin, rebiny = rebin)
# axx[1].set_title(r"$\delta_{y}$", fontsize = 24)
# axx[1].set_xlabel(r"$x$ (mm)", fontsize=25)
# axx[0].set_xticks(np.arange(-xlim, xlim + xlim/4, 2*xlim/5))
# axx[0].set_yticks(np.arange(-xlim, xlim + xlim/4, 2*xlim/4))
# axx[1].set_xticks(np.arange(-xlim, xlim + xlim/4, 2*xlim/4))
# axx[1].set_yticks(np.arange(-xlim, xlim + xlim/4, 2*xlim/4))

# ## add subplot for colorbar
# cbar = plt.colorbar(axx[1].images[0], ax=axx[1], ticks=np.arange(0, max+max/(max*10), max/(max*10)), label = r"$\delta$ (mm)")
# cbar.ax.set_yticklabels([f"{i:.1f}" if i != cbar.get_ticks()[-1] else r'$>$' + f"{i:.1f}" for i in cbar.get_ticks()])
# plt.savefig("Resolution_xy_secondstep_2025.pdf")
# plt.show()

"""
file = TFile("../Calibration_2025.root", "READ")
c = file.Get("c_Resolution_XY")
hxy, hx, hy = None, None, None
cmap = plt.get_cmap('Blues')

for p in c.GetListOfPrimitives():
    for h in p.GetListOfPrimitives():
            
        if h.GetName() == "H_Resolution_X_interpolated":
            hx = h
        if h.GetName() == "H_Resolution_Y_interpolated":
            hy = h
        if h.GetName() == "H_Resolution_XY_interpolated":
            hxy = h

xlim = 7
max = 0.6
rebin = 1

fig, ax = plt.subplots(figsize = (8, 6), constrained_layout=True)
DisplayTH2D(hxy, ax, color=cmap, xlabel=r"$x$ (mm)", ylabel=r"$y$ (mm)", title ="", xlim=(-xlim, xlim), ylim=(-xlim, xlim), vmax=max, rebinx = rebin, rebiny = rebin)
ax.set_xticks(np.arange(-xlim, xlim + xlim/4, 2*xlim/4))
ax.set_yticks(np.arange(-xlim, xlim + xlim/4, 2*xlim/4))
c = plt.colorbar(ax.images[0], ax=ax, ticks=np.arange(0, max+max/6, max/6), label = r"$\delta_{xy}$ (mm)")
#ticks labels of the color bar
labels = [str(round(i, 1)) for i in np.arange(0, max+max/6, max/6)]
labels[-1] = r'$>$ 0.6'
c.set_ticklabels(labels)
plt.savefig("Resolution_2025.pdf")
plt.show()

#################### 1T - 4T Scan #####################
max1T = 15
max4T = 15
rebin = 2
fig, ax = plt.subplots(figsize=(9, 4), constrained_layout=True, nrows=1, ncols=2, sharex=True, sharey=True)

f1T = TFile("/mnt/hgfs/shared-2/2025_DATA/MCP_DATA/CALIBRATED/run_006_MCP_Scan_calibrated.root", "READ")
c1T = f1T.Get("Measurement_2D_View")
for h in c1T.GetListOfPrimitives():
    if h.ClassName() == "TH2D":
        h1T = h

DisplayTH2D(h1T, ax[0], color="Blues", xlabel=r"$x$ (mm)", ylabel=r"$y$ (mm)", title="", xlim=(-8, 8), ylim=(-8, 8), vmax = max1T, rebinx = rebin, rebiny = rebin)
ax[0].text(0.85, 0.85, "1T", transform=ax[0].transAxes, ha='center', fontsize=30)

# f2T = TFile("/mnt/hgfs/shared-2/2025_DATA/MCP_DATA/CALIBRATED/run_034_MCP_32Ar_BeamScan_2T_calibrated.root", "READ")
# c2T = f2T.Get("Measurement_2D_View")
# for h in c2T.GetListOfPrimitives():
#     if h.ClassName() == "TH2D":
#         h2T = h

# DisplayTH2D(h2T, ax[0, 1], color="Blues", xlabel=r"", ylabel=r"", title="", xlim=(-8, 8), ylim=(-8, 8), vmax = max2T)

f3T = TFile("../run_035_MCP_32Ar_BeamScan_3T_calibrated_secondstep.root", "READ")
h3T = f3T.Get("H_reconstruction_interpolated_second_measurement")

DisplayTH2D(h3T, ax[1], color="Blues", xlabel=r"$x$ (mm)", ylabel=r"$y$ (mm)", title="", xlim=(-8, 8), ylim=(-8, 8), vmax = max1T)
ax[1].text(0.85, 0.85, "3T", transform=ax[1].transAxes, ha='center', fontsize=30)
# f4T = TFile("/mnt/hgfs/shared-2/2025_DATA/MCP_DATA/CALIBRATED/run_038_MCP_32Ar_BeamScan_4T_calibrated.root", "READ")
# c4T = f4T.Get("Measurement_2D_View")
# for h in c4T.GetListOfPrimitives():
#     if h.ClassName() == "TH2D":
#         h4T = h

# DisplayTH2D(h4T, ax[1], color="Blues", xlabel=r"$x$ (mm)", ylabel=r"", title="", xlim=(-8, 8), ylim=(-8, 8), vmax = max4T, rebinx = rebin, rebiny = rebin)
# ax[1].text(0.85, 0.85, "4T", transform=ax[1].transAxes, ha='center', fontsize=30)


ax[1].set_xticks(np.arange(-8, 9, 4))
ax[1].set_yticks(np.arange(-8, 9, 4))
ax[0].set_yticks(np.arange(-8, 9, 4))
cbar = plt.colorbar(ax[1].images[0], ax=ax[1], ticks=np.arange(0, max4T+max4T/5, max4T/5), label = "Counts")
cbar.ax.set_yticklabels([str(int(i)) if i != cbar.get_ticks()[-1] else r'$>$' + str(int(i)) for i in cbar.get_ticks()])
# plt.savefig("Fields_2D_Scan.pdf")
plt.show()
"""
#################### All Field Scan #####################
# max1T = 15
# max2T = 150
# max3T = 80
# max4T = 15
# rebin = 2
# fig, ax = plt.subplots(figsize=(9, 7), nrows=2, ncols=2, sharex=True, sharey=True)
# # // setiing spoace between plots
# fig.subplots_adjust(hspace=0.1, wspace=0.1, left=0.1, right=0.93, top=0.95, bottom=0.1)

# ##
# f1T = TFile("/mnt/hgfs/shared-2/2025_DATA/MCP_DATA/CALIBRATED/run_006_MCP_Scan_calibrated.root", "READ")
# c1T = f1T.Get("Measurement_2D_View")
# for h in c1T.GetListOfPrimitives():
#     if h.ClassName() == "TH2D":
#         h1T = h

# DisplayTH2D(h1T, ax[0, 0], color="Blues", xlabel=r"", ylabel=r"$y$ (mm)", title="", xlim=(-8, 8), ylim=(-8, 8), vmax = max1T, rebinx = rebin, rebiny = rebin)
# ax[0, 0].set_ylabel(r"$y$ (mm)", fontsize=25)
# ax[0,0].text(0.85, 0.85, "1T", transform=ax[0, 0].transAxes, ha='center', fontsize=30)
# colorbar = plt.colorbar(ax[0, 0].images[0], ax=ax[0, 0], ticks=np.arange(0, max1T+max1T/5, max1T/5), shrink=1.0)
# colorbar.ax.set_yticklabels([str(int(i)) if i != colorbar.get_ticks()[-1] else r'$>$' + str(int(i)) for i in colorbar.get_ticks()])
# # colorbar.set_label("Counts per bin", fontsize=15)
# ##
# # f2T = TFile("/mnt/hgfs/shared-2/2025_DATA/MCP_DATA/CALIBRATED/run_034_MCP_32Ar_BeamScan_2T_calibrated.root", "READ")
# # c2T = f2T.Get("Measurement_2D_View")
# # for h in c2T.GetListOfPrimitives():
# #     if h.ClassName() == "TH2D":
# #         h2T = h
# f2T = TFile("/mnt/hgfs/shared-2/2025_DATA/MCP_DATA/run_034_MCP_32Ar_BeamScan_2T_calibrated_secondstep.root", "READ")
# # f2T = TFile("../src/RateCorrectedPlot.root", "READ")
# h2T = f2T.Get("H_reconstruction_interpolated_second_measurement")
# DisplayTH2D(h2T, ax[0, 1], color="Blues", xlabel=r"", ylabel=r"", title="", xlim=(-8, 8), ylim=(-8, 8), rebinx = rebin, rebiny = rebin, vmax=   max2T)
# ax[0, 1].text(0.85, 0.85, "2T", transform=ax[0, 1].transAxes, ha='center', fontsize=30)
# colorbar = plt.colorbar(ax[0, 1].images[0], ax=ax[0, 1], ticks=np.arange(0, max2T+max2T/5, max2T/5), shrink=1.0)
# colorbar.set_label("Counts per bin", fontsize=15)
# colorbar.ax.set_yticklabels([str(int(i)) if i != colorbar.get_ticks()[-1] else r'$>$' + str(int(i)) for i in colorbar.get_ticks()])

# ##
# f3T = TFile("/mnt/hgfs/shared-2/2025_DATA/MCP_DATA/run_035_MCP_32Ar_BeamScan_3T_calibrated_secondstep.root", "READ")
# h3T = f3T.Get("H_reconstruction_interpolated_second_measurement")

# DisplayTH2D(h3T, ax[1,0], color="Blues", xlabel=r"$x$ (mm)", ylabel=r"$y$ (mm)", title="", xlim=(-8, 8), ylim=(-8, 8), vmax = max3T, vmin = 0, rebinx = rebin, rebiny = rebin)
# ax[1,0].set_xlabel(r"$x$ (mm)", fontsize=25)
# ax[1,0].set_ylabel(r"$y$ (mm)", fontsize=25)
# ax[1,0].text(0.85, 0.85, "3T", transform=ax[1,0].transAxes, ha='center', fontsize=30)
# colorbar = plt.colorbar(ax[1,0].images[0], ax=ax[1,0], ticks=np.arange(0, max3T+max3T/5, max3T/5), shrink=1.0)
# colorbar.ax.set_yticklabels([str(int(i)) if i != colorbar.get_ticks()[-1] else r'$>$' + str(int(i)) for i in colorbar.get_ticks()])
# # colorbar.set_label("Counts per bin", fontsize=10)
# ##
# f4T = TFile("/mnt/hgfs/shared-2/2025_DATA/MCP_DATA/CALIBRATED/run_038_MCP_32Ar_BeamScan_4T_calibrated.root", "READ")
# c4T = f4T.Get("Measurement_2D_View")
# for h in c4T.GetListOfPrimitives():
#     if h.ClassName() == "TH2D":
#         h4T = h

# DisplayTH2D(h4T, ax[1,1], color="Blues", xlabel=r"$x$ (mm)", ylabel=r"", title="", xlim=(-8, 8), ylim=(-8, 8), vmax = max4T, rebinx = rebin, rebiny = rebin)
# ax[1, 1].set_xlabel(r"$x$ (mm)", fontsize=25)
# ax[1,1].text(0.85, 0.85, "4T", transform=ax[1, 1].transAxes, ha='center', fontsize=30)
# colorbar = plt.colorbar(ax[1, 1].images[0], ax=ax[1, 1], ticks=np.arange(0, max4T+max4T/5, max4T/5), shrink=1.0)
# colorbar.set_label("Counts per bin", fontsize=15)
# colorbar.ax.set_yticklabels([str(int(i)) if i != colorbar.get_ticks()[-1] else r'$>$' + str(int(i)) for i in colorbar.get_ticks()])
# for i in range(2):
#     for j in range(2):
#         ax[i, j].set_xticks(np.arange(-8, 9, 4))
#         ax[i, j].set_yticks(np.arange(-8, 9, 4))
        

# # cbar = plt.colorbar(ax[1].images[0], ax=ax[1], ticks=np.arange(0, max4T+max4T/5, max4T/5), label = "Counts")
# # cbar.ax.set_yticklabels([str(int(i)) if i != cbar.get_ticks()[-1] else r'$>$' + str(int(i)) for i in cbar.get_ticks()])
# plt.savefig("Reviewed/AllFields_2D_Scan.pdf")
# plt.show()


##################### FUCNTION SKETCH #####################


# def linear(x, p1, p2):
#     a = (p1[1] - p2[1]) / (p1[0] - p2[0])
#     b = p1[1] - a * p1[0]
#     return a*x+b

# letter_color = "red"
# letter_size = 20
# line_color = "red"
# fig, ax = plt.subplots(figsize=(9, 8))
# #tight layout
# fig.tight_layout()
# x = np.linspace(-2, 2, 1000)

# A = (-1.5, -1)
# B = (1, -1.2)
# C = (1.1, 1)
# D = (-1, 1.3)

# xAB = np.linspace(A[0], B[0], 100)
# xBC = np.linspace(B[0], C[0], 100)
# xCD = np.linspace(C[0], D[0], 100)
# xDA = np.linspace(D[0], A[0], 100)

# xAB_out = np.append(np.linspace(x[0], A[0], 100), np.linspace(B[0], x[-1], 100))
# xBC_out = np.append(np.linspace(x[0], B[0], 100), np.linspace(C[0], x[-1], 100))
# xCD_out = np.append(np.linspace(x[0], C[0], 100), np.linspace(D[0], x[-1], 100))
# xDA_out = np.append(np.linspace(x[0], D[0], 100), np.linspace(A[0], x[-1], 100))

# # xAB, xBC, xCD, xDA = x, x, x, x

# y1 = linear(xAB, A, B)
# y2 = linear(xBC, B, C)
# y3 = linear(xCD, C, D)
# y4 = linear(xDA, D, A)

# y1_out = linear(xAB_out, A, B)
# y2_out = linear(xBC_out, B, C)
# y3_out = linear(xCD_out, C, D)
# y4_out = linear(xDA_out, D, A)

# ## contour
# ax.plot(xAB, y1, color=line_color, label = r"Cell $j$ contour")
# ax.text((A[0] + B[0]) / 2, (A[1] + B[1]) / 2 - 0.05, r'$f^{AB}_j(X)$', fontsize=letter_size, ha='center', va='top', color = letter_color)
# ax.plot(xBC, y2, color=line_color)
# ax.text((B[0] + C[0]) / 2 + 0.05, (B[1] + C[1]) / 2, r'$\;\;f^{BC}_j(Y)$', fontsize=letter_size, ha='left', va='top', color = letter_color)
# ax.plot(xCD, y3, color=line_color)
# ax.text((C[0] + D[0]) / 2, (C[1] + D[1]) / 2, r'$f^{CD}_j(X)$', fontsize=letter_size, ha='left', va='bottom', color = letter_color)
# ax.plot(xDA, y4, color=line_color)
# ax.text((D[0] + A[0]) / 2, (D[1] + A[1]) / 2, r'$f^{DA}_j(Y)$', fontsize=letter_size, ha='right', va='bottom', color = letter_color)

# ## extend 
# ax.plot(xAB_out, y1_out, color=line_color, ls = "dashed", label = r"$f_{j}$")
# ax.plot(xBC_out, y2_out, color=line_color, ls = "dashed")
# ax.plot(xCD_out, y3_out, color=line_color, ls = "dashed")
# ax.plot(xDA_out, y4_out, color=line_color, ls = "dashed")

# ax.text(*[A[0] - 0.05, A[1] + 0.05], 'A', fontsize=letter_size, ha='right', va='bottom', color = letter_color)
# ax.scatter(*A, color='red')
# ax.text(*[B[0] - 0.15, B[1] - 0.05], 'B', fontsize=letter_size, ha='left', va='top', color = letter_color)
# ax.scatter(*B, color='red')
# ax.text(*[C[0] + 0.05, C[1] - 0.05], 'C', fontsize=letter_size, ha='left', va='top', color = letter_color)
# ax.scatter(*C, color='red')
# ax.text(*[D[0] + 0.2, D[1]], 'D', fontsize=letter_size, ha='right', va='bottom', color = letter_color)
# ax.scatter(*D, color='red')

# x_fill = [A[0], B[0], C[0], D[0]]
# y_fill = [A[1], B[1], C[1], D[1]]
# ax.fill(x_fill, y_fill, color='#89BAD9', zorder = 0, label = r"$I_j(X,Y) = 1$")

# ax.set_xticks([])
# ax.set_yticks([])
# ax.set_xlabel(r'$X$', fontsize = 25)
# ax.set_ylabel(r'$Y$', fontsize = 25)
# ax.set_xlim(-2, 2)
# ax.set_ylim(-2, 2)

# ax.legend(loc='upper right', fontsize=15, framealpha=1.)
# plt.savefig("Function_Sketch.pdf")
# plt.show()

"""
##################### OSC #####################
fig, ax = plt.subplots(figsize=(9, 3), constrained_layout=True)
x = []
y = []
with open("osc_data.csv", "r") as f:
    for line in f:
        if line.startswith("x, y"):
            continue
        data = line.strip().split(",")
        x.append(float(data[0]))
        y.append(float(data[1]))
    
### sort x and chnage y in consequence
x, y = zip(*sorted(zip(x, y)))
        

plt.plot(x, y, color=color[0], label="Data", linewidth = 1.5)

plt.xlabel("t (ns)")
plt.ylabel("V (mV)")
plt.xlim(min(x), max(x))
plt.ylim(1.2*min(y), 1.1*max(y))
plt.savefig("osc_data.pdf")
"""


###################### 2T differents HV #####################
# cmap = plt.get_cmap('Blues')
# file = TFile("/mnt/hgfs/shared-2/2025_DATA/MCP_DATA/CALIBRATED/run_032_MCP_32Ar_Beam_2T_calibrated.root", "READ")
# fig, ax = plt.subplots(figsize = (10, 8), constrained_layout=True, nrows=1, ncols=2)
# max = 250
# min = 1
# zlog = True
# xlim = (-4, 4)
# ylim = (-4, 4)
# c = file.Get("Measurement_2D_View")
# h2 = None
# for h in c.GetListOfPrimitives():
#     if h.ClassName() == "TH2D":
#         h2 = h
# i = h2.GetEntries()
# DisplayTH2D(h2, ax[0], color=cmap, zlog=zlog, vmin = min, vmax=max, xlabel=r"$x$ (mm)", ylabel=r"$y$ (mm)", title="HV = 14 kV", xlim=xlim, ylim=ylim)

# file2 = TFile("/mnt/hgfs/shared-2/2025_DATA/MCP_DATA/CALIBRATED/run_033_MCP_32Ar_Beam_2T_calibrated.root", "READ")
# c2 = file2.Get("Measurement_2D_View")
# h22 = None
# for h in c2.GetListOfPrimitives():
#     if h.ClassName() == "TH2D":
#         h22 = h
# h22.Scale(i / h22.GetEntries())
# DisplayTH2D(h22, ax[1], color=cmap, zlog=zlog, vmin = min, vmax=max, xlabel=r"$x$ (mm)", ylabel=r"", title="HV = 15 kV", xlim=xlim, ylim=ylim)
# import matplotlib.colors as colors
# cbar = plt.colorbar(ax[1].images[0], ax=ax[1], ticks=np.arange(min, max+1, int(max/10)), label = "Counts", norm=colors.LogNorm(vmin=min, vmax=max))
# cbar.ax.set_yticklabels([str(int(i)) if i != cbar.get_ticks()[-1] else r'$>$' + str(int(i)) for i in cbar.get_ticks()])

# plt.savefig("2D_Scan_2T_different_HV.pdf")
# plt.show()
