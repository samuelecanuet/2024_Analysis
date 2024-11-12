from root2mpl import *
from ROOT import *
import numpy as np
dir = "../../../../../../../mnt/hgfs/shared-2/2024_DATA/DETECTOR_DATA/CALIBRATED/"

plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=-0.2)




################### 2024 ###################

def extranalbehavior(ax):
    ax[1].set_yticks([])
    ax[1].set_xlabel("Time [h]")
    ax[1].set_ylim(0, 1)
    #first target heating
    ax[1].fill_between([-53, -4.95], 0, 1, color='gray', alpha=0.5)
    ax[1].text(-29, 0.5, "First Target Heating", fontsize=11, ha='center')

    #DAQ
    ax[1].fill_between([-5.05, 0], 0, 1, color='gray', alpha=0.5)
    ax[1].text(-2.3, 0.4, "DAQ", fontsize=11, ha='center', rotation=90)

    #HRS
    ax[1].fill_between([9, 14], 0, 1, color='gray', alpha=0.5)
    ax[1].text(11.8, 0.05, "HRS tripped", fontsize=11, ha='center', rotation=90)

    #HRS
    ax[1].fill_between([46, 53], 0, 1, color='gray', alpha=0.5)
    ax[1].text(50.5, 0.05, "HRS tripped", fontsize=11, ha='center', rotation=90)

    #second target heating
    ax[1].fill_between([67, 137], 0, 1, color='gray', alpha=0.5)
    ax[1].text(100, 0.5, "Second Target Heating", fontsize=11, ha='center')

    return ax

file = TFile("../../../../2021_Analysis/New/rate2024.root")
plt.subplots_adjust(hspace=0.05)

fig, ax = plt.subplots(2, 1, figsize=(10, 5), sharex=True, gridspec_kw={'height_ratios': [3, 1]})
xlabel = ""
ylabel = r"Beam proportion [$\beta$/s]"
extranalbehavior(ax)
DisplayGraphError(file.Get("Betarate_N_ratio"), ax[0], color='blue', label=r"$^{18}\mathrm{N}$ decay chain", xlabel=xlabel, ylabel=ylabel, title="", ylim = (0, 0.51))
DisplayGraphError(file.Get("Betarate_Ar_ratio"), ax[0], color='red', label=r"$^{32}\mathrm{Ar}$ decay chain", xlabel=xlabel, ylabel=ylabel, title="")
ax[0].legend()
fig.savefig("2024_beam_proportion.svg")

##32Ar
xlabel = "Time [h]"
fig1, ax1 = plt.subplots(2, 1, figsize=(10, 5), sharex=True, gridspec_kw={'height_ratios': [3, 1]})
extranalbehavior(ax1)
DisplayGraphError(file.Get("Rate_Ar"), ax1[0], color='red', label=None, xlabel="", ylabel=r"$^{32}\mathrm{Ar}$/s", title="2024 Prodution rate", ylim=(0, 700))
fig1.savefig("2024_Ar_rate.svg")
##18N
xlabel = "Time [h]"
fig1, ax1 = plt.subplots(2, 1, figsize=(10, 5), sharex=True, gridspec_kw={'height_ratios': [3, 1]})
extranalbehavior(ax1)
DisplayGraphError(file.Get("Rate_N"), ax1[0], color='blue', label=None, xlabel="", ylabel=r"$^{18}\mathrm{N}$/s", title="2024 Prodution rate", ylim=(0, 85e3))
fig1.savefig("2024_N_rate.svg")
##beta
xlabel = "Time [h]"
fig1, ax1 = plt.subplots(2, 1, figsize=(10, 5), sharex=True, gridspec_kw={'height_ratios': [3, 1]})
extranalbehavior(ax1)
DisplayGraphError(file.Get("Beta_rate"), ax1[0], color='black', label=None, xlabel="", ylabel=r"$\beta$/s", title="2024 Prodution rate", ylim=(0, 20e4))
fig1.savefig("2024_beta_rate.svg")

##beta 32Ar 18N
xlabel = "Time [h]"
fig1, ax1 = plt.subplots(2, 1, figsize=(10, 5), sharex=True, gridspec_kw={'height_ratios': [3, 1]})
extranalbehavior(ax1)
DisplayGraphError(file.Get("Beta_rate"), ax1[0], color='black', label=r"All $\beta$", xlabel="", ylabel=r"$\beta$/s", title="2024 Prodution rate", ylim=(1e2-50, 30e4), ylog=True)
DisplayGraphError(file.Get("Rate_N"), ax1[0], color='blue', label=r"$\beta$ from $^{18}\mathrm{N}$ decay chain", xlabel="", ylabel=r"$\beta$/s", title="", ylim=(1e2-50, 30e4), ylog=True)
DisplayGraphError(file.Get("Rate_Ar"), ax1[0], color='red', label=r"$\beta$ from $^{32}\mathrm{Ar}$ decay chain", xlabel="", ylabel=r"$\beta$/s", title="2024 Prodution rate", ylim=(0, 20e4), ylog=True)
ax1[0].legend(loc = 'upper left', fontsize=10)
fig1.savefig("2024_beta_rate_Ar_N.svg")

##########################################

###### 2021 ######
file = TFile("../../../../2021_Analysis/New/rate2021.root")
plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=-0.2)
fig, ax = plt.subplots(2, 1, figsize=(10, 5), sharex=True, gridspec_kw={'height_ratios': [3, 1]})

xlabel = ""
ylabel = r"Beam proportion [$\beta$/s]"
##32Ar
xlabel = "Time [h]"
fig1, ax1 = plt.subplots(1, 1, figsize=(10, 3.5))
plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.2)
DisplayGraphError(file.Get("Rate_Ar"), ax1, color='red', label=None, xlabel=xlabel, ylabel=r"$^{32}\mathrm{Ar}$/s", title="2021 Prodution rate")
fig1.savefig("2021_Ar_rate.svg")
##18N
xlabel = "Time [h]"
fig1, ax1 = plt.subplots(1, 1, figsize=(10, 3.5))
plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.2)
DisplayGraphError(file.Get("Rate_N"), ax1, color='blue', label=None, xlabel=xlabel, ylabel=r"$^{18}\mathrm{N}$/s", title="2021 Prodution rate")
fig1.savefig("2021_N_rate.svg")

##32Ar 18N
xlabel = "Time [h]"
fig1, ax1 = plt.subplots(1, 1, figsize=(10, 3.5))
plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.2)
DisplayGraphError(file.Get("Rate_N"), ax1, color='blue', label=r"$^{18}\mathrm{N}$", xlabel=xlabel, ylabel=r"$^{18}\mathrm{N}$/s", title="2021 Prodution rate", ylog=True)
DisplayGraphError(file.Get("Rate_Ar"), ax1, color='red', label=r"$^{32}\mathrm{Ar}$", xlabel=xlabel, ylabel=r"$^{32}\mathrm{Ar}$/s", title="2021 Prodution rate")   
ax1.legend()
fig1.savefig("2021_Ar_N_rate.svg")
##### 18N spectrum on 32Ar #####

file = TFile("../../../../../../../../../mnt/hgfs/shared-2/2024_DATA/DETECTOR_DATA/CALIBRATED/Calibrated_Up.root")
fig, ax = plt.subplots(1, 1, figsize=(10, 5))
CANVAS =file.Get("32Ar_SUM")
hist = None
for primitive in CANVAS.GetListOfPrimitives():
    if primitive.GetName() == "H_Exp_All_32Ar":
        hist = primitive
        break



x = []
y = []

scale = 1300/5000
offset = 90
calib = 1.055
with open("../../../../../../../../../mnt/hgfs/shared-2/2024_DATA/R-MATRIX/18N_points.txt") as f:
    for line in f:
        x.append(float(line.split()[0])*calib-offset)
        y.append(float(line.split()[1])*scale)

#ordering the points over x
x, y = zip(*sorted(zip(x, y)))       

## spline interpolation for x and y
x = np.array(x)
y = np.array(y)
# remove all the data in x and y where x have deplicates
x_new = []
y_new = []
for i in range(len(x)):
    if x[i] not in x_new:
        x_new.append(x[i])
        y_new.append(y[i])
x_new = np.array(x_new)

xnew = np.linspace(x_new.min(), x_new.max(), 300)
from scipy.interpolate import make_interp_spline, BSpline
spl = make_interp_spline(x_new, y_new, k=2)  # type: BSpline
y = spl(xnew)


DisplayTH1D(hist, ax=ax, color='red', label=r"$^{32}\mathrm{Ar}$", xlabel="Energy [keV]", ylabel="Counts/keV", title="$^{32}\mathrm{Ar}$ spectrum", ylim=(1, 1.1e4), ylog=True, xlim=(800, 3500))
ax.plot(xnew, y, color="blue", lw =1, label=r"$^{18}\mathrm{N}$")

l = ax.legend(loc='upper right')
l.get_frame().set_alpha(1.0)
fig.savefig("18N_on_32Ar.svg")