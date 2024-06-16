from root2mpl import *
from ROOT import *
dir = "../../../../../../run/media/local1/Disque_Dur/Grouper/Merged/"

############silicium
# fig, ax = plt.subplots(figsize = (16, 4))
# file = TFile(dir + "runs_32Ar_merged.root", "READ")
# histo = file.Get("Silicon_Calibrated/HStrip_D5.1")
# Display(histo, ylog=True, xlim=(0, 6500), ylim=(1, 1e4), title = "", ax = ax, lw = 0.4)
# ax.set_ylim(1, 1e4)
# plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.2)
# plt.savefig("D5.1_all.png", dpi = 300)

# fig, ax = plt.subplots(figsize = (6, 4))
# file = TFile(dir + "runs_32Ar_merged.root", "READ")
# histo = file.Get("Silicon_Calibrated/Hist/D5.1")
# Display(histo, xlim=(3280, 3380), ylim=(1, 7e3), title = "", ax = ax, lw = 0.4)
# plt.subplots_adjust(left=0.25, right=0.9, top=0.9, bottom=0.2)
# plt.savefig("D5.1_IAS.png", dpi = 300)

# fig, ax = plt.subplots(figsize = (6, 4))
# file = TFile(dir + "runs_32Ar_merged.root", "READ")
# histo = file.Get("Silicon_Calibrated/Hist/D2.1")
# Display(histo, xlim=(3280, 3380), ylim=(1, 7e3), title = "", ax = ax, lw = 0.4)
# plt.subplots_adjust(left=0.25, right=0.9, top=0.9, bottom=0.2)
# plt.savefig("D2.1_IAS.png", dpi = 300)

xlim = (3280, 3380)
ylim = (1, 1e4)
dir = "../../../../../../mnt/hgfs/shared-2/plota/"

#a=-1
file = TFile(dir + "Ar_a-1_b0.root", "READ")
fig, ax = plt.subplots(figsize = (6, 4))
histo_m1 = file.Get("1Down_Strip_1_coinc")
histo_m1_2 = file.Get("2Down_Strip_1_coinc")
histo_m1.Add(histo_m1_2, 1)
histo_m1_3 = file.Get("3Down_Strip_1_coinc")
histo_m1.Add(histo_m1_3, 1)
histo_m1_4 = file.Get("4Down_Strip_1_coinc")
histo_m1.Add(histo_m1_4, 1)
histo_m1no = file.Get("1Down_Strip_1_nocoinc")
histo_m1no.Add(histo_m1, 1)
histo_m1_2no = file.Get("2Down_Strip_1_nocoinc")
histo_m1no.Add(histo_m1_2no, 1)
histo_m1_3no = file.Get("3Down_Strip_1_nocoinc")
histo_m1no.Add(histo_m1_3no, 1)
histo_m1_4no = file.Get("4Down_Strip_1_nocoinc")
histo_m1no.Add(histo_m1_4no, 1)
DisplayTH1D(histo_m1, xlim=xlim, ylim=ylim, title = "", ax = ax, lw = 0.4, color = "red", label = "1Down", rebin=5)
DisplayTH1D(histo_m1no, xlim=xlim, ylim=ylim, title = "", ax = ax, lw = 0.4, color = "black", label = "1Down", rebin=5)
plt.subplots_adjust(left=0.2, right=0.9, top=0.9, bottom=0.2)
ax.set_ylabel("Counts / 0.5 keV")
ax.set_xlabel("Energy [keV]")
plt.savefig("Ar_a-1_b0.png", dpi = 300)

#a=0
file = TFile(dir + "Ar_a0_b0.root", "READ")
fig, ax = plt.subplots(figsize = (6, 4))
histo_0 = file.Get("1Down_Strip_1_coinc")
histo_0_2 = file.Get("2Down_Strip_1_coinc")
histo_0.Add(histo_0_2, 1)
histo_0_3 = file.Get("3Down_Strip_1_coinc")
histo_0.Add(histo_0_3, 1)
histo_0_4 = file.Get("4Down_Strip_1_coinc")
histo_0.Add(histo_0_4, 1)
histo_0no = file.Get("1Down_Strip_1_nocoinc")
histo_0no.Add(histo_0, 1)
histo_0_2no = file.Get("2Down_Strip_1_nocoinc")
histo_0no.Add(histo_0_2no, 1)
histo_0_3no = file.Get("3Down_Strip_1_nocoinc")
histo_0no.Add(histo_0_3no, 1)
histo_0_4no = file.Get("4Down_Strip_1_nocoinc")
histo_0no.Add(histo_0_4no, 1)
DisplayTH1D(histo_0, xlim=xlim, ylim=ylim, title = "", ax = ax, lw = 0.4, color = "red", label = "1Down", rebin=5)
DisplayTH1D(histo_0no, xlim=xlim, ylim=ylim, title = "", ax = ax, lw = 0.4, color = "black", label = "1Down", rebin=5)
ax.set_ylabel("Counts / 0.5 keV")
plt.subplots_adjust(left=0.2, right=0.9, top=0.9, bottom=0.2)
plt.savefig("Ar_a0_b0.png", dpi = 300)

#a=1
file = TFile(dir + "Ar_a1_b0.root", "READ")
fig, ax = plt.subplots(figsize = (6, 4))
histo_1 = file.Get("1Down_Strip_1_coinc")
histo_1_2 = file.Get("2Down_Strip_1_coinc")
histo_1.Add(histo_1_2, 1)
histo_1_3 = file.Get("3Down_Strip_1_coinc")
histo_1.Add(histo_1_3, 1)
histo_1_4 = file.Get("4Down_Strip_1_coinc")
histo_1.Add(histo_1_4, 1)
histo_1no = file.Get("1Down_Strip_1_nocoinc")
histo_1no.Add(histo_1, 1)
histo_1_2no = file.Get("2Down_Strip_1_nocoinc")
histo_1no.Add(histo_1_2no, 1)
histo_1_3no = file.Get("3Down_Strip_1_nocoinc")
histo_1no.Add(histo_1_3no, 1)
histo_1_4no = file.Get("4Down_Strip_1_nocoinc")
histo_1no.Add(histo_1_4no, 1)
DisplayTH1D(histo_1, xlim=xlim, ylim=ylim, title = "", ax = ax, lw = 0.4, color = "red", label = "1Down", rebin=5)
DisplayTH1D(histo_1no, xlim=xlim, ylim=ylim, title = "", ax = ax, lw = 0.4, color = "black", label = "1Down", rebin=5)

ax.set_ylabel("Counts / 0.5 keV")
plt.subplots_adjust(left=0.2, right=0.9, top=0.9, bottom=0.2)
plt.savefig("Ar_a1_b0.png", dpi = 300)

plt.show()

############sipm
dir = "../../../2021_Analysis/Grouper/Merged/"

# file = TFile(dir + "runs_32Ar_merged.root", "READ")
# canvas = file.Get("SiPM_Multi_F0")


# for i, primitive in enumerate(canvas.GetListOfPrimitives()):
#     if primitive.GetName() == "SiPM_Multi_F0_1":
#         for i in primitive.GetListOfPrimitives():
#             if i.GetName() == "Fermi_M1":
#                 Exp_M1 = i
#             if i.GetName() == "TreeHist_Conv_F_1D1":
#                 Simu_M1 = i

#     if primitive.GetName() == "SiPM_Multi_F0_7":
#         for i in primitive.GetListOfPrimitives():
#             if i.GetName() == "Fermi_M7":
#                 Exp_M7 = i
#             if i.GetName() == "TreeHist_Conv_F_1D7":
#                 Simu_M7 = i
# fig, ax = plt.subplots(figsize = (12, 8))
# DisplayTH1D(Exp_M1, ax = ax, color = "black", label = "Experimental", lw = 0.8)
# # DisplayTH1D(Simu_M1, ax = ax, color = "red", label = "Simulation", lw = 0.8)
# ax.set_title("Multiplicity at least 1")  
# ax.set_xlim(0, 6500)  
# ax.set_xlabel("Energy [keV]")
# ax.set_ylabel("Counts / 10 keV")
# # plt.subplots_adjust(left=0.2, right=0.9, top=0.9, bottom=0.2)
# plt.savefig("SiPM1.png", dpi = 300)

# fig, ax = plt.subplots(figsize = (12, 8))
# DisplayTH1D(Exp_M7, ax = ax, color = "black", label = "Experimental", lw = 0.8)
# # DisplayTH1D(Simu_M7, ax = ax, color = "red", label = "Simulation", lw = 0.8)
# ax.set_title("Multiplicity at least 7")
# ax.set_xlim(0, 6500)
# ax.set_xlabel("Energy [keV]")
# ax.set_ylabel("Counts / 10 keV")
# # plt.subplots_adjust(left=0.2, right=0.9, top=0.9, bottom=0.2)
# plt.savefig("SiPM7.png", dpi = 300)        
        
                
# Display(histo, xlim=(3280, 3380), ylim=(1, 7e3), title = "", ax = ax, lw = 0.4)
# plt.subplots_adjust(left=0.25, right=0.9, top=0.9, bottom=0.2)
# plt.savefig("SiPM.png", dpi = 300)

# dir = "../../../2021_Analysis/Grouper/Merged/"
# fig, ax = plt.subplots(figsize = (16, 9))
# file = TFile(dir + "runs_32Ar_merged_all.root", "READ")
# canvas = file.Get("SiPM_Multi_F0")

# for i, primitive in enumerate(canvas.GetListOfPrimitives()):
#     if primitive.GetName() == "SiPM_Multi_F0_3":
#         for i in primitive.GetListOfPrimitives():
#             if i.GetName() == "Fermi_M3":
#                 Exp_M1_ALL = i

# file = TFile(dir + "runs_32Ar_merged.root", "READ")
# canvas = file.Get("SiPM_Multi_F0")

# for i, primitive in enumerate(canvas.GetListOfPrimitives()):
#     if primitive.GetName() == "SiPM_Multi_F0_3":
#         for i in primitive.GetListOfPrimitives():
#             if i.GetName() == "Fermi_M3":
#                 Exp_M1_F = i

# DisplayTH1D(Exp_M1_ALL, ax = ax, color = "black", label = "All silicon coincidences", lw = 0.8)
# DisplayTH1D(Exp_M1_F, ax = ax, color = "red", label = "IAS coincidence", lw = 0.8)
# ax.set_xlim(0, 12000)
# ax.set_xlabel("Energy [keV]")
# ax.set_ylabel("Counts / 10 keV")
# ax.set_yscale("log")
# ax.set_ylim(1, 3e4)
# ax.set_title("")
# # plt.subplots_adjust(left=0.2, right=0.9, top=0.9, bottom=0.2)
# plt.legend()
# plt.savefig('SiPM_M3.png', dpi = 300)

##################sources
dir = "../../../../../../run/media/local1/Disque_Dur/data/2024/"
####ALPHA
# fig, ax = plt.subplots(figsize = (16, 4))
# for run, scale, color in zip([127], [3427], ["black"]):
#     if run == 124:
#         file = TFile(dir + f"run_{run}_multifast_4alpha.fast/run_{run}_multifast_4alpha.root", "READ")
#     else:
#         file = TFile(dir + f"run_{run}_data_4alpha.root", "READ")

#     histo = file.Get("D7.1/D7.1_CCRC4")
#     DisplayTH1D(histo, xlim=(2800, 6000), ylim = (1, 1e4), title = "", ax = ax, lw = 0.8, scaler = scale/2527, color = color)
#     ax.set_yscale("log")
#     ax.set_xlabel("Energy [keV]")
#     ax.set_ylabel("Counts / 2 keV")

# plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.2)
# plt.savefig("4alpha.png", dpi = 300)

####90Sr
# fig, ax = plt.subplots(figsize = (10, 4))
# for run, scale, color, B in zip([125, 126, 129, 130, 133], [3001, 1504, 2013, 2035, 2794], ["orange", "black", "red", "blue", "green"], ["0", "1", "2", "3", "4"]):
#     file = TFile(dir + f"run_{run}_data_90Sr.root", "READ")
#     histo = file.Get("BetaHi4/QDC1/BetaHi4_QDC1")
#     ax.set_yscale("log")
#     DisplayTH1D(histo, xlim=(0, 3000), title = "", ax = ax, lw = 0.8, scaler = scale/3001, color = color, label = B+"T")
#     ax.set_xlabel("Energy [keV]")
#     ax.set_ylabel("Counts / 2 keV")
#     ax.legend()

# plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.2)
# plt.savefig("90Sr.png", dpi = 300)

####207Bi
# fig, ax = plt.subplots(figsize = (10, 4))
# for run, scale, color, B in zip([134, 135], [3493, 1753], ["orange", "green"], ["0", "4"]):
#     file = TFile(dir + f"run_{run}_data_207Bi.root", "READ")
#     histo = file.Get("BetaHi4/QDC1/BetaHi4_QDC1")
#     DisplayTH1D(histo, xlim=(0, 3000), title = "", ax = ax, lw = 0.8, scaler = scale/3493, color = color, label = B+"T")
#     ax.set_yscale("log")
#     ax.set_xlabel("Energy [keV]")
#     ax.set_ylabel("Counts / 2 keV")
#     ax.legend()

# plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.2)
# plt.savefig("207Bi.png", dpi = 300)

#############MCP

# dir = "../../../../../../run/media/local1/Disque_Dur/data/2024/data/"
# cmap = "Purples"
# xlim = (-7.5, 7.5)
# ylim = (-7.5, 7.5)
# fig, ax = plt.subplots(figsize = (5, 5))
# file = TFile(dir + "MCP_010_4T.root", "READ")
# histo = file.Get("anode_image_log_rec")
# DisplayTH2D(histo, title = "", ax = ax, label = "", xlabel="X [mm]", ylabel="Y [mm]", vmin=4, vmax=15, color=cmap, xlim=xlim, ylim=ylim)
# plt.subplots_adjust(left=0.15, right=0.9, top=0.9, bottom=0.15)
# plt.savefig("32Ar_4T.png", dpi = 300)

# fig, ax = plt.subplots(figsize = (5, 5))
# file = TFile(dir + "MCP_008_4T.root", "READ")
# histo = file.Get("anode_image_log_rec")
# DisplayTH2D(histo, title = "", ax = ax, label = "", xlabel="X [mm]", ylabel="Y [mm]", vmin=1, vmax=5, color=cmap, xlim=xlim, ylim=ylim)
# plt.subplots_adjust(left=0.15, right=0.9, top=0.9, bottom=0.15)
# plt.savefig("33Ar_4T.png", dpi = 300)

# fig, ax = plt.subplots(figsize = (5, 5))
# file = TFile(dir + "MCP_008_4T.root", "READ")
# histo = file.Get("anode_image")
# DisplayTH2D(histo, title = "", ax = ax, label = "", xlabel="X [au]", ylabel="Y [au]", vmin=1, vmax=5, color=cmap, xlim=(-0.6, 0.6), ylim=(-0.6, 0.6))
# plt.subplots_adjust(left=0.15, right=0.9, top=0.9, bottom=0.15)
# plt.savefig("33Ar_4T_raw.png", dpi = 300)


# GIF MAKER

# dir = "../../../../../../mnt/hgfs/shared-2/06_Tower/"
# images = []

# from PIL import Image

# def gen_frame(path):
#     im = Image.open(path)
#     alpha = im.getchannel('A')

#     # Convert the image into P mode but only use 255 colors in the palette out of 256
#     im = im.convert('RGB').convert('P', palette=Image.ADAPTIVE, colors=255)

#     # Set all pixel values below 128 to 255 , and the rest to 0
#     mask = Image.eval(alpha, lambda a: 255 if a <=128 else 0)

#     # Paste the color of index 255 and use alpha as a mask
#     im.paste(255, mask)

#     # The transparency index is 255
#     im.info['transparency'] = 255

#     return im


# im1 = gen_frame(dir+f'WISArD0685.png')
# for i in range(685, 900):
#     if (i < 10): i_str = f"000{i}"
#     elif (i < 100): i_str = f"00{i}"
#     elif (i < 1000): i_str = f"0{i}"
#     else: i_str = f"{i}"
#     images.append(gen_frame(dir+f'WISArD{i_str}.png'))
# im1.save(dir+'GIF.gif', save_all=True, append_images=images, duration=30, disposal = 2)


