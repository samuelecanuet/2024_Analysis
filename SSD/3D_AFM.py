import csv
import matplotlib.pyplot as plt 
import ROOT
import time
from root2mpl import *

filename = "/mnt/hgfs/shared-2/MICRON/MICRON_raw-values-image.txt"
data = []
    
    
    
H2D = ROOT.TH2D("H2D", "H2D", 300, 0, 25, 1024, 0, 85)
counter = 0
with open(filename, 'r') as file:
    reader = csv.reader(file, delimiter='\t')
    for row in reader:
        if (len(row) == 0 or row[0].startswith('#') or row[0].startswith(' ')):
            continue
    

        l = [float(value)*1e6 for value in row]
        for (i, value) in enumerate(l):
            H2D.Fill(counter, i*85/1024, value*820/227)
            
        counter += 25/300
        
 # display histopgram
# TCanvas = ROOT.TCanvas("TCanvas", "TCanvas", 800, 600)
H1D = H2D.ProjectionY("H1D")
H1D.SetTitle("X Projection")
H1D.GetXaxis().SetTitle("X Position (um)")
H1D.GetYaxis().SetTitle("Height (um)")

DisplayTH1D(H1D)
plt.text(15, 800, '10µm', horizontalalignment='center', verticalalignment='center')
plt.text(35, 800, '30µm', horizontalalignment='center', verticalalignment='center')
plt.text(55, 800, '10µm', horizontalalignment='center', verticalalignment='center')



plt.show()

    
