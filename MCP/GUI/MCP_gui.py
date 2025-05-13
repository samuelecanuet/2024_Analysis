import numpy as np
from matplotlib import pyplot as plt
from matplotlib.patches import Circle, Rectangle, Wedge, Polygon
from matplotlib.widgets import Button, Slider, RadioButtons
import csv
from ROOT import *
import subprocess
import argparse
import ctypes
import seaborn as sns
from matplotlib.backend_bases import MouseButton
list_corner=[]
custom_params = {
        "xtick.direction" : "out",
        "ytick.direction" : "out",
        "lines.markeredgecolor" : "k",
        "lines.markeredgewidth" : 0.5,
        "lines.linewidth" : 1,
        "lines.markersize" : 5,
        "figure.figsize" : (16,9),
        "font.family" : "serif",
        "ytick.labelsize" : 15,
        "xtick.labelsize" : 15,
        "axes.labelsize" : 20,
        "axes.titlesize" : 20,
        "legend.fontsize" : 15,
        "text.usetex" : True,
        # 'figure.subplot.left': 0.20, 
        # 'figure.subplot.bottom': 0.15, 
        # 'figure.subplot.right': 0.95, 
        # 'figure.subplot.top': 0.90
        }
#sns.set_theme(style = "ticks", rc=custom_params)

class MovePoint(object):
    def __init__(self, ax, graf, rec=None, index=None):
        self.ratio = 15
        self.ax = ax
        self.figcanvas = self.ax.figure.canvas
        self.graf = graf
        self.moved = None
        self.pointx = None
        self.pointy = None
        self.pressed = False
        self.start = False
        self.delete = False
        self.find = False
        self.initx = graf.get_center()[0]
        self.inity = graf.get_center()[1]
        self.label = graf.get_label()

        self.find = True
        # print(list_corner[index][1])
        # graf.set_center([1/25*self.initx, 1/25*self.inity])
        self.graf.set_center([list_corner[index][0], list_corner[index][1]])

        graf.set_radius(graf.get_radius()/self.ratio)

        self.figcanvas.mpl_connect('button_press_event', self.mouse_press)
        self.figcanvas.mpl_connect('button_release_event', self.mouse_release)
        self.figcanvas.mpl_connect('motion_notify_event', self.mouse_move)

    def mouse_release(self, event):
        if self.ax.get_navigate_mode()!= None: return
        if not event.inaxes: return
        if event.inaxes != self.ax: return

        if self.pressed: 
            self.pressed = False
            self.start = False
            self.pointx = None
            self.pointy = None
            if self.graf.get_facecolor()==(1.0, 0.0, 0.0, 1):
                self.graf.set_facecolor('green')
            return
        
    def mouse_press(self, event):
        if self.ax.get_navigate_mode()!= None: return
        if not event.inaxes: return
        if event.inaxes != self.ax: return
        if self.start: return
        if self.delete: return
        
        self.pointx = event.xdata
        self.pointy = event.ydata
        if self.graf.contains(event)[0]:
            if event.button is MouseButton.RIGHT:
                self.graf.remove() 
                self.delete = True
            self.pressed = True
        

    def mouse_move(self, event):
        if self.ax.get_navigate_mode()!= None: return
        if not event.inaxes: return
        if event.inaxes != self.ax: return
        if not self.pressed: return
        self.start = True
        self.graf.set_center((event.xdata, event.ydata))

    def GetCoordinates(self):
        return [self.initx, self.inity, self.graf.get_center()[0], self.graf.get_center()[1]]


class MoveRec(object):
    def __init__(self, ax, rec, points):
        self.ratio = 15
        self.ax = ax
        self.figcanvas = self.ax.figure.canvas
        self.moved = None
        self.pointx = None
        self.pointy = None
        self.pressed = False
        self.start = False
        self.delete = False
        self.rec = rec
        points_tmp = points[1]
        points[1] = points[2]
        points[2] = points_tmp
        self.points = points
        self.init1 = rec.get_xy()[0]
        self.init2 = rec.get_xy()[1]
        self.init3 = rec.get_xy()[2]
        self.init4 = rec.get_xy()[3]
        self.label = rec.get_label()

        index=False
        coordinate_change=None
        for i, point in enumerate(self.points):
                coo = self.rec.get_xy()
                coo[self.index_correction(i)] = point.graf.get_center()
        #self.rec.remove()
        self.rec = Polygon(coo[:4], fc='white', zorder=1, alpha=0.5)
        self.ax.add_patch(self.rec)
                
        self.figcanvas.mpl_connect('button_press_event', self.mouse_press)
        self.figcanvas.mpl_connect('button_release_event', self.mouse_release)
        self.figcanvas.mpl_connect('motion_notify_event', self.mouse_move)

    def mouse_release(self, event):
        if self.ax.get_navigate_mode()!= None: return
        if not event.inaxes: return
        if event.inaxes != self.ax: return
        if self.pressed: 
            self.pressed = False
            self.start = False
            self.pointx = None
            self.pointy = None
            return
        
    def mouse_press(self, event):
        if self.ax.get_navigate_mode()!= None: return
        if not event.inaxes: return
        if event.inaxes != self.ax: return
        if self.start: return
        self.pointx = event.xdata
        self.pointy = event.ydata
        if self.rec.contains(event)[0]:
            self.pressed = True
        

    def mouse_move(self, event):
        # if self.ax.get_navigate_mode()!= None: return
        # if not event.inaxes: return
        # if event.inaxes != self.ax: return
        if self.delete: return
        index=False
        coordinate_change=None
        for i, point in enumerate(self.points):
            if point.delete: 
                self.rec.remove() 
                self.delete = True
                break
            if point.pressed:
                coo = self.rec.get_xy()
                coo[self.index_correction(i)] = point.graf.get_center()
                self.rec.remove()
                self.rec = Polygon(coo[:4], fc='white', zorder=1, alpha=0.5)
                self.ax.add_patch(self.rec)
                break
                

    def index_correction(self, index ):
        if index == 1:
            return 2
        elif index == 2:
            return 1
        elif index == 3:
            return 3
        elif index == 0:
            return 0
            
    
class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

gROOT.ProcessLine("gROOT->SetBatch(kTRUE)")  
    
def DisplayTH2D(Hist, ax, color='plasma', label=None, title=None, xlabel=None, ylabel=None, xlim=None, ylim=None, xtick=None, ytick=None, ylog=None, xlog=None, zlog=None, rebinx=None, rebiny=None, vmax = None, vmin = None, view=None):
        if rebinx   != None: Hist.RebinX(rebinx)
        if rebiny   != None: Hist.RebinY(rebiny)

    
        nbins_x = Hist.GetNbinsX()
        nbins_y = Hist.GetNbinsY()

        hist_data = np.zeros((nbins_x, nbins_y))
        bin_centers_x = np.zeros(nbins_x)
        bin_centers_y = np.zeros(nbins_y)

        BD, BG, HG, HD = [0, 0], [0, 0], [0, 0], [0, 0]
        for i in range(1, nbins_x + 1):
            for j in range(1, nbins_y + 1):
                hist_data[i - 1, j - 1] = Hist.GetBinContent(i, j)
                bin_centers_x[i - 1] = Hist.GetXaxis().GetBinCenter(i)
                bin_centers_y[j - 1] = Hist.GetYaxis().GetBinCenter(j)

                if BD[0] < Hist.GetBinContent(i, j) and np.sqrt(Hist.GetXaxis().GetBinCenter(i)**2+Hist.GetYaxis().GetBinCenter(j)**2) > BD[1] and Hist.GetXaxis().GetBinCenter(i) > 0. and Hist.GetYaxis().GetBinCenter(j) < 0.:
                    BD = [Hist.GetBinContent(i, j), np.sqrt(Hist.GetXaxis().GetBinCenter(i)**2+Hist.GetYaxis().GetBinCenter(j)**2)]
                
                if BG[0] < Hist.GetBinContent(i, j) and np.sqrt(Hist.GetXaxis().GetBinCenter(i)**2+Hist.GetYaxis().GetBinCenter(j)**2) > BG[1] and Hist.GetXaxis().GetBinCenter(i) < 0. and Hist.GetYaxis().GetBinCenter(j) < 0.:
                    BG = [Hist.GetBinContent(i, j), np.sqrt(Hist.GetXaxis().GetBinCenter(i)**2+Hist.GetYaxis().GetBinCenter(j)**2)]
                
                if HG[0] < Hist.GetBinContent(i, j) and np.sqrt(Hist.GetXaxis().GetBinCenter(i)**2+Hist.GetYaxis().GetBinCenter(j)**2) > HG[1] and Hist.GetXaxis().GetBinCenter(i) < 0. and Hist.GetYaxis().GetBinCenter(j) > 0.:
                    HG = [Hist.GetBinContent(i, j), np.sqrt(Hist.GetXaxis().GetBinCenter(i)**2+Hist.GetYaxis().GetBinCenter(j)**2)]
                
                if HD[0] < Hist.GetBinContent(i, j) and np.sqrt(Hist.GetXaxis().GetBinCenter(i)**2+Hist.GetYaxis().GetBinCenter(j)**2) > HD[1] and Hist.GetXaxis().GetBinCenter(i) > 0. and Hist.GetYaxis().GetBinCenter(j) > 0.:
                    HD = [Hist.GetBinContent(i, j), np.sqrt(Hist.GetXaxis().GetBinCenter(i)**2+Hist.GetYaxis().GetBinCenter(j)**2)]

        lim = max([i[1] for i in [BD, BG, HG, HD]])/np.sqrt(2)*1.1


        if label  == None: label = Hist.GetTitle()
        if title  == None: title = Hist.GetTitle()
        if xlabel == None: xlabel = Hist.GetXaxis().GetTitle()
        if ylabel == None: ylabel = Hist.GetYaxis().GetTitle()
        if xlim   == None: xlim = ( bin_centers_x.min(), bin_centers_x.max() )
        if ylim   == None: ylim = (bin_centers_y.min(), bin_centers_y.max())
        if xlim   == 'auto': xlim = (-lim, lim)
        if ylim   == 'auto': ylim = (-lim, lim)
        if xtick  != None: ax.set_xticks(np.linspace(xlim[0], xlim[1], xtick))
        if ytick  != None: ax.set_yticks(np.linspace(ylim[0], ylim[1], ytick))
        if xlog   != None: ax.set_xscale('log')
        if ylog   != None: ax.set_yscale('log')

        cax = ax.imshow(hist_data.T, extent=(bin_centers_x.min(), bin_centers_x.max(), bin_centers_y.min(), bin_centers_y.max()), origin='lower', aspect='auto', cmap=color)

        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)

        if vmax      != None: cax.set_clim(vmax = vmax)
        if vmin      != None: cax.set_clim(vmin = vmin)

        return [cax, ax]

current=[]
current_patch=[]
sc = None

def on_move(event):
    global sc
    if event.inaxes is axs[0,0]:
        x, y = event.xdata, event.ydata
        del axs[1,0].get_children()[0]       
        axs[1,0].set_xlim(x - 0.3, x + 0.3)
        axs[1,0].set_ylim(y - 0.15, y + 0.15)
        
        if sc: sc.remove()
        sc = axs[1,0].scatter(x, y, marker = "+", s = 125, color="red")
        
        event.canvas.draw()
   
def writer(event):
    with open("guess_center_"+args.year+".txt", "w") as file:
        writer = csv.writer(file, delimiter='\t')
        counterp = 0
        listep=[]
        for i, point in enumerate(liste1):
            counterp+=1
            if not point.delete:
                values = point.GetCoordinates()
                listep.append(values[2])
                listep.append(values[3])
                if counterp == 4:
                    counterp = 0
                    writer.writerow([int((i+1)/4)] + listep)
                    listep=[]

def printer(nombre, error):
    if error == 0:
        return (nombre, error)
    else:
        str_er = list(str(error))
        if str_er[0] != "0": return (float(str(round(nombre))[0]), round(error))
        for i in range(len(str_er)):
            if str_er[i] != "0" and str_er[i] != ".":
                index = i
                break
        str_er.index(str_er[index])
        return float(str(round(nombre, index-1))), float(str(round(error, index-1)))

def update_max(val):
    HIST[0].set_clim(vmax = val)
    HIST_z[0].set_clim(vmax = val)
    if len(axs[0,1].get_images()) > 0: axs[0,1].get_images()[0].set_clim(vmax = val)

def update_min(val):
    HIST[0].set_clim(vmin = val)
    HIST_z[0].set_clim(vmin = val)
    if len(axs[0,1].get_images()) > 0: axs[0,1].get_images()[0].set_clim(vmin = val)


def loadpointcorner():
    l=[]
    with open("../guess_corner_"+args.year+"_test.txt", "r") as file:
        reader = csv.reader(file, delimiter=' ')
        for i, row in enumerate(reader):
            liste = row[1:9]
            liste = [float(i) for i in liste]
            l.append([liste[0], liste[1]])
            l.append([liste[2], liste[3]])
            l.append([liste[4], liste[5]])
            l.append([liste[6], liste[7]])

    return l


if __name__ == '__main__': 

    ### ARGUMENTS 
    parser = argparse.ArgumentParser(description="MCP reconstruction script",
                                    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-f", "--filename", action="store", help="Input filename", default=False)
    parser.add_argument("-y", "--year", action="store", help="Year of the data", default="2025")   
    args = parser.parse_args()
    config = vars(args)       

    ### INPUT MCP DATA  
    fig, ax = plt.subplots()
    
    if args.filename:
        input_file = args.filename

        if "root" in input_file:
            try : root_file = TFile(input_file); print(f"\n {bcolors.OKGREEN} Openning : {input_file} {bcolors.ENDC}")
            except OSError : print(f"\n {bcolors.FAIL} Wrong file name {bcolors.ENDC}\n"), exit(0)
            
            try : HIST = DisplayTH2D(root_file.Get("h_Image"), ax, title = "MCP data", xlim='auto' , ylim='auto'); print(f" {bcolors.OKGREEN} Openning : h_Image {bcolors.ENDC}"); name = True
            except AttributeError : name = None
           
            while name == None:
                name = input(f" {bcolors.WARNING} The TH2D 'h_Image' doesn't exist. Give a TH2D name : {bcolors.ENDC}")
                try : HIST = DisplayTH2D(root_file.Get(name), ax, title = "MCP data", xlim='auto' , ylim='auto'); print(f" {bcolors.OKGREEN} Openning : {name} {bcolors.ENDC}")
                except AttributeError : name = input(f" {bcolors.WARNING} The TH2D 'h_Image' doesn't exist. Give a TH2D name : {bcolors.ENDC}"); name = None
            


        elif "fast" in input_file:
            try : print(f"\n {bcolors.OKGREEN} Running  : fast2root {bcolors.ENDC}"); print(f" {bcolors.OKGREEN} Openning : {input_file} {bcolors.ENDC}"); resultat = subprocess.run(["fast2root", input_file[:-5]], stdout=subprocess.DEVNULL)
            except FileNotFoundError : print(f"\n {bcolors.FAIL} Wrong Executable Name File (Fast Reader Name group2tree_SL_MCP) {bcolors.ENDC}\n"), exit(0)

            input_file = input_file[:-5]+".root"
            try : root_file = TFile(input_file); print(f" {bcolors.OKGREEN} Openning : {input_file} {bcolors.ENDC}")
            except OSError : print(f"\n {bcolors.FAIL} Wrong file name {bcolors.ENDC}\n"), exit(0)
            
            try : HIST = DisplayTH2D(root_file.Get("h_Image"), ax, title = "MCP data", xlim='auto' , ylim='auto'); print(f" {bcolors.OKGREEN} Openning : h_Image {bcolors.ENDC}"); name = True
            except AttributeError : name = None
           
            while name == None:
                name = input(f" {bcolors.WARNING} The TH2D 'h_Image' doesn't exist. Give a TH2D name : {bcolors.ENDC}")
                try : HIST = DisplayTH2D(root_file.Get(name), ax, title = "MCP data", xlim='auto' , ylim='auto'); print(f" {bcolors.OKGREEN} Openning : {name} {bcolors.ENDC}")
                except AttributeError : name = input(f" {bcolors.WARNING} The TH2D 'h_Image' doesn't exist. Give a TH2D name : {bcolors.ENDC}"); name = None
            


        else:
            print(f"\n{bcolors.FAIL} Error : Please give a root or fast file as argument {bcolors.ENDC}\n")
            exit(0)
    else:
        print(f"\n{bcolors.FAIL} Error : Please give a root or fast file as argument {bcolors.ENDC}\n")
        exit(0)
    plt.close()    


    ### MODE SELECTION
    if True:
        ### CANVAS
        list_point = []
        global current_point
        for i in range(16*16):
            list_point.append([[None, None],[None, None]])


        fig, axs = plt.subplots(2, 2, figsize = (15, 9), gridspec_kw=dict(height_ratios=(2, 1)))
        axs[0,0].set_title("MCP data")
        axs[0,1].set_title("MCP reconstructed")
        axs[0,1].set_xlabel("X (mm)")
        axs[0,1].set_xlabel("Y (mm)")
        axs[1,1].remove()
        axs[1,1].set_xlabel("Distance from (0,0)")
        axs[1,1].set_ylabel("Distance from real point")

        HIST = DisplayTH2D(root_file.Get("h_Image"), axs[0,0], title = "MCP data", xlim='auto' , ylim='auto')
        HIST_z = DisplayTH2D(root_file.Get("h_Image"), axs[1,0], title = "MCP data", xlim='auto' , ylim='auto')

        ### CONSTRUCT MCP SKETCH
        pitch = 2
        size = 1.2
        n=  8
        if (args.year == "2025"): 
            n = 7
            size = 1.4
        radius = 0.15
        MCP = 10
        color_point='green'
        color_rec='black'

        liste = []
        liste1=[]
        list_corner = loadpointcorner()
        print(list_corner)
        for x in range(0,n):
            for y in range(0,n):
                allpoint=0
                x_corr = x*pitch-pitch*(n-1)/2 - size/2
                y_corr = y*pitch-pitch*(n-1)/2 - size/2

                # print(list_corner[len(liste)][0], list_corner[len(liste)][1])

                c1 = Circle((x_corr, y_corr), radius, fc=color_point, label=2*x+8*4*y, zorder=2)
                list_point[2*x+8*4*y][0] = [x_corr, y_corr]
                liste.append(axs[0,0].add_patch(c1))
                liste1.append(MovePoint(axs[0,0], liste[-1], index =len(liste1)))
                
                allpoint+=1
                c2 = Circle((x_corr+size, y_corr), radius, fc=color_point, label=2*x+8*4*y+1, zorder=2)
                list_point[2*x+8*4*y+1][0] = [round(x_corr+size, 5), round(y_corr, 5)]
                liste.append(axs[0,0].add_patch(c2))
                liste1.append(MovePoint(axs[0,0], liste[-1], index =len(liste1)))
                allpoint+=1
                
                c3 = Circle((x_corr+size, y_corr+size), radius, fc=color_point, label=2*x+8*4*y+2*8, zorder=2)
                list_point[2*x+8*4*y+2*8][0] = [round(x_corr+size, 5), y_corr+size]
                liste.append(axs[0,0].add_patch(c3))
                liste1.append(MovePoint(axs[0,0], liste[-1], index =len(liste1) ))
                allpoint+=1

                c4 = Circle((x_corr, y_corr+size), radius, fc=color_point, label=2*x+8*4*y+2*8+1, zorder=2)
                list_point[2*x+8*4*y+2*8+1][0] = [x_corr, round(y_corr+size, 5)]
                liste.append(axs[0,0].add_patch(c4))
                liste1.append(MovePoint(axs[0,0], liste[-1], index =len(liste1) ))
                allpoint+=1
                
                if allpoint == 4:
                    rec = Polygon(( (x_corr, y_corr), (x_corr+size, y_corr), (x_corr+size, y_corr+size), (x_corr, y_corr+size)), fc=color_rec, zorder=1, alpha=0.3)
                    liste.append(MoveRec(axs[0,0], rec, liste1[-4:]))
                    #axs[0,0].add_patch(rec)       


        ### WRITER BUTTON
        button_pos = plt.axes([0.6, 0.14, 0.2, 0.05])
        button_w = Button(button_pos, 'Write coordinates file')
        button_w.on_clicked(writer)

        ### VMAX SLIDER
        axfreq_max = fig.add_axes([0.6, 0.24, 0.2, 0.01])
        freq_slider_max = Slider(
        ax=axfreq_max,
        label='Max Value',
        valmin=1,
        valmax=150,
        valinit=10)
        # update_max(axfreq_max.valinit)   
        freq_slider_max.on_changed(update_max)

        ### VMIN SLIDER
        axfreq_min = fig.add_axes([0.6, 0.29, 0.2, 0.01])
        freq_slider_min = Slider(
        ax=axfreq_min,
        label='Min Value',
        valmin=1,
        valmax=100,
        valinit=1)   
        # update_min(axfreq_min.valinit)
        freq_slider_min.on_changed(update_min)

        ### LENS
        fig.canvas.mpl_connect('motion_notify_event', on_move)
        axs[1,0].set_xticks([])
        axs[1,0].set_yticks([])
        axs[1,0].set_title("")

        ###HIDER
        def hider(label):
            if label == 'hide':
                HIST[0].set_alpha(0)
            else:
                HIST[0].set_alpha(1)
            fig.canvas.draw()

        radio = RadioButtons(fig.add_axes([0.82, 0.1, 0.06, 0.10]), ('show', 'hide'))
        radio.on_clicked(hider)
        

        plt.subplots_adjust(bottom=0.03, top=0.95)
        plt.show()
