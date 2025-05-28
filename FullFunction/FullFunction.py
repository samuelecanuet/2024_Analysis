### computing recoil broadering on proton spectrum ###

import numpy as np
import matplotlib.pyplot as plt
from root2mpl import * 
from ROOT import *


## IMPUT PARAMETERS ##

M = 32.0 * 931.5 # MeV/c^2
Z = 18
A = 1
E_p = 3.356 #MeV
W0 = np.sqrt(6.046**2+0.511**2) #Me
m_p = 938.3 #MeV/c^2
m_e = 0.511 #MeV/c^2

a = 1
 
alpha=1./137 #fine structure constant
c= 3*10**8  #m/s

Phi_max = np.arccos(1)
Phi_min = np.arccos(-1)


## FUNCTIONS ##

def EnergyToMomentum (E, Mass):
    return np.sqrt(E**2 - Mass**2)

def MomentumToEnergy (p, Mass):
    return np.sqrt(p**2 + Mass**2)

def EnergyToVelocity (E, Mass):
    return E/Mass

def K (Ep):
    return np.sqrt(2*Ep/M * (m_p * (M-m_p) ) / M**2)


def Fermi (Z, E_b):
    v = EnergyToVelocity(E_b, M)
    nu = Z*alpha/v*c
    return 2*np.pi*nu/(1-np.exp(-2*np.pi*nu))

def C1(W0, W, t):
    value = 0
    p_b = EnergyToMomentum(W, m_e)
    if (p_b == 0):
        return 0
    return np.minimum(np.cos(Phi_max), np.maximum(np.cos(Phi_min), (W - W0 - t/k)/p_b))
    
def C2(W0, W, t):
    value = 0
    p_b = EnergyToMomentum(W, m_e)
    if (p_b == 0):
        return 0
    return np.maximum(np.cos(Phi_min), np.minimum(np.cos(Phi_max), (W0 - W - t/k)/p_b))

def fWmin(W0, k_, t_):
    if t_/k_ < -(W0-m_e):
        value = (m_e**2 + (W0+t_/k_)**2)/(2*(W0+t_/k_))
    if (abs(t_/k_) <= W0 - m_e):
        value = m_e
    if t_/k_ > W0 - m_e:
        value = (m_e**2 + (W0-t_/k_)**2)/(2*(W0-t_/k_))

    return value
    

def WeightKinematicShift(k_, W0_, W_, t_):
    c1 = C1(W0_, W_, t_)
    c2 = C2(W0_, W_, t_)
    p_b = EnergyToMomentum(W_, m_e)
    I1 = (c2 - c1) * W_ * (W0_ - W_)
    I2a = - a * (c2**2-c1**2)/2 *p_b*t_/k_
    I2b = - a * (c2**3-c1**3)/3 *p_b**2
    w = Fermi(Z, W_) * p_b/(2*k) * ( I1 + I2a + I2b )
    return w

def TotalWeightKinematicShift(W0_, k_, t_):
    wt = 0
    Wmin = fWmin(W0_, k_, t_)
    for Wb in np.linspace(Wmin, W0_, 1000):
        wt += WeightKinematicShift(k_, W0_, Wb, t_)    
    return wt


T = np.linspace(-0.02, 0.02, 2000)
k = K(E_p)
Wt = []
for t in T:
    if (abs(t/k) > np.sqrt(W0**2-m_e**2)):
        Wt.append(0)
    else:
        Wt.append(TotalWeightKinematicShift(W0, k, t))

#rescale
Wt = np.array(Wt)
sum = np.sum(Wt)
Wt = Wt/np.max(Wt)*12500

dir = "../../../../../../../mnt/hgfs/shared-2/2024_DATA/SIMULATED_DATA/fe/fe/test_analysed.root"
file = TFile(dir, "READ")
fig1, ax1 = plt.subplots()
h = file.Get("p/H_E0_2212")

h.GetXaxis().SetRangeUser(3300, 3400)
Display(h, ax=ax1, xlim = (3300, 3400), legend = None)
ax1.plot(E_p*1000+T*1000, Wt, color="red")
plt.show()