
import numpy as np

Cs = 1
Cv = 1


def Compute_a(Cs, Cv, Csp = None, Cvp = None, helicity = "None"):
    if helicity == "droit":
        Csp = -Cs
        Cvp = Cv
    elif "gauche":
        Csp = Cs
        Cvp = Cv
    else:
        if (Csp is None) or (Cvp is None) or (helicity == "None"):
            raise ValueError("For helicity other than 'droit' or 'gauche', Csp and Cvp must be provided.")
    return (-Cs**2-Csp**2+Cv**2+Cvp**2) / (Cs**2+Csp**2+Cv**2+Cvp**2)

def Compute_b(Cs, Cv, Csp = None, Cvp = None, helicity = "None"):
    if helicity == "droit":
        Csp = -Cs
        Cvp = Cv
    elif helicity == "gauche":
        Csp = Cs
        Cvp = Cv
    else :
        if (Csp is None) or (Cvp is None) or (helicity == "None"):
            raise ValueError("For helicity other than 'droit' or 'gauche', Csp and Cvp must be provided.")
    return 2*(Cs*Cv+Csp*Cvp) / (Cs**2+Csp**2+Cv**2+Cvp**2)


## gauche
helicity = "gauche" 
x = np.linspace(-1, 1, 1000)
y = np.linspace(-1, 1, 1000)

X, Y = np.meshgrid(x, y)
Z = Compute_a(X, Y, helicity=helicity)
Z1 = Compute_b(X, Y, helicity=helicity)

import matplotlib.pyplot as plt
fig, ax = plt.subplots(figsize=(14, 6), nrows=1, ncols=2, tight_layout=True)
ax[0].contourf(X, Y, Z, levels=50, cmap='viridis')
# ax[0].title(f'Contour plot of a for Cs and Cv with helicity {helicity}')
ax[0].set_xlabel('Cs')
ax[0].set_ylabel('Cv')
ax[1].contourf(X, Y, Z1, levels=50, cmap='viridis')
# ax[1].set_title(f'Contour plot of b for Cs and Cv with helicity {helicity}')
ax[1].set_xlabel('Cs')
ax[1].set_ylabel('Cv')


plt.colorbar(ax[0].contourf(X, Y, Z, levels=50, cmap='viridis'), ax=ax[0], label='a value', ticks=[-1, -0.5, 0, 0.5, 1])
plt.colorbar(ax[1].contourf(X, Y, Z1, levels=50, cmap='viridis'), ax=ax[1], label='b value', ticks=[-1, -0.5, 0, 0.5, 1])
# ax[0].axis("equal")
ax[0].set_xlim([-1, 1])
ax[0].set_ylim([-1, 1])
# ax[1].axis("equal")
ax[1].set_xlim([-1, 1])
ax[1].set_ylim([-1, 1])
plt.savefig('ab_CS_CV_gauche.png')


## droit
helicity = "droit"
Z = Compute_a(X, Y, helicity=helicity)
Z1 = Compute_b(X, Y, helicity=helicity)
fig, ax = plt.subplots(figsize=(8, 6), nrows=1, ncols=2, constrained_layout=True)
ax[0].contourf(X, Y, Z, levels=50, cmap='viridis')
# ax[0].title(f'Contour plot of a for Cs and Cv with helicity {helicity}')
ax[0].set_xlabel('Cs')
ax[0].set_ylabel('Cv')
ax[1].contourf(X, Y, Z1, levels=50, cmap='viridis')
# ax[1].set_title(f'Contour plot of b for Cs and Cv with helicity {helicity}')
ax[1].set_xlabel('Cs')
ax[1].set_ylabel('Cv')
cbar = plt.colorbar(ax[0].contourf(X, Y, Z, levels=50, cmap='viridis'), ax=ax[0], label='a value', ticks=[-1, -0.5, 0, 0.5, 1]) 
cbar = plt.colorbar(ax[1].contourf(X, Y, Z1, levels=50, cmap='viridis'), ax=ax[1], label='b value', ticks=[-1, -0.5, 0, 0.5, 1])
plt.savefig('ab_CS_CV_droit.png')

## CS/CV - CS'/CV plane

# helicity = "Mixed"
# Cs = np.linspace(-1, 1, 1000)
# Csp = np.linspace(-1, 1, 1000)
# Cv = np.linspace(-1, 1, 1000)

# rCsCv = Cs/Cv
# rCspCv = Csp/Cv
# X, Y = np.meshgrid(rCspCv, rCsCv)
# Z = Compute_a(Y*Cv, 1/Y*Cs, Csp=X*Cv, Cvp=Cv, helicity=helicity)
# Z1 = Compute_b(Y*Cv, 1/Y*Cs, Csp=X*Cv, Cvp=Cv, helicity=helicity)
# fig, ax = plt.subplots(figsize=(8, 6), nrows=1, ncols=2, constrained_layout=True)
# ax[0].contourf(X, Y, Z, levels=50, cmap='viridis')
# # ax[0].title(f'Contour plot of a for Cs/Cv and Csp/Cv with helicity {helicity}')
# ax[0].set_xlabel("Csp/Cv")
# ax[0].set_ylabel("Cs/Cv")
# ax[1].contourf(X, Y, Z1, levels=50, cmap='viridis')
# # ax[1].set_title(f'Contour plot of b for Cs/Cv and Csp/Cv with helicity {helicity}')
# ax[1].set_xlabel("Csp/Cv")
# ax[1].set_ylabel("Cs/Cv")
# plt.colorbar(ax[0].contourf(X, Y, Z, levels=50, cmap='viridis'), ax=ax[0], label='a value')
# plt.colorbar(ax[1].contourf(X, Y, Z1, levels=50, cmap='viridis'), ax=ax[1], label='b value')
# plt.savefig('ab_CsCv_CspCv_mixed.png')





