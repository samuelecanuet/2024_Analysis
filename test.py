import numpy as np
import matplotlib.pyplot as plt

# --- Physical constants ---
T_C = -10
T_K = T_C + 273.15        # Temperature in K
k_B = 8.617e-5            # Boltzmann constant in eV/K
E_a = 0.3                 # Activation energy in eV
D0 = 1e-6                 # Diffusion prefactor in cm^2/s
T_half = 0.098            # Half-life of Ar-32 in seconds
tau = T_half / np.log(2)  # Mean lifetime in s

# --- Compute diffusion coefficient at T ---
D = D0 * np.exp(-E_a / (k_B * T_K))  # cm^2/s

# --- Sampling parameters ---
N = 100000  # Number of Ar-32 ions
np.random.seed(0)

# Sample decay times from exponential distribution
t_decay = np.random.exponential(scale=tau, size=N)  # in seconds

# For each decay time, sample 3D displacements
# Each component is Gaussian with std = sqrt(2 D t)
sigma = np.sqrt(2 * D * t_decay)
x = np.random.normal(0, sigma)
y = np.random.normal(0, sigma)
z = np.random.normal(0, sigma)

# Compute radial distance
r = np.sqrt(x**2 + y**2 + z**2)  # in cm

# --- Plot histogram of distances ---
plt.figure(figsize=(8,5))
plt.hist(r*1e7, bins=150, density=True, color='steelblue', alpha=0.7)
plt.xlabel("Displacement before decay [nm]", fontsize=12)
plt.ylabel("Probability density", fontsize=12)
plt.title(rf"Diffusion of $^{{32}}$Ar in Mylar at {T_C}°C", fontsize=13)
plt.grid(True, ls="--", alpha=0.4)
plt.tight_layout()
plt.yscale('log')  # Log scale for better visibility
plt.show()