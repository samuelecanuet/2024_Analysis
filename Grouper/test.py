import numpy as np

# Given values
a12 =1./5
a23 = 15
a31 =1./3

sigma_proj_12 = 11.6e4
sigma_proj_23 = 80.8e4
sigma_proj_31 = 6.05e4

# Squaring the projected resolutions
sigma_proj_12_sq = sigma_proj_12**2
sigma_proj_23_sq = sigma_proj_23**2
sigma_proj_31_sq = sigma_proj_31**2

# Coefficient matrix for the system of equations
A = np.array([
    [1, a12**2, 0],
    [0, 1, a23**2],
    [a31**2, 0, 1]
])

# Right-hand side (RHS) vector
b = np.array([sigma_proj_12_sq, sigma_proj_23_sq, sigma_proj_31_sq])

# Solving for (sigma_1^2, sigma_2^2, sigma_3^2)
sigma_sq = np.linalg.solve(A, b)

# Taking square roots to get sigmas
print(sigma_sq)
sigma_1, sigma_2, sigma_3 = np.sqrt(sigma_sq)

print(f"Sigma 1: {sigma_1:.2e}")
print(f"Sigma 2: {sigma_2:.2e}")
print(f"Sigma 3: {sigma_3:.2e}")



import numpy as np
from scipy.optimize import curve_fit

# Given parameters
alpha1 = 30
alpha2 = 20
alpha3 = 10
beta1 = 10e3
beta2 = 20e3
beta3 = 30e3

# Given projected resolution equations
def sigma_proj12(E, gammaproj12, alpha12proj, beta12_proj):
    return gammaproj12 * E + alpha12proj * np.sqrt(E) + beta12_proj

def sigma_proj31(E, gammaproj31, alpha31proj, beta31_proj):
    return gammaproj31 * E + alpha31proj * np.sqrt(E) + beta31_proj

def sigma_proj23(E, gammaproj23, alpha23proj, beta23_proj):
    return gammaproj23 * E + alpha23proj * np.sqrt(E) + beta23_proj

# Experimental or given projected resolutions at different energies (E)
E_values = np.array([100, 200, 300, 400, 500])  # Example energy values, modify as needed
sigma12_proj_values = np.array([11.6e4, 12.3e4, 13.1e4, 13.8e4, 14.5e4])  # Example sigma12_proj values
sigma23_proj_values = np.array([80.8e4, 82.2e4, 83.7e4, 85.0e4, 86.3e4])  # Example sigma23_proj values
sigma31_proj_values = np.array([6.05e4, 6.15e4, 6.28e4, 6.40e4, 6.55e4])  # Example sigma31_proj values

# Fit the data to get the coefficients (gamma, alpha, beta) for each pair

# Fit function for sigma_proj12
params_12, _ = curve_fit(sigma_proj12, E_values, sigma12_proj_values, p0=[0.0, 0.0, 0.0])
gammaproj12, alpha12proj, beta12_proj = params_12

# Fit function for sigma_proj31
params_31, _ = curve_fit(sigma_proj31, E_values, sigma31_proj_values, p0=[0.0, 0.0, 0.0])
gammaproj31, alpha31proj, beta31_proj = params_31

# Fit function for sigma_proj23
params_23, _ = curve_fit(sigma_proj23, E_values, sigma23_proj_values, p0=[0.0, 0.0, 0.0])
gammaproj23, alpha23proj, beta23_proj = params_23

# Output the computed coefficients
print("Computed Coefficients:")
print(f"gammaproj12 = {gammaproj12}, alpha12proj = {alpha12proj}, beta12_proj = {beta12_proj}")
print(f"gammaproj31 = {gammaproj31}, alpha31proj = {alpha31proj}, beta31_proj = {beta31_proj}")
print(f"gammaproj23 = {gammaproj23}, alpha23proj = {alpha23proj}, beta23_proj = {beta23_proj}")






