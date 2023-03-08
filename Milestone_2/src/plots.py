import matplotlib.pyplot as plt
import numpy as np
import astropy.constants as const
import astropy.units as units
from scipy.optimize import curve_fit
#plt.rcParams['text.usetex'] = True

file = np.loadtxt("../cosmology.txt")

x           = file[:,0]
eta_of_x    = file[:,1]
Hp_of_x     = file[:,2]
dHpdx_of_x  = file[:,3]
OmegaB      = file[:,4]
OmegaCDM    = file[:,5]
OmegaLambda = file[:,6]
OmegaR      = file[:,7]
OmegaNu     = file[:,8]
OmegaK      = file[:,9]
dL_x        = file[:,10]
OmegaM = OmegaB + OmegaCDM
sum_omega = OmegaB + OmegaCDM + OmegaLambda + OmegaNu + OmegaK

"""
PLOT 1, H(x)
"""
fig, ax = plt.subplots(figsize=(6, 4), tight_layout=True)
ax.set_yscale('log')
ax.plot(x, Hp_of_x*const.pc*10**6/10**5)
ax.set_ylim(-1, 1000)
ax.set_xlim(-12, 0)
ax.set_xlabel('x = log(a)')
ax.set_title(r'$\mathcal{H}(x)\frac{100 km/s}{Mpc}$')
"""
PLOT 2, eta(x)
"""
fig, ax = plt.subplots(figsize=(6, 4), tight_layout=True)
ax.set_yscale('log')
ax.plot(x, eta_of_x/(const.pc*10**6))
ax.set_xlim(-12, 0)
ax.set_ylim(1, 10**5)
ax.set_xlabel('x = log(a)')
ax.set_title(r'$\eta(x)$ [Mpc]')
"""
PLOT 3, eta(x)*H(x)/c
"""
fig, ax = plt.subplots(figsize=(6, 4), tight_layout=True)
#ax.set_title(r"\frac{\eta(x)H}{c}")
ax.set_ylim(0.75, 3)
ax.set_xlim(-14, 0)
ax.plot(x, eta_of_x*Hp_of_x/const.c.value)
ax.set_xlabel('x = log(a)')
ax.set_title(r'$\frac{\eta(x)\mathcal{H}}{c}$')
"""
PLOT 4, Omegas(x)
"""
fig, ax = plt.subplots(figsize=(6, 4), tight_layout=True)
ax.plot(x, OmegaR+OmegaNu, label = r"$\Omega_{relativistic}$")
ax.plot(x, OmegaM, label = r"$\Omega_{matter}$")
ax.plot(x, OmegaLambda, label = r"$\Omega_{\Lambda}$")
ax.set_xlabel('x = log(a)')
ax.legend()
plt.ylim(0, 1)


"""
Calculating the age of the universe
(Done in BackgroundCosmology.cpp, run make hten ./cmb in ../src)
"""





"""
PLOT ????
"""
def x_to_z(x):
    return np.exp(-x) - 1

z = x_to_z(x)

data = np.loadtxt("../data/supernovadata.txt")
z_data = data[:,0]
dL_data =  data[:,1]
error_y = data[:,2]



fig, ax = plt.subplots(figsize=(6, 4), tight_layout=True)
ax.errorbar(z_data, dL_data/z_data, error_y/z_data, fmt=".k")
ax.fill_between(z_data, (dL_data - error_y)/z_data, (dL_data + error_y)/z_data, color="gray", alpha=0.2)
ax.plot(z, (dL_x*units.m).to('Gpc')/z, color="red")
ax.set_ylim(3.5, 8)
ax.set_xlim(0.0075, 1.4)
ax.set_xscale('log')
ax.set_xlabel("z")
ax.set_ylabel(r"dL/z [Gpc]")
#plt.show()

"""
PLOT 5
"""
supernova_chi2 = np.loadtxt("../results_supernovafitting.txt", skiprows=200)
chi2_fit = supernova_chi2[:,0]
h_fit = supernova_chi2[:,1]
OmegaM_fit = supernova_chi2[:,2]
OmegaK_fit = supernova_chi2[:,3]
least_chi2 = np.argmin(chi2_fit)
OmegaL_fit = 1 - (OmegaM_fit + OmegaK_fit)

best_OmegaM = OmegaM_fit[least_chi2]
best_OmegaK = OmegaK_fit[least_chi2]
best_OmegaL = 1 - (best_OmegaM + best_OmegaK)


"""
Plot Gaussian!!
"""
fig, ax = plt.subplots(figsize=(6, 4), tight_layout=True)
def gauss(x, C, mu, sigma):
    return C*np.exp(-(x-mu)**2/(2*sigma**2))
ax.hist(h_fit, density=True, bins=40, color='black', label='h')


h_linspace = np.linspace(np.min(h_fit), np.max(h_fit), len(h_fit))
mu = np.average(h_fit)
sigma  = np.sqrt(1/(len(h_fit) - 1)*sum((h_fit-mu)**2))
fit = gauss(h_linspace, 1/(sigma*np.sqrt(2*np.pi)), mu, sigma)
ax.plot(h_linspace, fit, color='red', label='Gaussian fit')
ax.legend()
"""
1 and 2 sigma confidence region
"""
fig, ax = plt.subplots(figsize=(6, 4), tight_layout=True)
ax.scatter(OmegaM_fit, OmegaL_fit, label='All points')
ax.scatter(OmegaM_fit[chi2_fit < chi2_fit[least_chi2] + 8.02], OmegaL_fit[chi2_fit < chi2_fit[least_chi2] + 8.02], label=r'2$\sigma$')
ax.scatter(OmegaM_fit[chi2_fit < chi2_fit[least_chi2] + 3.53], OmegaL_fit[chi2_fit < chi2_fit[least_chi2] + 3.53], label=r'1$\sigma$')
ax.scatter(best_OmegaM, best_OmegaL, label='Best fit')
ax.set_xlim(0, 1)
ax.set_ylim(0, 1.4)
ax.set_ylabel(r'$\Omega_\Lambda$')
ax.set_xlabel(r'$\Omega_M$')
ax.legend()
plt.show()