import matplotlib.pyplot as plt
import numpy as np
import astropy.constants as const
import astropy.units as units
from scipy.optimize import curve_fit
#plt.rcParams['text.usetex'] = True

file = np.loadtxt("../recombination.txt")

x                     = file[:,0]
Xe_of_x               = file[:,1]
ne_of_x               = file[:,2]
tau_of_x              = file[:,3]
dtaudx_of_x           = file[:,4]
ddtauddx_of_x         = file[:,5]
g_tilde_of_x          = file[:,6]
dgdx_tilde_of_x       = file[:,7]
ddgddx_tilde_of_x     = file[:,8]


fig, ax = plt.subplots()
ax.plot(x, tau_of_x, label=r"$\tau$(x)")
ax.plot(x, -dtaudx_of_x, label=r"$\tau^{\prime}$(x)")
ax.plot(x, ddtauddx_of_x, label=r"$\tau^{\prime\prime}$(x)")
ax.set_xlabel("ln(a)")
ax.set_xlim(-12, 0)
ax.set_yscale('log')
plt.legend()

fig, ax = plt.subplots()
ax.plot(x, Xe_of_x, label=r"$X_{e}$(x)")
ax.set_xlabel("ln(a)")
ax.set_xlim(-12, 0)
ax.set_yscale('log')
plt.legend()

fig, ax = plt.subplots(2, 1, sharex=True)
ax[0].plot(x, g_tilde_of_x, label=r"$\tilde{g}$(x)")
ax[0].set_xlim(-12, 0) 
ax[0].legend()
ax[1].plot(x, dgdx_tilde_of_x, label=r"$\tilde{g^{\prime}}$(x)")
ax[1].set_xlabel("ln(a)")
ax[1].set_xlim(-12, 0) 
ax[1].legend()
plt.show()