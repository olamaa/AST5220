import matplotlib.pyplot as plt
import numpy as np
import astropy.constants as const
import astropy.units as units
from scipy.optimize import curve_fit
#plt.rcParams['text.usetex'] = True

file_01 = np.loadtxt("../perturbations_k_0_1.txt")
file_001 = np.loadtxt("../perturbations_k_0_01.txt")
file_0001 = np.loadtxt("../perturbations_k_0_001.txt")

x                                                        = file_01[:,0]
Theta_0_k01,   Theta_0_k001,   Theta_0_k0001             = file_01[:,1],    file_001[:,1],     file_0001[:,1]
Theta_1_k01,   Theta_1_k001,   Theta_1_k0001             = file_01[:,2],    file_001[:,2],     file_0001[:,2]
Theta_2_k01,   Theta_2_k001,   Theta_2_k0001             = file_01[:,3],    file_001[:,3],     file_0001[:,3]
Phi_k01,       Phi_k001,       Phi_k0001                 = file_01[:,4],    file_001[:,4],     file_0001[:,4]
Psi_k01,       Psi_k001,       Psi_k0001                 = file_01[:,5],    file_001[:,5],     file_0001[:,5]
delta_cdm_k01, delta_cdm_k001, delta_cdm_k0001           = file_01[:,6],    file_001[:,6],     file_0001[:,6]
v_cdm_k01,     v_cdm_k001,     v_cdm_k0001               = file_01[:,7],    file_001[:,7],     file_0001[:,7]
delta_b_k01,   delta_b_k001,   delta_b_k0001             = file_01[:,8],    file_001[:,8],     file_0001[:,8]
v_b_k01,       v_b_k001,       v_b_k0001                 = file_01[:,9],    file_001[:,9],     file_0001[:,9]
eta_k01,       eta_k001,       eta_k0001                 = file_01[:,10],   file_001[:,10],    file_0001[:,10]

fig, ax = plt.subplots()
ax.set_title(r"$\delta_{CDM}, \delta_b$")
ax.plot(x, delta_cdm_k01, label=r"k = 0.1", color="green")
ax.plot(x, delta_cdm_k001, label=r"k = 0.01", color="red")
ax.plot(x, delta_cdm_k0001, label=r"k = 0.001", color="blue")

ax.plot(x, delta_b_k01, linestyle="dashed", color="green")
ax.plot(x, delta_b_k001, linestyle="dashed", color="red")
ax.plot(x, delta_b_k0001, linestyle="dashed", color="blue")

ax.set_xlabel(r"$x = ln(a)$")
ax.set_xlim(-20, 0)
ax.set_yscale('log')
plt.legend()


fig, ax = plt.subplots()
ax.set_title(r"$v_{CDM}, v_b$")
ax.plot(x, v_cdm_k01, label=r"k = 0.1", color="green")
ax.plot(x, v_cdm_k001, label=r"k = 0.01", color="red")
ax.plot(x, v_cdm_k0001, label=r"k = 0.001", color="blue")

ax.plot(x, v_b_k01, linestyle="dashed", color="green")
ax.plot(x, v_b_k001, linestyle="dashed", color="red")
ax.plot(x, v_b_k0001, linestyle="dashed", color="blue")

ax.set_xlabel("ln(a)")
ax.set_xlim(-20, 0)
ax.set_yscale('log')
plt.legend()


fig, ax = plt.subplots()
ax.set_title(r"$\Theta_0$")
ax.plot(x, Theta_0_k01, label=r"k = 0.1", color="green")
ax.plot(x, Theta_0_k001, label=r"k = 0.01", color="blue")
ax.plot(x, Theta_0_k0001, label=r"k = 0.001", color="red")

ax.set_xlabel("ln(a)")
ax.set_xlim(-20, 0)
plt.legend()


fig, ax = plt.subplots()
ax.set_title(r"$\Theta_1$")
ax.plot(x, Theta_1_k01, label=r"k = 0.1", color="green")
ax.plot(x, Theta_1_k001, label=r"k = 0.01", color="blue")
ax.plot(x, Theta_1_k0001, label=r"k = 0.001", color="red")

ax.set_xlabel("ln(a)")
ax.set_xlim(-20, 0)
plt.legend()


fig, ax = plt.subplots()
ax.set_title(r"$\Phi$")
ax.plot(x, Phi_k01, label=r"k = 0.1", color="green")
ax.plot(x, Phi_k001, label=r"k = 0.01", color="blue")
ax.plot(x, Phi_k0001, label=r"k = 0.001", color="red")
ax.set_xlabel("ln(a)")
ax.set_xlim(-20, 0)
#ax.set_yscale('log')
plt.legend()



fig, ax = plt.subplots()
ax.set_title(r"$\Phi + \Psi$")
ax.plot(x, Phi_k01+Psi_k01, label=r"k = 0.1", color="green")
ax.plot(x, Phi_k001+Psi_k001, label=r"k = 0.01", color="blue")
ax.plot(x, Phi_k0001+Psi_k0001, label=r"k = 0.001", color="red")
ax.set_xlabel("ln(a)")
ax.set_xlim(-20, 0)
#ax.set_yscale('log')
plt.legend()

"""
fig, ax = plt.subplots(2, 1, sharex=True)
ax[0].plot(x, g_tilde_of_x, label=r"$\tilde{g}$(x)")
ax[0].set_xlim(-12, 0) 
ax[0].legend()
ax[1].plot(x, dgdx_tilde_of_x, label=r"$\tilde{g^{\prime}}$(x)")
ax[1].set_xlabel("ln(a)")
ax[1].set_xlim(-12, 0) 
ax[1].legend()"""
plt.show()