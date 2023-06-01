import matplotlib.pyplot as plt
import numpy as np
import astropy.constants as const
import astropy.units as units
from scipy.optimize import curve_fit
#plt.rcParams['text.usetex'] = True
def Milestone_3():
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
    Source_k01,   Source_k001,   Source_k0001                = file_01[:,11],    file_001[:,11],     file_0001[:,11]


    fig, ax = plt.subplots(2, sharex=True)
    ax[0].set_title(r"$\delta_{CDM}, \delta_b$")
    ax[1].set_title(r"$\delta_\gamma$")

    ax[0].plot(x, delta_cdm_k01, label=r"k = 0.1", color="green")
    ax[0].plot(x, delta_cdm_k001, label=r"k = 0.01", color="red")
    ax[0].plot(x, delta_cdm_k0001, label=r"k = 0.001", color="blue")
    ax[0].legend()

    ax[1].plot(x, 4*Theta_0_k01, color="green")
    ax[1].plot(x, 4*Theta_0_k001, color="red")
    ax[1].plot(x, 4*Theta_0_k0001, color="blue")

    ax[0].plot(x, delta_b_k01, linestyle="dashed", color="green")
    ax[0].plot(x, delta_b_k001, linestyle="dashed", color="red")
    ax[0].plot(x, delta_b_k0001, linestyle="dashed", color="blue")

    ax[1].set_xlabel(r"$x = ln(a)$")
    ax[1].set_xlim(-20, 0)
    ax[0].set_yscale('log')


    fig, ax = plt.subplots(2, sharex=True)
    ax[0].set_title(r"$v_{CDM}, v_b$")
    ax[1].set_title(r"$v_\gamma$")
    ax[0].plot(x, v_cdm_k01, label=r"k = 0.1", color="green")
    ax[0].plot(x, v_cdm_k001, label=r"k = 0.01", color="red")
    ax[0].plot(x, v_cdm_k0001, label=r"k = 0.001", color="blue")
    ax[0].legend()

    ax[1].plot(x, -3*Theta_1_k01, color="green")
    ax[1].plot(x, -3*Theta_1_k001, color="red")
    ax[1].plot(x, -3*Theta_1_k0001, color="blue")

    ax[0].plot(x, v_b_k01, linestyle="dashed", color="green")
    ax[0].plot(x, v_b_k001, linestyle="dashed", color="red")
    ax[0].plot(x, v_b_k0001, linestyle="dashed", color="blue")

    ax[1].set_xlabel("ln(a)")
    ax[1].set_xlim(-20, 0)
    ax[0].set_yscale('log')


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
    ax.set_title(r"$\Theta_2$")
    ax.plot(x, Theta_2_k01, label=r"k = 0.1", color="green")
    ax.plot(x, Theta_2_k001, label=r"k = 0.01", color="blue")
    ax.plot(x, Theta_2_k0001, label=r"k = 0.001", color="red")

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
    plt.show()




## MILESTONE 4
## SOURCE FUNCTION
def Milestone_4():
    file_01 = np.loadtxt("../perturbations_k_0_1.txt")
    file_001 = np.loadtxt("../perturbations_k_0_01.txt")
    file_0001 = np.loadtxt("../perturbations_k_0_001.txt")
    x                                                        = file_01[:,0]
    Source_k01,   Source_k001,   Source_k0001                = file_01[:,11],    file_001[:,11],     file_0001[:,11]
    fig, ax = plt.subplots()
    ax.set_title(r"SOURCE FUNCTION")
    ax.plot(x, Source_k01, label=r"k = 0.1", color="green")
    ax.plot(x, Source_k001, label=r"k = 0.01", color="blue")
    ax.plot(x, Source_k0001, label=r"k = 0.001", color="red")
    ax.set_xlabel("ln(a)")
    ax.set_xlim(-10, 0)
    #ax.set_yscale('log')
    plt.legend()
    plt.show()

    thetas = np.loadtxt('../thetak_ells.txt')
    cells = np.loadtxt("../cells.txt")
    source_data = np.loadtxt("../source.txt")
    Matter_power_spectrum = np.loadtxt("../matter_power_spectrum_result.txt") 
    SDSSDR7LRG_data = np.loadtxt("SDSSDR7LRG.txt")
    WMAP_data = np.loadtxt("WMAP_ACT.txt")
    lowellTTdata = np.loadtxt("../lowelltt.txt")
    highellTTdata = np.loadtxt("../highellTT.txt")




    kc_H0, k_norm, theta_ell_min1, theta_ell_min2, theta_ell_min3, theta_ell_min4, theta_ell_min5 = thetas[:,0], thetas[:,1], thetas[:,9], thetas[:,21], thetas[:,26], thetas[:,34], thetas[:,44]
    ells, values = cells[:,0], cells[:,1]
    k, P_k = Matter_power_spectrum[:,0], Matter_power_spectrum[:,1]
    SDSSDR7LRG_k, SDSSDR7LRG_Pk, SDSSDR7LRG_error = SDSSDR7LRG_data[1:,0], SDSSDR7LRG_data[1:,1], SDSSDR7LRG_data[1:,2]
    WMAP_k, WMAP_Pk, WMAP_upper = WMAP_data[1:,0], WMAP_data[1:,1], WMAP_data[1:,2]


    x_source, sourcejell, integrand = source_data[:,0], source_data[:,1], source_data[:,2]
    fig, ax = plt.subplots()
    ax.set_title(r"Source function")
    ax.plot(x_source, sourcejell)
    ax.set_xlabel("x = ln(a)")
    ax.set_ylabel(r"$\tilde{S}(k, x)j_\ell[\eta_0 - \eta(x)] / 10^{-3}$")
    #ax.set_xlim(-20, 0)
    ax.set_ylim(-1, 2)





    theta_ell_min = [theta_ell_min1, theta_ell_min2, theta_ell_min3, theta_ell_min4, theta_ell_min5]
    labels = [r"$\ell=10$", r"$\ell=100$", r"$\ell=200$", r"$\ell=500$", r"$\ell=1000$"]
    fig, ax = plt.subplots()

    for i in range(len(theta_ell_min)):
        ax.plot(kc_H0, theta_ell_min[i], label=labels[i])
        #theta_min_index = np.argwhere(theta_ell_min[i] == np.min(theta_ell_min[i]))[0][0]
        #ax.scatter(kc_H0[theta_min_index], theta_ell_min[i][theta_min_index])
    ax.set_ylim(-0.001, 0.01)
    #ax.set_xlim(0, 500)
    ax.set_xlabel(r"$k\eta_0$")
    ax.set_title(r"$\Theta_\ell$")
    plt.legend()

    fig, ax = plt.subplots()
    ell = [6, 20, 200, 400, 800, 1200]
    for i in range(len(theta_ell_min)):
        ax.plot(kc_H0, ell[i]*(ell[i]+1)*theta_ell_min[i]**2/kc_H0, label=labels[i])
    ax.set_ylim(10**(-13), 5*10**(-3))
    ax.set_xlim(0, 500)
    ax.set_xlabel(r"$k\eta_0$")
    ax.set_title(r"$\frac{\ell(\ell+1)\Theta_\ell^2}{10^6k\eta^{-1}}$")
    plt.legend()



    

    fig, ax = plt.subplots()
    ax.plot(k, P_k, label="Theoretical prediction")
    upper_error = WMAP_upper-WMAP_Pk
    ax.errorbar(WMAP_k, WMAP_upper, yerr = upper_error, uplims = True, label="WMAP + ACT", linestyle='dotted', color='orange', capsize=0)
    ax.errorbar(SDSSDR7LRG_k, SDSSDR7LRG_Pk, SDSSDR7LRG_error, label="SDSS DR7 LRG",linestyle='dotted', color='black')
    ax.set_title(r"P(k)")
    ax.set_xlabel(r"k (h/Mpc)$^3$")
    ax.set_ylabel(r"(Mpc/h)$^3$")
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim(2*10**(-3), 10**0)
    ax.legend()



    lowellTT_ell, lowellTT_best, elllolims, elluplims = lowellTTdata[:,0], lowellTTdata[:,1], lowellTTdata[:,2], lowellTTdata[:,3]
    highellTT_ell, highellTT_best = highellTTdata[:,0], highellTTdata[:,-1]
    fig, ax = plt.subplots()
    ax.fill_between(lowellTT_ell, lowellTT_best+elllolims, lowellTT_best-elluplims, alpha=0.2)
    ax.plot(lowellTT_ell, lowellTT_best, label=r"Low $\ell$ - TT")
    ax.plot(highellTT_ell, highellTT_best, label=r"High $\ell$ - TT")
    C_ell = values*2*np.pi/(ells*(ells+1))
    ax.plot(ells, values, label="Our model")
    ax.set_title(r"$\frac{\ell(\ell+1)}{2\pi} C_\ell$")
    ax.set_ylabel(r"[$\mu$K$^2$]")
    ax.set_xlabel(r"$\ell$")
    #ax.set_xlim(0, 50)
    ax.set_xscale('log')
    plt.legend()



    plt.show()


    
Milestone_4()
