import numpy as np 
import matplotlib.pyplot as plt
from scipy.constants import c, hbar, e, m_e, alpha, epsilon_0, pi, k as k_b, eV
from scipy.optimize import fsolve
from scipy.special import lambertw
# constants 
r_e = e**2 / (4*pi*epsilon_0*m_e*c**2)
r_c = hbar / (m_e * c) 

def g4_SB(E1 = 5*1e6*eV, Z = 29):

    fname = 'data/g4_br'+str(Z)
    # photon energy / electron energy - 32 values 
    k_T1 = np.genfromtxt(fname, skip_header=1, skip_footer=58)
    k_T1 = np.reshape(k_T1, (1,32))

    # electron kinetic energy - 57 values
    T1 = np.genfromtxt(fname, skip_header=2, skip_footer=57)
    T1 = np.exp(T1)  # MeV
    T1 = np.reshape(T1, (57,1))
    gamma1 = T1 / 0.511 + 1.
    beta1 = np.sqrt(gamma1**2-1)/gamma1
    
    k_mat = np.matmul(T1,k_T1) / 0.511 # photon energy [MeV] / mc2

    # differential cross section = dsigma / dk * beta**2 / Z**2 * k - 57 x 32 values
    tab = np.genfromtxt(fname, skip_header=3, skip_footer=0)

    k_dsigma_dk = Z**2 * (tab / beta1**2) # millibarn 
    k_gamma = (k_mat / (gamma1-1))

    index=np.argmin(np.abs(T1 - E1 / (1e6*eV))**2)
    print('************************', T1[index])


    return k_gamma[index,:], k_dsigma_dk[index,:]*1e-28*1e-3    


