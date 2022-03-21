import numpy as np 
import matplotlib.pyplot as plt
from scipy.constants import c, hbar, e, m_e, alpha, epsilon_0, pi, k as k_b, eV

from non_rel import non_rel_sdcs
from mod_rel import mod_rel_sdcs
from ult_rel import ult_rel_sdcs

from threebn import three_bn_a, three_bn_b
from g4 import g4_SB
amc = 1.6605390666e-27 # kg - atomic mass constant

#_________________________________________________________________________________
# copper 
Z = 29 
mass_dens = 8.94e3 # kg / m**3
atom_weigh = 63.546 * amc # kg 
n_at = mass_dens / atom_weigh 


#_________________________________________________________________________________
# non relativistic 
fig, axs = plt.subplots(1,2,figsize=(2000./300., 1000./300.), dpi = 300.) 

# incident electron 
E1 = 100*1e3*eV 

# neutral medium - no correction 
x, y = non_rel_sdcs(E1 = E1, Z = Z, Z_star = 0, n_at = n_at, T = 0, corr = False)
axs[0].plot(x, y, lw = 2, label = 'Z*=0, T=0') 

# neutral medium - with correction 
x, y = non_rel_sdcs(E1 = E1, Z = Z, Z_star = 0, n_at = n_at, T = 0, corr = True)
axs[0].plot(x, y, lw = 2, label = 'Z*=0, T=0, elwert')


x, y = three_bn_a(E1 = E1, Z = Z)
axs[0].plot(x, y, lw = 2, label = '3BNa')


x, y = g4_SB(E1 = E1, Z = Z)
axs[0].plot(x, y, lw = 2,  label='geant4')

#mart=np.loadtxt('martinez_nr.csv', delimiter=',')
#axs[0].plot(mart[:,0], mart[:,1]*1e-27, lw = 2, label = 'Martinez')

axs[0].set_xlabel(r'k / ($\gamma_1$ -1)')
axs[0].set_ylabel(r'k d$\sigma$/dk')
axs[0].set_title(r'E$_1$' + str(E1/eV/1e3) + 'keV') 

axs[0].set_ylim(0,6e-27)
axs[0].legend()

# ionized mediums - no correction 
for T in (0.1*1e3*eV, ): # 10*1e3*eV): 
    for Z_star in range(Z):
        x, y = non_rel_sdcs(E1 = E1, Z = Z, Z_star = Z_star, n_at = n_at, T = T, corr = False)
        axs[1].plot(x, y, lw = 2, label = 'Z* = %01d, T = %.2f keV' % (Z_star, T/1e3/eV))

#axs[1].legend()
plt.tight_layout()
fig.savefig('nr.png') 

#_________________________________________________________________________________
# mod rel electrons
fig, axs = plt.subplots(1,2,figsize=(2000./300., 1000./300.), dpi = 300.) 

E1 = 5*1e6*eV
Z_star = 0
T = 0 

x, y = mod_rel_sdcs(E1 = E1, Z = Z, Z_star = Z_star, n_at = n_at, T = T)
axs[0].plot(x, y, lw = 2, label = 'TFD')


x, y = three_bn_b(E1 = E1, Z = Z)
axs[0].plot(x, y, lw = 2, label = '3BNb')



x, y = g4_SB(E1 = E1, Z = Z)
axs[0].plot(x, y, lw = 2,  label='geant4')



axs[0].set_title('Z* = %01d, T = %.2f keV' % (Z_star, T/1e3/eV))
axs[0].legend()
axs[0].set_xlabel(r'k / ($\gamma_1$ -1)')
axs[0].set_ylabel(r'k d$\sigma$/dk')


plt.tight_layout()
fig.savefig('mr.png') 





#_________________________________________________________________________________
# mod rel electrons
fig, axs = plt.subplots(1,2,figsize=(2000./300., 1000./300.), dpi = 300.) 

E1 = 100*1e6*eV
Z_star = 0
T = 0 

x, y = ult_rel_sdcs(E1 = E1, Z = Z, Z_star = Z_star, n_at = n_at, T = T)
axs[0].plot(x, y, lw = 2, label = 'TFD')


x, y = three_bn_b(E1 = E1, Z = Z)
axs[0].plot(x, y, lw = 2, label = '3BNb')



x, y = g4_SB(E1 = E1, Z = Z)
axs[0].plot(x, y, lw = 2,  label='geant4')



axs[0].set_title('Z* = %01d, T = %.2f keV' % (Z_star, T/1e3/eV))
axs[0].legend()
axs[0].set_xlabel(r'k / ($\gamma_1$ -1)')
axs[0].set_ylabel(r'k d$\sigma$/dk')


plt.tight_layout()
fig.savefig('ur.png') 



