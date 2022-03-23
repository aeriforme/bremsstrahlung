import numpy as np 
import matplotlib.pyplot as plt
from scipy.constants import c, hbar, e, m_e, alpha, epsilon_0, pi, k as k_b, eV

from non_rel import non_rel_sdcs
from mod_rel import mod_rel_sdcs
from ult_rel import ult_rel_sdcs

from three_b import three_bn_a, three_bn_b, three_bn
from g4 import g4_SB
amc = 1.6605390666e-27 # kg - atomic mass constant

#_________________________________________________________________________________
# copper 
Z = 29 
mass_dens = 8.94e3 # kg / m**3
atom_weigh = 63.546 * amc # kg 
n_at = mass_dens / atom_weigh 


# all in 1 

ele_kin_en = (100*1e3*eV, 5*1e6*eV, 100*1e6*eV) # 3 kinetic energies for 3 regimes
labels = ('non-relativistic', 'moderately-relativistic', 'ultra-relativistic')
func = (non_rel_sdcs, mod_rel_sdcs, ult_rel_sdcs)

for i in range(3):
    E1 = ele_kin_en[i]

    fig, axs = plt.subplots(1,3,figsize=(3200./300., 1200./300.), dpi = 300.) 

    Z_star = 0 
    T = 0 

    # neutral medium - with correction 
    x, y = func[i](E1 = E1, Z = Z, Z_star = Z_star, n_at = n_at, T = T)
    axs[0].plot(x, y, lw = 2, label = 'TFD')

    # 3bn
    x, y = three_bn(E1 = E1, Z = Z)
    axs[0].plot(x, y, lw = 2, label = '3BN')

    # seltzer-berger
    x, y = g4_SB(E1 = E1, Z = Z)
    axs[0].plot(x, y, lw = 2,  label='SB')

    axs[0].set_xlabel(r'k / ($\gamma_1$ -1)')   
    axs[0].set_ylabel(r'k d$\sigma$/dk')
    axs[0].set_title(r'E$_1$=%.1f MeV, Z*=%d, T=%d keV' % (E1/(1e6*eV), Z_star, T), pad=15) 
    axs[0].set_xlim(0,1)
    axs[0].legend()


    # fully ionized medium  
    Z_star = Z
    for T in (0.01*1e3*eV, 0.1*1e3*eV, 10*1e3*eV, 100*1e3*eV): 
        x, y = func[i](E1 = E1, Z = Z, Z_star = Z_star, n_at = n_at, T = T)
        axs[1].plot(x, y, lw = 2, label = 'T=%.2f keV' % (T/1e3/eV))

    axs[1].set_xlabel(r'k / ($\gamma_1$ -1)')
    axs[1].set_ylabel(r'k d$\sigma$/dk')
    axs[1].set_title(r'E$_1$=%.1f MeV, Z*=%d' % (E1/(1e6*eV), Z_star), pad=15) 
    axs[1].set_xlim(0,1)
    axs[1].legend()


    # hot medium 
    T = 10.*1e3*eV 
    for Z_star in (int(0.5*Z),int(0.75*Z),Z): 
        x, y = func[i](E1 = E1, Z = Z, Z_star = Z_star, n_at = n_at, T = T)
        axs[2].plot(x, y, lw = 2, label = 'Z*=%01d' % (Z_star))

    axs[2].set_xlabel(r'k / ($\gamma_1$ -1)')
    axs[2].set_ylabel(r'k d$\sigma$/dk')
    axs[2].set_title(r'E$_1$=%.1f MeV, T=%.0f keV' % (E1/(1e6*eV),T/(1e3*eV)), pad=15) 
    axs[2].set_xlim(0,1)
    axs[2].legend()


    fig.subplots_adjust(top=0.8)
    fig.suptitle('particle = %s electron | target = copper' % labels[i])
    plt.tight_layout()
    fig.savefig('%s.png' % labels[i]) 
    plt.close('all')









#_________________________________________________________________________________
# non relativistic 
fig, axs = plt.subplots(1,3,figsize=(3200./300., 1200./300.), dpi = 300.) 

# incident electron 
E1 = 100*1e3*eV # kinetic energy 
Z_star = 0 
T = 0 

# neutral medium - no correction 
x, y = non_rel_sdcs(E1 = E1, Z = Z, Z_star = Z_star, n_at = n_at, T = T, corr = False)
axs[0].plot(x, y, lw = 2, label = 'TFD') 

# neutral medium - with correction 
x, y = non_rel_sdcs(E1 = E1, Z = Z, Z_star = Z_star, n_at = n_at, T = T, corr = True)
axs[0].plot(x, y, lw = 2, label = 'TFD + elwert')

# 3bn
x, y = three_bn(E1 = E1, Z = Z)
axs[0].plot(x, y, lw = 2, label = '3BN')

# seltzer-berger
x, y = g4_SB(E1 = E1, Z = Z)
axs[0].plot(x, y, lw = 2,  label='SB')

axs[0].set_xlabel(r'k / ($\gamma_1$ -1)')
axs[0].set_ylabel(r'k d$\sigma$/dk')
axs[0].set_title(r'E$_1$=%.1f MeV, Z*=%d, T=%d keV' % (E1/(1e6*eV), Z_star, T), pad=15) 
axs[0].set_xlim(0,1)
#axs[0].set_ylim(0,6e-27)
axs[0].legend()


# fully ionized medium  
Z_star = 29 
for T in (0.01*1e3*eV, 0.1*1e3*eV, 10*1e3*eV, 100*1e3*eV): 
    x, y = non_rel_sdcs(E1 = E1, Z = Z, Z_star = Z_star, n_at = n_at, T = T, corr = True)
    axs[1].plot(x, y, lw = 2, label = 'T=%.2f keV' % (T/1e3/eV))

axs[1].set_xlabel(r'k / ($\gamma_1$ -1)')
axs[1].set_ylabel(r'k d$\sigma$/dk')
axs[1].set_title(r'E$_1$=%.1f MeV, Z*=%d' % (E1/(1e6*eV), Z_star), pad=15) 
axs[1].set_xlim(0,1)
#axs[1].set_ylim(0,6e-27)
axs[1].legend()


# hot medium 
T = 10.*1e3*eV 
for Z_star in (int(0.5*Z),int(0.75*Z),Z): 
    x, y = non_rel_sdcs(E1 = E1, Z = Z, Z_star = Z_star, n_at = n_at, T = T, corr = True)
    axs[2].plot(x, y, lw = 2, label = 'Z*=%01d' % (Z_star))

axs[2].set_xlabel(r'k / ($\gamma_1$ -1)')
axs[2].set_ylabel(r'k d$\sigma$/dk')
axs[2].set_title(r'E$_1$=%.1f MeV, T=%.0f keV' % (E1/(1e6*eV),T/(1e3*eV)), pad=15) 
axs[2].set_xlim(0,1)
#axs[2].set_ylim(0,6e-27)
axs[2].legend()


fig.subplots_adjust(top=0.8)
fig.suptitle('particle = non relativistic electron | target = copper')
plt.tight_layout()
fig.savefig('nr.png') 

#_________________________________________________________________________________
# mod rel electrons
fig, axs = plt.subplots(1,3,figsize=(3200./300., 1200./300.), dpi = 300.) 

E1 = 5*1e6*eV
Z_star = 0
T = 0 

# tfd
x, y = mod_rel_sdcs(E1 = E1, Z = Z, Z_star = Z_star, n_at = n_at, T = T)
axs[0].plot(x, y, lw = 2, label = 'TFD')

# 3bn
x, y = three_bn(E1 = E1, Z = Z)
axs[0].plot(x, y, lw = 2, label = '3BN')

x, y = three_bn_b(E1 = E1, Z = Z)
axs[0].plot(x, y, lw = 2, label = '3BNb')

# seltzer-berger
x, y = g4_SB(E1 = E1, Z = Z)
axs[0].plot(x, y, lw = 2,  label='geant4')


axs[0].set_xlabel(r'k / ($\gamma_1$ -1)')
axs[0].set_ylabel(r'k d$\sigma$/dk')
axs[0].set_title(r'E$_1$=%.0f MeV, Z*=%d, T=%d keV' % (E1/(1e6*eV), Z_star, T), pad=15) 
axs[0].set_xlim(0,1)
#axs[0].set_ylim(0,2.5e-27)
axs[0].legend()



# fully ionized medium 
Z_star = Z 
for T in (0.01*1e3*eV, 0.1*1e3*eV, 10*1e3*eV, 100*1e3*eV): 
    x, y = mod_rel_sdcs(E1 = E1, Z = Z, Z_star = Z_star, n_at = n_at, T = T)
    axs[1].plot(x, y, lw = 2, label = 'T = %.0f keV' % (T/1e3/eV))

axs[1].set_xlabel(r'k / ($\gamma_1$ -1)')
axs[1].set_ylabel(r'k d$\sigma$/dk')
axs[1].set_title(r'E$_1$=%.0f MeV, Z*=%d' % (E1/(1e6*eV), Z_star), pad=15) 
axs[1].set_xlim(0,1)
#axs[1].set_ylim(0,2.5e-27)
axs[1].legend()



# hot medium 
T = 10.*1e3*eV 
for Z_star in (int(0.5*Z),int(0.75*Z),Z): 
    x, y = mod_rel_sdcs(E1 = E1, Z = Z, Z_star = Z_star, n_at = n_at, T = T)
    axs[2].plot(x, y, lw = 2, label = 'Z*=%01d' % (Z_star))

axs[2].set_xlabel(r'k / ($\gamma_1$ -1)')
axs[2].set_ylabel(r'k d$\sigma$/dk')
axs[2].set_title(r'E$_1$=%.0f MeV, T=%.0f keV' % (E1/(1e6*eV), T/(1e3*eV)), pad=15) 
axs[2].set_xlim(0,1)
#axs[2].set_ylim(0,2.5e-27)
axs[2].legend()


fig.subplots_adjust(top=0.8)
fig.suptitle('particle = moderately relativistic electron | target = copper')
plt.tight_layout()
fig.savefig('mr.png') 





#_________________________________________________________________________________
# ult rel electrons
fig, axs = plt.subplots(1,3,figsize=(3200./300., 1200./300.), dpi = 300.) 

E1 = 100*1e6*eV
Z_star = 0
T = 0 

x, y = ult_rel_sdcs(E1 = E1, Z = Z, Z_star = Z_star, n_at = n_at, T = T)
axs[0].plot(x, y, lw = 2, label = 'TFD')

x, y = three_bn(E1 = E1, Z = Z)
axs[0].plot(x, y, lw = 2, label = '3BN')

x, y = three_bn_b(E1 = E1, Z = Z)
axs[0].plot(x, y, lw = 2, label = '3BNb')

x, y = g4_SB(E1 = E1, Z = Z)
axs[0].plot(x, y, lw = 2,  label='geant4')


axs[0].set_title('Z* = %01d, T = %.0f keV' % (Z_star, T/1e3/eV), pad=15)
axs[0].legend()
axs[0].set_xlabel(r'k / ($\gamma_1$ -1)')
axs[0].set_ylabel(r'k d$\sigma$/dk')
#axs[0].set_ylim(0, 3.5e-27)


# fully ionized medium 
Z_star = Z 
for T in (0.01*1e3*eV, 0.1*1e3*eV, 10*1e3*eV, 100*1e3*eV): 
    x, y = ult_rel_sdcs(E1 = E1, Z = Z, Z_star = Z_star, n_at = n_at, T = T)
    axs[1].plot(x, y, lw = 2, label = 'T = %.2f keV' % (T/1e3/eV))

axs[1].set_xlabel(r'k / ($\gamma_1$ -1)')
axs[1].set_ylabel(r'k d$\sigma$/dk')
axs[1].set_title(r'E$_1$=%.0f MeV, Z*=%d' % (E1/(1e6*eV), Z_star), pad=15) 
axs[1].set_xlim(0,1)
#axs[1].set_ylim(0,3.5e-27)
axs[1].legend()


# hot medium 
T = 10.*1e3*eV 
for Z_star in (int(0.5*Z),int(0.75*Z),Z): 
    x, y = ult_rel_sdcs(E1 = E1, Z = Z, Z_star = Z_star, n_at = n_at, T = T)
    axs[2].plot(x, y, lw = 2, label = 'Z*=%01d' % (Z_star))

axs[2].set_xlabel(r'k / ($\gamma_1$ -1)')
axs[2].set_ylabel(r'k d$\sigma$/dk')
axs[2].set_title(r'E$_1$=%.0f MeV, T=%.0f keV' % (E1/(1e6*eV), T/(1e3*eV)), pad=15) 
axs[2].set_xlim(0,1)
#axs[2].set_ylim(0,3.5e-27)
axs[2].legend()

fig.subplots_adjust(top=0.8)
fig.suptitle('particle = ultra relativistic electron | target = copper')
plt.tight_layout()
fig.savefig('ur.png') 



