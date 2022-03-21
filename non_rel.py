import numpy as np 
import matplotlib.pyplot as plt
from scipy.constants import c, hbar, e, m_e, alpha, epsilon_0, pi, k as k_b, eV

# constants 
r_e = e**2 / (4*pi*epsilon_0*m_e*c**2)
r_c = hbar / (m_e * c) 

def non_rel_sdcs(E1 = 100*1e3*eV, Z = 29, Z_star = 0, n_at = 8.47e+28, T = 0, corr = True):

    # electron in
    gamma1 = E1 / (m_e*c**2) + 1.  
    p1 = np.sqrt(gamma1**2-1) 
    beta1 = p1/gamma1

    # photon 
    k = np.linspace(0.0,gamma1-1,1000, endpoint=False)  # energy normalized to m_e*c**2


    # electron out 
    gamma2 = gamma1 - k   
    p2 = np.sqrt(gamma2**2-1)
    beta2 = np.sqrt(gamma2**2-1)/gamma2


    # max min mom transfer to nucleus
    dp_plus = p1 + p2 # 2*p1 - k 
    dp_minus = p1 - p2 #k 



    # material
    n_i = n_at
    n_e = Z_star * n_i  

    # thomas fermi length 
    L_tf = 4*pi*epsilon_0*hbar**2/m_e/e**2*Z**(-1./3)

    # debye length 
    if Z_star ==0 or T==0 : 
        L_d = 0.
    else: L_d = np.sqrt(epsilon_0*T / (e**2*n_i*Z_star*(Z_star+1)))

    eta_tf = 1 - Z_star / Z
    eta_d = Z_star / Z 

    # differential cross-section 

    def g(dp_plus, dp_minus, eta, LL): 
        num_log = dp_plus**2 * LL**2 + 1
        den_log = dp_minus**2 * LL**2 + 1
        return 0.5*eta**2 * ( np.log( num_log / den_log ) + (dp_plus**2 * LL**2 + 1)**(-1) - (dp_minus**2 * LL**2 +1)**(-1) ) 


    LL_d = L_d / r_c
    LL_tf = L_tf / r_c 
    arg1 = ( dp_minus**2 * LL_d**2 + 1 ) / ( dp_plus**2 * LL_d**2 + 1 )
    arg2 = ( dp_plus**2 * LL_tf**2 + 1 ) / ( dp_minus**2 * LL_tf**2 + 1 )

    Gamma_c = eta_tf * eta_d / (LL_d**2 - LL_tf**2) * ( LL_tf**2 * np.log(arg1) + LL_d**2 * np.log(arg2) )

    term1 = g(dp_plus, dp_minus, eta_tf, LL_tf)
    term2 = g(dp_plus, dp_minus, eta_d, LL_d)
    term3 = Gamma_c
    term = term1 + term2 + term3


    k_dsigma_dk = 16. * r_e**2 * Z**2 * alpha / (3. * p1**2) * term
    
    if corr==True:
        f_E = beta1 / beta2 * ( 1. - np.exp(-2*pi*Z*alpha/beta1) ) / ( 1. - np.exp(-2*pi*Z*alpha/beta2) )
    else: f_E = 1.
    
    
    return k/(gamma1-1), f_E*k_dsigma_dk 


