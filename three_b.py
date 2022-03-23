import numpy as np 
import matplotlib.pyplot as plt
from scipy.constants import c, hbar, e, m_e, alpha, epsilon_0, pi, k as k_b, eV
from scipy.optimize import fsolve
from scipy.special import lambertw
# constants 
r_e = e**2 / (4*pi*epsilon_0*m_e*c**2)
r_c = hbar / (m_e * c) 

def three_bn(E1 = 5*1e6*eV, Z = 29):

    # in koch motz E1 is the total energy 
    E1 = E1 + m_e*c**2 
    
    # electron in
    gamma1 = E1 / (m_e*c**2) 
    p1 = np.sqrt(gamma1**2-1) 

    # photon 
    k = np.linspace(0.00001,gamma1-1,1000, endpoint=False)  # energy normalized to m_e*c**2

    # electron out 
    gamma2 = gamma1 - k   
    p2 = np.sqrt(gamma2**2-1)
    
    E1 = gamma1
    E2 = gamma2 
    
    L = 2.*np.log((E1*E2 + p1*p2 -1)/k)  
    eps0 = np.log((E1+p1)/(E1-p1)) 
    eps  = np.log((E2+p2)/(E2-p2))
    
    t1 = 8.*E1*E2/(3.*p1*p2)
    t2 = k**2 *((E1*E2)**2 + (p1*p2)**2)/(p1*p2)**3
    
    tonde = eps0 * (E1*E2+p1**2)/p1**3 - eps * (E1*E2+p2**2)/p2**3 + 2*k*E1*E2/(p1*p2)**2

    t3 = k / (2*p1*p2) * tonde
    squared =  t1 + t2 + t3  
    
    curly = 4./3. - 2*E1*E2* (p1**2+p2**2)/(p1*p2)**2 + eps0*E2/p1**3 + eps*E1/p2**3 + eps*eps0/(p1*p2) + L*squared
    k_dsigma_dk = Z**2 * r_e**2 * alpha * p2/p1 * curly 
    
    return k/(gamma1-1), k_dsigma_dk 
    
    

def three_bn_a(E1 = 5*1e6*eV, Z = 29):
    # in koch motz E1 is the total energy 
    E1 = E1 + m_e*c**2 
    
    # electron in
    gamma1 = E1 / (m_e*c**2)   
    p1 = np.sqrt(gamma1**2-1) 

    # photon 
    k = np.linspace(0.00001,gamma1-1,1000, endpoint=False)  # energy normalized to m_e*c**2

    # electron out 
    gamma2 = gamma1 - k   
    p2 = np.sqrt(gamma2**2-1)

    k_dsigma_dk = Z**2 * r_e**2 * 16. / 3. * alpha / p1**2 * np.log( (p1+p2) / (p1-p2) ) 
    
    return k/(gamma1-1), k_dsigma_dk 


def three_bn_b(E1 = 5*1e6*eV, Z = 29):

    # in koch motz E1 is the total energy 
    E1 = E1 + m_e*c**2 
    
    # electron in
    gamma1 = E1 / (m_e*c**2)  
    p1 = np.sqrt(gamma1**2-1) 

    # photon 
    k = np.linspace(0.00001,gamma1-1,1000, endpoint=False)  # energy normalized to m_e*c**2

    # electron out 
    gamma2 = gamma1 - k   
    p2 = np.sqrt(gamma2**2-1)



    term1 = 1 + (gamma2 / gamma1)**2 - 2./3. * gamma2 / gamma1
    term2 = np.log(2*gamma1*gamma2 / k ) - 0.5 


    k_dsigma_dk = 4 * Z**2 * r_e**2 * alpha * term1 * term2
    
    return k/(gamma1-1), k_dsigma_dk 

