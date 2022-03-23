import numpy as np 
import matplotlib.pyplot as plt
from scipy.constants import c, hbar, e, m_e, alpha, epsilon_0, pi, k as k_b, eV
from scipy.optimize import fsolve
from scipy.special import lambertw, zeta
# constants 
r_e = e**2 / (4*pi*epsilon_0*m_e*c**2)
r_c = hbar / (m_e * c) 

def ult_rel_sdcs(E1 = 100*1e6*eV, Z = 29, Z_star = 0, n_at = 8.47e+28, T = 0):

    # electron in
    gamma1 = E1 / (m_e*c**2) + 1.  
    p1 = np.sqrt(gamma1**2-1) 
    beta1 = p1/gamma1

    # photon 
    k = np.linspace(0.00001,gamma1-1,1000, endpoint=False)  # energy normalized to m_e*c**2

    # electron out 
    gamma2 = gamma1 - k   
    beta2 = np.sqrt(gamma2**2-1)/gamma2

    # material
    n_i = n_at
    n_e = Z_star * n_i  

    # thomas fermi length 
    L_tf = 4*pi*epsilon_0*hbar**2/m_e/e**2*Z**(-1./3)

    # debye length 
    if Z_star ==0 or T==0 : 
        L_d = 0.
    else: L_d = np.sqrt(epsilon_0*T / (e**2*n_i*Z_star*(Z_star+1)))

    LL_d = L_d / r_c
    LL_tf = L_tf / r_c 


    q = -e         
    q_tf = q * (1 - Z_star / Z) 
    q_d = q * Z_star / Z 

    
    
    def I1(d, L): 
        term1 = L * d * (np.arctan(d*L) - np.arctan(L))
        term2 = -0.5*L**2 * (1.-d)**2 / (1.+L**2) 
        term3 = 0.5*np.log( (1.+L**2) / (1.+(L*d)**2) )
        return term1 + term2 + term3

    def I2(d, L):
        term1 = 4*L**3 * d**3 * (np.arctan(L*d) - np.arctan(L))
        term2 = (1+3*L**2 * d**2)*np.log( (1+L**2) / (1+L**2*d**2) )
        term3 = 6* L**4 * d**2 *np.log(d) / (1+L**2)
        term4 = L**2 * (d-1)*(d+1-4*L**2 * d**2)/(1+L**2)
        return 0.5*(term1+term2+term3+term4) 
        
        
    delta = 0.5 * k / gamma1 / (gamma1-k)    
    #print('................................................', delta) 
    def I_tfd():
        num = (1+LL_tf**2)*np.log(1+LL_tf**2)-LL_tf**2
        den = 1 + LL_tf**2
        term1 = 0.5*q_tf**2/q**2 * num / den 

        num = (1+LL_d**2)*np.log(1+LL_d**2)-LL_d**2
        den = 1 + LL_d**2
        term2 = 0.5*q_d**2/q**2 * num / den   

        num = LL_d**2*np.log(1+LL_tf**2)-LL_tf**2*np.log(1+LL_d**2) 
        den = LL_d**2 - LL_tf**2
        term3 = q_tf*q_d/q**2 * num / den       
       
        return term1 + term2 + term3
        
    
    def LL_r():
        a = I_tfd()
        W = lambertw(-np.exp(-1-2*a))
        if np.imag(W)==0:
            W = np.real(W)
        else: print('errorrrrrrrrrrrrrrrrrrrr') 
        return np.sqrt( np.exp(W + 1 + 2*a) -1)
    

    def f_c(Z): 
        
        #return 0. #(Z*alpha)**2
        s = 0
        for i in range(1,5): 
            #s = s + ( i * (i**2+alpha**2) )**(-1)
            s = s + (-alpha**2 * Z**2)**i * ( zeta(2*i+1) - 1.)    
            #print(i, '-----', zeta(2*i+1, q=1))
        return alpha**2*Z**2/(1+alpha**2*Z**2) * s 
    
    term1 = 1. + ( (gamma1-k)/gamma1 )**2
    term2 = I1(delta, LL_r()) + 1. - f_c(Z)
    
    term3 = -2./3.*(gamma1-k)/gamma1
    term4 = I2(delta, LL_r()) + 5./6. - f_c(Z) 
   

    #print(term2)
    #plt.plot(term1*term2)
    #plt.plot(term3*term4)
    #plt.show()
    
    term = term1*term2 + term3*term4 
    k_dsigma_dk = 4 * Z**2 * r_e**2 * alpha * term 
    
    
    
    return k/(gamma1-1), k_dsigma_dk 


