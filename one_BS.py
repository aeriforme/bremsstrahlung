import numpy as np 
import matplotlib.pyplot as plt
from scipy.constants import e, m_e, c, hbar, pi, eV, epsilon_0, alpha
from numpy import cos, sin, sqrt
from scipy.integrate import tplquad 

Z = 29

r_e = e**2 / (4*pi*epsilon_0*m_e*c**2)
r_c = hbar / (m_e * c) 
a_0 = 4.*pi*epsilon_0*hbar**2/(m_e*e**2)


from time import perf_counter






def one_BS_num(E1 = 100*1e3*eV, Z = 29, Z_star = 0, n_at = 8.47e+28, T = 0, N = 50):

    # material
    n_i = n_at
    n_e = Z_star * n_i  

    # thomas fermi length 
    L_tf = a_0 * Z**(-1./3)

    # debye length 
    if Z_star ==0 or T==0 : 
        L_d = 0.
    else: L_d = np.sqrt(epsilon_0*T / (e**2*n_i*Z_star*(Z_star+1)))


    LL_tf = L_tf / r_c
    LL_d  = L_d  / r_c


    g0 = E1 / (m_e*c**2) + 1.
    p0 = sqrt(g0**2-1.) 

    k_vec = np.linspace(0, g0-1, N, endpoint=False)
    x = k_vec / (g0-1) 
    


    def k_dsigma_dk(t, t0, f, g0, k): 
        g =  g0 - k 
        p = sqrt(g**2-1.) 
    
        q2 = p**2 + p0**2 + k**2 - 2*p0*k*cos(t0) + 2*p*k*cos(t) - 2*p*p0*(cos(t)*cos(t0) + sin(t)*sin(t0)*cos(f))
    

        term1 = p**2  * sin(t)**2  / (g  - p*cos(t)  )**2   * (4*g0**2 - q2) 
        term2 = p0**2 * sin(t0)**2 / (g0 - p0*cos(t0))**2 * (4*g**2  - q2)
        term3 = -2*p*p0*sin(t)*sin(t0)*cos(f)*(4*g*g0-q2+2*k**2) / (g - p*cos(t)) / (g0 - p0*cos(t0))
        term4 = 2*k**2*(p**2*sin(t)**2 + p0**2*sin(t0)**2) / (g - p*cos(t)) / (g0 - p0*cos(t0))

        curly = term1 + term2 + term3 + term4
    
        if L_d == 0: 
            FT_V = e*Z / (4*pi*epsilon_0*r_c) *   (1. - Z_star / Z) * LL_tf**2 / ( q2*LL_tf**2 + 1 ) 
        else: 
            FT_V = e*Z / (4*pi*epsilon_0*r_c) * ( (1. - Z_star / Z) * LL_tf**2 / ( q2*LL_tf**2 + 1 ) + Z_star/Z * LL_d**2 / ( q2*LL_d**2 + 1 ) )
    
        one_minus_F = FT_V * q2 * 4*pi*epsilon_0*r_c/(e*Z) 
        return Z**2 * alpha * r_e**2 /(2*pi) * one_minus_F**2 * curly * p / p0 / q2**2 * sin(t) * sin(t0)  

    t  = np.linspace(0, pi,   N)
    t0 = np.linspace(0, pi,   N)
    f  = np.linspace(0, 2*pi, 2*N)


    T,T0,F = np.meshgrid(t,t0,f)



    y1 = []
    y2 = []
    
    start = perf_counter()

    for k in k_vec:    
        y1.append( np.trapz(np.trapz(np.trapz(k_dsigma_dk(T,T0,F,g0,k), t, axis=0), t0, axis=0), f, axis=0) ) 

    end = perf_counter() 
    execution_time = (end - start)
    print('trapz time=', execution_time)
    
    
    start = perf_counter()
    for k in k_vec:
        y2.append( tplquad(lambda z,y,x: k_dsigma_dk(x,y,z,g0,k), 0, pi, lambda x: 0, lambda x:  pi, lambda x,y: 0, lambda x,y: 2*pi)[0] )

    end = perf_counter() 
    execution_time = (end - start)
    print('tplquad time=', execution_time)


    return x, y1, y2 


    






