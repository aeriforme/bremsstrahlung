import numpy as np 
import matplotlib.pyplot as plt
from scipy.constants import e, m_e, c, alpha, hbar, pi, cos, sin, sqrt


Z = 29
r0 = e**2 / (m_e * c**2) 
r_c = hbar / (m_e * c) 


def k_dsigma_dk(t, t0, f, g0, k):

    p0 = sqrt(g0**2-1.) 
        
    g =  k + g0 
    p = sqrt(g**2-1.) 
    
    q = p**2 + p0**2 + k**2 - 2*p0*k*cos(t0) + 2*p*k*cos(t) - 2*p*p0*(cos(t)*cos(t0) + sin(t)*sin(t0)*cos(f))
    
    term1 = p**2 * sin(t)**2 / (g - p*cos(t))**2 *(4*g0**2 - q**2) 
    term2 = p0**2 * sin(t0)**2 / (g0 - p0*cos(t0))**2 * (4*g**2 - q**2)
    term3 = -2*p*p0*sin(t)*sin(t0)*cos(f)*(4*g*g0-q**2) / (g - p*cos(t)) / (g0 - p0*cos(t0))
    term4 = 2*k**2*(p**2*sin(t)**2 + p0**2*sin(t0)**2 - 2*p*p0*sin(t)*sin(t0)*cos(f)) / (g - p*cos(t)) / (g0 - p0*cos(t0))

    curly = term1 + term2 + term3 + term4
 
    FT_V = 4*pi*Z*e*(hbar*c)**2 / q**2 

    one_minus_F = FT_V * q**2 / (4*pi*e*Z*h**2*c**2)

    return Z**2 * alpha * r0**2 / 4 / pi**2 * (1. - F)**2 * p / p0 / q**4 

#    q = np.sqrt(p0**2+p**2-2*p*p0*(np.cos(t0)*np.cos(t)+np.sin(t0)*np.sin(t)*np.cos(f)))
#    curly = 4*p**2*np.sin(t)**2 + 4*p0**2*np.sin(t0)**2 - 8*p*p0*np.sin(t)*np.sin(t0)*np.cos(f)
#    dvol = np.sin(t)*np.sin(t0)
#    const1 = (4.*pi*e*Z*hbar**2*c**2) 
#    const2 =  4*pi*Z*e*hbar**2*c**2
#    FT_V = const2 / q**2
#    one_minus_F = FT_V * q**2 / const1
#    const3 = Z**2*alpha*e**4/(2.*pi)
#    return one_minus_F**2 * p / (p0 * q**4) * curly * dvol
    
    
t = np.linspace(0,pi,100)
t0 = np.linspace(0,pi,100)
f = np.linspace(0,2*pi,100)

T,T0,F = np.meshgrid(t,t0,f)

g0 = 1e6*eV / (m_e*c**2) 

x = []
y = []

for k in np.linspace(0.0,g0-1,100, endpoint=False):    
    x.append( k / (g0-1) )
    y.append( np.trapz(np.trapz(np.trapz(k_dsigma_dk(T,T0,F,g0,k), axis=0),axis=0),axis=0) ) 

plt.plot(x,y)
plt.show()






