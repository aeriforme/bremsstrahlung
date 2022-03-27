import numpy as np 
import matplotlib.pyplot as plt
from scipy.constants import e, m_e, c, alpha, hbar, pi


Z = 1
p0 = 0.1
p = 0.08

eps = 1e-6

def integrand(t, t0, q):
    cos_phi = (q**2-p0**2-p**2+2*p0*p*np.cos(t)*np.cos(t0)) / (-2*p*p0*np.sin(t)*np.sin(t0))
    curly = 4*p**2*np.sin(t)**2 + 4*p0**2*np.sin(t0)**2 - 8*p*p0*np.sin(t)*np.sin(t0)*cos_phi
    dvol = np.sin(t)*np.sin(t0)
    const1 = (4.*pi*e*Z*hbar**2*c**2) 
    const2 =  4*pi*Z*e*hbar**2*c**2
    FT_V = const2 / q**2
    one_minus_F = FT_V * q**2 / const1
    const3 = Z**2*alpha*e**4/(2.*pi)
    return one_minus_F**2 * p / (p0 * q**4) * curly * dvol
    
    
t = np.linspace(0+eps,pi-eps,400)
t0 = np.linspace(0+eps,pi-eps,400)
q = np.linspace(p0-p+eps, p0+p-eps, 400)

T,T0,Q = np.meshgrid(t,t0,q)

I=np.trapz(np.trapz(integrand(T,T0,Q), axis=0),axis=0)


plt.plot(q,I)
plt.show()






