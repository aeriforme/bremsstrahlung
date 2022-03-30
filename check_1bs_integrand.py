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




E1 = 100*1e6*eV
Z = 29
Z_star = 0
n_at = 8.47e+28
T = 0
N = 50

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
    print(g, g0)
    
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



        
#k = (g0-1)*0.5
#t = pi*0.5
#f = pi

#for k in ((g0-1)*0.1, (g0-1)*0.5, (g0-1)*0.7):
#    for t0 in (0.5*pi, ): 
#        plt.plot(t,k_dsigma_dk(t,t0,f,g0,k))

#for k in ((g0-1)*0.1, (g0-1)*0.5, (g0-1)*0.7):
#    for t in (0.5*pi, ): 
#        plt.plot(t0,k_dsigma_dk(t,t0,f,g0,k))


N = 400

fig, axs = plt.subplots(1,3,figsize=(3200./300., 1200./300.), dpi = 300.) 

# varying phi 
f  = np.linspace(0, 2*pi, 2*N)
for k in ((g0-1)*0.5, (g0-1)*0.7):
    for t in (0.5*pi, ):
        for t0 in (0.8*pi, 0.5*pi): 
            axs[0].plot(f,k_dsigma_dk(t,t0,f,g0,k), label='%.2f, %.2f, %.2f' % (k, t0, t))

axs[0].legend(title = r'k, $\theta_0$, $\theta$')
axs[0].set_xlabel(r'$\phi$')


# varying theta0
t0  = np.linspace(0, pi, N)
for k in ((g0-1)*0.5,(g0-1)*0.9):
    for t in (0.5*pi, ):
        for f in (0.,  ): 
            axs[1].plot(t0,k_dsigma_dk(t,t0,f,g0,k), label='%.2f, %.2f, %.2f' % (k, t, f))

axs[1].legend(title = r'k, $\theta$, $\phi$')
axs[1].set_xlabel(r'$\theta_0$')



# varying theta
t  = np.linspace(0, pi, N)
for k in ((g0-1)*0.5, (g0-1)*0.7):
    for t0 in (0.5*pi, ):
        for f in (0., ): 
            axs[2].plot(t,k_dsigma_dk(t,t0,f,g0,k), label='%.2f, %.2f, %.2f' % (k, t0, f))

axs[2].legend(title = r'k, $\theta_0$, $\phi$')
axs[2].set_xlabel(r'$\theta$')

plt.tight_layout()
fig.savefig('integrand.png') 
plt.close('all')



#    T,T0,F = np.meshgrid(t,t0,f)

#    y1 = []
#    y2 = []
#    for k in k_vec:    
#        y1.append( np.trapz(np.trapz(np.trapz(k_dsigma_dk(T,T0,F,g0,k), t, axis=0), t0, axis=0), f, axis=0) ) 

#    end = perf_counter() 
#    execution_time = (end - start)
#    print('trapz time=', execution_time)
    
 

from one_BS import one_BS_num_limit

x, y = one_BS_num_limit(E1 = E1, Z = Z, Z_star = Z_star, n_at = n_at, T = T, N = 100)
plt.plot(x,y)



from ult_rel import ult_rel_sdcs

x, y = ult_rel_sdcs(E1 = E1, Z = Z, Z_star = Z_star, n_at = n_at, T = T)
plt.plot(x, y, lw = 2, label = 'TFD')

plt.show()




