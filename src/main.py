from rlc import RLC
from photon import Photon
from matplotlib import pyplot as plt
import numpy as np

m = 10 # iteration
dt = 1e-11 #time step
n = 100 #steps per time step
I = 1 #intensity in photons per time step
T = 9e-2 #effective temperature of light source


r = .4e6
l = 10e-6
c = 1e-15
a = 500e-9

energy_l = []

j=0
while m >= j:
    x_arr = np.array([])
    step = 0
    t_c = [0,5e-10]
    ic = [1,0]
    systemp = RLC(r,l,c,t_c,ic)
    for step in np.arange(t_c[0],t_c[1]-dt,dt):
        p = Photon(T, I)

        omega = systemp.omega
        wshift= np.sum(p.random_frequency()*p.coupler(a,a/3)) # photon that we sent
            
        systemp.omega = omega+wshift
        systemp.tint = [step, step+dt]
        
        t, x, v = systemp.solve_ivp(n)

        x_arr = np.concatenate([x_arr, x], axis=None)

        plt.plot(t, x, 'b')
        systemp.ic=[x[-1],v[-1]]

    t_c = [0,5e-10]    
    ic = [1,0]    
    sys = RLC(r,l,c,t_c,ic)
    t_i,x_i,v_i = sys.solve_ivp(5000)
    initial_pow =sys.power(r, x_i)
    final_pow = sys.power(r, x_arr)
    delta_e = abs(final_pow*t[-1] - initial_pow*t_i[-1])
    energy_l.append(delta_e)
    plt.plot(t_i,x_i,'r')
    plt.xlabel('time')
    plt.ylabel('current')
    plt.show()
    j = j + 1


energy_arr = np.array(energy_l)
std = np.std(energy_arr)
median = np.median(energy_arr)
mean = np.mean(energy_arr)

print('std:', std, 'median:', median, 'mean:', mean)