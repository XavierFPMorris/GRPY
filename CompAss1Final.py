#%%
from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt
#%%
initial = [0,0,1,0]
time_points_0 = np.linspace(-1,0,4)
time_points_1 = np.linspace(0,25,100)
time_points_2 = np.linspace(25,75,200)
time_points_3 = np.linspace(75,100,100)
time_points_4 = np.linspace(100,101,4)


def rocket(state,tau):
    t,x,T,X= state
    dt = T
    dT = a*X
    dx = X
    dX = a*T

    return [dt,dx, dT, dX]
a = 0
tx_0 = odeint(rocket, initial, time_points_0, rtol = 1e-10, atol = 1e-10)
a = 0.2
tx_1 = odeint(rocket, tx_0[-1,:], time_points_1, rtol = 1e-10, atol = 1e-10)
a = -0.2
tx_2 = odeint(rocket, tx_1[-1,:], time_points_2, rtol = 1e-10, atol = 1e-10)
a = 0.2
tx_3 = odeint(rocket, tx_2[-1,:], time_points_3, rtol = 1e-10, atol = 1e-10)
a = 0
tx_4 = odeint(rocket, tx_3[-1,:], time_points_4, rtol = 1e-10, atol = 1e-10)
#%%
tx = np.vstack((tx_0,tx_1,tx_2,tx_3,tx_4))
tp = np.hstack((time_points_0, time_points_1, time_points_2, time_points_3, time_points_4))
#%% (x,t)
fig, ax = plt.subplots()
ax.plot(tx[:,1],tx[:,0])
ax.set(xlabel = 'x', ylabel = 't', title = 'Path of rocket, spacetime diagram')
plt.show()
# %% (t,dx/dt)
fig, ax = plt.subplots()
ax.plot(tx[:,0],tx[:,3]/tx[:,2])
ax.set(xlabel = 't', ylabel = 'dx/dt', title = 'Speed as a function of time')
plt.show()
# %% (tau,u.u) zoom out
fig, ax = plt.subplots()
ax.plot(tp, -1*(tx[:,2])**2 + (tx[:,3])**2)
ax.set_ylim([-2,0])
ax.set(xlabel = 'tau', ylabel = 'u.u', title = 'u.u as a function of tau, large range')
plt.show()
# %% (tau,u.u) 
fig, ax = plt.subplots()
ax.plot(tp, -1*(tx[:,2])**2 + (tx[:,3])**2)
ax.set(xlabel = 'tau', ylabel = 'u.u', title = 'u.u as a function of tau')
plt.show()
# %% (tau,t)
fig, ax = plt.subplots()
ax.plot(tp,tx[:,0])
ax.set(xlabel = 'tau', ylabel = 't', title = 'Time as a function of proper time')
plt.show()
# %%
