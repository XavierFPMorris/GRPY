#%%
from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt
#%%
initial = [0,0,1,0]
time_points_1 = np.linspace(0,25,100)
time_points_2 = np.linspace(25,75,200)
time_points_3 = np.linspace(75,100,100)

def rocket(state,tau):
    t,x,T,X= state
    dt = T
    dT = a*X
    dx = X
    dX = a*T

    return [dt,dx, dT, dX]

a = 0.2
tx_1 = odeint(rocket, initial, time_points_1)
a = -0.2
tx_2 = odeint(rocket, tx_1[-1,:], time_points_2)
a = 0.2
tx_3 = odeint(rocket, tx_2[-1,:], time_points_3)

#%%
# (x,t)
fig, ax = plt.subplots()
ax.plot(tx_1[:,1],tx_1[:,0])
ax.plot(tx_2[:,1],tx_2[:,0])
ax.plot(tx_3[:,1],tx_3[:,0])
ax.set(xlabel = 'x', ylabel = 't', title = 'Path of rocket, spacetime diagram')
plt.show()
# %%
# (t,dx/dt)
fig, ax = plt.subplots()
ax.plot(tx_1[:,0],tx_1[:,3]/tx_1[:,2])
ax.plot(tx_2[:,0],tx_2[:,3]/tx_2[:,2])
ax.plot(tx_3[:,0],tx_3[:,3]/tx_3[:,2])
ax.set(xlabel = 't', ylabel = 'dx/dt', title = 'Speed as a function of time')
plt.show()
# %%
# (t,dx/dt)
# fig, ax = plt.subplots()
# ax.plot(tx_1[:,0],tx_1[:,3])
# ax.plot(tx_2[:,0],tx_2[:,3])
# ax.plot(tx_3[:,0],tx_3[:,3])
# ax.set(xlabel = 't', ylabel = 'dx/dtau', title = 'Speed as a function of time')
# plt.show()
# %%
# (tau,u.u)
fig, ax = plt.subplots()
ax.plot(time_points_1, -1*(tx_1[:,2])**2 + (tx_1[:,3])**2)
ax.plot(time_points_2, -1*(tx_2[:,2])**2 + (tx_2[:,3])**2)
ax.plot(time_points_3, -1*(tx_3[:,2])**2 + (tx_3[:,3])**2)
ax.set_ylim([-2,0])
ax.set(xlabel = 'tau', ylabel = 'u.u', title = 'u.u as a function of tau')
plt.show()
# %%
# (tau,t)
fig, ax = plt.subplots()
ax.plot(time_points_1,tx_1[:,0])
ax.plot(time_points_2,tx_2[:,0])
ax.plot(time_points_3,tx_3[:,0])
ax.set(xlabel = 'tau', ylabel = 't', title = 'Time as a function of proper time')
plt.show()
# %%
