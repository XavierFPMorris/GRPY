#%%
from scipy.integrate import odeint
import numpy as np 
import matplotlib.pyplot as plt

#%%

R = 10
C = 1/(10 *np.sqrt(7))
m = 1
uu = -1
initial = [0,R,np.pi/2,0,np.sqrt(((R**2*C**2-uu)*R)/(R-2*m)),0,0,C]

#%%

def orbit(state,tau):
    t,r,th,p,T,R,TH,P = state

    dt = T
    dr = R
    dth = 0#TH
    dp = P
    dT = -2*m/(r*(r-2*m))*T*R
    dR = m/(r*(r-2*m))*R**2 - m*(r-2*m)/(r**3)*T**2 + (r-2*m)*P**2
    dTH = 0 
    dP = -2/r*R*P

    return [dt,dr,dth,dp,dT,dR,dTH,dP]

#%%
time_points = np.linspace(0,1000,4000)

int_ret = odeint(orbit, initial, time_points)

#%%
t = int_ret[:,0]
r = int_ret[:,1]
phi = int_ret[:,3]
# %%
# phi func of tau
plt.plot(time_points, int_ret[:,3])
# %%
# r func of tau
plt.plot(time_points, int_ret[:,1])
# %%
# circ show 
fig, ax = plt.subplots()
ax.plot(r*np.cos(phi), r*np.sin(phi))
ax.set_aspect(1)
# %%
#try again for bigger c
C = 1/(10 *np.sqrt(7))*1.1
R = 10
m = 1
uu = -1
initial = [0,R,np.pi/2,0,np.sqrt(((R**2*C**2-uu)*R)/(R-2*m)),0,0,C]
int_ret_2 = odeint(orbit, initial, time_points)
t_2 = int_ret_2[:,0]
r_2 = int_ret_2[:,1]
phi_2 = int_ret_2[:,3]
# %%
# circ show 2
fig, ax = plt.subplots()
ax.plot(r_2*np.cos(phi_2), r_2*np.sin(phi_2))
ax.set_aspect(1)
# %%
