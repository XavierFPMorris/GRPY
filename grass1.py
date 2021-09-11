#%%
from scipy.integrate import odeint
import numpy as np 
import matplotlib.pyplot as plt

#%%
# initial conditions
R = 10
C = 1/(10 *np.sqrt(7))
m = 1
uu = -1
initial = [0,R,np.pi/2,0,np.sqrt(((R**2*C**2-uu)*R)/(R-2*m)),0,0,C]

#%%
# integrator function
def orbit(state,tau):
    t,r,th,p,T,R,TH,P = state

    dt = T
    dr = R
    dth = 0
    dp = P
    dT = -2*m/(r*(r-2*m))*T*R
    dR = m/(r*(r-2*m))*R**2 - m*(r-2*m)/(r**3)*T**2 + (r-2*m)*P**2
    dTH = 0 
    dP = -2/r*R*P

    return [dt,dr,dth,dp,dT,dR,dTH,dP]

#%%
#integrate
time_points = np.linspace(0,1000,4000)
int_ret = odeint(orbit, initial, time_points, rtol = 1e-10, atol = 1e-10)

#%%
t = int_ret[:,0]
r = int_ret[:,1]
phi = int_ret[:,3]
#%%
# t func of tau
plt.plot(time_points,t)
plt.xlabel(r'$\tau$')
plt.ylabel('t')
plt.savefig('ass1/C0_t.svg')
plt.close()
# %%
# phi func of tau
plt.plot(time_points, int_ret[:,3])
plt.xlabel(r'$\tau$')
plt.ylabel(r'$\phi$')
plt.savefig('ass1/C0_phi.svg')
plt.close()
# %%
# r func of tau
plt.plot(time_points, int_ret[:,1])
plt.xlabel(r'$\tau$')
plt.ylabel('r')
plt.savefig('ass1/C0_r.svg')
plt.close()
# %%
# circ show 
fig, ax = plt.subplots()
ax.plot(r*np.cos(phi), r*np.sin(phi))
ax.set_aspect(1)
plt.xlabel(r'$x$')
plt.ylabel(r'$y$')
plt.savefig('ass1/C0_xy.svg')
plt.close()
#%%
# u.u vs tau
dt = int_ret[:,4]
dr = int_ret[:,5]
dp = int_ret[:,7]
gtt = -1*(1-2*m/r)
gpp = r**2
grr = (1-2*m/r)**(-1)
plt.plot(time_points, gtt*dt**2 +grr*dr**2 + gpp*dp**2)
plt.xlabel(r'$\tau$')
plt.ylabel(r'$u\cdot u$')
plt.savefig('ass1/C0_uu.svg')
plt.close()
# %%
#try again for bigger c
C = 1/(10 *np.sqrt(7))*1.1
R = 10
m = 1
uu = -1
initial = [0,R,np.pi/2,0,np.sqrt(((R**2*C**2-uu)*R)/(R-2*m)),0,0,C]
int_ret_2 = odeint(orbit, initial, time_points,rtol = 1e-10, atol = 1e-10)
t_2 = int_ret_2[:,0]
r_2 = int_ret_2[:,1]
phi_2 = int_ret_2[:,3]

#%%
# t v tau
plt.plot(time_points,t_2)
plt.xlabel(r'$\tau$')
plt.ylabel('t')
plt.savefig('ass1/C1_t.svg')
plt.close()
#%%
# r v tau
plt.plot(time_points,r_2)
plt.xlabel(r'$\tau$')
plt.ylabel('r')
plt.savefig('ass1/C1_r.svg')
plt.close()
#%%
# phi v tau
plt.plot(time_points,phi_2)
plt.xlabel(r'$\tau$')
plt.ylabel(r'$\phi$')
plt.savefig('ass1/C1_phi.svg')
plt.close()
# %%
# circ show 2
fig, ax = plt.subplots()
ax.plot(r_2*np.cos(phi_2), r_2*np.sin(phi_2))
ax.set_aspect(1)
plt.xlabel(r'$x$')
plt.ylabel(r'$y$')
plt.savefig('ass1/C1_xy.svg')
plt.close()
#%%
# u.u vs tau
dt_2 = int_ret_2[:,4]
dr_2 = int_ret_2[:,5]
dp_2 = int_ret_2[:,7]
gtt = -1*(1-2*m/r_2)
gpp = r_2**2
grr = (1-2*m/r_2)**(-1)
plt.plot(time_points, gtt*dt_2**2 +grr*dr_2**2 + gpp*dp_2**2)
plt.xlabel(r'$\tau$')
plt.ylabel(r'$u\cdot u$')
plt.savefig('ass1/C1_uu.svg')
plt.close()
# %%
# now for massless 

R = 10
m = 1
uu = 0
lambdas = np.linspace(0,10,100)
# c -> 0.5
C = 0.5
initial = [0,R,np.pi/2,0,np.sqrt(((R**2*C**2-uu)*R)/(R-2*m)),0,0,C]
int_ret_p1 = odeint(orbit, initial, lambdas,rtol = 1e-10, atol = 1e-10)
t_p1 = int_ret_p1[:,0]
r_p1 = int_ret_p1[:,1]
phi_p1 = int_ret_p1[:,3]
# %%

# t v tau
plt.plot(lambdas,t_p1)
plt.xlabel(r'$\lambda$')
plt.ylabel('t')
plt.savefig('ass1/p1_t.svg')
plt.close()
#%%
# r v tau
plt.plot(lambdas,r_p1)
plt.xlabel(r'$\lambda$')
plt.ylabel('r')
plt.savefig('ass1/p1_r.svg')
plt.close()
#%%
# phi v tau
plt.plot(lambdas,phi_p1)
plt.xlabel(r'$\lambda$')
plt.ylabel(r'$\phi$')
plt.savefig('ass1/p1_phi.svg')
plt.close()
# %%
# circ show 2
fig, ax = plt.subplots()
ax.plot(r_p1*np.cos(phi_p1), r_p1*np.sin(phi_p1))
plt.xlabel(r'$x$')
plt.ylabel(r'$y$')
plt.savefig('ass1/p1_xy.svg')
plt.close()
#%%
# u.u vs tau
dt_p1 = int_ret_p1[:,4]
dr_p1 = int_ret_p1[:,5]
dp_p1 = int_ret_p1[:,7]
gtt = -1*(1-2*m/r_p1)
gpp = r_p1**2
grr = (1-2*m/r_p1)**(-1)
plt.plot(lambdas, gtt*dt_p1**2 +grr*dr_p1**2 + gpp*dp_p1**2)
plt.xlabel(r'$\lambda$')
plt.ylabel(r'$u\cdot u$')
plt.savefig('ass1/p1_uu.svg')
plt.close()
# %%
# c -> 1
C = 1
initial = [0,R,np.pi/2,0,np.sqrt(((R**2*C**2-uu)*R)/(R-2*m)),0,0,C]
int_ret_p2 = odeint(orbit, initial, lambdas,rtol = 1e-10, atol = 1e-10)
t_p2 = int_ret_p2[:,0]
r_p2 = int_ret_p2[:,1]
phi_p2 = int_ret_p2[:,3]
# %%

# t v tau
plt.plot(lambdas,t_p2)
plt.xlabel(r'$\lambda$')
plt.ylabel('t')
plt.savefig('ass1/p2_t.svg')
plt.close()
#%%
# r v tau
plt.plot(lambdas,r_p2)
plt.xlabel(r'$\lambda$')
plt.ylabel('r')
plt.savefig('ass1/p2_r.svg')
plt.close()
#%%
# phi v tau
plt.plot(lambdas,phi_p2)
plt.xlabel(r'$\lambda$')
plt.ylabel(r'$\phi$')
plt.savefig('ass1/p2_phi.svg')
plt.close()
# %%
# circ show 2
fig, ax = plt.subplots()
ax.plot(r_p2*np.cos(phi_p2), r_p2*np.sin(phi_p2))
plt.xlabel(r'$x$')
plt.ylabel(r'$y$')
plt.savefig('ass1/p2_xy.svg')
plt.close()
#%%
# u.u vs tau
dt_p2 = int_ret_p2[:,4]
dr_p2 = int_ret_p2[:,5]
dp_p2 = int_ret_p2[:,7]
gtt = -1*(1-2*m/r_p2)
gpp = r_p2**2
grr = (1-2*m/r_p2)**(-1)
plt.plot(lambdas, gtt*dt_p2**2 +grr*dr_p2**2 + gpp*dp_p2**2)
plt.xlabel(r'$\lambda$')
plt.ylabel(r'$u\cdot u$')
plt.savefig('ass1/p2_uu.svg')
plt.close()
# %%
# %%
# c -> 2
C = 2
initial = [0,R,np.pi/2,0,np.sqrt(((R**2*C**2-uu)*R)/(R-2*m)),0,0,C]
int_ret_p3 = odeint(orbit, initial, lambdas,rtol = 1e-10, atol = 1e-10)
t_p3 = int_ret_p3[:,0]
r_p3 = int_ret_p3[:,1]
phi_p3 = int_ret_p3[:,3]
# %%

# t v tau
plt.plot(lambdas,t_p3)
plt.xlabel(r'$\lambda$')
plt.ylabel('t')
plt.savefig('ass1/p3_t.svg')
plt.close()
#%%
# r v tau
plt.plot(lambdas,r_p3)
plt.xlabel(r'$\lambda$')
plt.ylabel('r')
plt.savefig('ass1/p3_r.svg')
plt.close()
#%%
# phi v tau
plt.plot(lambdas,phi_p3)
plt.xlabel(r'$\lambda$')
plt.ylabel(r'$\phi$')
plt.savefig('ass1/p3_phi.svg')
plt.close()
# %%
# circ show 2
fig, ax = plt.subplots()
ax.plot(r_p3*np.cos(phi_p3), r_p3*np.sin(phi_p3))
plt.xlabel(r'$x$')
plt.ylabel(r'$y$')
plt.savefig('ass1/p3_xy.svg')
plt.close()
#%%
# u.u vs tau
dt_p3 = int_ret_p3[:,4]
dr_p3 = int_ret_p3[:,5]
dp_p3 = int_ret_p3[:,7]
gtt = -1*(1-2*m/r_p3)
gpp = r_p3**2
grr = (1-2*m/r_p3)**(-1)
plt.plot(lambdas, gtt*dt_p3**2 +grr*dr_p3**2 + gpp*dp_p3**2)
plt.xlabel(r'$\lambda$')
plt.ylabel(r'$u\cdot u$')
plt.savefig('ass1/p3_uu.svg')
plt.close()
# %%
