#%%
# from scipy.integrate import odeint
import numpy as np 
import matplotlib.pyplot as plt
#%%
# Function for all the differential Equations
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
# intial conditions for first massive particle 
R = 10
C = 1/(10 *np.sqrt(7))
m = 1
uu = -1
initial = [0,R,np.pi/2,0,np.sqrt(((R**2*C**2-uu)*R)/(R-2*m)),0,0,C]
