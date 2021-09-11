#%%
import matplotlib.pyplot as plt 
import numpy as np 

a = 0.2

tau = np.linspace(0,12,1000)
t = 1/a*np.sinh(a*tau)
x= 1/a*(np.cosh(a*tau) -1 )
plt.plot(x,t,'k',label = 'Rocket')


plt.plot(np.linspace(0,np.max([x])),np.linspace(0,np.max([x]))+1.5,'--b', label = 'Observed Light')
plt.plot(np.linspace(0,np.max([x])),np.linspace(0,np.max([x]))+3,'--b')
x1,y1 = np.linspace(0,np.max([x])),np.linspace(0,np.max([x]))+4.63
plt.plot(x1,y1,'--r', label = 'Rindler Horizon')
plt.text(x1[20]-5, y1[20]+1,'Rindler Horizon',c ='r')
plt.xlabel('x')
plt.ylabel('t')
plt.title('Spacetime path of the rocket and the Rindler Horizon')
plt.legend()
plt.savefig('Rindler.svg')
# %%
