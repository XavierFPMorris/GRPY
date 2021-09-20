#%%
import matplotlib.pyplot as plt
import numpy as np
#%%

N = 1000
r = np.linspace(0,4,N)

#%%
# Q3 Plot
plt.rcParams["figure.figsize"] = (16,10)
plt.rcParams.update({'font.size': 12})

plt.figure()

# constant v 
for v in np.linspace(0.5,10,6):
    v_vec = v*np.ones(N)
    tt = v_vec-r
    plt.plot(r,tt,'k')

for v in np.linspace(0.5,10,6):
    v_vec = v + 2*( r + 2*np.log(abs(r/2 - 1)))
    tt = v_vec-r
    plt.plot(r,tt,'b')

plt.plot([2 ,2],[0 ,8],'r',label = 'Event Horizon')
plt.plot([ -1 -2],[-1 -2], 'k',label = 'Inwards Rays')
plt.plot([ -1 -2],[-1 -2], 'b',label = 'Outwards Rays')
plt.xlim([0,4])
plt.ylim([0,8])
plt.xlabel('r')
plt.ylabel(r'$\tilde{t}$')
plt.legend()
plt.show()
# %%
# Q4 Plot
plt.rcParams["figure.figsize"] = (16,10)
plt.rcParams.update({'font.size': 12})

plt.figure()

# constant v 
for v in np.linspace(0.5,10,6):
    v_vec = v*np.ones(N)
    tt = v_vec-r - 2*np.log(abs(r/2 - 1))
    plt.plot(r,tt,'k')

for v in np.linspace(0.5,10,6):
    v_vec = v + 2*( r + 2*np.log(abs(r/2 - 1)))
    tt = v_vec-r- 2*np.log(abs(r/2 - 1))
    plt.plot(r,tt,'b')

plt.plot([2 ,2],[0 ,8],'r',label = 'Event Horizon')
plt.plot([ -1 -2],[-1 -2], 'k',label = 'Inwards Rays')
plt.plot([ -1 -2],[-1 -2], 'b',label = 'Outwards Rays')
plt.xlim([0,4])
plt.ylim([0,8])
plt.xlabel('r')
plt.ylabel('t')
plt.legend()
plt.show()
# %%
