#%%
from scipy.integrate import quad
import numpy as np
from scipy.optimize import fsolve
import matplotlib.pyplot as plt

#%%
print('''
Question 1 b
''')
r_1 = 10

f = lambda r,r0: 1/(r**2 * np.sqrt( 1/r0**2*(1-2/r0) -  1/r**2*(1-2/r)))


def findZ(input_r):
    y,err = quad(f, input_r, r_1, args=(input_r,))
    return y - np.pi/2


r_0 = fsolve(findZ, 5)[0]

print('r_0 = {:.2f}'.format(r_0))

#%%
print('''
Question 1 c
''')

ib2 = 1/r_0**2*(1-2/r_0)

ib = np.sqrt(ib2)

tfunc = lambda r: ib*(1-2/r)**(-1)*( ib2 - 1/r**2*(1-2/r))**(-1/2)

t_0, e = quad(tfunc, r_1, r_0)
t_1, e = quad(tfunc, r_0, r_1)
t_2, e = quad(tfunc, r_1, r_0)
t_3, e = quad(tfunc, r_0, r_1) 


print('Travel Time Schwarzschild: {:.2f}'.format( - t_0 + t_1 - t_2 + t_3))

#%%
print('''
Question 1 d
''')
t_flat_segment = np.sqrt(r_0**2 + r_1**2)
print('Travel Time Flat: {:.2f}'.format(t_flat_segment*4))
print('Flat < Schwarzschild: {}'.format(t_flat_segment*4 < - t_0 + t_1 - t_2 + t_3))

#%%
print('''
Question 3 d
''')

R = 10
m = 1
w_s = 1
r = np.arange(100,10.1, -0.1)
w_r = (1 + np.sqrt(2*m/r))*(r/(r-2*m))*np.sqrt(1 - 2*m/R)

plt.plot(r, w_r)
plt.xlabel('r')
plt.ylabel(r'$\omega_r$')
# %%
