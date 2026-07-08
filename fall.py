#!/usr/bin/env python3

### Calculate free fall with newtonian air resistance

import numpy as np
import matplotlib.pyplot as plt

# input values:
cw = float(input("c_w-Wert: ")) # dimensionless, c_w value, sphere=0.45, eichel=0.3, zeppelin(with l:d=6) = 0.05
d = float(input("d (in m): ")) # in m, diameter of body, eichel=0.016
m = float(input("m (in kg): ")) # in kg, mass of body, eichel=0.006
A = (d/2)**2*np.pi
#v0 = 0 # initial velocity, in m/s
x0=0
# constants:
rho = 1.20 # in kg/m^3, density of air, normal pressure and temperature
g = 9.805 # in m/s^2, stand. grav. acc.

# calc terminal velocity (from F_g = m*g = F_r = cw*rho*A/(2*m)*v^2)
vt = np.sqrt(2*m*g/(cw*rho*A))

# time to reach terminal velocity
tt = (vt*np.arctanh(0.95))/g
tt1 = tt + 1

#initialize array of time
t=np.arange(0,tt1,0.01)

## calc functions
# Diff eq: x" = -g + k*v^2 with k=cw*rho*A/(2m)
# Ansatz for diff eq: v' = -g + k*v^2 -> v=p*tanh(qt), because d/dx tanh(x) = 1-tanh^2(x) -> v=-vt*tanh(g*t/vt)+v0
v = -vt*np.tanh(g*t/vt)
# integrate tanh via tanh=sinh/cosh-> substitute u=cosh
x=x0-(vt)**2/g*np.log(np.cosh(g/vt*t))
# derivative
a=-g+cw*rho*A/(2*m)*v**2

# plotting
plt.figure(figsize=(9,6))
plt.xlim(0,tt1)
plt.plot(t,x/10, color='blue', label=r'$x$ in $10\cdot\mathrm{m}$')
plt.plot(t,v, color='green', label=r'$v$ in $\mathrm{m}/\mathrm{s}$')
plt.plot(t,a, color='red', label=r'$a$ in $\mathrm{m}/\mathrm{s}^2$')
plt.plot([0.,tt1],[-vt,-vt], color='gray', ls='--', label=r'$v_t$ in m/s')
plt.plot([], [], ' ', label=f"v_t = {vt:.2f} m/s")
plt.plot([], [], ' ', label=f"v = .95 v_t after {tt:.2f} s")
plt.xlabel(r'$t$ in s')
plt.title(r'Free fall for object with $c_w= %.2f , d= %.3f \mathrm{ m} , m= %.3f \mathrm{ kg}$, in norm atmosphere' % (cw, d, m))
plt.grid()
plt.legend()
plt.show()
