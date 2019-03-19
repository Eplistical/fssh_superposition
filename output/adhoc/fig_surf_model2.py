import numpy as np
from matplotlib import pyplot as plt
from scipy.special import erf
import sys

plt.rcParams.update({'font.size': 15})

A = 0.1
B = 3.0
C = 3.0

c0, c1 = 0.6, 0.8
sigma = 1.0
xinit = -8

def cal_theta(x):
    return np.pi / 2 * (erf(B * (x - C)) - erf(B * (x + C)))

def cal_der_theta(x):
    return np.sqrt(np.pi) * B * np.exp(-B * B * (x - C)**2) - np.sqrt(np.pi) * B * np.exp(-B * B * (x + C)**2)


Nx = 500
xx = np.linspace(-10,10,Nx)

Psi0_2 = c0**2 * np.exp(-2 * (xx - xinit)**2 / sigma**2) / 10
Psi1_2 = c1**2 * np.exp(-2 * (xx - xinit)**2 / sigma**2) / 10 

theta = cal_theta(xx)
der_theta = cal_der_theta(xx)

H00 = -A * np.cos(theta)
H11 =  A * np.cos(theta)
E0 = np.ones(Nx) * -A
E1 = np.ones(Nx) * A

d01 = der_theta / 2

fig= plt.figure(figsize=(6,6))
ax = fig.gca()

ax.fill_between(xx, E0, Psi0_2 + E0, color='#eec881')
ax.fill_between(xx, E1, Psi1_2 + E1, color='#81d0ee')

ax.arrow(-8,0.13,1, 0, width=0.001, head_width=0.01, head_length=0.2, color='k')
ax.arrow(-8,-0.08,1, 0, width=0.001, head_width=0.01, head_length=0.2, color='k')
ax.annotate('$p_{init}$', xy=(-7.2, 0.14))
ax.annotate('$p_{init}$', xy=(-7.2, -0.07))


ax.plot(xx, E0, marker='', color='red', linestyle='solid', linewidth=2, label='$E_0$')
ax.plot(xx, E1, marker='', color='green', linestyle='solid', linewidth=2, label='$E_1$')
ax.plot(xx, H00, marker='', color='blue', linestyle='dashed', linewidth=2, label='$H_{00}$')
ax.plot(xx, H11, marker='', color='black', linestyle='dashed', linewidth=2, label='$H_{11}$')
ax.plot(xx, d01 / 10, marker='', color='brown', linestyle='dotted', linewidth=2, label='$d_{01} / 10$')

#ax.annotate('diabat 0', xy=(0.03, 0.14), xycoords='axes fraction')
#ax.annotate('adiabat 0', xy=(0.02, 0.22), xycoords='axes fraction')
#ax.annotate('diabat 1', xy=(0.03, 0.54), xycoords='axes fraction')
#ax.annotate('adiabat 1', xy=(0.02, 0.62), xycoords='axes fraction')

ax.annotate('diabat 1', xy=(0.79, 0.28), xycoords='axes fraction')
ax.annotate('adiabat 0', xy=(0.78, 0.35), xycoords='axes fraction')
ax.annotate('diabat 0', xy=(0.79, 0.62), xycoords='axes fraction')
ax.annotate('adiabat 1', xy=(0.78, 0.69), xycoords='axes fraction')

ax.set_title('Energy Surfaces')
ax.set_xlim([-10,10])
ax.set_xlabel('$x$')
ax.set_ylim([-0.3,0.3])
#ax.set_ylabel('energy')

ax.legend(loc='upper right', bbox_to_anchor=(1.1,1.15))
fig.savefig('fig_surf_model2.eps')
#plt.show()

