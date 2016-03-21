## Test to check accuracy of lensing lookup table

import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import os, sys

PANGLOSS_DIR = os.path.expandvars("$PANGLOSS_DIR")
sys.path.append(PANGLOSS_DIR)
import pangloss

table = pangloss.lensing.LensingTable(x_lim=[1e-6,500],t_lim=[50.0,200.0],Nx=1000,Nt=100)

X = table.x_grid
T = table.t_grid
Z = table.BMO1F_grid

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

ax.plot_surface(X, T, Z)

ax.set_xlabel('x')
ax.set_ylabel('t')
ax.set_zlabel('BMO1F')

plt.grid()
plt.show()

fig2 = plt.figure()
ax2 = fig2.add_subplot(111, projection='3d')

x2 = np.arange(1e-6,500,1)
t2 = np.arange(50,200)

X2,T2 = np.meshgrid(x2,t2)
X2 = np.transpose(X2)
T2 = np.transpose(T2)
Z2 = table.BMO1F_spline.ev(X2,T2)
print 'Z2 = '+str(Z2)
ax2.plot_surface(X2, T2, Z2)

ax2.set_xlabel('x')
ax2.set_ylabel('t')
ax2.set_zlabel('BMO1F')

plt.grid()
plt.show()
