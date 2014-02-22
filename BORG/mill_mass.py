#!/usr/bin/env python
# ===========================================================================
# Plots Mstell as fn of z for Millenium Sim
# ===========================================================================

__author__ = 'cmason'
# ===========================================================================
# Import useful packages
#from matplotlib import rc_file
#rc_file('matplotlibrc')
import os
import mpmath
import numpy as np
import matplotlib.pyplot as plt
import scipy.special as sp
from scipy import integrate
from matplotlib import rc_file
rc_file('matplotlibrc')

# ===========================================================================
# Define useful functions
# ---------------------------------------------------------------------------


# ===========================================================================    
# Differential Optical Depth - probability for lensing for source at z=8
# ---------------------------------------------------------------------------

# Input from Millenium Simulation
catalog = np.genfromtxt('../calib/Millennium/catalog_example.txt', skip_header=1)

z = catalog[:,5]
Mhalo = catalog[:,9]
Mstell = catalog[:,11]

# Plot
plt.figure(1)
plt.subplot(111)
#plt.semilogy(z, Mhalo, 'bo', label='Halo Mass')
#plt.semilogy(z, Mstell, 'bo', alpha=0.5, label='Stellar Mass')



def find_index(zl,zu):
    ind = np.transpose(np.nonzero((z >= zl)))# & (z <= zu)))
    return ind
        
maxMstell, minMstell, zmid = [], [], []

for i in range(7):
    print 'finding max and min stellar mass in range ' + str(i/2.0) + ' < z < ' + str((i+1)/2.0)
    bini = find_index(i/2.0,(i+1.0)/2.0)
    zmidi = 0.5*i+0.25
    Mstelli = Mstell[bini[:,0]]
    maxMstelli = np.amax(Mstelli)
    minMstelli = np.amin(Mstelli)
    maxMstell.append(maxMstelli)
    minMstell.append(minMstelli)
    zmid.append(zmidi)

plt.fill_between(zmid, maxMstell, minMstell,
    alpha=0.5, edgecolor='#CC4F1B', facecolor='#FF9848')
    
plt.semilogy(zmid, maxMstell, 'r-')
plt.semilogy(zmid, minMstell, 'r-')

plt.xlabel(r'Redshift, $z$', fontsize=16)
plt.ylabel(r'Mass $(M_{\odot}/h)$',fontsize=16)
plt.title('Stellar Mass in Millenium Simulation', fontsize=16)

plt.tight_layout()
#plt.legend()

savedfile = "mill_mass_z.pdf"
plt.savefig(savedfile,dpi=300)
print "Plot saved as "+os.getcwd()+"/"+savedfile

plt.show()