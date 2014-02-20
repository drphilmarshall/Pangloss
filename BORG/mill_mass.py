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
plt.semilogy(z, Mstell, 'bo', alpha=0.5, label='Stellar Mass')

nbins=10

histz, edges = np.histogram(z, bins=nbins)

bins = np.linspace(0,4,10)

izbin = np.digitize(z,edges)

minMstell, maxMstell = [], []
for n in range(nbins):
    ind = np.nonzero(izbin == n+1.0)
    ind = np.array(ind)
    
    iMstell=[]
    
    for i in range(len(ind.T)):
        iMstell.append(Mstell[i])
    
    maxiMstell=np.amax(iMstell)
    maxMstell.append(maxiMstell)
    miniMstell=np.amin(iMstell)
    minMstell.append(miniMstell)

print len(edges[0:10]), len(minMstell)
plt.plot(edges[0:10], maxMstell, 'r-')
plt.plot(edges[0:10], minMstell, 'r-')

plt.xlabel(r'Redshift, $z$', fontsize=16)
plt.ylabel(r'Mass $(M_{\odot}/h)$',fontsize=16)
plt.title('Stellar Mass in Millenium Simulation', fontsize=16)

plt.tight_layout()
#plt.legend()

savedfile = "mill_mass_z.pdf"
plt.savefig(savedfile,dpi=300)
print "Plot saved as "+os.getcwd()+"/"+savedfile

plt.show()