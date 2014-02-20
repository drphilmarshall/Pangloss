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
plt.semilogy(z, Mhalo, 'bo', label='Halo Mass')
#plt.semilogy(z, Mstell, 'ro', label='Stellar Mass')

plt.xlabel(r'Redshift, $z$', fontsize=16)
plt.ylabel(r'Mass $(M_{\odot}/h)$',fontsize=16)
plt.title('Halo Mass in Millenium Simulation', fontsize=16)

plt.tight_layout()
#plt.legend()

savedfile = "mill_mass_z.pdf"
#plt.savefig(savedfile,dpi=300)
print "Plot saved as "+os.getcwd()+"/"+savedfile

plt.show()