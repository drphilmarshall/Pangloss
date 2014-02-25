#!/usr/bin/env python
# ===========================================================================
# This file plots modified luminosity function
# ===========================================================================

__author__ = 'cmason'
# ===========================================================================
# Import useful packages
from matplotlib import rc_file
rc_file('matplotlibrc')
import os
import mpmath
import numpy as np
import matplotlib.pyplot as plt
import scipy.special as sp
from scipy.stats.kde import gaussian_kde
from scipy import integrate

# ===========================================================================
# Define useful functions
# ---------------------------------------------------------------------------

def schechterSMF(M, Phistar, Mstar, alpha):
    """
    Schechter SMF function in terms of log (Mstar/Msun).
    Number volume density of galaxies in range (M, M + dM) is phi(M) * dM

    """
    Mdiff = 10 ** (M - Mstar)
    phi_M = Phistar * np.log(10.0) * Mdiff ** (alpha + 1.0) * np.exp(-Mdiff)
    return phi_M   

def find_index(redshift,zl,zu):
    """
    Finds the indices of all elements in array that are within given
    redshift range
    """
    ind = np.transpose(np.nonzero((redshift >= zl) & (redshift <= zu)))
    return ind

# Comoving Distance
WM = 0.3
WV = 0.7
H0 = 70 # km/s/Mpc
c = 299792.458 # km/s

def Hubble(x):
    """c / Hubble parameter in Mpc"""
    return c/(H0 * np.sqrt(WV + WM*(1 + x)**3))
    
def CMD(xmin,xmax):
    """Calculate comoving distance"""
    cmd = integrate.quad(Hubble, xmin, xmax)
    return cmd[0]
    
def proper_size(z):
    """Calculates proper size of 1 degree length at redshift z"""
    theta = 0.01745 # 1 degree in radians
    return theta * CMD(0.0, z)/(1. + z)
        
# ===========================================================================
# Data from Muzzin et al. 2013 best fit Schechter params
#smf = np.genfromtxt('Muzzin2013_MaxLik_Schechter_QUIESCENT.dat', skip_header=26)
smf = np.genfromtxt('Muzzin2013_MaxLik_Schechter_ALL_alpha1.2.dat', skip_header=26)

z_low = smf[:,0]
z_high = smf[:,1]
z = (z_high + z_low)/2.0

Mstar = smf[:,4]
Mstar_u = smf[:,7]
Mstar_l = smf[:,8]

Phistar = (smf[:,9])/10
Phistar_u = (smf[:,12])/10
Phistar_l = (smf[:,13])/10

alpha = smf[:,14]
alpha_u = smf[:,17]
alpha_l = smf[:,18]

# ===========================================================================
# Millennium Simulation Catalog

catalog = np.genfromtxt('../calib/Millennium/catalog_example.txt', skip_header=1)

z_ms = catalog[:,5]
Mhalo_ms = catalog[:,9]
Mstell_ms = catalog[:,11]

plt.figure(1, figsize=(4,16))

phi_M = []

for i in range(len(z)):
    plt.subplot(7,1,i+1)    
    zrange = find_index(z_ms, z_low[i], z_high[i])
    zmid = 0.5 * (z_low[i] + z_high[i])
    
    Mstelli = np.log10(Mstell_ms[zrange[:,0]])

    # proper volume at zmid
    volume = proper_size(zmid) ** 3.0
    
    print len(Mstelli),' objects at z = ' ,zmid
    print 'in a proper volume of ' ,volume, 'Mpc^3'
    
    # Kernel density
    flens_kde = gaussian_kde(Mstelli)
    x = np.linspace(6,12,100)
    plt.plot(x, flens_kde(x), 'r') # distribution function
    
    n, bins, patches = plt.hist(Mstelli, 8, facecolor='0.2', alpha=0.4, normed=True)
    plt.setp(patches, 'edgecolor', 'none')    
    
    Mstell_mid =[]
    for j in range(len(bins)-1):
        bin_mid = (bins[j]+bins[j+1]) * 0.5
        Mstell_mid.append(bin_mid)
   
    Mstell_mid = np.array(Mstell_mid)
    phi_mill = np.array(n*len(Mstelli)/volume)

    phi_mstell = np.column_stack((phi_mill, Mstell_mid))
    phi_M.append(phi_mstell)
                            
    plt.ylabel(str(z_low[i])+r' < z < '+str(z_high[i]))
    plt.tight_layout()
phi_M = np.array(phi_M)
print phi_M    
plt.xlabel(r'$\log(M_{stell}/M_{\odot})$', fontsize=16)

savedfile = "mill_mass_histograms.pdf"
plt.savefig(savedfile,dpi=300,bbox_inches='tight')
print "Plot saved as "+os.getcwd()+"/"+savedfile


# ===========================================================================
# Plot SMF
# ---------------------------------------------------------------------------
Mstell = np.arange(7,12,0.1)

plt.figure(2)
colors=['k','r','b','g','orange','lightblue','lightgreen']
plt.subplot(111)
for i in range(len(z)):
    phi_mill = phi_M[i,:,0]
    M_mill = phi_M[i,:,1]
    
    plt.scatter(M_mill,phi_mill, color=colors[i])
    plt.plot(M_mill,phi_mill, color=colors[i], linestyle='dashed')

    plt.semilogy(Mstell, schechterSMF(Mstell, Phistar[i], Mstar[i], alpha[i]), 
                label=str(z_low[i])+r'$< z <$'+str(z_high[i]),color=colors[i])
    print 'Plotting SMF in range ',z_low[i],' < z < ',z_high[i]

plt.xlabel(r'$\log(M_{stell}/M_{\odot})$', fontsize=16)
plt.ylabel(r'$\log$ ${\Phi}$ $($Mpc$^{-3})$', fontsize=16)
plt.title(r'Stellar Mass Function')
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)  

savedfile = "mill_smf.pdf"
plt.savefig(savedfile,dpi=300,bbox_inches='tight')
print "Plot saved as "+os.getcwd()+"/"+savedfile

plt.show()
