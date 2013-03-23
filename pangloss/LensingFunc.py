'''
This file is part of the Pangloss project.
Copyright 2012 Tom Collett (IoA) and Phil Marshall (Oxford).

description:
------------
Functions that are useful for lensing (and necessary to the Pangloss Project)

to-do:
------


issues:
-------


'''

import distances
import pylab
import matplotlib.pyplot as plt
import numpy, numpy.random as rnd, atpy
from time import clock

D = distances.Distance()

# ----------------------------------------------------------------------------
def SigmaCrit_Da(Da,Da_s,Da_tosource): 
    # NOTE zl here is the lensing object NOT necessarily the primary lens
    return (1.663*10**18)*(Da_s/(Da*Da_tosource)) 
    # numerical factor is c^2/(4 pi G) in Solarmasses per megaparsec
def SigmaCrit_z(zd,zs): #As above but with redshifts (much slower)
    return (1.663*10**18)*(D.Da(zs)/(D.Da()*Da_tosource)) 

# ----------------------------------------------------------------------------

# Keeton Beta parameter for a perturber at j:  
def beta(i,j,k):  
    if j>k:
        print "z_pert > z_source? you shouldn't be asking for this"
    if j>i:
        R1a = D.Da(i,j)/D.Da(j)
        R2a = D.Da(i,k)/D.Da(k)
        return R1a/R2a
    if i>j:
        R1b = D.Da(j,i)/D.Da(i)
        R2b = D.Da(j,k)/D.Da(k)
        return R1b/R2b
    if i == j:
        return 1.0
    if j==k:
        return 0.0

# Same but for known angular diameter distances (much faster)
def beta_Da(D12,Ds,D2,D1s): 
    return (D12*Ds)/(D2*D1s)


# ----------------------------------------------------------------------------

# Kappa Keeton, following Keeton (2003) and Momcheva et al. (2006)
def KappaKeeton(zl,zd,zs,kappa,shear):
    output = numpy.zeros(len(zd))
    for i in range(len(zd)):
        if zd[i] < zs:
            B=self.beta(zl,zd[i],zs)
            K=kappa[i]
            G=shear[i]
            D= K**2-G**2
            output[i] = (1.-B) * (K- B*(D)) /  ( (1-B*K)**2   - (B*G)**2   )
        else: 
            output[i]= 0.0
            print "dont call me! - kappakeeton"
    return output
# ----------------------------------------------------------------------------
# Kappa Keeton, following Keeton (2003) and Momcheva et al. (2006)
def KappaScale_beta(beta,kappa,shear,scaling="tom"):
            B=beta
            K=kappa
            G=shear
            D= K**2-G**2

            if scaling == "keeton" or scaling == "Keeton":
                output = (1.-B) * (K- B*(D)) /  ( (1-B*K)**2   - (B*G)**2   )

            if scaling == "tom" or scaling =="Tom" or scaling =="TC" or scaling =="Collett" or scaling =="collett" or scaling == "prime":
                output=(1.-B) * (K)

            if scaling == "additive" or scaling =="addition" or scaling =="add":
                #print "boo"
                output = (K)#K*1.0
            return output
# ----------------------------------------------------------------------------
