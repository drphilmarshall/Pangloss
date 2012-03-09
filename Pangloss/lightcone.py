'''
This file is part of the Pangloss project.
Copyright 2012 Tom Collett (IoA) and Phil Marshall (Oxford).

to-do
-----
- finish basic functions
- carry out curve of growth test
- decide whether to continue
'''

import distances,pylab
import matplotlib.pyplot as plt
import numpy
import numpy.random as rnd
#import time
#t0=time.clock()    
D = distances.Distance()
zl=0.3
zs=1.5


# ============================================================================
    
class lightcone:

   def __init__(self, position, radius, zd, zs):
        self.name = 'Light cone containing a gravitational lens'
        self.position = np.array(position)
        self.radius = radius
        self.zd = zd
        self.zs = zs
        return None

    def __str__(self):
        return '%.2f %.2f ' % (position)

# ----------------------------------------------------------------------------
    
    def selectlens()

    def beta(i,j,k):  
        if j>i:
            R1 = D.Da(i,j)/D.Da(j)
            R2 = D.Da(i,k)/D.Da(k)
            return R2/R1
        if i>j:
            R1 = D.Da(j,i)/D.Da(i)
            R2 = D.Da(j,k)/D.Da(k)
            return R1/R2

    def Kappaindiv()

    def Shearindiv()



    def Kappa_keeton(i=zl,j,k=zs,):
        B=beta(i,j,k)
        K=Kappaindiv()
        G=Shearindiv()
        return (1.-B)*(K-B(K^2-G^2))/((1-BK)^2-(G)^2)

    def Kappa_sum(i,j,k):


# ============================================================================

# TESTS:

def cog_test():

    print "Plotting curve of growth around a random lens"
    
    pos = [0.0, 0.0]
    r = 0.05 # radians
    lc = lightcone(pos, r, zd, zs)

    return
   
# ============================================================================

if __name__ == '__main__':

    cog_test()
    
# ============================================================================





"""
kappaN=numpy.ones(n)
z_phot=numpy.ones(n)



n=numpy.size(kappaN)
KappaN=numpy.zeros(n)
KappaNtrue=numpy.zeros(n)

#{ Test 1: compare one line of sight with randomised kappas and redshifts
for i in range(n):
    kappaI=kappaN[i]*(rnd.randn()) #*?
    zN=z_phot[i]*(1+rnd.randn()*.05) #?

    KappaN[i]=kappaI*beta(zl,zN,zs)
    KappaNtrue[i]=kappaN[i]*beta(zl,z_phot[i],zs)



Kappa=KappaN.sum()
KappaTrue=KappaNtrue.sum()

Kappasum=KappaN.cumsum()
KappaTruesum=KappaNtrue.cumsum()
y=range(n)
plt.plot(KappaTruesum,y,'b',linestyle='--')
plt.plot(KappaTruesum,y,'b')
plt.show()
#}



#{ Test2: calculate lots of kappas, and see how we do on average:

kappaTrue=0.0
kappa=0.0
m=1000
Kappa=numpy.zeros(m)
for i in range(n):
    KappaNtrue[i]=kappaN[i]*beta(zl,z_phot[i],zs)
    for j in range(m):
        print m
        zN=z_phot[i]*rnd.randn()*.05 #?
        kappaI=kappaN[i]*rnd.randn()*1 #?
        Kappa[m]+=kappaI*beta(zl,zN,zs)

Kappatrue=KappaNtrue.sum()    

plt.hist(Kappa)
plt.axvline(x=Kappatrue, ymin=0, ymax=1,color='black', ls='dotted')
plt.show()
#}
"""
