'''
This file is part of the Pangloss project.
Copyright 2012 Tom Collett (IoA) and Phil Marshall (Oxford).

to-do
-----
- finish basic functions
- carry out curve of growth test
- decide whether to continue
'''

# import distances
import pylab, atpy
import matplotlib.pyplot as plt
import numpy
import numpy.random as rnd

#import time
#t0=time.clock()    

# D = distances.Distance()
# zl=0.3
# zs=1.5

arcmin2rad = (1.0/60.0)*numpy.pi/180.0
rad2arcmin = 1.0/arcmin2rad

# ============================================================================
    
class catalog:

   def __init__(self, datafile):
        self.name = 'Catalog of objects'
        self.datafile = datafile
#         self.catalog = atpy.Table(datafile, type='ascii')
#         print "Read in table, length",len(self.catalog)
        return None

   def __str__(self):
#         return 'Lightcone extracted from %s, of radius %.2f arcmin, centred on (%.2f,%.2f) deg' % (datafile,radius,position)
        return 'Lightcone extracted from %s' % (datafile)

# ============================================================================
    
# Drill out the narrow pencil beam defined by position vector xc and radius.
class lightcone:

   def __init__(self, catalog, position, radius):

        self.name = 'Lightcone through the observed Universe'
        self.master = catalog
        self.xc = position
        self.rmax = radius
        
        # Drill out galaxies in a box centred on xc:
        dx = self.rmax*arcmin2rad
        index = numpy.where(self.catalog.x > xc[0]-dx and
                            self.catalog.x < xc[0]+dx and
                            self.catalog.y > xc[1]-dx and
                            self.catalog.y < xc[1]+dx)
        galaxies = self.catalog[index]
        # Recentre the coordinate system on the cone centroid, and 
        # convert to arcmin:
        galaxies.add_empty_column('r')
        galaxies.x = (galaxies.x - xc[0])*rad2arcmin
        galaxies.y = (galaxies.y - xc[1])*rad2arcmin
        galaxies.r = numpy.sqrt(galaxies.x*galaxies.x + galaxies.y*galaxies.y)
        index = numpy.where(galaxies.r < self.rmax)
        self.galaxies = objects[index]
        
        return None


   def __str__(self):
        return 'Lightcone of radius %.2f arcmin, centred on (%.2f,%.2f) deg' % (datafile,radius,position)


   def plot(self):
       plt.scatter(self.galaxies.x, self.galaxies.y, c='k', marker='.')
       plt.xlabel('x (arcmin)')
       plt.ylabel('y (arcmin)')
       plt.axes().set_aspect('equal')
       plt.title('%s' % self)
       return None

# ----------------------------------------------------------------------------
#     
#     def selectlens()
# 
# # ----------------------------------------------------------------------------
#     
#     def beta(i,j,k):  
#         if j>i:
#             R1 = D.Da(i,j)/D.Da(j)
#             R2 = D.Da(i,k)/D.Da(k)
#             return R2/R1
#         if i>j:
#             R1 = D.Da(j,i)/D.Da(i)
#             R2 = D.Da(j,k)/D.Da(k)
#             return R1/R2
# 
# # ----------------------------------------------------------------------------
#     
#     def Kappaindiv()
# 
# # ----------------------------------------------------------------------------
#     
#     def Shearindiv()
# 
# 
# 
# # ----------------------------------------------------------------------------
#     
#     def Kappa_keeton(i=zl,j,k=zs,):
#         B=beta(i,j,k)
#         K=Kappaindiv()
#         G=Shearindiv()
#         return (1.-B)*(K-B(K^2-G^2))/((1-BK)^2-(G)^2)
# 
#     def Kappa_sum(i,j,k):
# 
# 
# ============================================================================

# TESTS:

def test(catalog):

    print "Initialising lightcone data..."
    
    xc = [0.0, 0.0]
    rmax = 2 # arcmin
    lc = lightcone(catalog,xc,rmax)

    print "Making halos..."
    lc.make_halos()

    print "Making stellar masses..."
    lc.make_galaxy_stellar_masses()

    print "Making stellar masses..."

    print "Plotting positions..."
    lc.plot()
    
    return
   
# ============================================================================

if __name__ == '__main__':

    datafile = "../../data/GGL_los_8_0_0_1_1_N_4096_ang_4_STARS_SA_galaxies_ANALYTIC_SA_galaxies_on_plane_27_to_63.images.txt"
    master = atpy.Table(datafile, type='ascii')
    print "Read in master table, length",len(master)
    
    test(master)
    
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
