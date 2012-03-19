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
import distances
from mpl_toolkits.axes_grid1 import ImageGrid

#import time
#t0=time.clock()    

D = distances.Distance()


arcmin2rad = (1.0/60.0)*numpy.pi/180.0
rad2arcmin = 1.0/arcmin2rad

# ============================================================================
    
# Given a catalog, drill out a narrow pencil beam defined by position vector 
# xc and some radius.

class lightcone:

   def __init__(self, catalog, position, radius, zlens,zsource):

        self.name = 'Lightcone through the observed Universe'
        self.catalog = catalog
        self.xc = position 
        self.rmax = radius
        self.zl = zlens
        self.zs = zsource

       
        # Drill out galaxies in a box centred on xc:
        dx = self.rmax*arcmin2rad
        self.galaxies = self.catalog.where((self.catalog['pos_0[rad]'] > (self.xc[0]-dx)) & \
                                           (self.catalog['pos_0[rad]'] < (self.xc[0]+dx)) & \
                                           (self.catalog['pos_1[rad]'] > (self.xc[1]-dx)) & \
                                           (self.catalog['pos_1[rad]'] < (self.xc[1]+dx)))

        # Recentre the coordinate system on the cone centroid, and 
        # convert to arcmin:
        x = (self.galaxies['pos_0[rad]'] - self.xc[0])*rad2arcmin
        y = (self.galaxies['pos_1[rad]'] - self.xc[1])*rad2arcmin
        r = numpy.sqrt(x*x + y*y)
        self.galaxies.add_column('x',x)
        self.galaxies.add_column('y',y)
        self.galaxies.add_column('r',r)
        self.galaxies = self.galaxies.where(self.galaxies.r < self.rmax)
        
        # Note that these are placeholder galaxies: they may not have any 
        # observable properties yet.
        
        return None

# ----------------------------------------------------------------------------

   def __str__(self):
        return 'Lightcone of radius %.2f arcmin, centred on (%.3f,%.3f) rad' % (self.rmax,self.xc[0],self.xc[1])

# ----------------------------------------------------------------------------

   def plot(self):
 
       #ax1.ylabel('y / arcmin')


       # Galaxy positions:
       ax1=plt.subplot(2,1,1, aspect ='equal')
       ax1.scatter(self.galaxies.x, self.galaxies.y, c='k', marker='o',s=((numpy.log(self.galaxies['M_Halo[M_sol/h]']))/3) )
       
       ax1.scatter(self.galaxies.x, self.galaxies.y, c='y', marker='o',s=((numpy.log(self.galaxies['M_Stellar[M_sol/h]']))/2) )     
 
   
       # Lightcone boundary and centroid:
       circ=pylab.Circle(self.xc,radius=self.rmax,fill=False,linestyle='dotted')
       #circ2=pylab.Circle(self.xc,radius=self.rmax/2.,fill=False,linestyle='dotted')
       ax1.add_patch(circ)
       #ax.add_patch(circ2)
       ax1.plot([self.xc[0]],[self.xc[1]], c='k', marker='+',markersize=10)

       # Labels:

       plt.xlabel('x / arcmin')
       plt.ylabel('y / arcmin')



       ax1.axis([self.xc[0]-self.rmax-0.1,self.xc[0]+self.rmax+0.1,self.xc[1]-self.rmax-0.1,self.xc[1]+self.rmax+0.1])
       #plt.title('%s' % self)


       #subplot for view along redshift axis
       ax2=plt.subplot(2,1,2)

       ax2.scatter(self.galaxies['z_spec'],(self.galaxies.x**2+self.galaxies.y**2)**.5, c='k', marker='o',s=1)
       plt.xlabel('Redshift')
       plt.ylabel('LoS distance / arcmin')
       zmax = self.galaxies['z_spec'].max()

       if self.zs > zmax:
          ax2.axis([0,self.zs+0.1,0,self.rmax+0.1])
       else:
          ax2.axis([0,zmax+0.1,0,self.rmax+0.1])      
       #add lines fof source and lens plane
       ax2.axvline(x=self.zl, ymin=0, ymax=1,color='black', ls='dashed',label='bla')
       ax2.axvline(x=self.zs, ymin=0, ymax=1,color='black', ls='solid')




       return None

# ----------------------------------------------------------------------------
#     
#     def selectlens()
# 
# # ----------------------------------------------------------------------------
   #beta parameter for a perturber at j:  
   def beta(i,j,k):  
      if j>k:
         print "z_pert > z_source?, wtf?"
         df
      if j>i:
         R1 = D.Da(i,j)/D.Da(j)
         R2 = D.Da(i,k)/D.Da(k)
         return R2/R1
      if i>j:
         R1 = D.Da(j,i)/D.Da(i)
         R2 = D.Da(j,k)/D.Da(k)
         return R2/R1
 
# # ----------------------------------------------------------------------------
    # SIS model for halo 
#    def Kappaindiv(self):
#       return
 
# # ----------------------------------------------------------------------------
    # SIS model for halo
#    def Shearindiv(self):
# 
#       return
 
# # ----------------------------------------------------------------------------
#     
#     def Kappa_keeton(i=zl,j,k=zs,):
#       if j<k:
#         B=beta(i,j,k)
#         K=Kappaindiv()
#         G=Shearindiv()
#         return (1.-B)*(K-B(K^2-G^2))/((1-BK)^2-(G)^2)
#       else return 0.0
#
#     def Kappa_sum(i,j,k):
# 
# 
# ============================================================================

# TESTS:

def test(catalog):

    plt.clf()
    print "Initialising lightcone data..."
    zl=0.3
    zs=2.0
    #xc = [-1*arcmin2rad,-1*arcmin2rad]
    xc = [-.004, -0.010] # radians
    rmax = 2 # arcmin
    lc = lightcone(catalog,xc,rmax,zl,zs)

    print "Plotting positions of objects..."
    lc.plot()

#     print "Distributing dark matter in halos..."
#     lc.make_dmhalos()
#     lc.plot(dmhalos=True)
# 
#     print "Distributing stars in galaxies..."
#     lc.make_galaxies()
#     lc.plot(galaxies=True)
# 
#     print "Computing Keeton (2003) convergence at optical axis..."
#     lc.keeton_kappa()


    pngfile = 'test.png'
    plt.savefig(pngfile)
    print "Plot saved in",pngfile
    
    return
   
# ============================================================================

if __name__ == '__main__':

# Simulated (Millenium) lightcone data from Hilbert:
#     datafile = "../../data/GGL_los_8_0_0_1_1_N_4096_ang_4_STARS_SA_galaxies_ANALYTIC_SA_galaxies_on_plane_27_to_63.images.txt"

# First 99 lines of the same, for testing:
    datafile = "test.txt"
    
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
