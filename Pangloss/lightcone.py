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

# PJM: I am not sure that lightcone objects should necessarily have the 
# properties zl and zs. lens_lightcones, yes, but lightcones no. The logical
# way to code this would be to define a new class, lens_lightcone, and have it
# *inherit* the properties (parameters and methods) of the lightcone object.
# We should try and do this, in order to learn how inheritance works!

#TC: Good idea. But I have no idea how to do it! 

#Note to Phil. Anything marked '###' is still unfinished/preliminary and probably wrong!

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
 
# ----------------------------------------------------------------------------
   #Function needed to calculate kappa and gamma for an NFW halo.
   def Ffunc(x):
       if x>1:
          y=(x**2-1)**.5
          return (1./y)*numpy.arctan(y)
       else: 
          print "WARNING You are very close to a halo"
          if x==1:
             return 1.
          else:
             y=(-x**2+1)**.5
             return (1./y)*numpy.arctanh(y)

# ----------------------------------------------------------------------------
   def SigmaCrit(zl,zs): #NOTE zl here is the lensing object NOT necessarily the primary lens
       return (1.663*10**18)*(D.Da(zs)/(D.Da(zl)*D.Da(zl,zs))) # numerical factor is c^2/(4 pi G) in Solarmasses per megaparsec
# ----------------------------------------------------------------------------
   def MCrelation(M_200):
       c_200 = 4.67*(M_200/(10**14))**0.11 #Neto et al. equation 5
       return c_200
# ----------------------------------------------------------------------------
   def delta_c(c):
       return (200./3)*c^3/(numpy.log(1+c)-c/(1+c))
# ----------------------------------------------------------------------------
   def Hsquared(z):
       H0 =D.h*3.241*10**-18
       Hsq=(H0**2)*(D.Omega_M*(1+z)**3+(1-D.Omega_M)) #Lambda CDM only at this stage
       return Hsq
    
# ----------------------------------------------------------------------------
   def rho_crit_univ(z):   #critical density of the universe at z
       ro= 2.642*10**46**Hsquared(z) #units of solar mass per cubic megaparsec, H(z) must be in units of persecond.
       return ro
# ----------------------------------------------------------------------------
   # NFW model for halo ###Not finished yet
   def make_kappa_contributions(self):

#        X=(self.galaxies.x**2+self.galaxies.y**2)**.5
# PJM: this is already stored as self.galaxies.r!

       zd = self.galaxies['z_spec']
       DA = D.Da(zd)
       self.galaxies.add_column('DA',DA)
       
       rphys=self.galaxies.r*DA  # Mpc
       ###Is this the right conversion to angular units?###
       # PJM: I think so! Note units. Check H0 is included...
       
       M200 = self.galaxies['M_Halo[M_sol/h]']
       
       # Compute NFW quantities, and store for later:
       c200 = MCrelation(M200)
       self.galaxies.add_column('c200',c)       
       rs = r_200/c200
       rhos = delta_c(c200)*rho_crit_univ(zd)       
       self.galaxies.add_column('rs',rs)
       self.galaxies.add_column('rhos',rhos)

       x = rphys/rs
       sigmacrit = SigmaCrit(zd,self.zs)
       kappas = rhos*rs/sigma_crit
       kappa = 2*kappas(1-Ffunc(x))/(x*x-1)
       
       # Contributions to simple weighted sum (Keeton 2003):
       self.galaxies.add_column('kappa',kappa)
       
       return None

# ----------------------------------------------------------------------------
   # 2-panel plot, showing view from Earth and also line of sight section:
   def plot(self,starlight=False,dmglow=False):
 

       #ax1.ylabel('y / arcmin')


       # Galaxy positions:
       ax1=plt.subplot(2,1,1, aspect ='equal')
       ax1.scatter(self.galaxies.x, self.galaxies.y, c='k', marker='o',s=((numpy.log(self.galaxies['M_Halo[M_sol/h]']))/3) )
       
       ax1.scatter(self.galaxies.x, self.galaxies.y, c='y', marker='o',s=((numpy.log(self.galaxies['M_Stellar[M_sol/h]']))/2) )     

       # View from Earth (needs to be larger!):
       plt.subplot(2,1,1,aspect='equal')
       empty = True
       if dmglow:
         plt.scatter(self.galaxies.x, self.galaxies.y, c='k', marker='o',s=((numpy.log(self.galaxies['M_Halo[M_sol/h]']))/3) )
         empty = False
       if starlight:
         plt.scatter(self.galaxies.x, self.galaxies.y, c='y', marker='o',s=((numpy.log(self.galaxies['M_Stellar[M_sol/h]']))/2) )     
         empty = False
       if empty:
         plt.scatter(self.galaxies.x, self.galaxies.y, c='k', marker='o',s=1)
       
#        plt.scatter(self.galaxies.x, self.galaxies.y, c='r', marker='o',s=  ( Kappaindiv(self) ))###This line is broken: why?###
       plt.scatter(self.galaxies.x, self.galaxies.y, c='r', marker='o',s=(self.galaxies.kappa) )
       # PJM: This plot script is getting messy: can you try and keep things separate,
       # and organised by options please? if this then that etc 
 
 
       # Lightcone boundary and centroid:
       circ=pylab.Circle(self.xc,radius=self.rmax,fill=False,linestyle='dotted')

       #circ2=pylab.Circle(self.xc,radius=self.rmax/2.,fill=False,linestyle='dotted')
       ax1.add_patch(circ)
       #ax.add_patch(circ2)
       ax1.plot([self.xc[0]],[self.xc[1]], c='k', marker='+',markersize=10)

       # circ2=pylab.Circle(self.xc,radius=self.rmax/2.,fill=False,linestyle='dotted')
       ax=pylab.gca()
       ax.add_patch(circ)
       # ax.add_patch(circ2)
       ax.plot([self.xc[0]],[self.xc[1]], c='k', marker='+',markersize=10)


       # Labels:
       plt.xlabel('x / arcmin')
       plt.ylabel('y / arcmin')




       ax1.axis([self.xc[0]-self.rmax-0.1,self.xc[0]+self.rmax+0.1,self.xc[1]-self.rmax-0.1,self.xc[1]+self.rmax+0.1])
       #plt.title('%s' % self)


       #Tom Collett's subplot for view along redshift axis
       ax2=plt.subplot(2,1,2)

       ax2.scatter(self.galaxies['z_spec'],(self.galaxies.x**2+self.galaxies.y**2)**.5, c='k', marker='o',s=1)
       plt.xlabel('Redshift')
       plt.ylabel('LoS distance / arcmin')
       zmax = self.galaxies['z_spec'].max()

       if self.zs > zmax:
          ax2.axis([0,self.zs+0.1,0,self.rmax+0.1])
       else:
          ax2.axis([0,zmax+0.1,0,self.rmax+0.1])      
       #add lines for source and lens plane
       ax2.axvline(x=self.zl, ymin=0, ymax=1,color='black', ls='dashed',label='bla')
       ax2.axvline(x=self.zs, ymin=0, ymax=1,color='black', ls='solid')




       plt.axis([self.xc[0]-self.rmax-0.1,self.xc[0]+self.rmax+0.1,self.xc[1]-self.rmax-0.1,self.xc[1]+self.rmax+0.1])
       # plt.title('%s' % self)


       # Phil Marshall's View along redshift axis: 

#TC: I prefered the old plot! y position gives the confusing impression that some objects 
#    are much closer to the optical axis than they really are.

       plt.subplot(2,1,2)
       empty = True
       if dmglow:
         plt.scatter(self.galaxies['z_spec'],self.galaxies.y, c='k', marker='o',s=((numpy.log(self.galaxies['M_Halo[M_sol/h]']))/3) )
         empty = False
       if starlight:
         plt.scatter(self.galaxies['z_spec'],self.galaxies.y, c='y', marker='o',s=((numpy.log(self.galaxies['M_Stellar[M_sol/h]']))/2) )     
         empty = False
       if empty:
         plt.scatter(self.galaxies['z_spec'],self.galaxies.y, c='k', marker='o',s=1)

       plt.xlabel('redshift z')
#        plt.ylabel('distance from LoS / arcmin')
       plt.ylabel('y / arcmin')
       
       zmax = max(self.galaxies['z_spec'].max(),self.zs)
#        print self.zs,self.galaxies['z_spec'].max()
#        print zmax
       plt.axis([0,zmax+0.1,-self.rmax-0.1,self.rmax+0.1])

       # Add lines marking source and lens plane, and optical axis:
       plt.axvline(x=self.zl, ymin=0, ymax=1,color='black', ls='dotted',label='bla')
       plt.axvline(x=self.zs, ymin=0, ymax=1,color='black', ls='dotted')
       plt.axhline(y=0.0, xmin=0.0, xmax=zmax, color='black', ls='dashed')


       return None

# ----------------------------------------------------------------------------
#     
#     def selectlens()
# 
# # ----------------------------------------------------------------------------
# # ----------------------------------------------------------------------------
    # NFW model for halo
#    def Shearindiv(self):
# 
#       return
 
# # ----------------------------------------------------------------------------
#     #Kappa Keeton, following Keeton (200?)and Momcheva et al. (200?)
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

#     print "Distributing dark matter in halos..."
#     lc.make_dmhalos()
#     lc.plot(dmhalos=True)
# 
#     print "Distributing stars in galaxies..."
#     lc.make_galaxies()
#     lc.plot(galaxies=True)
# 
    print "Computing Keeton (2003) convergence at optical axis due to each halo..."
    lc.make_kappa_contributions()
    
    print "Total external convergence =",numpy.sum(lc.galaxies.kappa)
    
    # Now make an illustrative plot:
    
    print "Plotting objects in lightcone..."
    lc.plot(starlight=True,dmglow=True)

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
