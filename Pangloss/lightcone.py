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
D.h = 0.7

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
        self.Da_l = D.Da(zlens)
        self.Da_s = D.Da(zsource)
       
        # Drill out galaxies in a box centred on xc:
        dx = self.rmax*arcmin2rad
        self.galaxies = self.catalog.where((self.catalog['pos_0[rad]'] > (self.xc[0]-dx)) & \
                                           (self.catalog['pos_0[rad]'] < (self.xc[0]+dx)) & \
                                           (self.catalog['pos_1[rad]'] > (self.xc[1]-dx)) & \
                                           (self.catalog['pos_1[rad]'] < (self.xc[1]+dx))) # & \
                                           # zspec< z_source

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
#   # Function needed to calculate kappa and gamma for an NFW halo.
#   def Ffunc(self,x):
#       z=numpy.zeros(len(x))
#       for i in range(len(x)):
#          if x[i]>1:
#             y=(x[i]**2-1)**.5
#             z[i]= (1./y)*numpy.arctan(y)
#          else: 
#             print "WARNING You are very close to a halo"
#             if x[i]==1:
#                z[i] =1
#             else:
#                y=(-x[i]**2+1)**.5
#                z[i] = (1./y)*numpy.arctanh(y)
#       #print x
#       #print z
#       return z
#
# ----------------------------------------------------------------------------
   # Function needed to calculate kappa for an NFW halo. 
   # Note that Ffunc and Ffunc2 are identical.
   def Ffunc2(self,x):
       z=numpy.zeros(len(x))
       for i in range(len(x)):
          if x[i]>1:
             z[i]= 1-(2./(x[i]**2-1)**.5)*numpy.arctan(((x[i]-1)/(x[i]+1))**.5)
          else: 
             print "WARNING You are very close to a halo"
             if x[i]==1:
                z[i] =1
             else:
                y=(-x[i]**2+1)**.5
                z[i] = 1-(2./(1-x[i]**2)**.5)*numpy.arctan(((1-x[i])/(x[i]+1))**.5)
       #print x
       #print z
       return z

# ----------------------------------------------------------------------------
   # Function needed to calculate gamma for an NFW halo. 
   # Form is ridiculously long, but follows http://arxiv.org/pdf/astro-ph/9908213v1.pdf
   def Gfunc(self,x):
       z=numpy.zeros(len(x))
       for i in range(len(x)):
          X=x[i]
          if x[i]>1:
             y=(((X-1)/(X+1))**.5)
             z[i]= (8* numpy.arctan(y) / (X**2*(X**2-1)**0.5)) +\
                 (4/X**2)*numpy.log(X/2) - \
                 2/(X**2-1) +\
                 4*numpy.arctan(y)/((X**2-1)**(3/2))

          else: 
             #print "WARNING You are very close to a halo"
             if x[i]==1:
                z[i] =(10./3+4*numpy.log(0.5))
             else:
                y=(((1-X)/(X+1))**.5)
                z[i]= (8* numpy.arctanh(y) / (X**2*(1-X**2)**0.5)) +\
                    (4/X**2)*numpy.log(X/2) - \
                    2/(X**2-1) +\
                    4*numpy.arctanh(y)/((X**2-1)*(1-X**2)**(1/2))
       #print x
       #print z
       return z

# ----------------------------------------------------------------------------

   def SigmaCrit(self,zl,zs): 
   # NOTE zl here is the lensing object NOT necessarily the primary lens
       return (1.663*10**18)*(self.Da_s/(self.galaxies.Da*self.galaxies.Da_tosource)) # numerical factor is c^2/(4 pi G) in Solarmasses per megaparsec

# ----------------------------------------------------------------------------

   def MCrelation(self,M200):
       c_200 = 4.67*(M200/(10**14))**0.11 #Neto et al. equation 5
       return c_200

# ----------------------------------------------------------------------------

   def delta_c(self,c):
       return (200./3)*(c**3)/(numpy.log(1+c)-c/(1+c))

# ----------------------------------------------------------------------------

   def Hsquared(self,z):
       H0 =D.h*3.241*10**-18
       Hsq=(H0**2)*(D.OMEGA_M*(1+z)**3+(1-D.OMEGA_M)) #Lambda CDM only at this stage
       return Hsq
    
# ----------------------------------------------------------------------------

   def rho_crit_univ(self,z):   #critical density of the universe at z
       ro= 2.642*10**46**self.Hsquared(z) #units of solar mass per cubic megaparsec, H(z) must be in units of persecond.
       return ro
 
# ----------------------------------------------------------------------------

   # Beta parameter for a perturber at j:  
   def beta(self,i,j,k):  
      if j>k:
         print "z_pert > z_source? you shouldn't be asking for this"
      if j>i:
         R1 = D.Da(i,j)/D.Da(j)
         R2 = D.Da(i,k)/D.Da(k)
         return R2/R1
      if i>j:
         R1 = D.Da(j,i)/D.Da(i)
         R2 = D.Da(j,k)/D.Da(k)
         return R2/R1

# ----------------------------------------------------------------------------

   # Kappa Keeton, following Keeton (2003)and Momcheva et al. (2006)
   def KappaKeeton(self,zl,zd,zs,kappa,shear):
      output = numpy.zeros(len(zd))
      for i in range(len(zd)):
       if zd[i] < zs:
         B=self.beta(zl,zd[i],zs)
         K=kappa[i]
         G=shear[i]
         output[i] = (1.-B)*(K-B*(K**2-G**2))/((1-B*K)**2-(B*G)**2)
       else: 
         output[i]= 0.0
      return output

# ----------------------------------------------------------------------------

   # NFW model for halo
   def make_kappa_contributions(self):

       # Compute distance to each galaxy 
       # (note Da is not a function of an array):
       zd = self.galaxies['z_spec']
       Da = numpy.zeros(len(zd))
       Da_tosource = numpy.zeros(len(zd))
       for i in range(len(zd)):
         Da[i] = D.Da(zd[i])
         Da_tosource[i] = D.Da(zd[i])
       self.galaxies.add_column('Da',Da)
       self.galaxies.add_column('Da_tosource',Da_tosource)
       rphys=self.galaxies.r*Da  # Mpc
       M200 = self.galaxies['M_Halo[M_sol/h]']
       
       # Compute NFW quantities, and store for later:
       c200 = self.MCrelation(M200)
       self.galaxies.add_column('c200',c200)     

       
       #rs = (M200/(4*3.14159*self.rho_crit_univ(zd))*(numpy.log(1+c200)-c200/(1+c200)))**(1./3) #units: megaparsecs
       
       r200 = (3*M200/(800*3.14159*self.rho_crit_univ(zd)))**(1./3) #units: megaparsecs         #http://arxiv.org/pdf/astro-ph/9908213v1.pdf
       rs = r200/c200                                                                           #(wright and brainerd)
       ### Is this right - I'm not 100% clear on what the concentration 
       ### parameter does... PJM: yes, this is right! c200 is about 4 for 10^15Msun clusters

       rhos = self.delta_c(c200)*self.rho_crit_univ(zd)  # units: solar mass per cubic megaparsec
       
       self.galaxies.add_column('rs',rs)
       self.galaxies.add_column('rhos',rhos)

       x = rphys/rs
       sigmacrit = self.SigmaCrit(zd,self.zs)  # units: Solarmasses per megaparsec^2
       self.galaxies.add_column('SigmaCrit',sigmacrit)

       kappas = rhos*rs/sigmacrit
       # kappa = 2*kappas*(1-self.Ffunc(x))/(x**2-1)
           #following keeton's catalog of mass functions.
       kappa = 2*kappas*(1-self.Ffunc2(x))/(x**2-1)  # Note that Ffunc and Ffunc2 are identical. (but written differently in different references)
           #fololowing http://arxiv.org/pdf/astro-ph/9908213v1.pdf

       shear = kappas * self.Gfunc(x)

       # kappa or shear is zero if behind the source!
       for i in range(len(zd)):
          if zd[i] > self.zs:
             kappa[i] = 0.0
             shear[i] = 0.0

       kappa_keeton=self.KappaKeeton(self.zl,zd,self.zs,kappa,shear)

       # Contributions to simple weighted sum (Keeton 2003):
       self.galaxies.add_column('kappa',kappa)
       self.galaxies.add_column('gamma',shear)
       self.galaxies.add_column('kappa_keeton',kappa_keeton)
      
       return None

# ----------------------------------------------------------------------------
   # 2-panel plot, showing view from Earth and also line of sight section:
   def plot(self,starlight=False,dmglow=False,kappa_indiv=False,kappa_keeton=False):

       # Galaxy positions:
       ax1=plt.subplot(2,1,1, aspect ='equal')
       plt.subplot(2,1,1,aspect='equal')
       empty = True

       # plot different properties depending on options:
       if dmglow==True:
         plt.scatter(self.galaxies.x, self.galaxies.y, c='k', marker='o',s=((numpy.log(self.galaxies['M_Halo[M_sol/h]']))/3) )
         empty = False
       if starlight==True:
         plt.scatter(self.galaxies.x, self.galaxies.y, c='y', marker='o',s=((numpy.log(self.galaxies['M_Stellar[M_sol/h]']))/2) )     
         empty = False
       if kappa_indiv:
         plt.scatter(self.galaxies.x, self.galaxies.y, c='r', marker='o',s=(self.galaxies.kappa)*10**12, edgecolor = 'none')     
         empty = False
       if kappa_keeton:
         plt.scatter(self.galaxies.x, self.galaxies.y, c='r', marker='o',s=(self.galaxies.kappa_keeton)*10**12, edgecolor = 'none')     
         empty = False
       if empty:
         plt.scatter(self.galaxies.x, self.galaxies.y, c='k', marker='o',s=1)      

       # Lightcone boundary and centroid:
       circ=pylab.Circle(self.xc,radius=self.rmax,fill=False,linestyle='dotted')
       ax1.add_patch(circ)
       ax1.plot([self.xc[0]],[self.xc[1]], c='k', marker='+',markersize=10)


       # Labels:
       plt.xlabel('x / arcmin')
       plt.ylabel('y / arcmin')

       #axis limits
       ax1.axis([self.xc[0]-self.rmax-0.1,self.xc[0]+self.rmax+0.1,self.xc[1]-self.rmax-0.1,self.xc[1]+self.rmax+0.1])

       #plt.title('%s' % self)


#      -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.


#      #Tom Collett's subplot for view along redshift axis
#      ax2=plt.subplot(2,1,2)
#
#       ax2.scatter(self.galaxies['z_spec'],(self.galaxies.x**2+self.galaxies.y**2)**.5, c='k', marker='o',s=1)
#       plt.xlabel('Redshift')
#       plt.ylabel('LoS distance / arcmin')
#       zmax = self.galaxies['z_spec'].max()
#
#       if self.zs > zmax:
#          ax2.axis([0,self.zs+0.1,0,self.rmax+0.1])
#       else:
#          ax2.axis([0,zmax+0.1,0,self.rmax+0.1])      
#       #add lines for source and lens plane
#       ax2.axvline(x=self.zl, ymin=0, ymax=1,color='black', ls='dashed',label='bla')
#       ax2.axvline(x=self.zs, ymin=0, ymax=1,color='black', ls='solid')
#
#
#
#
#      plt.axis([self.xc[0]-self.rmax-0.1,self.xc[0]+self.rmax+0.1,self.xc[1]-self.rmax-0.1,self.xc[1]+self.rmax+0.1])
#       # plt.title('%s' % self)
#

#      -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.


       # Phil Marshall's subplot  View along redshift axis: 

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
    lc.plot(starlight=False,dmglow=False, kappa_indiv=True)

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
