'''
This file is part of the Pangloss project.
Copyright 2012 Tom Collett (IoA) and Phil Marshall (Oxford).

to-do
-----
- finish basic functions
- carry out curve of growth test
- decide whether to continue
'''
###Things that I think are flawed:###
# Starlight lensing
# --> Lines of sight close to halos.







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
c = 299792458.
G = 4.3e-6


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


#### there is a bug in here somewhere (I _always_ get very small negative kappa_keetons, this can't be right can it?). Functions I haven't checked are marked "####"


class lightcone:

   def __init__(self, catalog,radius,zsource, position=[], zlens=[], lensindex=-1):

        self.name = 'Lightcone through the observed Universe'
        self.catalog = catalog
        # selects a lens, if not already given a lens index.
        if position==[]: #1-4) inside the light cone. 5-6) zspec sensible. 7-8) halo mass. 9) r band colour cut (currently disabled)
           if lensindex <0: print "Choosing a lens..."
           else: print "evaluating for lens %i" % lensindex
           xmax=self.catalog['pos_0[rad]'].max()-radius*arcmin2rad
           xmin=self.catalog['pos_0[rad]'].min()+radius*arcmin2rad
           ymax=self.catalog['pos_0[rad]'].max()-radius*arcmin2rad
           ymin=self.catalog['pos_0[rad]'].min()+radius*arcmin2rad
           self.potential_lenses = self.catalog.where((self.catalog['pos_0[rad]'] > xmin) & \
                                                         (self.catalog['pos_0[rad]'] < xmax) & \
                                                         (self.catalog['pos_1[rad]'] > ymin) & \
                                                         (self.catalog['pos_1[rad]'] < ymax) & \
                                                         (self.catalog['z_spec']<0.6)& \
                                                         (self.catalog['z_spec']>0.3)& \
                                                         (self.catalog['M_Stellar[M_sol/h]']>10**(10))& \
                                                         (self.catalog['M_Stellar[M_sol/h]']<10**(12)))#& (self.catalog['mag_SDSS_r'] < 21.5))


           #print len(self.potential_lenses['z_spec'])

           rndnum=rnd.randint(0,len(self.potential_lenses['z_spec']))
           if lensindex >-1: rndnum=lensindex
           self.lens=self.potential_lenses[rndnum]
           x=self.lens['pos_0[rad]']
           y=self.lens['pos_1[rad]']
           zlens=self.lens['z_spec']
           position =[x,y]
           if lensindex <0:
              print "Lens %i at [%.4f,%.4f], redshift of z=%.2f and log([stellar mass, halo mass]/M_sun) = [%.2f, %.2f]" % \
                  (rndnum, self.lens['pos_0[rad]'],self.lens['pos_1[rad]'],self.lens['z_spec'],self.lens['M_Stellar[M_sol/h]'],self.lens['M_Halo[M_sol/h]'])



        # Calculate/save some necessary parameters
        self.rmax = radius
        self.Da_l = D.Da(zlens)
        self.Da_s = D.Da(zsource)
        self.xc = position 
        self.zl = zlens
        self.zs = zsource
        self.Da_l = D.Da(zlens)
        self.Da_s = D.Da(zsource)

    
        # Drill out galaxies in a box centred on the lens, insist they are lower redshift than the source, and remove the lens itself:
        dx = self.rmax*arcmin2rad
        self.galaxies = self.catalog.where((self.catalog['pos_0[rad]'] > (self.xc[0]-dx)) & \
                                           (self.catalog['pos_0[rad]'] < (self.xc[0]+dx)) & \
                                           (self.catalog['pos_1[rad]'] > (self.xc[1]-dx)) & \
                                           (self.catalog['pos_1[rad]'] < (self.xc[1]+dx)) & \
                                           (self.catalog['z_spec'] < zsource ) & \
                                           (self.catalog['z_spec'] != self.lens['z_spec']) & \
                                           (self.catalog['pos_0[rad]'] != self.lens['pos_0[rad]']) & \
                                           (self.catalog['pos_1[rad]'] != self.lens['pos_1[rad]']))

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



# ----------------------------------------------------------------------------
   # Function needed to calculate kappa for an NFW halo. 

   def Ffunc(self,x):
       z=numpy.zeros(len(x))
       for i in range(len(x)):
          if x[i]>1:
             z[i]= (1-(2./(x[i]**2-1)**.5)*numpy.arctan(((x[i]-1.)/(x[i]+1))**.5))/(x[i]**2-1.)
          else: 
             #print "WARNING You are very close to a halo"
             if x[i]==1:
                z[i] =1./3
             else:
                y=(-x[i]**2+1)**.5
                z[i] = (1.-(2./(1-x[i]**2)**.5)*numpy.arctanh(((1.-x[i])/(x[i]+1))**.5))/(x[i]**2-1)


          if z[i] < 0: print 'warning Ffunc'
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
                 4*numpy.arctan(y)/(((X**2)-1)**(3./2))
          else: 
             #print "WARNING You are very close to a halo"
             if x[i]==1:
                z[i] =(10./3+4*numpy.log(0.5))
             else:
                y=(((1-X)/(X+1))**.5)
                z[i]= (8* numpy.arctanh(y) / (X**2*(1-X**2)**0.5)) +\
                    (4/X**2)*numpy.log(X/2) - \
                    2/(X**2-1) +\
                    4*numpy.arctanh(y)/((X**2-1)*(1-X**2)**(1./2))
          if z[i]<0: print 'warning Gfunc'; print x[i]; print z[i]

       return z

# ----------------------------------------------------------------------------

   def SigmaCrit(self,zl,zs): 
   # NOTE zl here is the lensing object NOT necessarily the primary lens
       return (1.663*10**18)*(self.Da_s/(self.galaxies.Da*self.galaxies.Da_tosource)) 
               # ^ numerical factor is c^2/(4 pi G) in Solarmasses per megaparsec

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
       Hsq=(H0**2.)*(D.OMEGA_M*(1.+z)**3.+(1.-D.OMEGA_M)) #Flat LambdaCDM only at this stage
       return Hsq
    
# ----------------------------------------------------------------------------

   def rho_crit_univ(self,z):   #critical density of the universe at z
       ro= (2.642*10**46)*self.Hsquared(z) #units of solar mass per cubic megaparsec, H(z) must be in units of per second.
       return ro 
 
# ----------------------------------------------------------------------------

   # Beta parameter for a perturber at j:  
   def beta(self,i,j,k):  
      if j>k:
         print "z_pert > z_source? you shouldn't be asking for this"
      if j>i:
         R1a = D.Da(i,j)/D.Da(j)
         R2a = D.Da(i,k)/D.Da(k)
         return R2a/R1a
      if i>j:
         R1b = D.Da(j,i)/D.Da(i)
         R2b = D.Da(j,k)/D.Da(k)
         return R2b/R1b

# ----------------------------------------------------------------------------

   # Kappa Keeton, following Keeton (2003)and Momcheva et al. (2006)
   def KappaKeeton(self,zl,zd,zs,kappa,shear):
      output = numpy.zeros(len(zd))
      for i in range(len(zd)):
       if zd[i] < zs:
         B=self.beta(zl,zd[i],zs)
         K=kappa[i]
         G=shear[i]
         output[i] = (1.-B)   * (K-    B*(K**2-G**2)   )  /    ( (1-B*K)**2   - (B*G)**2   )
       else: 
         output[i]= 0.0
      return output

# ----------------------------------------------------------------------------

   # NFW model for halo
   def make_kappa_contributions(self): ####

       # Compute distance to each galaxy 
       zd = self.galaxies['z_spec']
       Da = numpy.zeros(len(zd))
       Da_tosource = numpy.zeros(len(zd))
       for i in range(len(zd)):
         Da[i] = D.Da(zd[i])
         Da_tosource[i] = D.Da(zd[i])
       self.galaxies.add_column('Da',Da)
       self.galaxies.add_column('Da_tosource',Da_tosource)
       rphys=self.galaxies.r*Da  # Mpc

       # Compute NFW quantities, and store for later:
       M200 = self.galaxies['M_Halo[M_sol/h]']
       c200 = self.MCrelation(M200)
       self.galaxies.add_column('c200',c200)     
       rho=self.rho_crit_univ(zd)
       self.galaxies.add_column('rho_crit',rho)
       r200 = (3*M200/(800*3.14159*self.galaxies.rho_crit))**(1./3) #units: megaparsecs         #http://arxiv.org/pdf/astro-ph/9908213v1.pdf
       rs = r200/c200                                                                           #(wright and brainerd)
       rhos = self.delta_c(c200)*self.galaxies.rho_crit # units: solar mass per cubic megaparsec

       self.galaxies.add_column('rs',rs)
       self.galaxies.add_column('rhos',rhos)

       x = rphys/(200*rs)
       sigmacrit = self.SigmaCrit(zd,self.zs)  # units: Solarmasses per megaparsec^2
       self.galaxies.add_column('SigmaCrit',sigmacrit)

       kappas = rhos*rs/sigmacrit

       kappaNFW = 2.*kappas*(self.Ffunc(x)) #following http://arxiv.org/pdf/astro-ph/9908213v1.pdf

       shearNFW = kappas * self.Gfunc(x)
       
       kappa = kappaNFW   #these should include a component from starlight too
       shear = shearNFW

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
       #print self.galaxies.gamma
       #print self.galaxies.kappa
       #print self.galaxies.kappa_keeton
       return None

# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------

   # 2-panel plot, showing view from Earth and also line of sight section:
   def plot(self,starlight=False,dmglow=False,kappa_indiv=False,kappa_keeton=False,observed_light=True):
       scale=numpy.max([ numpy.absolute((numpy.min(self.galaxies.kappa_keeton))), numpy.max(self.galaxies.kappa_keeton)])/200
       scale2= (numpy.max(self.galaxies.kappa))/200
       # Galaxy positions:
       ax1=plt.subplot(2,1,1, aspect ='equal')
       plt.subplot(2,1,1,aspect='equal')
       empty = True

       # plot different properties depending on options:
       if dmglow==True:
         plt.scatter(self.galaxies.x, self.galaxies.y, c='k', marker='o',s=((numpy.log(self.galaxies['M_Halo[M_sol/h]']))/3),edgecolor = 'none' )
         empty = False
       if starlight==True:
         plt.scatter(self.galaxies.x, self.galaxies.y, c='y', marker='o',s=((numpy.log(self.galaxies['M_Stellar[M_sol/h]']))/2),edgecolor = 'none' )     
         empty = False
       if kappa_keeton:
         for galaxy in self.galaxies:
            kappa = galaxy['kappa_keeton']
            if kappa <0:
               plt.scatter(galaxy['x'], galaxy['y'], c='b', marker='o',s=-kappa/scale , edgecolor = 'none')    
            else:
               plt.scatter(galaxy['x'], galaxy['y'],  c='r', marker='o',s=kappa/scale , edgecolor = 'none')
         empty = False
       if kappa_indiv: 
         plt.scatter(self.galaxies.x, self.galaxies.y, c='none', marker='o',s=(self.galaxies.kappa)/scale2)     
         empty = False
       if observed_light==True:
         plt.scatter(self.galaxies.x, self.galaxies.y, c='y', marker='o',s=2**(25-(self.galaxies['mag_SDSS_r'])),edgecolor = 'none' )     
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
       ""

       plt.subplot(2,1,2)
       empty = True
       

       if dmglow==True:
         plt.scatter(self.galaxies.z_spec, self.galaxies.r, c='k', marker='o',s=((numpy.log(self.galaxies['M_Halo[M_sol/h]']))/3),edgecolor = 'none' )
         empty = False
       if starlight==True:
         plt.scatter(self.galaxies.z_spec, self.galaxies.r, c='y', marker='o',s=((numpy.log(self.galaxies['M_Stellar[M_sol/h]']))/2),edgecolor = 'none' )     
         empty = False
       if kappa_keeton:
         for galaxy in self.galaxies:
            kappa = galaxy['kappa_keeton']
            if kappa <0:
               plt.scatter(galaxy['z_spec'], galaxy['r'], c='b', marker='o',s=-kappa/scale , edgecolor = 'none')    
            else:
               plt.scatter(galaxy['z_spec'], galaxy['r'],  c='r', marker='o',s=kappa/scale , edgecolor = 'none')
         empty = False
       if kappa_indiv: 
         plt.scatter(self.galaxies.z_spec, self.galaxies.r, c='none', marker='o',s=(self.galaxies.kappa)/scale2)     
         empty = False
       if observed_light==True:
         plt.scatter(self.galaxies.z_spec, self.galaxies.r, c='y', marker='o',s=2**(25-(self.galaxies['mag_SDSS_r'])),edgecolor = 'none' )     
         empty = False


       if empty:
         plt.scatter(self.galaxies.z_spec, self.galaxies.r, c='k', marker='o',s=1)


       plt.xlabel('Redshift')
       plt.ylabel('LoS distance / arcmin')
       zmax = self.galaxies['z_spec'].max()


       plt.axis([0,self.zs+0.1,0,self.rmax+0.1])

       #add lines for source and lens plane
       plt.axvline(x=self.zl, ymin=0, ymax=1,color='black', ls='dashed',label='bla')
       plt.axvline(x=self.zs, ymin=0, ymax=1,color='black', ls='solid')       # plt.title('%s' % self)

       ""

#      -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.


       # Phil Marshall's subplot  View along redshift axis: 

       #TC: I prefered the old plot! y position gives the confusing impression that some objects 
       #    are much closer to the optical axis than they really are.
       """
       plt.subplot(2,1,2)
       empty = True
       if dmglow:
         plt.scatter(self.galaxies['z_spec'],self.galaxies.y, c='k', marker='o',s=((numpy.log(self.galaxies['M_Halo[M_sol/h]']))/3) )
         empty = False
       if starlight:
         plt.scatter(self.galaxies['z_spec'],self.galaxies.y, c='y', marker='o',s=((numpy.log(self.galaxies['M_Stellar[M_sol/h]']))/2) )     
         empty = False
       if kappa_indiv: 
         plt.scatter(self.galaxies.x, self.galaxies.y, c='r', marker='o',s=(self.galaxies.kappa)*10**11, edgecolor = 'none')     
         empty = False
       if kappa_keeton:
         for galaxy in self.galaxies:
            kappa = galaxy['kappa_keeton']*10**11
            if kappa <0:
               plt.scatter(galaxy['x'], galaxy['y'], c='b', marker='o',s=-kappa, edgecolor = 'none')     
            else:
               plt.scatter(galaxy['x'], galaxy['y'], c='r', marker='o',s=kappa, edgecolor = 'none')
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

       """
       return None


# # ----------------------------------------------------------------------------

   def curve_of_growth(self):
      #set up things to plot:
      args=numpy.argsort(-numpy.absolute(self.galaxies.kappa_keeton))
      ordered_kappas=numpy.take(self.galaxies.kappa_keeton,args)
      zero=numpy.zeros(1)
      cumtot=numpy.concatenate((zero,numpy.cumsum(ordered_kappas)))

      corresponding_r=numpy.concatenate((zero,numpy.take(self.galaxies.r,args)))

      n=range(len(args)+1)

      #work out limits for plot
      #maximum=ordered_kappas.max()      
      #minimum=ordered_kappas.min()
      #if numpy.absolute(maximum) < numpy.absolute(minimum):
      #   max = minimum
      #   min = numpy.max([maximum,0])
      #else: 
      #   max = maximum
      #   min = numpy.min([minimum,0])

      #choose first m points
      n1=n[:]
      cumtot1=cumtot[:]
      corresponding_r1=corresponding_r[:]

      #start plotting:
      plt.clf()
      plt.subplot(1,1,1)
      plt.plot(n,cumtot)
      #plt.xscale('log')
      plt.scatter(n1,cumtot1,c=corresponding_r1,s=15,edgecolor='none')
      plt.colorbar()
      plt.xlim([0,20])
      plt.axhline(y=cumtot[-1], xmin=0, xmax=1000,color='black', ls='dotted')
      plt.axhline(y=0, xmin=0, xmax=1000,color='black', ls='solid')
      #plt.ylim([1.1*min,1.1*max])

# ============================================================================

# TESTS:

def test1(catalog):

    plt.clf()
    rmax = 2

    zs=2.0 #should be selecting a source redshift (possibly use forecaster from Collett et al. 2012)


    #xc = [] # radians, leave as[] to select a lens at random
    #zl=[] # leave as [] to select a lens at random
    
    print "Initialising lightcone data..."
    lc = lightcone(catalog,rmax,zs,lensindex=-1)


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
    
    print "Total external convergence =",numpy.sum(lc.galaxies.kappa_keeton)
    # Now make an illustrative plot:
    
    print "Plotting objects in lightcone..."
    lc.plot(starlight=False,dmglow=False, kappa_indiv=False, kappa_keeton=True,observed_light=True)
    

    pngfile = 'test.png'
    plt.savefig(pngfile)
    print "Plot saved in",pngfile

    print "Plotting curve of growth..."
    plt.clf()
    lc.curve_of_growth()
    pngfile = 'curve_of_growth.png'
    plt.savefig(pngfile)
    print "Plot saved in",pngfile
    
    #print lc.galaxies.kappa/lc.galaxies.gamma
    
    
    return

#-------------------------------------------------------------------------

def test2(catalog): #plots a kappa distribution:
###not finished###
#we could consider feeding in our lens selection parameters here)


    plt.clf()
    rmax = 2

    zs=1.39 #should be selecting a source redshift (possibly use forecaster from Collett et al. 2012)
    xc = [] # radians, leave as[] to select a lens at random
    zl=[] # leave as [] to select a lens at random


    iterations=85
    K=numpy.zeros(iterations)
    for j in range(iterations):
       i = j+0 #systematically evaluate for all j
       #i = rnd.randint(0,13594) # or pick random sample
       lc = lightcone(catalog,rmax,zs,xc,zl,lensindex=i)
       lc.make_kappa_contributions()
       K[j]=numpy.sum(lc.galaxies.kappa_keeton)
       #if K[j]>0.01: print j
    #for j in range(iterations):      
       #if K[j]>0.01: print j
    plt.hist(K,bins=30)
    #plt.yscale('log')
    pngfile = 'Kappa_keeton_distribution.png'
    plt.savefig(pngfile)
    plt.ylim([0,500])
    print "Plot saved in",pngfile

    
    return

#-------------------------------------------------------------------------

#def test3(catalogue):


# ============================================================================

if __name__ == '__main__':

# Simulated (Millenium) lightcone data from Hilbert:
    datafile = "../../data/GGL_los_8_0_0_1_1_N_4096_ang_4_STARS_SA_galaxies_ANALYTIC_SA_galaxies_on_plane_27_to_63.images.txt"

# First 99 lines of the same, for testing:
#    datafile = "test.txt"
    
    master = atpy.Table(datafile, type='ascii')
    print "Read in master table, length",len(master)
    
    test1(master)
    
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
