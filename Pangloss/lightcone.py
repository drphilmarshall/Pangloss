'''
This file is part of the Pangloss project.
Copyright 2012 Tom Collett (IoA) and Phil Marshall (Oxford).

description:
------------

Given a catalog, drill out a narrow pencil beam defined by position
vector  xc and some radius.

# NB. Anything marked '###' is still unfinished/preliminary and 
probably wrong! #### are finished functions that haven't been 
bug-checked.

to-do:
------
- finish basic functions
- carry out curve of growth test
- decide whether to continue

issues:
-------
- stellar mass is not included yet
- so lines of sight passing close to halos will not be accurate
- Born approximation will break down in this case anyway
- Keeton approximation may well be way off, need better approx from 
    Roger, Sherry and Stefan.
- lightcone and lenslightcone should be separate classes, where llc 
    inherits from lc
'''

# ======================================================================

import pylab
import matplotlib.pyplot as plt
import numpy, numpy.random as rnd, atpy
import distances
from mpl_toolkits.axes_grid1 import ImageGrid
from time import clock
import LensingProfiles as LP

#import time
#t0=time.clock()    

D = distances.Distance()
D.h = 0.7


arcmin2rad = (1.0/60.0)*numpy.pi/180.0
rad2arcmin = 1.0/arcmin2rad

vb = False

# ============================================================================
    

class lightcone:

   def __init__(self,catalog,radius,zsource,position=[],lensindex=-1,deterministic=True):

        self.name = 'Lightcone through the observed Universe'
        self.catalog = catalog
        self.deterministic = deterministic
        self.xmax = self.catalog['pos_0[rad]'].max()
        self.xmin = self.catalog['pos_0[rad]'].min()
        self.ymax = self.catalog['pos_1[rad]'].max()
        self.ymin = self.catalog['pos_1[rad]'].min() 

        self.Bin_MstarsRUN = False

        if position ==[]: 
            flag=1
        else:
            flag=0

        # selects a lens, if not already given a lens position/index.
        if position==[]: #1-4) inside the light cone. 5-6) zspec sensible. 7-8) halo mass. 9) r band colour cut (currently disabled)
           if lensindex <0: print "Choosing a lens..."
           else: print "evaluating for lens %i" % lensindex
           xmax=self.xmax-radius*arcmin2rad
           xmin=self.xmin+radius*arcmin2rad
           ymax=self.ymax-radius*arcmin2rad
           ymin=self.ymin+radius*arcmin2rad
           self.potential_lenses = self.catalog.where((self.catalog['pos_0[rad]'] > xmin) & \
                                                         (self.catalog['pos_0[rad]'] < xmax) & \
                                                         (self.catalog['pos_1[rad]'] > ymin) & \
                                                         (self.catalog['pos_1[rad]'] < ymax) & \
                                                         (self.catalog['z_spec']<0.62)& \
                                                         (self.catalog['z_spec']>0.58)& \
                                                         (self.catalog['M_Stellar[M_sol/h]']>10**(0))& \
                                                         (self.catalog['M_Stellar[M_sol/h]']<10**(15)))#& \
                                                         #(self.catalog['mag_SDSS_r'] < 21.5))


           #print len(self.potential_lenses['z_spec'])

           rndnum=rnd.randint(0,len(self.potential_lenses['z_spec']))
           if lensindex >-1: rndnum=lensindex
           self.lens=self.potential_lenses[rndnum]
           x=self.lens['pos_0[rad]']
           y=self.lens['pos_1[rad]']
           zlens=self.lens['z_spec']
           position = [x,y,zlens]
           if lensindex <0: # or lensindex >0:
              print "Lens %i at [%.4f,%.4f], redshift of z=%.2f and log([stellar mass, halo mass]/M_sun) = [%2.2f, %2.2f]" % \
                  (rndnum, self.lens['pos_0[rad]'],self.lens['pos_1[rad]'],self.lens['z_spec'],numpy.log10(self.lens['M_Stellar[M_sol/h]']),numpy.log10(self.lens['M_Halo[M_sol/h]']))



          



        # Calculate/save some necessary parameters
        zlens=position[2]
        self.rmax = radius
        self.xc = [position[0],position[1]] 
        self.zl = position[2]
        self.zs = zsource
        self.Da_l = D.Da(zlens)
        self.Da_s = D.Da(zsource)

    
        # Drill out galaxies in a box centred on the lens, insist they are lower redshift than the source and (if a lens was selected, i.e. flag == 1) remove the lens itself:
        dx = self.rmax*arcmin2rad
        if flag == 1:
           self.galaxies = self.catalog.where((self.catalog['pos_0[rad]'] > (self.xc[0]-dx)) & \
                                           (self.catalog['pos_0[rad]'] < (self.xc[0]+dx)) & \
                                           (self.catalog['pos_1[rad]'] > (self.xc[1]-dx)) & \
                                           (self.catalog['pos_1[rad]'] < (self.xc[1]+dx)) & \
                                           (self.catalog['z_spec'] < zsource ) & \
                                           (self.catalog['z_spec'] != self.lens['z_spec']) & \
                                           (self.catalog['pos_0[rad]'] != self.lens['pos_0[rad]']) & \
                                           (self.catalog['pos_1[rad]'] != self.lens['pos_1[rad]']))
        else:
           self.galaxies = self.catalog.where((self.catalog['pos_0[rad]'] > (self.xc[0]-dx)) & \
                                           (self.catalog['pos_0[rad]'] < (self.xc[0]+dx)) & \
                                           (self.catalog['pos_1[rad]'] > (self.xc[1]-dx)) & \
                                           (self.catalog['pos_1[rad]'] < (self.xc[1]+dx)) & \
                                           (self.catalog['z_spec'] < zsource ))

        # Recentre the coordinate system on the cone centroid, and 
        # convert to arcmin:
        x = (self.galaxies['pos_0[rad]'] - self.xc[0])*rad2arcmin
        y = (self.galaxies['pos_1[rad]'] - self.xc[1])*rad2arcmin
        r = numpy.sqrt(x*x + y*y)
        self.galaxies.add_column('x',x)
        self.galaxies.add_column('y',y)
        self.galaxies.add_column('r',r)
        self.galaxies = self.galaxies.where(self.galaxies.r < self.rmax)


        F814 = (self.galaxies.mag_SDSS_i + self.galaxies.mag_SDSS_z)/2. #approximate F814 colour.
        self.galaxies.add_column("mag_F814W",F814)

        # Note that these are placeholder galaxies: they may not have any 
        # observable properties yet.
        
        return None

# ----------------------------------------------------------------------------

   def __str__(self):
        return 'Lightcone of radius %.2f arcmin, centred on (%.3f,%.3f) rad' % (self.rmax,self.xc[0],self.xc[1])

# ----------------------------------------------------------------------------

   # function to calculate what kappa to subtract due to an empty cone.
   def kappa_expected(self): ####
      #print clock()

      nplanes=200
      zp,dz=numpy.linspace(0,self.zs,nplanes,endpoint=False,retstep=True)
      D_p =numpy.zeros(len(zp))
      #D_pl=numpy.zeros(len(zp))
      D_ps=numpy.zeros(len(zp))
      box_vol_p=numpy.zeros(len(zp))
      sigma_crit_p=numpy.zeros(len(zp))
      kappa_p=numpy.zeros(len(zp))
      gamma_p=numpy.zeros(len(zp)) #dummy variable - gamma is trivially zero for a sheet.

      theta=self.rmax*arcmin2rad
      for i in range(len(zp)):
         D_p[i] =D.Da(zp[i])
         #D_pl[i]=D.Da(zp[i],self.zl)
         D_ps[i]=D.Da(zp[i],self.zs)

         #comoving volume of element, is total comoving volume * solid angle/4pi: omega = 2*pi(1-cos(theta))
         #box_vol_p[i]=(D.comoving_volume((zp[i]-dz/2.),(zp[i]+dz/2.)))*((1-numpy.cos(theta)))/2
         box_vol_p[i]=(D.comoving_distance((zp[i]-dz/2.),(zp[i]+dz/2.)))

         sigma_crit_p[i]=(1.663*10**18)*(self.Da_s/(D_p[i]*D_ps[i])) 
         kappa_p[i]=(self.rho_crit_univ(zp[i])*D.OMEGA_M*box_vol_p[i])/sigma_crit_p[i] #All times Omega_halos!!!

      kappa_keeton_p=self.KappaKeeton(self.zl,zp,self.zs,kappa_p,gamma_p) 
      self.kappa_empty= numpy.sum(kappa_keeton_p)
      return numpy.sum(kappa_keeton_p)

      #print clock()

# ----------------------------------------------------------------------------
   #Tell me the number of galaxies within a certain radius, that pass a certain magnitude cut.
   def N_radius_cat(self,radius,cut=[18.5,24.5], band="F814W", radius_unit="arcsec"):
       if band == "u" or band ==  "g" or band == "r" or band ==  "i" or band == "z":
           col = "mag_SDSS_%s" % band
       elif band == "F814" or band == "F814W" or band == "814" or band == 814:
           col = "mag_F814W"
       else:
           col = "mag_%s" % band
       if radius < 10: print "Warning: Default units for N_radius are arcsec!"
       if radius_unit == "arcsec":
           radius = radius/60.
       if col != "warning":
           self.N_cut=self.galaxies.where((self.galaxies.r < radius)  & \
                                              (self.galaxies["%s"%col] < cut[1])& \
                                          (self.galaxies["%s"%col] > cut[0]))

           return self.galaxies.where((self.galaxies.r < radius) & \
                                          (self.galaxies["%s"%col] < cut[1])& \
                                          (self.galaxies["%s"%col] > cut[0]))

   def N_radius(self,radius,cut=[18.5,24.5],band="F814W", radius_unit="arcsec"):
       Ntable=self.N_radius_cat(radius,cut,band, radius_unit)
       #print len(Ntable.r)
       #plt.hist(Ntable.r)
       #plt.show()
       return len(Ntable.r)

# ----------------------------------------------------------------------------

   def SigmaCrit(self, deterministic=True): 
   # NOTE zl here is the lensing object NOT necessarily the primary lens
      if deterministic ==False:
         return (1.663*10**18)*(self.Da_s/(self.galaxies.b_Da*self.galaxies.b_Da_tosource)) 
      else: 
         return (1.663*10**18)*(self.Da_s/(self.galaxies.Da*self.galaxies.Da_tosource)) 
   #              ^ numerical factor is c^2/(4 pi G) in Solarmasses per megaparsec

# ----------------------------------------------------------------------------
   def logerr(self,l,m,s):
      c=10*rnd.normal(m,s)
      return c

   def MCrelation(self,M200,MCerror=False):
      if MCerror==False:
         c_200 = 4.67*(M200/(10**14))**0.11 #Neto et al. equation 5
      if MCerror==True:
       c_200=4.67*(M200/(10**14))**0.11
       logc_200=numpy.log10(c_200)
       lM200 = numpy.log10(M200)
       for i in range(len(M200)):    #best fit scatter parameters of neto et al (double log normal)
         if lM200[i]<11.875: ### FILLER
            f=0.205 ###FILLER
            if rnd.rand()<f:
               c_200[i]=self.logerr(logc_200[i],0.683,0.147) ###FILLER
            else:
               c_200[i]=self.logerr(logc_200[i],0.920,0.106) ###FILLER
         if lM200[i]<12.125:
            f=0.205
            if rnd.rand()<f:
               c_200[i]=self.logerr(logc_200[i],0.683,0.147)
            else:
               c_200[i]=self.logerr(logc_200[i],0.920,0.106)
         elif lM200[i]<12.375:
            f=0.171
            if rnd.rand()<f:
               c_200[i]=self.logerr(logc_200[i],0.658,0.150)
            else:
               c_200[i]=self.logerr(logc_200[i],0.903,0.108)
         elif lM200[i]<12.625:
            f=0.199
            if rnd.rand()<f:
               c_200[i]=self.logerr(logc_200[i],0.646,0.139)
            else:
               c_200[i]=self.logerr(logc_200[i],0.881,0.099)        
         elif lM200[i]<12.875:
            f=0.229
            if rnd.rand()<f:
               c_200[i]=self.logerr(logc_200[i],0.605,0.158)
            else:
               c_200[i]=self.logerr(logc_200[i],0.838,0.101)          
         elif lM200[i]<13.125:
            f=0.263
            if rnd.rand()<f:
               c_200[i]=self.logerr(logc_200[i],0.603,0.136)
            else:
               c_200[i]=self.logerr(logc_200[i],0.810,0.100)          
         elif lM200[i]<13.375:
            f=0.253
            if rnd.rand()<f:
               c_200[i]=self.logerr(logc_200[i],0.586,0.140)
            else:
               c_200[i]=self.logerr(logc_200[i],0.793,0.099)         
         elif lM200[i]<13.625:
            f=0.275
            if rnd.rand()<f:
               c_200[i]=self.logerr(logc_200[i],0.566,0.142)
            else:
               c_200[i]=self.logerr(logc_200[i],0.763,0.095)            
         elif lM200[i]<13.875:
            f=0.318
            if rnd.rand()<f:
               c_200[i]=self.logerr(logc_200[i],0.543,0.140)
            else:
               c_200[i]=self.logerr(logc_200[i],0.744,0.094)
         elif lM200[i]<14.125:
            f=0.361
            if rnd.rand()<f:
               c_200[i]=self.logerr(logc_200[i],0.531,0.131)
            else:
               c_200[i]=self.logerr(logc_200[i],0.716,0.088)            
         elif lM200[i]<14.375:
            f=0.383
            if rnd.rand()<f:
               c_200[i]=self.logerr(logc_200[i],0.510,0.121)
            else:
               c_200[i]=self.logerr(logc_200[i],0.689,0.095)       
         elif lM200[i]<14.625:
            f=0.370
            if rnd.rand()<f:
               c_200[i]=self.logerr(logc_200[i],0.490,0.133)
            else:
               c_200[i]= self.logerr(logc_200[i],0.670,0.094)           
         elif lM200[i]<14.875:
            f=0.484
            if rnd.rand()<f:
               c_200[i]=self.logerr(logc_200[i],0.519,0.121)
            else:
               c_200[i]=self.logerr(logc_200[i],0.635,0.091)    
         elif lM200[i]<15.125:
            f=0.578
            if rnd.rand()<f:
               c_200[i]=self.logerr(logc_200[i],0.493,0.094)
            else:
               c_200[i]=self.logerr(logc_200[i],0.661,0.061)
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
         return R1a/R2a
      if i>j:
         R1b = D.Da(j,i)/D.Da(i)
         R2b = D.Da(j,k)/D.Da(k)
         return R1b/R2b
      if i == j:
         return 1.0
      if j==k:
          return 0.0

# ----------------------------------------------------------------------------

   # Kappa Keeton, following Keeton (2003)and Momcheva et al. (2006)
   def KappaKeeton(self,zl,zd,zs,kappa,shear):
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
      return output

# ----------------------------------------------------------------------------
   #Function using the Behroozi M*-Mhalo relationship to recreate MHalos from M*s

   def Mstar_to_M200_Behroozi(self,M_Star,redshift,scatter=True):
      #Following Behroozi et al. 2010.
      M_200=numpy.zeros(len(M_Star))
      #parameters:
      for i in range(len(M_Star)):
         z=redshift[i]
         if z<0.9:
            Mstar00 = 10.72
            Mstar0a = 0.55
            Mstar0aa=0.0
            M_10 = 12.35
            M_1a = 0.28
            beta0 = 0.44
            betaa = 0.18
            delta0 = 0.57
            deltaa = 0.17
            gamma0 = 1.56
            gammaa = 2.51
         else:
            Mstar00 = 11.09
            Mstar0a = 0.56
            Mstar0aa= 6.99
            M_10 = 12.27
            M_1a = -0.84
            beta0 = 0.65
            betaa = 0.31
            delta0 = 0.56
            deltaa = -0.12
            gamma0 = 1.12
            gammaa = -0.53


      #scaled parameters:
         a=1./(1.+z)
         M_1=10**(M_10+M_1a*(a-1))
         beta=beta0+betaa*(a-1)
         Mstar0=10**(Mstar00+Mstar0a*(a-1)+Mstar0aa*(a-0.5)**2)
         delta=delta0+deltaa*(a-1)
         gamma=gamma0+gammaa*(a-1)

      #reltationship ****NO SCATTER****

         M_200[i] =10.0**(numpy.log10(M_1)+beta*numpy.log10(M_Star[i]/Mstar0)+((M_Star[i]/Mstar0)**delta)/(1.+(M_Star[i]/Mstar0)**-gamma)-0.5)

         if scatter==True:
             M_200[i]=10.0**(numpy.log10(M_200[i])*rnd.normal(0,0.2))
      return M_200
# ----------------------------------------------------------------------------

   # NFW model for halo
   def make_kappa_contributions(self,BehrooziHalos=False,hardcut=False,truncation=False): 
       # Compute distance to each galaxy 
       zd = self.galaxies['z_spec']
       Da = numpy.zeros(len(zd))
       Da_tosource = numpy.zeros(len(zd))
       for i in range(len(zd)):
         Da[i] = D.Da(zd[i])
         Da_tosource[i] = D.Da(zd[i],self.zs)
       self.galaxies.add_column('Da',Da)
       self.galaxies.add_column('Da_tosource',Da_tosource)
       rphys=self.galaxies.r*arcmin2rad*Da  # Mpc
       self.galaxies.add_column('rphys',rphys)     

       # ---------------------------------------------------------------------

       # Compute NFW quantities, and store for later:
       if BehrooziHalos==False:
           M200 = self.galaxies['M_Subhalo[M_sol/h]']
       if BehrooziHalos==True:
           M_star=self.galaxies['M_Stellar[M_sol/h]']
           M200 = self.Mstar_to_M200_Behroozi(M_star,zd,scatter=False)
       c200 = self.MCrelation(M200,MCerror=False)
       cscatter=self.MCrelation(M200,MCerror=True)
       self.galaxies.add_column('c200',c200)     
       rho=self.rho_crit_univ(zd)
       self.galaxies.add_column('rho_crit',rho)
       r200 = (3*M200/(800*3.14159*self.galaxies.rho_crit))**(1./3) #units: megaparsecs         #http://arxiv.org/pdf/astro-ph/9908213v1.pdf
       rs = r200/c200                                                                           #(wright and brainerd)
       rhos = self.delta_c(c200)*self.galaxies.rho_crit # units: solar mass per cubic megaparsec

       self.galaxies.add_column('rs',rs)
       self.galaxies.add_column('rhos',rhos)

       sigmacrit = self.SigmaCrit()  # units: Solarmasses per megaparsec^2
       self.galaxies.add_column('SigmaCrit',sigmacrit)
       kappas = rhos*rs/sigmacrit


       x = rphys/rs
       if hardcut!=False:
           for i in range(len(x)):
               if rphys[i] > hardcut*r200[i]: x[i]=-1 # Flag for hard cutoff.
           kappaNFW = kappas*(LP.Ffunc(x)) #following http://arxiv.org/pdf/astro-ph/9908213v1.pdf
           shearNFW = kappas * LP.Gfunc(x)
       if truncation=="BMO2":
           kappaNFW = kappas*(LP.BMO2Ffunc(x)) #BMO profile
           shearNFW = kappas * LP.BMO2Gfunc(x)    


       #---------------------------------------------------------------------------------

       # Compute starlight lensing component.

       #---------------------------------------------------------------------------------


       kappa = kappaNFW   #these could include a component from starlight too
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

       return None


# ----------------------------------------------------------------------------
   #Function binning MStar by MStar and redshift.
   def Bin_Mstars(self):

      linbin,delta=numpy.linspace(7,13,12,retstep=True)
      Mstarbins=numpy.power(10,linbin)
      self.Mstarbins=Mstarbins
      zbins=[0,1,2,99]
      self.zbins=zbins

      halomean = numpy.zeros((len(Mstarbins)-1,len(zbins)-1))
      halodist = numpy.zeros((len(Mstarbins)-1,len(zbins)-1))

      for i in range(len(Mstarbins)-1):
         for j in range(len(zbins)-1):
            M=self.catalog.where((self.catalog['M_Stellar[M_sol/h]']>Mstarbins[i]) & \
                                   (self.catalog['M_Stellar[M_sol/h]']<Mstarbins[i+1]) & \
                                   (self.catalog['z_spec']>zbins[j])&\
                                   (self.catalog['z_spec']<zbins[j+1]))
            #if i==6 and j == 1: print M['M_Stellar[M_sol/h]'][9], M['z_spec'][9]
            halomean[i,j]=numpy.mean(M['M_Subhalo[M_sol/h]'])
            if numpy.isnan(halomean[i,j])!=True:
               halodist[i,j]=numpy.std(numpy.log10(M['M_Subhalo[M_sol/h]']))
               #print halodist[i,j]

      #for i in range(len(Mstarbins)):
      #   for j in range(len(zbins-1)):
      #      if halodist[i,j]==0:
      #         halodist[i,j]=0.2 ###FILLER

      self.binnedhalomeans=halomean
      self.binnedhalodists=halodist
      self.Bin_MstarsRUN=True

   #Function taking Mstar (binned relationship from catalog) and pick a hlao mass. 
   def Mstar_to_M200(self,M_Star,redshift):
      #Using catalog's M*-Mh relation:
      M_200=numpy.zeros(len(M_Star))
      i=numpy.digitize(M_Star,self.Mstarbins)
      j=numpy.digitize(redshift,self.zbins)


      for k in range(len(M_Star)):
         halomean=self.binnedhalomeans[i[k]-1,j[k]-1]
         halodist=self.binnedhalodists[i[k]-1,j[k]-1]
         M_200[k]=10**(numpy.log10(halomean)+rnd.normal(0,halodist))
      return M_200
# ----------------------------------------------------------------------------

   #function to reconstruct the lightcone, starting with minimal information (currently M* and z_spec, but should ideally just be colours)
   def reconstruct_lightcone(self,photozerr=0.05,Behroozi=False):
      self.reconstruct = self.galaxies.where(self.galaxies['z_spec'] >0.0 )     
      self.reconstruct.keep_columns(['pos_0[rad]','pos_1[rad]','z_spec','M_Stellar[M_sol/h]', 'r'])

      z=self.reconstruct['z_spec']
      if photozerr>0:
         z_phot=rnd.normal(z,photozerr*(1+z))
      else:
         z_phot=z
      Da = numpy.zeros(len(z_phot))
      Da_tosource = numpy.zeros(len(z_phot))
      for i in range(len(z_phot)):
         Da[i] = D.Da(z_phot[i])
         Da_tosource[i] = D.Da(z_phot[i],self.zs)
      rphys=self.reconstruct.r*Da*arcmin2rad

      self.reconstruct.add_column('z_phot',z_phot)

      M_StarTrue=self.reconstruct['M_Stellar[M_sol/h]']
      M_Star=M_StarTrue*0.0
      for i in range(len(M_StarTrue)):
         M_Star[i]=10**(numpy.log10(M_StarTrue[i])+0*rnd.normal(0,0.2)) ### Fake M* scatter. (0.2 error in the log)

      #print self.galaxies['M_Stellar[M_sol/h]']
      #print M_Star.max()
      #print M_Star.min()
      
      self.reconstruct.add_column('M_Star',M_Star)

      if Behroozi==True:
         M_200=self.Mstar_to_M200_Behroozi(M_Star,z_phot,scatter=True) ###needs a scatter
      else:
         if self.Bin_MstarsRUN==False:
            self.Bin_Mstars()        
         M_200=self.Mstar_to_M200(M_Star,z_phot)

      self.reconstruct.add_column('M_200',M_200)

      c_200= self.MCrelation(M_200, MCerror=True) ###needs correct scatter for low mass halos.
      
      self.reconstruct.add_column('c_200',c_200)      

      rho=self.rho_crit_univ(z_phot)
      """
      l= M_200[:]/self.galaxies['M_Subhalo[M_sol/h]'][:]
      #plt.scatter((l),self.reconstruct.z_phot,s=1, edgecolor='none',c='k')
      plt.scatter((l),numpy.log10(self.galaxies['M_Subhalo[M_sol/h]']),s=1, edgecolor='none',c='k')
      plt.xscale('log')
      #plt.yscale('log')
      plt.xlim([0.1,10**3])
      #plt.ylim([0,3])
      #plt.ylim([0,3])
      plt.xlabel('Overestimate of halo mass from Behroozi compared to M_halo_MS')
      plt.ylabel('log(M_subhalo/M_sol)')   
      #plt.ylabel('redshift')
      plt.axvline(x=1,ymin=0,ymax=1)
      plt.savefig('M_200.png')
      plt.show()
      """

      r_200 = (3*M_200/(800*3.14159*rho))**(1./3) #units: megaparsecs         #http://arxiv.org/pdf/astro-ph/9908213v1.pdf
      rs = r_200/c_200                                                                           #(wright and brainerd)
      rhos = self.delta_c(c_200)*rho # units: solar mass per cubic megaparsec

      x = rphys/rs
      sigmacrit=(1.663*10**18)*(self.Da_s/(Da_tosource*Da)) 

      kappas = rhos*rs/sigmacrit
     

      #plt.hist(numpy.log10(M_200))
      #plt.show()


       # need to cut out sub halos - can be done in several ways.
       #for i in range(kappas.size()):
       #   if self.galaxies

      kappaNFW = 2.*kappas*(LP.Ffunc(x)) #following http://arxiv.org/pdf/astro-ph/9908213v1.pdf

      shearNFW = kappas * LP.Gfunc(x)
       
       #---------------------------------------------------------------------------------

       # Compute starlight lensing component.

       #---------------------------------------------------------------------------------


      kappa = kappaNFW   #these could include a component from starlight too
      shear = shearNFW

      # kappa or shear is zero if behind the source!
      for i in range(len(z_phot)):
          if z_phot[i] > self.zs:
             kappa[i] = 0.0
             shear[i] = 0.0

      kappa_keeton=self.KappaKeeton(self.zl,z_phot,self.zs,kappa,shear)

       # Contributions to simple weighted sum (Keeton 2003):
      self.reconstruct.add_column('kappa',kappa)
      self.reconstruct.add_column('gamma',shear)
      self.reconstruct.add_column('kappa_keeton',kappa_keeton)

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
         #plt.scatter(self.galaxies.x, self.galaxies.y, c='k', marker='o',s=((numpy.log(self.galaxies['M_Subhalo[M_sol/h]']))/3),edgecolor = 'none' )
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


       # Phil Marshall's subplot  View along redshift axis: 

       plt.subplot(2,1,2)
       empty = True
       

       if dmglow==True:
         plt.scatter(self.galaxies.z_spec, self.galaxies.y, c='k', marker='o',s=(((self.galaxies['M_Subhalo[M_sol/h]']))*3e-12),edgecolor = 'none' )
         empty = False
       if starlight==True:
         plt.scatter(self.galaxies.z_spec, self.galaxies.y, c='y', marker='o',s=((numpy.log(self.galaxies['M_Stellar[M_sol/h]']))/2),edgecolor = 'none' )     
         empty = False
       if kappa_keeton:
         for galaxy in self.galaxies:
            kappa = galaxy['kappa_keeton']
            if kappa <0:
               plt.scatter(galaxy['z_spec'], galaxy['y'], c='b', marker='o',s=-kappa/scale , edgecolor = 'none')    
            else:
               plt.scatter(galaxy['z_spec'], galaxy['y'],  c='r', marker='o',s=kappa/scale , edgecolor = 'none')
         empty = False
       if kappa_indiv: 
         plt.scatter(self.galaxies.z_spec, self.galaxies.y, c='none', marker='o',s=(self.galaxies.kappa)/scale2)     
         empty = False
       if observed_light==True:
         plt.scatter(self.galaxies.z_spec, self.galaxies.y, c='y', marker='o',s=2**(25-(self.galaxies['mag_SDSS_r'])),edgecolor = 'none' )     
         empty = False


       if empty:
         plt.scatter(self.galaxies.z_spec, self.galaxies.y, c='k', marker='o',s=1)


       plt.xlabel('Redshift')
       plt.ylabel('y / arcmin')
       zmax = self.galaxies['z_spec'].max()


       plt.axis([0,self.zs+0.1,0,self.rmax+0.1])

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


# # ----------------------------------------------------------------------------

   def curve_of_growth(self,ordering="contribution",starlight=False,dmglow=False,kappa_indiv=False,kappa_keeton=False,observed_light=True):
       plt.clf()
       scale=numpy.max([ numpy.absolute((numpy.min(self.galaxies.kappa_keeton))), \
                            numpy.max(self.galaxies.kappa_keeton)])/200
       scale2= (numpy.max(self.galaxies.kappa))/200

       zeropoint=-self.kappa_expected()
       zero=numpy.ones(1)*zeropoint


       plt.subplot(2,1,1)
       #set up things to plot:
       if ordering=="contribution":
          args=numpy.argsort(-numpy.absolute(self.galaxies.kappa_keeton))
       if ordering=="distance":
          args=numpy.argsort((self.galaxies.r))
          dist=numpy.concatenate((zero,numpy.take(self.galaxies.r,args)))
       if ordering=="r_mag":
          args=numpy.argsort((self.galaxies['mag_SDSS_r']))
          rmag=numpy.concatenate((zero,numpy.take(self.galaxies['mag_SDSS_r'],args)))
       ordered_kappas=numpy.take(self.galaxies.kappa_keeton,args)


       cat=numpy.concatenate((zero,(ordered_kappas)))
       cumtot=numpy.cumsum(cat)
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
       #n1=n[:20]
       #cumtot1=cumtot[:20]
       #corresponding_r1=corresponding_r[:20]

       if ordering=="contribution":
          plt.plot(n,cumtot)
          #plt.xlim([0,200])
          plt.axhline(y=cumtot[-1], xmin=0, xmax=1000,color='black', ls='dotted')
          plt.axhline(y=0, xmin=0, xmax=1000,color='black', ls='solid')
       if ordering=="distance":
          plt.plot(dist,cumtot)
          #plt.xlim([0,200])
          plt.axhline(y=cumtot[-1], xmin=0, xmax=1000,color='black', ls='dotted')
          plt.axhline(y=0, xmin=0, xmax=1000,color='black', ls='solid')
          plt.xlabel('LOS distance (arcmin)')
       if ordering=="r_mag":
          plt.plot(rmag,cumtot)
          #plt.xlim([0,200])
          plt.axhline(y=cumtot[-1], xmin=0, xmax=1000,color='black', ls='dotted')
          plt.axhline(y=0, xmin=0, xmax=1000,color='black', ls='solid')
          plt.xlim([16,26])
       plt.ylabel('$\kappa_{ext}$ (cumulative)')


#      -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.


#      #Redshift axis vs LOS distance


       plt.subplot(2,1,2)
       empty = True

       self.galaxies=self.galaxies.where(self.galaxies.r<8)
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

       if dmglow==True:
         plt.scatter(self.galaxies.z_spec, self.galaxies.r, c='k', marker='o',s=((numpy.log(self.galaxies['M_Subhalo[M_sol/h]']))/3),edgecolor = 'none' )
         empty = False
       if starlight==True:
         plt.scatter(self.galaxies.z_spec, self.galaxies.r, c='g', marker='o',s=2**(numpy.log10(self.galaxies['M_Stellar[M_sol/h]'])-7),edgecolor = 'none' )     
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

# ============================================================================

# TESTS:

def test1(catalog):

    plt.clf()
    rmax = 20
    print rmax

    zs=2.0 #should be selecting a source redshift (possibly use forecaster from Collett et al. 2012)

#-0.00633207	-0.0173126
    xpos = -0.008
    ypos = -0.008
    zl =  1.0 
    xc = []
    #xc = [xpos,ypos,zl] #leave as [] to select a lens at random

    print "Initialising lightcone data..."
    lc = lightcone(catalog,rmax,zs,lensindex=-1,position=xc)


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

    print "Total external convergence =",numpy.sum(lc.galaxies.kappa_keeton)-lc.kappa_expected()
    # Now make an illustrative plot:
    
    if rmax < 8:
       print "Plotting objects in lightcone..."
       lc.plot(starlight=False,dmglow=False, kappa_indiv=False, kappa_keeton=True,observed_light=True)
       

       pngfile = 'test.png'
       plt.savefig(pngfile)
       print "Plot saved in",pngfile

    print "Plotting curve of growth..."
    plt.clf()
    lc.curve_of_growth(ordering="distance")
    pngfile = 'curve_of_growth.png'
    plt.savefig(pngfile)
    print "Plot saved in",pngfile
    

    return

#-------------------------------------------------------------------------

def test2(catalog): #plots a kappa distribution:
###not finished###
#we could consider feeding in our lens selection parameters here)


    plt.clf()
    rmax = 5

    zs=1.4 #should be selecting a source redshift (possibly use forecaster from Collett et al. 2012)
    position = [] # radians, leave as [] to select a lens at random
    zl=0.6 #
    xc = []

    xmax = catalog['pos_0[rad]'].max()
    xmin = catalog['pos_0[rad]'].min()
    ymax = catalog['pos_1[rad]'].max()
    ymin = catalog['pos_1[rad]'].min()

    iterations=1000
    K=numpy.zeros(iterations)
    h=3
    for j in range(iterations):
       x = rnd.uniform(xmin+rmax*arcmin2rad,xmax-rmax*arcmin2rad)
       y = rnd.uniform(ymin+rmax*arcmin2rad,ymax-rmax*arcmin2rad)
       xc=[x,y,zl]
       lc = lightcone(catalog,rmax,zs,position=xc)
       if j==0: kappa_empty=lc.kappa_expected()
       lc.make_kappa_contributions(hardcut=h)
       K[j]=numpy.sum(lc.galaxies.kappa_keeton)-kappa_empty
       if j % 100 ==0: print K[j],j

    bins=numpy.arange(-0.05,0.1,0.005)
    plt.hist(K,bins,normed=True)
    plt.xlabel("$\kappa_{ext}$")
    plt.ylabel("pdf($\kappa_{ext}$)")
    #plt.yscale('log')
    pngfile = 'Kappa_keeton_distribution_%.0f_Rvir_cut.png'%h
    #plt.ylim([0,500])
    #plt.xlim([-0.2,.3])
 

    plt.savefig(pngfile)


    print "Plot saved in",pngfile

    return

#-------------------------------------------------------------------------

def test3(catalog):
    plt.clf()
    rmax = 2
    xc = []
    zs=1.4
    print "Initialising lightcone data..."
    lc = lightcone(catalog,rmax,zs,lensindex=-1,position=xc)

    print "Computing Keeton (2003) convergence at optical axis due to each halo..."
    lc.make_kappa_contributions()

    print "Total external convergence =",numpy.sum(lc.galaxies.kappa_keeton)
    # Now make an illustrative plot:
    
    lc.reconstruct_lightcone()

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
