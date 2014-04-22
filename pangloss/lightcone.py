# ===========================================================================

import pangloss

import cPickle
import numpy
import pylab as plt
from math import pi


# ======================================================================

class Lightcone(object):
    """
    NAME
        Lightcone

    PURPOSE
        Define a conical region of space containing many galaxies, and 
        enable mass to be assigned to those galaxies and consequent 
        useful quantities to be calculated.

    COMMENTS
        The masses of galaxies in a lightcone are always uncertain. 
        Methods are provided to characterise this uncertainty by drawing
        sample masses for each galaxy, from various specified relations.

    INITIALISATION
        catalog       Filename of the parent galaxy catalog
        flavor        Is the catalog 'real' or 'simulated'?
        position      The sky position (J2000 deg) of the cone centre
        radius        The radius of the lightcone field of view (arcmin)
        maglimit      The depth of the galaxy selection (magnitudes)
        band          The band in which the selection is made
    
    METHODS
        galaxiesWithin(self,radius,cut=[18.5,24.5],band="F814W",radius_unit="arcsec"):
        
        numberWithin(self,radius,cut=[18.5,24.5],band="F814W",radius_unit="arcsec"):
        
        define_system(self,zl,zs,cosmo=[0.25,0.75,0.73]):
        
        loadGrid(self, Grid):
        
        mimicPhotozError(self,sigma=0.1):
        
        writeColumn(self,string,values):
        
        snapToGrid(self, Grid):
        
        drawMstars(self,model): Needs updating to take SHMR object
        
        drawMhalos(self,modelT):
        
        drawConcentrations(self,errors=False):
        
        makeKappas(self,errors=False,truncationscale=5,profile="BMO1"):
        
        combineKappas(self):

    BUGS

    AUTHORS
      This file is part of the Pangloss project, distributed under the
      GPL v2, by Tom Collett (IoA) and  Phil Marshall (Oxford). 
      Please cite: Collett et al 2013, http://arxiv.org/abs/1303.6564

    HISTORY
      2013-03-23  Collett & Marshall (Cambridge)
    """

# ----------------------------------------------------------------------------

    def __init__(self,catalog,flavor,position,radius,maglimit=99,band="r"):
        
        self.name = 'Lightcone through the Universe'
        self.flavor = flavor   # 'real' or 'simulated'
        self.catalog = catalog
        
        # Simulated lightcones have "true" (ray-traced) convergence:
        self.kappa_hilbert = None # until set!
        
        # Catalog limits:
        self.xmax = self.catalog['nRA'].max()
        self.xmin = self.catalog['nRA'].min()
        self.ymax = self.catalog['Dec'].max()
        self.ymin = self.catalog['Dec'].min() 

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

        # Cut out a square cone:
        self.rmax = radius
        self.xc = [position[0],position[1]]

        dx = self.rmax*pangloss.arcmin2rad
        self.galaxies = self.catalog.where((self.catalog.nRA > (self.xc[0]-dx)) & \
                                           (self.catalog.nRA < (self.xc[0]+dx)) & \
                                           (self.catalog.Dec > (self.xc[1]-dx)) & \
                                           (self.catalog.Dec < (self.xc[1]+dx))   )

        # Trim it to a circle:
        x = (self.galaxies.nRA - self.xc[0])*pangloss.rad2arcmin
        y = (self.galaxies.Dec - self.xc[1])*pangloss.rad2arcmin
        r = numpy.sqrt(x*x + y*y)
        phi=numpy.arctan(y/x)
        self.galaxies.add_column('x',x)
        self.galaxies.add_column('y',y)
        self.galaxies.add_column('r',r)
        self.galaxies.add_column('phi',phi)
        self.galaxies = self.galaxies.where(self.galaxies.r < self.rmax)

        try: 
            self.galaxies = self.galaxies.where(self.galaxies.Type != 2) 
        except AttributeError: pass

        # Now we have a small catalog, just for the requested lightcone.
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

        # Halo parameters, to be varied during sampling analysis:
        self.galaxies.add_column('z',self.galaxies.z_obs*1.0)

        if self.flavor == 'simulated':
            # Take the log of the halo mass, and set up the parameter array:
            self.galaxies.add_column('Mh_obs',numpy.log10(self.galaxies.Mhalo_obs))
            self.galaxies.add_column('Mh',self.galaxies.Mh_obs*1.0)
            # Stellar masses will be added by drawMstars
            # Halo masses will be replaced by drawMhalos

        elif self.flavor == 'real':
            # Mstar is already given as log M...
            self.galaxies.add_column('Mstar',self.galaxies.Mstar_obs*1.0)
            # Halo masses will be added by drawMhalos
          
        if len(self.galaxies) == 0: 
            print "Lightcone: WARNING: no galaxies here!"

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

        # Save memory! 
        del self.catalog
        del catalog
        
        return None

# ----------------------------------------------------------------------------

    def __str__(self):
        return 'Lightcone of radius %.2f arcmin, centred on (%.3f,%.3f) rad' % (self.rmax,self.xc[0],self.xc[1])

# ----------------------------------------------------------------------------
# Tell me the number of galaxies within a certain radius, that pass a 
# certain magnitude cut.

    def galaxiesWithin(self,radius,cut=[18.5,24.5], band="F814W", radius_unit="arcsec"):

        if band == "u" or band ==  "g" or band == "r" or band ==  "i" or band == "z":
            col = "mag_SDSS_%s" % band
        elif band == "F814" or band == "F814W" or band == "814" or band == 814:
            col = "mag_F814W"
        else:
            col = "mag_%s" % band
        if radius < 10: 
            print "Warning: Default units for radius are arcsec!"
        if radius_unit == "arcsec":
            radius = radius/60.
        if col != "warning":
            self.N_cut=self.galaxies.where((self.galaxies.r < radius)  & \
                                              (self.galaxies["%s"%col] < cut[1])& \
                                          (self.galaxies["%s"%col] > cut[0]))

            return self.galaxies.where((self.galaxies.r < radius) & \
                                          (self.galaxies["%s"%col] < cut[1])& \
                                          (self.galaxies["%s"%col] > cut[0]))

    def numberWithin(self,radius,cut=[18.5,24.5],band="F814W",units="arcsec"):
        Ntable = self.galaxiesWithin(radius,cut,band,units)
        return len(Ntable.r)

# ----------------------------------------------------------------------------

    def defineSystem(self,zl,zs,cosmo=[0.25,0.75,0.73]):
        self.zl = zl
        self.zs = zs
        self.cosmo = cosmo
        self.galaxies = self.galaxies.where(self.galaxies.z_obs<zs+0.2)
        return

# ----------------------------------------------------------------------------

    def loadGrid(self, Grid):
        if numpy.abs(self.zl-Grid.zltrue) > 0.05: print "Grid zl != lens zl" 
        if numpy.abs(self.zs-Grid.zs)     > 0.05: print "Grid zs != lens zs" 
        self.redshifts,self.dz = Grid.redshifts,Grid.dz
        self.Da_l,self.Da_s,self.Da_ls = Grid.Da_l,Grid.Da_s,Grid.Da_ls
        return

# ----------------------------------------------------------------------------
# Simple model for spectroscopic coverage - to be applied to the 
# calibration lightcones to match the observed catalog.

    def configureForSurvey(self, experiment):
    
        PR = experiment.parameters['PhotometricRadius']
        PD = experiment.parameters['PhotometricDepth']
        assert len(PR) == len(PD)

        SR = experiment.parameters['SpectroscopicRadius']
        SD = experiment.parameters['SpectroscopicDepth']
        assert len(SR)==len(SD)

        band = experiment.parameters['LightconeDepthBand']

        if band == "u" or band ==  "g" or band == "r" or band ==  "i" or band == "z":
            col = "mag_SDSS_%s" % band
        elif band == "F814" or band == "F814W" or band == "814" or band == 814:
            col = "mag_F814W" #note that this isn't included atm
        else:
            col = "mag_%s" % band
            
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
        # Only include galaxies observed by photometry:
        
        self.writeColumn('photo_flag',False)
        self.writeColumn('identifier',range(len(self.galaxies.x)))
        if PR != ['']:
            for i in range(len(PR)):
                R=PR[i]*60 # positions are stored in arcseconds
                D=PD[i]

                goodset = set(self.galaxies.where((self.galaxies.r < R) & \
                   (self.galaxies["%s"%col] < D)).identifier)

                self.galaxies.photo_flag[numpy.array(\
                   [_ in goodset for _ in self.galaxies.identifier])]=True
        
        self.galaxies = self.galaxies.where(self.galaxies.photo_flag==True)

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
        # Set spectroscopic flag of any galaxy that should have 
        # spectroscopy, according to the given radius and depth:
        
        self.writeColumn('spec_flag',False)

        if SR!=['']:
            for i in range(len(SR)):
                R=SR[i]*60 # positions are stored in arcseconds
                D=SD[i]
                goodset=set(self.galaxies.where((self.galaxies.r < R) & \
                   (self.galaxies["%s"%col] < D)).identifier)

                self.galaxies.spec_flag[numpy.array(\
                   [_ in goodset for _ in self.galaxies.identifier])]=True

        return

# ----------------------------------------------------------------------------
# The following methods are designed to be run multiple times
# (previous ones are single use per lightcone)
# ----------------------------------------------------------------------------

# Add galaxy property column, overwriting any values that already exist:

    def writeColumn(self,string,values):
        try:
            self.galaxies.add_column('%s'%string,values)
        except ValueError:
            self.galaxies["%s"%string]=values

# ----------------------------------------------------------------------------
 
    def mimicPhotozError(self,sigma=0.1):

        # Start with the original, constant z_obs in the catalog:
        z_obs = self.galaxies.z_obs.copy()
        
        # RMS error on this z is either sigma or 0.0:
        e = ( self.galaxies.spec_flag == False )*sigma
        
        # Add Gaussian noise to get this realisation's z:
        z = z_obs + e*(1+z_obs)*numpy.random.randn(len(z_obs))
        
        # Over-write the z parameter array:
        self.writeColumn('z',z)

        return

# ----------------------------------------------------------------------------
# Snap the parameters z onto the grid, to speed up calculations:

    def snapToGrid(self, Grid):
        z = self.galaxies.z
        sz,p = Grid.snap(z)
        self.writeColumn('Da_p',Grid.Da_p[p])
        self.writeColumn('rho_crit',Grid.rho_crit[p])
        self.writeColumn('sigma_crit',Grid.sigma_crit[p])
        self.writeColumn('beta',Grid.beta[p])
        rphys = self.galaxies.r*pangloss.arcmin2rad*self.galaxies.Da_p
        self.writeColumn('rphys',rphys)

# ----------------------------------------------------------------------------
# Given Mhalo and z, draw an Mstar, and an identical Mstar_obs:

    def drawMstars(self,model):
        Mstar = model.drawMstars(self.galaxies.Mh,self.galaxies.z)
        self.writeColumn('Mstar',Mstar)
        # Copy this to Mstar_obs - mimicMstarError will deal with noise
        self.writeColumn('Mstar_obs',Mstar)
        return
   
# ----------------------------------------------------------------------------
# Given an observed Mstar_obs, what could the parameter Mstar be?

    def mimicMstarError(self,sigmaP,sigmaS):
        
        Mstar = self.galaxies.Mstar_obs.copy()

        # Add offset due to photometric or spectroscopic error:
        Mstar[self.galaxies.spec_flag==False] += numpy.random.randn(Mstar[self.galaxies.spec_flag==False].size)*sigmaP
        Mstar[self.galaxies.spec_flag==True] += numpy.random.randn(Mstar[self.galaxies.spec_flag==True].size)*sigmaS

        # Over-write the Mstar parameter array:
        self.writeColumn('Mstar',Mstar)
        return
           

# ----------------------------------------------------------------------------
# Given an Mstar and z, what could the parameter Mh be?

    def drawMhalos(self,model):
        Mh = model.drawMhalos(self.galaxies.Mstar,self.galaxies.z)
        self.writeColumn('Mh',Mh)
        return

# ----------------------------------------------------------------------------
# Given an Mh, what could the halo concentration be?

    def drawConcentrations(self,errors=False):
        M200 = 10**self.galaxies.Mh        
        r200 = (3*M200/(800*3.14159*self.galaxies.rho_crit))**(1./3)
        self.writeColumn("r200",r200)
        
        c200 = pangloss.MCrelation(M200,scatter=errors)
        
        self.writeColumn("c200",c200)
        r_s = r200/c200        
        self.writeColumn('rs',r_s)
        x = self.galaxies.rphys/r_s
        self.writeColumn('X',x)
        return

# ----------------------------------------------------------------------------
# Compute halos' contributions to the convergence:

    def makeKappas(self,errors=False,truncationscale=5,profile="BMO1"):
        c200 = self.galaxies.c200
        r200 = self.galaxies.r200
        x = self.galaxies.X
        r_s = self.galaxies.rs
        rho_s = pangloss.delta_c(c200)*self.galaxies.rho_crit
        kappa_s = rho_s * r_s /self.galaxies.sigma_crit  #kappa slice for each lightcone
        
        r_trunc = truncationscale*r200
        xtrunc = r_trunc/r_s
        kappaHalo = kappa_s*1.0
        gammaHalo = kappa_s*1.0
        
        if profile=="BMO1":
            F=pangloss.BMO1Ffunc(x,xtrunc)
            G=pangloss.BMO1Gfunc(x,xtrunc)
        
        if profile=="BMO2":
            F=pangloss.BMO2Ffunc(x,xtrunc)
            G=pangloss.BMO2Gfunc(x,xtrunc)
        
        kappaHalo *= F
        gammaHalo *= (G-F)

        phi = self.galaxies.phi
        
        kappa = kappaHalo 
        gamma = gammaHalo
        gamma1 = gamma*numpy.cos(2*phi)
        gamma2 = gamma*numpy.sin(2*phi)
        
        mu = 1.0/(((1.0 - kappa)**2.0) - (gamma**2.0))

        self.writeColumn('kappa',kappa)
        self.writeColumn('gamma',gamma)
        self.writeColumn('gamma1',-gamma1)
        self.writeColumn('gamma2',-gamma2)
        self.writeColumn('mu',mu)
        
        return
        
# ----------------------------------------------------------------------------

    def combineKappas(self):
    
        B=self.galaxies.beta
        K=self.galaxies.kappa
        G=self.galaxies.gamma
        G1=self.galaxies.gamma1
        G2=self.galaxies.gamma2
        D= K**2-G**2

        kappa_keeton  = (1.-B) * (K- B*(D)) /  ( (1-B*K)**2   - (B*G)**2   )    
        gamma1_keeton = (1.-B) * (G1) /  ( (1-B*K)**2   - (B*G)**2   )  
        gamma2_keeton = (1.-B) * (G2) /  ( (1-B*K)**2   - (B*G)**2   )  

        kappa_tom  = (1.-B) * K
        gamma1_tom = (1.-B) * G1
        gamma2_tom = (1.-B) * G2

        self.writeColumn('kappa_keeton',kappa_keeton)
        self.writeColumn('gamma1_keeton',gamma1_keeton)
        self.writeColumn('gamma2_keeton',gamma2_keeton)

        self.writeColumn('kappa_tom',kappa_tom)
        self.writeColumn('gamma1_tom',gamma1_tom)
        self.writeColumn('gamma2_tom',gamma2_tom)

        self.writeColumn('kappa_add',K)
        self.writeColumn('gamma1_add',G1)
        self.writeColumn('gamma2_add',G2)


        self.kappa_add_total=numpy.sum(self.galaxies.kappa)
        self.kappa_keeton_total=numpy.sum(self.galaxies.kappa_keeton)
        self.kappa_tom_total=numpy.sum(self.galaxies.kappa_tom)

        self.gamma1_add_total=numpy.sum(self.galaxies.gamma1)
        self.gamma1_keeton_total=numpy.sum(self.galaxies.gamma1_keeton)
        self.gamma1_tom_total=numpy.sum(self.galaxies.gamma1_tom)

        self.gamma2_add_total=numpy.sum(self.galaxies.gamma2)
        self.gamma2_keeton_total=numpy.sum(self.galaxies.gamma2_keeton)
        self.gamma2_tom_total=numpy.sum(self.galaxies.gamma2_tom)
        

        #print self.galaxies.Mstar.max()
        #print self.galaxies.Mh.max()
        #print self.galaxies.kappa.max()
        #print self.kappa_add_total

        return self.kappa_add_total
        
# ----------------------------------------------------------------------------

    def combineMus(self):
    
        M=self.galaxies.mu
        K=self.galaxies.kappa
        G=self.galaxies.gamma
        G1=self.galaxies.gamma1
        G2=self.galaxies.gamma2
        
        self.writeColumn('mu_add',M)
        
        self.kappa_add_total=numpy.sum(K)
        self.gamma1_add_total=numpy.sum(G1)
        self.gamma2_add_total=numpy.sum(G2)
        
        Ksum = self.kappa_add_total
        G1sum = self.gamma1_add_total
        G2sum = self.gamma2_add_total
        Gsum = numpy.sqrt(G1sum**2 + G2sum**2)
        
        Msum = 1.0/(((1.0 - Ksum)**2.0) - (Gsum**2.0))
     
        self.mu_add_total=Msum

        return self.mu_add_total        

# ----------------------------------------------------------------------------

    def allMass(self, tau=5):
        
        Mvir = self.galaxies.Mh
        self.M0 = 4.*Mvir/(pi - 2.)
        self.Mhalo = self.M0 * (tau**2/(tau**2 + 1.)**2) * ((tau**2 - 1.)*numpy.log(tau) + tau*pi - tau**2 - 1)
        self.writeColumn('truncNFW_Mhalo',self.Mhalo)

        self.totalMhalo = numpy.sum(self.Mhalo)
        
        return self.totalMhalo
# ----------------------------------------------------------------------------
# Plotting
# ----------------------------------------------------------------------------

    def plotFieldOfView(self,quantity,AX):
       
       slicehalfwidth = self.rmax / 6.0
       
       if quantity == 'mass':
           self.writeColumn('rtrunc_arcmin',(self.galaxies.r200/ self.galaxies.Da_p) * pangloss.rad2arcmin)
           self.writeColumn('rscale_arcmin',(self.galaxies.rs/ self.galaxies.Da_p) * pangloss.rad2arcmin)
           for i in range(len(self.galaxies.x)):
               trunc = plt.Circle([self.galaxies.x[i], self.galaxies.y[i]],radius=self.galaxies.rtrunc_arcmin[i],fill=True,fc="b",alpha=0.05)
               AX.add_patch(trunc)
           for i in range(len(self.galaxies.x)):
               core = plt.Circle([self.galaxies.x[i], self.galaxies.y[i]],radius=self.galaxies.rscale_arcmin[i],fill=True,fc="r",alpha=0.5)
               AX.add_patch(core)
           plt.title('Halo Mass')

       elif quantity == 'kappa':
           plt.scatter(self.galaxies.x, self.galaxies.y, c='r', marker='o', s=(self.galaxies.kappa)*30000)    
           plt.title('Convergence')
       
       elif quantity == 'mu':
           plt.scatter(self.galaxies.x, self.galaxies.y, c='g', marker='o', s=((self.galaxies.mu-1.0)*3E4))    
           plt.title('Magnification')
       
       elif quantity == 'stellarmass':
           plt.scatter(self.galaxies.x, self.galaxies.y, c='y', marker='o', s=(numpy.log(self.galaxies.Mstar)/2),edgecolor = 'none' )     
           plt.title('Stellar Mass')
       
       elif quantity == 'light':
           plt.scatter(self.galaxies.x, self.galaxies.y, c='y', marker='o', s=(2**(25-self.galaxies.mag)),edgecolor = 'none' )     
           plt.title('Galaxy Light')

       else:
           raise "Lightcone plotting error: unknown quantity "+quantity
           
       # Lightcone boundary and centroid:
       circ = plt.Circle(self.xc,radius=self.rmax,fill=False,linestyle='dotted')
       AX.add_patch(circ)
       AX.plot([self.xc[0]],[self.xc[1]], c='k', marker='+',markersize=20,markeredgewidth=1)
       axlimits = [self.xc[0]-self.rmax-0.1,self.xc[0]+self.rmax+0.1,self.xc[1]-self.rmax-0.1,self.xc[1]+self.rmax+0.1]
       AX.axis(axlimits)

       # Show slice:
       plt.axvline(x=slicehalfwidth, ymin=axlimits[2], ymax=axlimits[3],color='black', ls='dotted')
       plt.axvline(x=-slicehalfwidth, ymin=axlimits[2], ymax=axlimits[3],color='black', ls='dotted')
       
       # Labels:
       plt.xlabel('x / arcmin')
       # plt.ylabel('y / arcmin')
       
       return

# ----------------------------------------------------------------------------

    def plotLineOfSight(self,quantity,AX):
       
       
       # Only plot a subset of points, in a slice down the middle of 
       # the light cone:
       subset = numpy.abs(self.galaxies.x)<0.3
       #print self.galaxies.mu[subset]
 
       # Point positions:
       z = self.galaxies.z[subset]
       y = self.galaxies.y[subset]
       
       # Plot the points:
       if quantity == 'mass':
           size = (10.0**(self.galaxies.Mh[subset]-11.0))
           plt.scatter(z, y, c='k', marker='o', s=size, edgecolor='none' )
           plt.title('Line-of-sight Halo Mass')
      
       elif quantity == 'kappa':
           size = ((self.galaxies.kappa[subset])*30000.0)
           plt.scatter(z, y, c='r', marker='o', s=size, edgecolor='k' )
           plt.title('Line-of-sight Convergence')
       
       elif quantity == 'mu':
           size = ((self.galaxies.mu[subset]-1.0)*3E4)
           plt.scatter(z, y, c='g', marker='o', s=size, edgecolor='k' )
           plt.title(r'Line-of-sight Magnification $(\mu - 1)$')  

       elif quantity == 'stellarmass':
           size = ((numpy.log(self.galaxies.Mstar[subset]))/2.0)
           plt.scatter(z, y, c='y', marker='o', s=size, edgecolor='none' )
           plt.title('Line-of-sight Stellar Mass')

       elif quantity == 'light':
           size = (2**(25-(self.galaxies.mag[subset])))     
           plt.scatter(z, y, c='y', marker='o', s=size, edgecolor='none' )
           plt.title('Line-of-sight Galaxy Light')

       else:
           raise "Lightcone plotting error: unknown quantity "+quantity

       # Axis limits:
       zmax = max(self.galaxies.z.max(),self.zs)
       AX.axis([0,zmax+0.1,-self.rmax-0.1,self.rmax+0.1])

       # Labels:
       plt.xlabel('redshift z')
       plt.ylabel('y / arcmin')
      
       # Add lines marking source and lens plane, and optical axis:
       plt.axvline(x=self.zl, ymin=0, ymax=1,color='black', ls='dotted',label='bla')
       plt.axvline(x=self.zs, ymin=0, ymax=1,color='black', ls='dotted')
       plt.axhline(y=0.0, xmin=0.0, xmax=zmax, color='black', ls='dashed')

       return

# ----------------------------------------------------------------------------

    def plot(self,var,output=None):

       plt.clf()

       # Panel 1: Galaxy positions:
       ax1 = plt.subplot(3,3,(1,4), aspect ='equal')
       self.plotFieldOfView('light',ax1)
       
       # Panel 2: Halo mass distributions:
       ax2 = plt.subplot(3,3,(2,5), aspect ='equal')
       self.plotFieldOfView('mass',ax2)

       # Panel 3: Kappa contributions:
       ax3 = plt.subplot(3,3,(3,6), aspect ='equal')
       self.plotFieldOfView(var,ax3)
       
       # Lower panel: View along redshift axis
       ax4 = plt.subplot(3,3,(7,9))
       self.plotLineOfSight(var,ax4)
       
       if output != None:
           pangloss.rm(output)
           plt.savefig(output,dpi=300)
       
       return None
       
# ----------------------------------------------------------------------------
    

#=============================================================================
