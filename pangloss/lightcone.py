# ===========================================================================

import pangloss

#import pylab
#import matplotlib.pyplot as plt
#from mpl_toolkits.axes_grid1 import ImageGrid
import cPickle
import numpy
# import Relations as Rel
# import LensingProfiles as LP
# import LensingFunc as LF
# import pylab as plt

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
        
        tryColumn(self,string,values):
        
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
      Please cite: Collett et al 2013, arxiv/###

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
        
        self.xmax = self.catalog['pos_0[rad]'].max()
        self.xmin = self.catalog['pos_0[rad]'].min()
        self.ymax = self.catalog['pos_1[rad]'].max()
        self.ymin = self.catalog['pos_1[rad]'].min() 
        self.rmax = radius
        self.xc = [position[0],position[1]]

        dx = self.rmax*pangloss.arcmin2rad
        self.galaxies = self.catalog.where((self.catalog['pos_0[rad]'] > (self.xc[0]-dx)) & \
                                           (self.catalog['pos_0[rad]'] < (self.xc[0]+dx)) & \
                                           (self.catalog['pos_1[rad]'] > (self.xc[1]-dx)) & \
                                           (self.catalog['pos_1[rad]'] < (self.xc[1]+dx))   )

 
        x = (self.galaxies['pos_0[rad]'] - self.xc[0])*pangloss.rad2arcmin
        y = (self.galaxies['pos_1[rad]'] - self.xc[1])*pangloss.rad2arcmin
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

        #log the mass:
        try: self.galaxies.add_column('Mh',numpy.log10(self.galaxies['M_Subhalo[M_sol/h]']))
        except ValueError: pass


        self.galaxies.add_column('z_obs',self.galaxies.z_spec)
        self.galaxies.add_column('spec_flag',False)
        
        if len(self.galaxies) == 0: 
            print "Lightcone: WARNING: no galaxies here!"

        del self.catalog
        del catalog
        
        """
        F814 = (self.galaxies.mag_SDSS_i + self.galaxies.mag_SDSS_z)/2. #approximate F814 colour.
        self.galaxies.add_column("mag_F814W",F814)
        if band == "u" or band ==  "g" or band == "r" or band ==  "i" or band == "z":
            col = "mag_SDSS_%s" % band
        elif band == "F814" or band == "F814W" or band == "814" or band == 814:
            col = "mag_F814W"
        else:
            col = "mag_%s" % band
        self.galaxies=self.galaxies.where(self.galaxies["%s"%col] < maglimit)
        """
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
        self.zl=zl
        self.zs=zs
        self.cosmo=cosmo
        self.galaxies=self.galaxies.where(self.galaxies.z_spec<zs)

# ----------------------------------------------------------------------------

    def loadGrid(self, Grid):
        if numpy.abs(self.zl-Grid.zltrue)>0.05: print "Grid zl != lens zl" 
        if numpy.abs(self.zs-Grid.zs)    >0.05: print "Grid zs != lens zs" 
        self.redshifts,self.dz=Grid.redshifts,Grid.dz
        self.Da_l,self.Da_s,self.Da_ls=Grid.Da_l,Grid.Da_s,Grid.Da_ls

# ----------------------------------------------------------------------------

    def configureForSurvey(self, experiment):
        PR=experiment.parameters['PhotometricRadius']
        PD=experiment.parameters['PhotometricDepth']
        assert len(PR) == len(PD)

        SR=experiment.parameters['SpectroscopicRadius']
        SD=experiment.parameters['SpectroscopicDepth']
        assert len(SR)==len(SD)

        #SC=experiment.parameters['SpectroscopicCompleteness']
        #assert len(SR) == len (SC)
        #(this isn't currently included in the analysis. 
        # would need a mimik spectroscopic incompletness function)

        band=experiment.parameters['LightconeDepthBand']
        if band == "u" or band ==  "g" or band == "r" or band ==  "i" or band == "z":
            col = "mag_SDSS_%s" % band
        elif band == "F814" or band == "F814W" or band == "814" or band == 814:
            col = "mag_F814W" #note that this isn't included atm
        else:
            col = "mag_%s" % band
        # -------------------------------------------------

        #exclude galaxies not observed by photometry:
        self.tryColumn('photo_flag',False)
        self.tryColumn('identifier',range(len(self.galaxies.x)))

        if PR!=['']:
          for i in range(len(PR)):
            R=PR[i]*60 #positions are stored in arcseconds
            D=PD[i]
           
            goodset=set(self.galaxies.where((self.galaxies.r < R) & \
               (self.galaxies["%s"%col] < D)).identifier)

            self.galaxies.photo_flag[numpy.array(\
               [_ in goodset for _ in self.galaxies.identifier])]=True
        
        self.galaxies=self.galaxies.where(self.galaxies.photo_flag==True)

        # ------------------------------------------------- 
  
        #change spectrscopicflag of any galaxy that should have spectroscopy:
        self.tryColumn('spec_flag',False)

        if SR!=['']:
          for i in range(len(SR)):
            R=SR[i]*60 #positions are stored in arcseconds
            D=SD[i]
            goodset=set(self.galaxies.where((self.galaxies.r < R) & \
               (self.galaxies["%s"%col] < D)).identifier)

            self.galaxies.spec_flag[numpy.array(\
               [_ in goodset for _ in self.galaxies.identifier])]=True

        return

# ----------------------------------------------------------------------------
# Following functions are designed to be run multiple times
# (Previous functions are single use/lightcone)
# ----------------------------------------------------------------------------

    def mimicPhotozError(self,sigma=0.1):
        #this code is not written well at the moment.

        e=sigma
        flag=self.galaxies.spec_flag==False ## Will give true for objects that have no spectrum
        z_obs=self.galaxies.z_spec*1.0
        for i in range(len(z_obs)):
            z=z_obs[i]
            if flag[i]==True: z=numpy.random.normal(z,e*(1+z))
            z_obs[i]=z
        self.galaxies.z_obs=z_obs

# ----------------------------------------------------------------------------

    # shortcut function for adding columns that might already exist
    def tryColumn(self,string,values):
        try:
            self.galaxies.add_column('%s'%string,values)
        except ValueError:
            self.galaxies["%s"%string]=values

# ----------------------------------------------------------------------------

    def snapToGrid(self, Grid):
        z=self.galaxies.z_obs
        sz,p=Grid.snap(z)
        self.tryColumn('Da_p',Grid.Da_p[p])
        self.tryColumn('rho_crit',Grid.rho_crit[p])
        self.tryColumn('sigma_crit',Grid.sigma_crit[p])
        self.tryColumn('beta',Grid.beta[p])
        rphys=self.galaxies.r*pangloss.arcmin2rad*self.galaxies.Da_p
        self.tryColumn('rphys',rphys)

# ----------------------------------------------------------------------------
#  One line description here!

    def drawMstars(self,model):
        Mhlist=self.galaxies
        redshiftList=self.galaxies.z_obs
        Ms=model.drawMstars(Mhlist,redshiftList)
        self.tryColumn('Mstar',Ms)

    def mimicMstarError(self,sigmaP,sigmaS):
        Ms=self.galaxies.Mstar.copy()
        Ms[self.galaxies.spec_flag==False]+=numpy.random.randn(Ms[self.galaxies.spec_flag==False].size)*sigmaP
        Ms[self.galaxies.spec_flag==True]+=numpy.random.randn(Ms[self.galaxies.spec_flag==True].size)*sigmaS
        try:
            self.galaxies.add_column('Ms_obs',Ms)
        except ValueError:
            self.galaxies.Ms_obs=Ms
            
# ----------------------------------------------------------------------------
#  One line description here!

    def drawMhalos(self,model):
        Mslist=self.galaxies.Ms_obs
        redshiftList=self.galaxies.z_obs
        Mhlist=model.drawMhalos(Mslist,redshiftList)
        self.tryColumn("Mh_obs",Mhlist)

# ----------------------------------------------------------------------------

    def drawConcentrations(self,errors=False):
        M200=10**self.galaxies.Mh_obs        
        r200 = (3*M200/(800*3.14159*self.galaxies.rho_crit))**(1./3)
        self.tryColumn("r200",r200)
        c200 = pangloss.MCrelation(M200,scatter=errors)
        self.tryColumn("c200",c200)
        r_s = r200/c200        
        self.tryColumn('rs',r_s)
        x=self.galaxies.rphys/r_s
        self.tryColumn('X',x)

# -----------------------------------------------------------------------------

    def makeKappas(self,errors=False,truncationscale=5,profile="BMO1"):
        c200=self.galaxies.c200
        r200=self.galaxies.r200
        x=self.galaxies.X
        r_s=self.galaxies.rs
        rho_s = pangloss.delta_c(c200)*self.galaxies.rho_crit
        kappa_s = rho_s * r_s /self.galaxies.sigma_crit
        r_trunc=truncationscale*r200
        xtrunc=r_trunc/r_s
        kappaHalo=kappa_s*1.0
        gammaHalo=kappa_s*1.0
        if profile=="BMO1":
            F=pangloss.BMO1Ffunc(x,xtrunc)
            G=pangloss.BMO1Gfunc(x,xtrunc)
        if profile=="BMO2":
            F=pangloss.BMO2Ffunc(x,xtrunc)
            G=pangloss.BMO2Gfunc(x,xtrunc)
        kappaHalo*=F
        gammaHalo*=(G-F)

        phi=self.galaxies.phi

        kappa = kappaHalo 
        gamma = gammaHalo
        gamma1 = gamma*numpy.cos(2*phi)
        gamma2 = gamma*numpy.sin(2*phi)

        self.tryColumn('kappa',kappa)
        self.tryColumn('gamma',gamma)
        self.tryColumn('gamma1',-gamma1)
        self.tryColumn('gamma2',-gamma2)
        
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

        self.tryColumn('kappa_keeton',kappa_keeton)
        self.tryColumn('gamma1_keeton',gamma1_keeton)
        self.tryColumn('gamma2_keeton',gamma2_keeton)

        self.tryColumn('kappa_tom',kappa_tom)
        self.tryColumn('gamma1_tom',gamma1_tom)
        self.tryColumn('gamma2_tom',gamma2_tom)

        self.tryColumn('kappa_add',K)
        self.tryColumn('gamma1_add',G1)
        self.tryColumn('gamma2_add',G2)


        self.kappa_add_total=numpy.sum(self.galaxies.kappa_add)
        self.kappa_keeton_total=numpy.sum(self.galaxies.kappa_keeton)
        self.kappa_tom_total=numpy.sum(self.galaxies.kappa_tom)

        self.gamma1_add_total=numpy.sum(self.galaxies.gamma1_add)
        self.gamma1_keeton_total=numpy.sum(self.galaxies.gamma1_keeton)
        self.gamma1_tom_total=numpy.sum(self.galaxies.gamma1_tom)

        self.gamma2_add_total=numpy.sum(self.galaxies.gamma2_add)
        self.gamma2_keeton_total=numpy.sum(self.galaxies.gamma2_keeton)
        self.gamma2_tom_total=numpy.sum(self.galaxies.gamma2_tom)

#=============================================================================
