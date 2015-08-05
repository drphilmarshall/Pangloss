import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import os, random, math, cmath, sys, timeit
import cPickle as pickle
from astropy.table import Table, Column
from matplotlib.patches import Ellipse
import treecorr
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar

# Import Pangloss:
PANGLOSS_DIR = os.path.expandvars("$PANGLOSS_DIR")
sys.path.append(PANGLOSS_DIR)
import pangloss

# Verbose
vb = True

# ============================================================================

class BackgroundCatalog(pangloss.Catalog):
    """
    NAME
        BackgroundCatalog

    PURPOSE
        Generate a catalog of source galaxies, and be able to plot or write out
        catalog data.

    COMMENTS
        Inherits from the super class Catalog in catalog.py

    INITIALISATION
        ???

    METHODS
        generate()
        write(output)

    BUGS

    AUTHORS
      This file is part of the Pangloss project, distributed under the
      GPL v2, by Tom Collett (IoA) and  Phil Marshall (Oxford).
      Please cite: Collett et al 2013, http://arxiv.org/abs/1303.6564

    HISTORY
      2015-06-29  Started Everett (SLAC)
    """
    def __init__(self,domain=None,field=None,N=10,mag_lim=[24.0,0.0],mass_lim=[10.0**6,10.0**12],z_lim=[1.3857,1.3857],sigma_e=0.2):
        self.type = 'background'

        # Only a subplot OR a field is allowed, not both!
        assert not (domain is not None and field is None)

        # The catalog can be created over a field corresponding to a foreground catalog (1 deg^2) or over an inputted subplot
        if domain is not None and field is not None:
            self.map_x = field[0]
            self.map_y = field[1]
            self.field_i = field[2]
            self.field_j = field[3]
            # NB: the domain and field may still be incompatible! Should really check for this.

        elif domain is None and field is not None:
            # Set domain based upon inputted field
            self.map_x = field[0]
            self.map_y = field[1]
            self.field_i = field[2]
            self.field_j = field[3]

            # Set ra and dec limits based upon field (x,y,i,j)
            ra_i = 2.0-self.map_x*4.0-self.field_i*1.0
            ra_f = 1.0-self.map_x*4.0-self.field_i*1.0
            dec_i = -2.0+self.map_y*4.0+self.field_j*1.0
            dec_f = -1.0+self.map_y*4.0+self.field_j*1.0

            # Set the domain to the inputted field
            domain = [ra_i,ra_f,dec_i,dec_f]

        elif domain is None and field is None:
            # If neither are inputted, use the field x=y=i=j=0:
            domain = [2.0,1.0,-2.0,-1.0]
            self.map_x = 0
            self.map_y = 0
            self.field_i = 0
            self.field_j = 0

        # Generate the background catalog
        self.generate(domain,N,mag_lim,mass_lim,z_lim,sigma_e)

        # The catalog keeps track of the number of excluded strongly-lensed galaxies, and add strong lensing flag
        self.strong_lensed_removed = 0
        self.galaxies['strong_flag'] = 0

        # Set source and strong-lens redshifts
        self.zl = 0       # There is no strong-lens present
        self.zs = 1.3857  # All source galaxies are at redshift 1.3857

        # Needed for lensing by halos
        self.planes = 100
        self.grid = None
        self.lightcones = None

        # Initialize ellipticity-ellipticity correlation attributes
        self.gg_none = None
        self.gg_map = None
        self.gg_halo = None

        # Initialize galaxy-galaxy correlation attributes
        self.ng_none = None
        self.ng_map = None
        self.ng_halo = None

        # Calls the superclass initialization for useful catalog attributes
        pangloss.Catalog.__init__(self)

        return

    def __str__(self):
        # *!Need to fix with new attributes!*
        return 'Background catalog with {} galaxies, '+ \
               'with redshifts ranging from {} to {}'\
                .format(self.galaxyCount,self.minZ,self.maxZ)

    def write(self,output=os.getcwd()):
        # Writes catalog data to current directory unless otherwise specified
        self.galaxies.write(output,format = 'ascii')
        return

# ----------------------------------------------------------------------------

    def generate(self,domain=None,N=10,mag_lim=[24.0,0.0],mass_lim=[10.0**6,10.0**12],z_lim=[1.3857,1.3857],sigma_e=0.2):
        '''
        Draw N-generated world-coordinate positions of galaxies in the sky per
        square arcminute inside a given domain of the form
        domain=[ra_init,ra_final,dec_init,dec_final]. The other optional inputs
        are value limits; any generated galaxy will have attributes within these
        values. Will make a scatter plot of the generated catalog only if
        plot = True.
        '''

        if domain == None:
            # Make a default domain (shouldn't be used except for testing or demo purposes)
            ra_init = np.deg2rad(2)    # initial value is larger as ra is left-handed
            dec_init = np.deg2rad(-2)
            ra_final = np.deg2rad(-2)
            dec_final = np.deg2rad(2)

        else:
            # Set ra and dec limits from domain. domain = [ra_init,ra_final,dec_init,dec_final]
            ra_init = np.deg2rad(domain[0])
            ra_final = np.deg2rad(domain[1])
            dec_init = np.deg2rad(domain[2])
            dec_final = np.deg2rad(domain[3])

        # Determine area of domain and the number of generated galaxies contained in it
        # (expecting wcs in degrees)
        self.Lx, self.Ly = abs(np.rad2deg(ra_final)-np.rad2deg(ra_init)), abs(np.rad2deg(dec_final)-np.rad2deg(dec_init))
        area = 3600*self.Lx*self.Ly # square arcminutes
        self.galaxy_count = int(N*area) # N galaxies per square arcminute

        # Initialize generated variables
        ra = np.zeros(self.galaxy_count)
        dec = np.zeros(self.galaxy_count)
        mag = np.zeros(self.galaxy_count)
        mass = np.zeros(self.galaxy_count)
        z = np.zeros(self.galaxy_count)
        e1_int = np.zeros(self.galaxy_count)
        e2_int = np.zeros(self.galaxy_count)
        eMod_int = np.zeros(self.galaxy_count)
        ePhi_int = np.zeros(self.galaxy_count)

        # Populate the generated variables
        ID = np.arange(self.galaxy_count)
        ra = np.random.uniform(ra_init,ra_final,self.galaxy_count)
        dec = np.random.uniform(dec_init,dec_final,self.galaxy_count)
        mag = np.random.uniform(mag_lim[0],mag_lim[1],self.galaxy_count)
        mass = np.random.uniform(mass_lim[0],mass_lim[1],self.galaxy_count)
        z = np.random.uniform(z_lim[0],z_lim[1],self.galaxy_count)
        e1_int = np.random.normal(0.0,sigma_e,self.galaxy_count)
        e2_int = np.random.normal(0.0,sigma_e,self.galaxy_count)

        # Change any |e|> 1 ellipticity components
        while (e1_int>1.0).any() or (e2_int>1.0).any():

            for i in [j for j in range(len(e1_int)) if e1_int[j]>1.0]:
                e1_int[i] = np.random.normal(0.0,sigma_e)

            for i in [j for j in range(len(e2_int)) if e2_int[j]>1.0]:
                e2_int[i] = np.random.normal(0.0,sigma_e)

        # Calculate Cartesian components of intrinsic complex ellipticity
        e_int = e1_int+1.0j*e2_int
        eMod_int = abs(e_int)
        ePhi_int = np.rad2deg(np.arctan2(e2_int,e1_int))/2.0

        # Save generated catalog as an astropy table
        self.galaxies = Table([ID,ra,dec,mag,mass,z,eMod_int,ePhi_int,e1_int,e2_int],names=['ID','RA','Dec','mag','Mstar_obs','z_obs','eMod_int','ePhi_int','e1_int','e2_int'], \
                              meta={'name':'background catalog','size':N,'mag_lim':mag_lim, \
                                    'mass_lim':mass_lim,'z_lim':z_lim,'sigma_e':sigma_e})

        return

    def lens_by_map(self,kappamap,shearmap,mag_lim=[0,24],mass_lim=[0,10**20],z_lim=[0,1.3857]):
        '''
        Lense background galaxies by the shear and convergence in their respective Kappamaps and Shearmaps.
        '''

        # Exctract needed data from catalog galaxies
        #galaxies = pangloss.Catalog.return_galaxies(self,mag_lim,mass_lim,z_lim,ra_lim,dec_lim)
        ra = np.rad2deg(self.galaxies['RA'])
        dec = np.rad2deg(self.galaxies['Dec'])
        e1_int = self.galaxies['e1_int']
        e2_int = self.galaxies['e2_int']

        # Initialize new variables
        kappa = np.zeros(self.galaxy_count)
        gamma1 = np.zeros(self.galaxy_count)
        gamma2 = np.zeros(self.galaxy_count)
        g = np.zeros(self.galaxy_count)
        e1 = np.zeros(self.galaxy_count)
        e2 = np.zeros(self.galaxy_count)
        eMod = np.zeros(self.galaxy_count)
        ePhi = np.zeros(self.galaxy_count)

        # Extract convergence and shear values at each galaxy location from maps
        for i in range(self.galaxy_count):
            kappa[i] = kappamap.at(ra[i],dec[i],mapfile=0)
            gamma1[i] = shearmap.at(ra[i],dec[i],mapfile=0)
            gamma2[i] = shearmap.at(ra[i],dec[i],mapfile=1)

        # Calculate the reduced shear g and its conjugate g_conj
        g = (gamma1 + 1j*gamma2)/(1.0-kappa)
        g_conj = np.array([val.conjugate() for val in g])

        # Flag any galaxy that has been strongly lensed
        self.galaxies['strong_flag'][abs(g)>1.0] = 1
        #assert abs(g).all() < 1.0, 'error: strong lensing for {} galaxies. k = {}, gamma = {}, g = {}. Locations at ra={}, dec={}.'.format(len(g[abs(g)>1.0]),kappa[abs(g)>1.0],gamma1[abs(g)>1.0]+1j*gamma2[abs(g)>1.0],g[abs(g)>1.0],ra[abs(g)>1.0],dec[abs(g)>1.0])

        # Calculate the observed ellipticity
        e = ((e1_int + 1j*e2_int) + g)/(1.0+g_conj * (e1_int + 1j*e2_int))
        e1, e2 = e.real, e.imag
        eMod = np.abs(e)
        ePhi = np.rad2deg(np.arctan2(e2,e1))/2.0

        # Add convergence and shear values to catalog
        self.galaxies['kappa'] = kappa
        self.galaxies['gamma1'] = gamma1
        self.galaxies['gamma2'] = gamma2
        self.galaxies['g'] = g
        self.galaxies['e1'] = e1
        self.galaxies['e2'] = e2
        self.galaxies['eMod'] = eMod
        self.galaxies['ePhi'] = ePhi

        return

    def lens_by_halos(self,save=False,methods=['add'],use_method='add'):
        '''
        Lens background galaxies by the combined shear and convergence in their respective lightcones.
        '''

        # Grid should already be setup, but set if not
        if self.grid is None:
            self.setup_grid()

        # Lightcones should already be drilled, but drill if not
        if self.lightcones is None:
            self.drill_lightcones()

        # Initialize new variables
        assert self.galaxy_count == len(self.lightcones)
        kappa = np.zeros(self.galaxy_count)
        gamma1 = np.zeros(self.galaxy_count)
        gamma2 = np.zeros(self.galaxy_count)
        g = np.zeros(self.galaxy_count)
        e1 = np.zeros(self.galaxy_count)
        e2 = np.zeros(self.galaxy_count)
        eMod = np.zeros(self.galaxy_count)
        ePhi = np.zeros(self.galaxy_count)

        # Keep track of how long each lightcone takes to process
        runtimes = np.zeros(self.galaxy_count)

        # Set the counter to be 10%
        counter = np.ceil(len(self.lightcones)/10.0)

        # Calculate lensing in each lightcone:
        for lightcone in self.lightcones:
            start_time = timeit.default_timer()
            if lightcone.ID%counter == 0 and vb is True:
                print lightcone.ID,' ',np.ceil(100*lightcone.ID/self.galaxy_count),'%'

            '''
            # Set the stellar mass - halo mass relation
            shmr = pangloss.SHMR(method='Behroozi')
            HMFfile = PANGLOSS_DIR+'/calib/SHMR/HaloMassRedshiftCatalog.pickle'
            shmr.makeHaloMassFunction(HMFfile)
            shmr.makeCDFs()

            # Simulated lightcones need mock observed Mstar_obs values
            # drawing from their Mhalos:
            lightcone.drawMstars(shmr)

            # Draw Mhalo from Mstar, and then c from Mhalo:
            lightcone.drawMhalos(shmr)
            '''

            lightcone.drawConcentrations(errors=True)

            # Compute each halo's contribution to the convergence and shear:
            lightcone.makeKappas(truncationscale=10)

            # Combine all contributions into a single kappa and gamma for the lightcone
            lightcone.combineKappas(methods=methods)

            # Populate the kappa and gamma values using the 'use_method'
            if use_method == 'add' :
                # Use the add_total combination
                kappa[lightcone.ID] = lightcone.kappa_add_total
                gamma1[lightcone.ID] = lightcone.gamma1_add_total
                gamma2[lightcone.ID] = lightcone.gamma2_add_total

            elif use_method == 'keeton':
                # Use the keeton_total combination
                kappa[lightcone.ID] = lightcone.kappa_keeton_total
                gamma1[lightcone.ID] = lightcone.gamma1_keeton_total
                gamma2[lightcone.ID] = lightcone.gamma2_keeton_total

            elif use_method == 'tom':
                # Use the tom_total combination
                kappa[lightcone.ID] = lightcone.kappa_tom_total
                gamma1[lightcone.ID] = lightcone.gamma1_tom_total
                gamma2[lightcone.ID] = lightcone.gamma2_tom_total

            elapsed = timeit.default_timer() - start_time
            runtimes[lightcone.ID] = elapsed

        if vb is True:
            print 'average CPU time per background galaxy: ',np.mean(runtimes),'+/-',np.std(runtimes)

        #-------------------------------------------------------------------------------------
        # Use the halo model's kappa and gamma values to compute the new galaxy ellipticities

        # Extract background galaxy intrinsic ellipticites
        e1_int = self.galaxies['e1_int']
        e2_int = self.galaxies['e2_int']

        # Calculate the reduced shear g and its conjugate g_conj
        g = (gamma1 + 1j*gamma2)/(1.0-kappa)
        g_conj = np.array([val.conjugate() for val in g])

        # Flag any galaxy that has been strongly lensed
        self.galaxies['strong_flag'][abs(g)>1.0] = 1
        #assert abs(g).all() < 1.0, 'error: strong lensing for {} galaxies. k = {}, gamma = {}, g = {}. Locations at ra={}, dec={}.'.format(len(g[abs(g)>1.0]),kappa[abs(g)>1.0],gamma1[abs(g)>1.0]+1j*gamma2[abs(g)>1.0],g[abs(g)>1.0],ra[abs(g)>1.0],dec[abs(g)>1.0])

        # Calculate the observed ellipticity
        e = ((e1_int + 1j*e2_int) + g)/(1.0+g_conj * (e1_int + 1j*e2_int))
        e1, e2 = e.real, e.imag
        eMod = np.abs(e)
        ePhi = np.rad2deg(np.arctan2(e2,e1))/2.0

        # Add convergence and shear values to catalog
        self.galaxies['kappa_halo'] = kappa
        self.galaxies['gamma1_halo'] = gamma1
        self.galaxies['gamma2_halo'] = gamma2
        self.galaxies['g_halo'] = g
        self.galaxies['e1_halo'] = e1
        self.galaxies['e2_halo'] = e2
        self.galaxies['eMod_halo'] = eMod
        self.galaxies['ePhi_halo'] = ePhi

        # Save catalog with new lensing values to enable a restart
        if save == True:
            self.save()

        return


    def add_noise(self,M=1,sigma_obs=0.1):
        '''
        Add measurement and shape noise to the background galaxy intrinsic shapes.
        '''

        # Extract data that is to have noise added
        e1 = self.galaxies['e1']
        e2 = self.galaxies['e2']

        # Multiplicative shear calibration error:
        # We tend to systematically underestimate the ellipticity of background galaxies.
        # Multiplying by M < 1 accounts for this error.
        e1 = M*e1
        e2 = M*e2

        # Measurement noise:
        e1 += np.random.normal(0.0,sigma_obs,self.galaxy_count)
        e2 += np.random.normal(0.0,sigma_obs,self.galaxy_count)

        # Change any |e|> 1 ellipticity components
        while (e1>1.0).any() or (e2>1.0).any():

            for i in [j for j in range(len(e1)) if e1[j]>1.0]:
                e1[i] = np.random.normal(0.0,sigma_obs)

            for i in [j for j in range(len(e2)) if e2[j]>1.0]:
                e2[i] = np.random.normal(0.0,sigma_obs)

        # Calculate noisy modulus
        eMod = np.sqrt(e1**2+e2**2)

        # Save new noisy ellipticities
        self.galaxies['e1'] = e1
        self.galaxies['e2'] = e2
        self.galaxies['eMod'] = eMod

        return

    def setup_grid(self):
        '''
        Create the distance grid to simplify distance calculations.
        '''

        # Make redshift grid:
        self.grid = pangloss.Grid(self.zl,self.zs,nplanes=self.planes)

        return

    def drill_lightcones(self,radius=2.0,foreground=None,save=False):
        '''
        Drill a lightcone at each background source with radius in arcmin. Will
        write the lightcones to pangloss/data/lightcones. The method can only be
        used with BackgroundCatalogs that were created as a field, not a subplot.
        '''

        # Retrieve background galaxy data and initialize the lightcones
        self.lightcones = np.zeros(self.galaxy_count,dtype='object')

        # Set lightcone parameters
        flavor = 'simulated'

        # Setup the grid
        self.setup_grid()

        # If a foreground catalog object is not passed, load in the appropriate catalog.
        if foreground is None:
            # Load in the corresponding foreground catalog
            config = pangloss.Configuration(PANGLOSS_DIR+'/example/example.config')
            foreground = pangloss.ForegroundCatalog(PANGLOSS_DIR+'/data/GGL_los_8_'+str(self.map_x)+'_'+str(self.map_y)+'_'+str(self.field_i)+'_'+str(self.field_j)+'_N_4096_ang_4_Guo_galaxies_on_plane_27_to_63.images.txt',config)

        # Only take the columns from foreground that are needed for a lightcone
        lc_catalog = Table(names=['RA','Dec','z_obs','Mhalo_obs','Type'],data=[foreground.galaxies['RA'],foreground.galaxies['Dec'],foreground.galaxies['z_obs'],foreground.galaxies['Mhalo_obs'],foreground.galaxies['Type']])
        del foreground

        # Set the counter to be 10%
        counter = np.ceil(self.galaxy_count/10.0)

        # Drill a lightcone at each galaxy location
        for i in range(self.galaxy_count):
            if i%counter == 0 and vb is True:
                print i,' ',np.ceil(100*i/self.galaxy_count),'%'
            # Set galaxy positions
            ra0 = self.galaxies['RA'][i]
            dec0 = self.galaxies['Dec'][i]
            position = [ra0,dec0]

            # Create the lightcone for background galaxy i
            self.lightcones[i] = pangloss.Lightcone(lc_catalog,flavor,position,radius,ID=i)

            # Create the redshift scaffolding for lightcone i:
            self.lightcones[i].defineSystem(self.zl,self.zs)
            self.lightcones[i].loadGrid(self.grid)
            self.lightcones[i].snapToGrid(self.grid)

        if save == True:
            # Write out this entire catalog to the 'data' directory
            self.save()

        return

# ----------------------------------------------------------------------------

    def save(self,filename=None):
        '''
        Save the background galaxy catalog in '/data'.
        '''

        if vb == True:
            print 'Pickling this background catalog to disk...'

        # If no filename is given, a default is set using the field numbers (x,y,i,j).
        if filename is None:
            filename = PANGLOSS_DIR+'/data/background_catalog_'+str(self.map_x)+'_'+str(self.map_y)+'_'+str(self.field_i)+'_'+str(self.field_j)+'.pkl'

        self.filename = filename

        pickle_file = open(self.filename, 'wb')
        pickle.dump(self.__dict__,self.filename,2)
        pickle_file.close()

        return

    def load(self,filename=None):
        '''
        Load in an old background galaxy catalog. Note: Not sure when we will use this yet.
        '''
        if filename is None:
            try: filename = self.filename
            except:
                raise Exception('No file to load from! Exiting.')

        pickle_file = open(self.filename,'rb')
        tmp_dict = pickle.load(pickle_file)
        pickle_file.close()

        self.__dict__.update(tmp_dict)

        return

    def bin_to_map(self,data=['kappa','gamma1','gamma2'],binsize=1.0,bin_units='arcmin'):
        '''
        Bin the background galaxies into WLMaps with data set in the list 'data'. Binsize is in units of arcmin.
        '''
        # If a single string is passed for data, turn into list
        if type(data) == str:
            data = [data]

        # Finish!!

# ----------------------------------------------------------------------------

    def calculate_corr(self,corr_type='gg',min_sep=0.1,max_sep=30.0,sep_units='arcmin',binsize=None,N=15.0,lensed='map',foreground=None):
        '''
        Calculate the inputted correlation function type from min_sep<dtheta<max_sep. If no binsize or
        number of bins (N) are inputted, the binsize is automatically calculated using 15 bins. The 'lensed'
        argument is only used for shear-shear correlation (gg).
        '''

        galaxies = self.galaxies

        # If none is given, calculate (log) binsize based upon separation limit values
        if binsize == None:
            binsize = np.log10(1.0*max_sep/min_sep)/(1.0*N)

        # Calculate the shear-shear (or ellipticity-ellipticity) correlation function
        if corr_type == 'gg':
            # Create catalog of the pre or post-lensed background galaxies and their ellipticities
            if lensed == 'map':
                corr_cat = treecorr.Catalog(ra=galaxies['RA'], dec=galaxies['Dec'], g1=galaxies['e1'], g2=galaxies['e2'], ra_units='rad', dec_units='rad')
            elif lensed == 'halo':
                corr_cat = treecorr.Catalog(ra=galaxies['RA'], dec=galaxies['Dec'], g1=galaxies['e1_halo'], g2=galaxies['e2_halo'], ra_units='rad', dec_units='rad')
            elif lensed == 'none':
                corr_cat = treecorr.Catalog(ra=galaxies['RA'], dec=galaxies['Dec'], g1=galaxies['e1_int'], g2=galaxies['e2_int'], ra_units='rad', dec_units='rad')

            # Set g-g correlation parameters
            gg = treecorr.GGCorrelation(bin_size=binsize, min_sep=min_sep, max_sep=max_sep, sep_units=sep_units, bin_slop=0.05/binsize)

            # Calculate g-g correlation function
            gg.process(corr_cat)

            # Check to make sure none of the values are Nan's (Fix in fugure using 0 weights for galaxies not in K/S maps)
            assert not np.isnan(gg.xip).any()

            # Save correlation object to catalog
            if lensed == 'none': self.gg_none = gg
            elif lensed == 'map': self.gg_map = gg
            elif lensed == 'halo': self.gg_halo = gg

            return gg

        # Calculate the galaxy-galaxy correlation function
        elif corr_type == 'ng':
            # Load in foreground catalog if none is passed
            if foreground is None:
                config = pangloss.Configuration(PANGLOSS_DIR+'/example/example.config')
                foreground = pangloss.ForegroundCatalog(PANGLOSS_DIR+'/data/GGL_los_8_'+str(self.map_x)+'_'+str(self.map_y)+'_'+str(self.field_i)+'_'+str(self.field_j)+'_N_4096_ang_4_Guo_galaxies_on_plane_27_to_63.images.txt',config)

            # Create catalog of the foreground galaxy locations
            corr_cat1 = treecorr.Catalog(ra=foreground.galaxies['RA'], dec=foreground.galaxies['Dec'], ra_units='rad', dec_units='rad')

            # Create catalog of the background galaxy ellipticities
            if lensed == 'map':
                corr_cat2 = treecorr.Catalog(ra=galaxies['RA'], dec=galaxies['Dec'], g1=galaxies['e1'], g2=galaxies['e2'], ra_units='rad', dec_units='rad')
            elif lensed == 'halo':
                corr_cat2 = treecorr.Catalog(ra=galaxies['RA'], dec=galaxies['Dec'], g1=galaxies['e1_halo'], g2=galaxies['e2_halo'], ra_units='rad', dec_units='rad')
            elif lensed == 'none':
                corr_cat2 = treecorr.Catalog(ra=galaxies['RA'], dec=galaxies['Dec'], g1=galaxies['e1_int'], g2=galaxies['e2_int'], ra_units='rad', dec_units='rad')

            # Set n-g correlation parameters
            ng = treecorr.NGCorrelation(bin_size=binsize, min_sep=min_sep, max_sep=max_sep, sep_units=sep_units, bin_slop=0.05/binsize)

            # Calculate n-g correlation function
            ng.process(corr_cat1,corr_cat2)

            # Check to make sure none of the values are Nan's (Fix in fugure using 0 weights for galaxies not in K/S maps)
            assert not np.isnan(ng.xi).any()

            # Save correlation object to catalog
            if lensed == 'none': self.ng_none = ng
            elif lensed == 'map': self.ng_map = ng
            elif lensed == 'halo': self.ng_halo = ng

            return ng

        else:
            # Add other correlation types later if necessary
            pass

    def compare_corr(self,corr1,corr2,corr_type='gg',corr_comp='plus',percent_type='error'):
        '''
        Compares two correlation function components (with the same binning and separation values) using
        a chi-squared approximation and a mean percent difference (if both correlations are predicted quantities)
        or a mean percent error (if one of the correlations is an observed quantity). Note that by default the
        method computes a percent error rather than a percent difference, and that the observed (which is our theoretical)
        value must be inputted as corr2.
        '''

        # For ellipticity-ellipticity correlation functions:
        if corr_type == 'gg':
            # Extract the correlation values
            if corr_comp == 'plus':
                y1, y2 = corr1.xip, corr2.xip
                var1, var2 = corr1.varxi, corr2.varxi

            elif corr_comp == 'minus':
                y1, y2 = corr1.xim, corr2.xim
                var1, var2 = corr1.varxi, corr2.varxi

            elif corr_comp == 'cross':
                y1, y2 = 0.5*(corr1.xim_im-corr1.xip_im), 0.5*(corr2.xim_im-corr2.xip_im)
                var1, var2 = 0.5*np.sqrt(2*corr1.varxi), 0.5*np.sqrt(2*corr2.varxi)

            elif corr_comp == 'cross_prime':
                y1, y2 = corr1.xip_im, corr2.xip_im
                var1, var2 = corr1.varxi, corr2,varxi

        # For galaxy-galaxy correlation functions:
        elif corr_type == 'ng':
            # Extract the correlation values
            if corr_comp == 'real':
                y1, y2 = corr1.xi, corr2.xi

            elif corr_comp == 'imag':
                y1, y2 = corr.xi_im, corr2.xi_im

            # Extract the variance
            var1, var2 = corr1.varxi, corr2.varxi

        else:
            # Can add more correlation types here if needed
            pass

        # Check that each correlation object is the same size, and set the degrees of freedom
        assert np.size(y1) == np.size(y2)
        N = np.size(y1)

        # Calculate the chi-squared value
        chi2 = np.sum( (y1-y2)**2 / (1.0*var1 + 1.0*var2) )

        # Determine the significance
        n_sigma = np.sqrt(2.0*chi2) - np.sqrt(2.0*N)

        # Calculate the weight of each correlation value pair
        w = 1.0/(np.sqrt( 1.0*var1 + 1.0*var2 ))

        if percent_type == 'error':
            # Calculate mean percent error between the predicted and observed correlation function
            # values at different separations
            percent_err = abs( (y1-y2) / (1.0*y2) ) * 100.0
            mean_err = np.average(percent_err,weights=w)

        elif percent_type == 'difference':
            # Calculate mean percent difference between correlation function values at different separations
            percent_diff = abs( (y1-y2) / (0.5*y1+0.5*y2) ) * 100.0
            mean_err = np.average(percent_diff,weights=w)

        return chi2, n_sigma, mean_err


    def plot(self,subplot=None,mag_lim=[0,24],mass_lim=[0,10**20],z_lim=[0,1.3857],fig_size=10,graph='scatter',lensed='map'):
        '''
        Make scatter plot of generated galaxies.
        '''

        # Get current figure (or make one if it doesn't exist)
        fig = plt.gcf()

        # If there is a Pangloss map open:
        if fig._label == 'Pangloss Map':
            # Adopt axes from the open Kappamap:
            imshow = fig.axes[0]
            world = fig.axes[1]

            # If the Kappamap subplot was not passed to this Shearmap:
            if subplot == None:
                # Adopt subplot from the open Kappamap:
                fig.sca(world)
                subplot = plt.axis()

            # Adopt figure size from open Kappamap:
            fig_size = plt.gcf().get_size_inches()[0]

        # Otherwise:
        else:
            if subplot is None:
                # Default subplot is entire catalog
                ai, di = self.ra_max, self.dec_min
                af, df = self.ra_min, self.dec_max
                subplot = [ai,af,di,df]

            # Adjust the subplot in wcs by half a pixel
            #subplot = [subplot[0]-self.PIXSCALE[0]/2.0,subplot[1]-self.PIXSCALE[0]/2.0,subplot[2]-self.PIXSCALE[0]/2.0,subplot[3]-self.PIXSCALE[0]/2.0]

            # Create new imshow and world axes
            imshow, world = pangloss.make_axes(fig,subplot)

        ai, af = subplot[0], subplot[1]    # RA limits for subplot
        di, df = subplot[2], subplot[3]    # DEC limits for subplot
        Lx, Ly = abs(ai-af), abs(di-df)    # Length of axes in wcs
        L = np.mean([Lx,Ly])

        # Find the galaxies that are within the limits, and extract the useful data from them
        ra_lim, dec_lim = [ai, af], [di, df]
        galaxies = self.return_galaxies(mag_lim,mass_lim,z_lim,ra_lim,dec_lim)
        ra = np.rad2deg(galaxies['RA'])
        dec = np.rad2deg(galaxies['Dec'])
        mass = galaxies['Mstar_obs']

        # The angles are flipped as we are using a left-handed coordinate reference axis for plotting
        if graph == 'ellipse' or graph == 'stick':
            if lensed == 'none':
                # Extract intrinsic ellipticity
                eMod_int = galaxies['eMod_int']
                ePhi_int = -galaxies['ePhi_int']

            elif lensed == 'map':
                # Extract lensed-by-map ellipticity
                eMod = galaxies['eMod']
                ePhi = -galaxies['ePhi']

            elif lensed == 'halo':
                # Extract lensed-by-halos ellipticity
                eMod_halo = galaxies['eMod_halo']
                ePhi_halo = -galaxies['ePhi_halo']

            elif lensed == 'both':
                # Extract both the intrinsic and lensed-by-map ellipticity
                eMod_int = galaxies['eMod_int']
                ePhi_int = -galaxies['ePhi_int']
                eMod = galaxies['eMod']
                ePhi = -galaxies['ePhi']

            elif lensed == 'all':
                # Extract both the intrinsic and both lensed ellipticities
                eMod_int = galaxies['eMod_int']
                ePhi_int = -galaxies['ePhi_int']
                eMod = galaxies['eMod']
                ePhi = -galaxies['ePhi']
                eMod_halo = galaxies['eMod_halo']
                ePhi_halo = -galaxies['ePhi_halo']

        # Set current axis to world coordinates and set the limits
        fig.sca(world)
        world.set_xlim(subplot[0],subplot[1])
        world.set_ylim(subplot[2],subplot[3])

        if graph == 'scatter':
            # Scale size of point by the galaxy mass
            s = [math.log(mass[i]) for i in range(0,len(mass))]
            plt.scatter(ra,dec,s,alpha=0.5,edgecolor=None,color='blue')

        elif graph == 'ellipse':
            # Scale galaxy plot size by its mass?
            # scale = ((np.log10(mass)-9.0)/(12.0-9.0))
            scale = 0.5
            floor = 0.01
            size = 0.01*(scale*(scale > 0) + floor)

            # Plot each galaxy as an ellipse
            for i in range(np.shape(galaxies)[0]):
                if lensed == 'none':
                    # Plot intrinsic ellipticities
                    alpha = 0.25
                    pangloss.plotting.plot_ellipse(ra[i],dec[i],size,eMod_int[i],ePhi_int[i],world,'blue',alpha)

                elif lensed == 'map':
                    # Plot lensed-by-map ellipticities
                    alpha = 0.3
                    pangloss.plotting.plot_ellipse(ra[i],dec[i],size,eMod[i],ePhi[i],world,'green',alpha)

                elif lensed == 'halo':
                    # Plot lensed-by-halo ellipticities
                    alpha = 0.3
                    pangloss.plotting.plot_ellipse(ra[i],dec[i],size,eMod_halo[i],ePhi_halo[i],world,'purple',alpha)

                elif lensed == 'both':
                    # Plot both lensed-by-map and intrinsic ellipticities
                    alpha1 = 0.25
                    alpha2 = 0.3
                    pangloss.plotting.plot_ellipse(ra[i],dec[i],size,eMod_int[i],ePhi_int[i],world,'blue',alpha1)
                    pangloss.plotting.plot_ellipse(ra[i],dec[i],size,eMod[i],ePhi[i],world,'green',alpha2)

                elif lensed == 'all':
                    # Plot both types of lensed and intrinsic ellipticities
                    alpha1 = 0.25
                    alpha2 = 0.3
                    alpha3 = 0.3
                    pangloss.plotting.plot_ellipse(ra[i],dec[i],size,eMod_int[i],ePhi_int[i],world,'blue',alpha1)
                    pangloss.plotting.plot_ellipse(ra[i],dec[i],size,eMod[i],ePhi[i],world,'green',alpha2)
                    pangloss.plotting.plot_ellipse(ra[i],dec[i],size,eMod_halo[i],ePhi_halo[i],world,'purple',alpha3)

        elif graph == 'stick':
            if lensed == 'none':
                # Plot intrinsic ellipticity sticks
                pangloss.plotting.plot_sticks(ra,dec,eMod_int,ePhi_int,world,'blue')

            elif lensed == 'map':
                # Plot lensed-by-map ellipticity sticks
                pangloss.plotting.plot_sticks(ra,dec,eMod,ePhi,world,'green')

            elif lensed == 'halo':
                # Plot lensed-by-halos ellipticity sticks
                pangloss.plotting.plot_sticks(ra,dec,eMod_halo,ePhi_halo,world,'purple')

            elif lensed == 'both':
                # Plot both lensed and intrinsic ellipticity sticks
                pangloss.plotting.plot_sticks(ra,dec,eMod_int,ePhi_int,world,'blue')
                pangloss.plotting.plot_sticks(ra,dec,eMod,ePhi,world,'green')

            elif lensed == 'all':
                # Plot both types of lensed and intrinsic ellipticity sticks
                pangloss.plotting.plot_sticks(ra,dec,eMod_int,ePhi_int,world,'blue')
                pangloss.plotting.plot_sticks(ra,dec,eMod,ePhi,world,'green')
                pangloss.plotting.plot_sticks(ra,dec,eMod_halo,ePhi_halo,world,'purple')

            # Add scale bar
            if lensed == 'map':
                # Plot as green
                color = (0,0.6,0,1)
            if lensed == 'halo':
                # Plot as purple
                color = (1,0,1,1)
            if lensed == 'none':
                # Plot as blue
                color = (0,0,1,1)
            else:
                # Plot as black
                color = (0,0,0,1)

            bar = AnchoredSizeBar(world.transData,L/10.0,'10% Ellipticity',pad=0.5,loc=4,sep=5,borderpad=0.25,frameon=True)
            bar.size_bar._children[0]._linewidth = 2
            bar.size_bar._children[0]._edgecolor = color
            world.add_artist(bar)

        # Label axes and set the correct figure size
        plt.xlabel('Right Ascension / deg')
        plt.ylabel('Declination / deg')
        pangloss.set_figure_size(fig,fig_size,Lx,Ly)

        return
