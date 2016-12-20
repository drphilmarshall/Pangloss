import cPickle as pickle
import math
import os
import timeit
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
import pangloss
from pandas import DataFrame

# Fast correlation functions:
try:
    import treecorr
except ImportError:
    import pangloss.nocorr as treecorr


# Verbose
vb = False

# Record CPU time per lightcone?
time = False

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
        generate
        write
        lens_by_map
        lens_by_halos
        add_noise
        setup_grid
        set_importance
        drill_lightcones
        save
        load
        bin_to_map
        calculate_corr
        compare_corr
        plot

    BUGS

    HISTORY
      2015-06-29  Started Everett (SLAC)
    """

# ============================================================================

    def __init__(self,domain=None,field=None,N=10,mag_lim=[24.0,0.0],mass_lim=[10.0**6,10.0**12],z_lim=[1.3857,1.3857],sigma_e=0.2,spacing=None):
        self.type = 'background'

        # A domain must come with a corresponding field!
        assert not (domain is not None and field is None)

        # The catalog can be created over a field corresponding to a foreground catalog (1 deg^2) or over an inputted subplot
        if domain is not None and field is not None:
            self.domain = domain
            self.field = field
            self.map_x = field[0]
            self.map_y = field[1]
            self.field_i = field[2]
            self.field_j = field[3]
            # NB: the domain and field may still be incompatible! Should really check for this.

        elif domain is None and field is not None:
            # Set domain based upon inputted field
            self.field = field
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
            self.domain = [ra_i,ra_f,dec_i,dec_f]

        elif domain is None and field is None:
            # If neither are inputted, use the field x=y=i=j=0:
            self.domain = [2.0,1.0,-2.0,-1.0]
            self.field = [0,0,0,0]
            self.map_x = 0
            self.map_y = 0
            self.field_i = 0
            self.field_j = 0

        # Generate the background catalog
        self.generate(self.domain,N,mag_lim,mass_lim,z_lim,sigma_e,spacing)

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

# ----------------------------------------------------------------------------

    def __str__(self):
        # *!Need to fix with new attributes!*
        return 'Background catalog with {} galaxies, with redshifts ranging from {} to {}'.format(self.galaxy_count,self.minZ,self.maxZ)

# ----------------------------------------------------------------------------

    def write(self,output=os.getcwd()):
        # Writes catalog data to current directory unless otherwise specified
        self.galaxies.write(output,format = 'ascii')
        return

# ----------------------------------------------------------------------------

    def generate(self,domain=None,N=10,mag_lim=[24.0,0.0],mass_lim=[10.0**6,10.0**12],z_lim=[1.3857,1.3857],sigma_e=0.2,spacing=None):
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

        # Determine area of domain and the number of generated galaxies contained in it,
        # as well as fenerate positions based upon random locations or a uniform spacing
        # (expecting wcs in degrees)
        self.Lx, self.Ly = abs(np.rad2deg(ra_final)-np.rad2deg(ra_init)), abs(np.rad2deg(dec_final)-np.rad2deg(dec_init))
        if spacing is None:
            # Determine number of galaxies at uniformly distributed positions
            # by the passed galaxy density N
            area = 3600*self.Lx*self.Ly # square arcminutes
            self.galaxy_count = int(N*area) # N galaxies per square arcminute
            ra = np.random.uniform(ra_init,ra_final,self.galaxy_count)
            dec = np.random.uniform(dec_init,dec_final,self.galaxy_count)
        else:
            # If a spacing is passed, overrule N to make a uniformly spaced
            # grid of background galaxies with separation distance 'spacing'
            # in rad (spacing = 1.705e-5 rad/pix for kappa/shear map density)
            r = np.arange(ra_final,ra_init,spacing) # swapped b/c left-handed
            d = np.arange(dec_init,dec_final,spacing)
            self.galaxy_count = np.size(r)*np.size(d)
            ra, dec = np.meshgrid(r,d)
            ra, dec = ra.reshape(self.galaxy_count), dec.reshape(self.galaxy_count)

        # Populate the generated variables
        ID = np.arange(self.galaxy_count)
        mag = np.random.uniform(mag_lim[0],mag_lim[1],self.galaxy_count)
        mass = np.random.uniform(mass_lim[0],mass_lim[1],self.galaxy_count)
        z = np.random.uniform(z_lim[0],z_lim[1],self.galaxy_count)
        e1_int = np.random.normal(0.0,sigma_e,self.galaxy_count)
        e2_int = np.random.normal(0.0,sigma_e,self.galaxy_count)

        # Save intrinsic ellipticity std
        self.std_int = sigma_e

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

        # Save generated catalog as a pandas dataframe
        columns = ['ID','RA','Dec','mag','Mstar_obs','z_obs','eMod_int','ePhi_int','e1_int','e2_int']
        data = np.matrix([ID,ra,dec,mag,mass,z,eMod_int,ePhi_int,e1_int,e2_int]).transpose()
        self.galaxies = DataFrame(columns=columns, data=data)
        return

    def lens_by_map(self,kappamap,shearmap):
        '''
        Lense background galaxies by the shear and convergence in their respective Kappamaps and Shearmaps.
        '''

        # Exctract needed data from catalog galaxies
        #galaxies = pangloss.Catalog.return_galaxies(self,mag_lim,mass_lim,z_lim,ra_lim,dec_lim)
        ra = np.rad2deg(self.galaxies['RA'])
        dec = np.rad2deg(self.galaxies['Dec'])
        e1_int = self.galaxies['e1_int']
        e2_int = self.galaxies['e2_int']

        # Initialize new variables (note: e and g have to be initialized as complex for memory allocation)
        kappa = np.zeros(self.galaxy_count)
        gamma1 = np.zeros(self.galaxy_count)
        gamma2 = np.zeros(self.galaxy_count)
        g = (gamma1 + 1j*gamma2)/(1.0-kappa)
        e1 = np.zeros(self.galaxy_count)
        e2 = np.zeros(self.galaxy_count)
        e = e1+1j*e2
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

        # Flag any galaxy that has been strongly (or near-strongly) lensed
        self.galaxies['strong_flag'] = np.abs(g)>0.5

        # Calculate the observed ellipticity for weak lensing events
        index = np.abs(g) < 1.0
        e[index] = ( (e1_int[index] + 1j*e2_int[index]) + g[index]) / (1.0+g_conj[index] * (e1_int[index] + 1j*e2_int[index]) )

        # Calculate the observed ellipticity for strong lensing events
        index = ~index
        e1_int_conj = np.array([val.conjugate() for val in e1_int])
        e2_int_conj = np.array([val.conjugate() for val in e2_int])
        e[index] = (1.0 + g[index]*(e1_int_conj[index]+1j*e2_int_conj[index]) ) / ( (e1_int_conj[index]+1j*e2_int_conj[index]) + g_conj[index])

        # Calculate Cartesian and polar components
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

        # Note: For now, we are removing any background galaxy that is strongly lensed by a map.
        self.galaxies = self.galaxies[self.galaxies['strong_flag']!=1]
        self.strong_lensed_removed = np.sum(self.galaxies['strong_flag'])
        self.galaxy_count -= self.strong_lensed_removed

        return

    def lens_by_halos(self, save=False, methods=['add'], use_method='add', relevance_lim=10e-6,
                      lookup_table=False, smooth_corr=False, foreground_corr=False):
        '''
        Lens background galaxies by the combined shear and convergence in their respective lightcones using
        the method given by `use_method`. By default all foreground objects in a lightcone are used in the
        calculation, but this can be changed by setting the 'relevance_lim' higher.
        '''

        # Grid should already be setup, but set if not
        if self.grid is None:
            self.setup_grid()

        # Lightcones should already be drilled, but drill if not
        if self.lightcones is None:
            self.drill_lightcones()

        # Create lensing lookup table
        if lookup_table is True: lens_table = pangloss.lensing.LensingTable()

        # Initialize new variables (note: e and g have to be initialized as complex for memory allocation)
        assert self.galaxy_count == len(self.lightcones)
        kappa = np.zeros(self.galaxy_count)
        gamma1 = np.zeros(self.galaxy_count)
        gamma2 = np.zeros(self.galaxy_count)
        g = (gamma1 + 1j*gamma2)/(1.0-kappa)
        e1 = np.zeros(self.galaxy_count)
        e2 = np.zeros(self.galaxy_count)
        e = e1+1j*e2
        eMod = np.zeros(self.galaxy_count)
        ePhi = np.zeros(self.galaxy_count)

        # Keep track of how long each lightcone takes to process
        runtimes = np.zeros(self.galaxy_count)

        # Set the counter to be 10%
        counter = np.ceil(len(self.lightcones)/10.0)

        # Calculate lensing in each lightcone:
        for lightcone in self.lightcones:
            if time is True: start_time = timeit.default_timer()
            if lightcone.ID%counter == 0 and vb is True:
                print lightcone.ID,' ',np.ceil(100*lightcone.ID/self.galaxy_count),'%'

            # Remove galaxies in lightcone that do not meet the relevance limit
            lightcone.galaxies = lightcone.galaxies[lightcone.galaxies['relevance'] >= relevance_lim]
            lightcone.galaxy_count = len(lightcone.galaxies)

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

            lightcone.drawConcentrations(errors=False)

            # Compute each halo's contribution to the convergence and shear:
            if lookup_table is True:
                # Use lookup table to speed up kappa calculations
                lightcone.makeKappas(truncationscale=10,lensing_table=lens_table)
            else:
                # Calculate kappa values explicitly
                lightcone.makeKappas(truncationscale=10)

            # Combine all contributions into a single kappa and gamma for the lightcone
            if foreground_corr is True:
                # Implement foreground correction using calcualted mean foreground kappas
                lightcone.combineKappas(methods=methods,smooth_corr=smooth_corr,foreground_kappas=self.foreground_kappas)
            else:
                # No foreground correction - simply add all kappas
                lightcone.combineKappas(methods=methods,smooth_corr=smooth_corr)

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

            if time is True:
                elapsed = timeit.default_timer() - start_time
                runtimes[lightcone.ID] = elapsed

            assert lightcone.galaxy_count == len(lightcone.galaxies)

        if vb and time is True:
            print 'average CPU time per background galaxy: ',np.mean(runtimes),'+/-',np.std(runtimes)

        # Calculate mean/std relevant galaxies per lightcone
        self.relevant_counts = [lightcone.galaxy_count for lightcone in self.lightcones]
        self.mean_relevant_halos = np.mean(self.relevant_counts)
        self.std_relevant_halos = np.std(self.relevant_counts)

        #-------------------------------------------------------------------------------------
        # Use the halo model's kappa and gamma values to compute the new galaxy ellipticities

        # Extract background galaxy intrinsic ellipticites
        e1_int = self.galaxies['e1_int']
        e2_int = self.galaxies['e2_int']

        # Calculate the reduced shear g and its conjugate g_conj
        g = (gamma1 + 1j*gamma2)/(1.0-kappa)
        g_conj = np.array([val.conjugate() for val in g])

        # Flag any galaxy that has been strongly (or near-strongly) lensed
        self.galaxies['strong_flag'] = abs(g)>0.5

        # Calculate the observed ellipticity for weak lensing events
        index = np.abs(g) < 1.0
        e[index] = ( (e1_int[index] + 1j*e2_int[index]) + g[index]) / (1.0+g_conj[index] * (e1_int[index] + 1j*e2_int[index]) )

        # Calculate the observed ellipticity for strong lensing events
        index = ~index
        e1_int_conj = np.array([val.conjugate() for val in e1_int])
        e2_int_conj = np.array([val.conjugate() for val in e2_int])
        e[index] = (1.0 + g[index]*(e1_int_conj[index]+1j*e2_int_conj[index]) ) / ( (e1_int_conj[index]+1j*e2_int_conj[index]) + g_conj[index])

        # Calculate Cartesian and polar components
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

# ----------------------------------------------------------------------------

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

# ----------------------------------------------------------------------------

    def setup_grid(self):
        '''
        Create the distance grid to simplify distance calculations.
        '''

        # Make redshift grid:
        self.grid = pangloss.Grid(self.zl,self.zs,nplanes=self.planes)

        return

    def set_relevance(self,lightcone,metric='curtis'):
        '''
        Give each foreground galaxy in a lightcone an relevance based upon the
        inputted metric.
        '''

        if metric == 'curtis':
            # Set mass and radius thresholds
            M = 10.0**12 # Solar masses
            R = 0.01 # Mpc

            # Calculate relevance according to scheme given by McCully et al. in http://arxiv.org/abs/1601.05417
            lightcone.galaxies['relevance'] = (lightcone.galaxies['Mh']/M)*(R/lightcone.galaxies['rphys'])**3

        elif metric == 'linear':
            # Find the min and max physical distance and
            r_min = np.min(lightcone.galaxies['rphys'])
            r_max = np.max(lightcone.galaxies['rphys'])
            Mh_min = np.min(lightcone.galaxies['Mh'])
            Mh_max = np.max(lightcone.galaxies['Mh'])

            # Compute the relevance using a linear metric from 0 to 1
            relevance_r = (r_max - self.galaxies['rphys']) / r_max
            relevance_m = (self.galaxies['Mh'] - Mh_min) / (Mh_max - Mh_min)
            #relevance = relevance_r*relevance_m
            relevance = np.sqrt(relevance_r**2+relevance_m**2)

            # Set the relevance of each galaxy normalized by the maximum relevance
            lightcone.galaxies['relevance'] = relevance/np.max(relevance)

        # Make sure all relevance values are positive
        assert (lightcone.galaxies['relevance'] > 0).all()

        return

# ----------------------------------------------------------------------------

    def drill_lightcones(self,radius=2.0,foreground=None,save=False,smooth_corr=False):
        '''
        Drill a lightcone at each background source with radius in arcmin. Will
        write the lightcones to pangloss/data/lightcones. The method can only be
        used with BackgroundCatalogs that were created as a field, not a subplot.
        If smooth_corr is True, will pre-calculate the necessary volumes to
        speed up the correction.
        '''

        # Retrieve background galaxy data and initialize the lightcones
        self.lightcones = np.zeros(self.galaxy_count,dtype='object')

        # Set lightcone parameters
        flavor = 'sampled'

        # Setup the grid
        self.setup_grid()

        # If a foreground catalog object is not passed, load in the appropriate catalog.
        if foreground is None:
            # Load in the corresponding foreground catalog
            config = pangloss.Configuration.example()
            foreground = pangloss.ForegroundCatalog(pangloss.GUO_FILE, config)

        # Only take the columns from foreground that are needed for a lightcone
        lc_catalog = foreground.galaxies[['RA','Dec','z_obs','Mhalo_obs','Type','GalID',
                                          'mag_SDSS_i', 'Mh']]
        lc_catalog.rename(columns={'mag_SDSS_i': 'mag'})

        # Save mean kappas in foreground redshift slices for foreground correction
        self.foreground_kappas = foreground.mean_kappas
        del foreground

        # If smooth_corr is True, pre-calculate volumes to be stored in self.grid
        if smooth_corr is True: self.grid.calculate_bin_volumes(lc_radius=radius)

        # Set the counter to be 10%
        counter = np.ceil(self.galaxy_count/10.0)

        # Drill a lightcone at each galaxy location
        for i in range(self.galaxy_count):
            if i % counter == 0 and vb is True:
                print i, ' ', np.ceil(100*i/self.galaxy_count), '%'
            # Set galaxy positions
            ra0 = self.galaxies['RA'][i]
            dec0 = self.galaxies['Dec'][i]
            position = [ra0, dec0]

            # Create the lightcone for background galaxy i
            self.lightcones[i] = pangloss.Lightcone(lc_catalog, flavor, position, radius, ID=i)

            # Create the redshift scaffolding for lightcone i:
            self.lightcones[i].defineSystem(self.zl,self.zs)
            self.lightcones[i].loadGrid(self.grid)
            self.lightcones[i].snapToGrid(self.grid,self.foreground_kappas)

            # Set relevance of each foreground object in the lightcone for lensing
            self.set_relevance(self.lightcones[i])

        if save == True:
            # Write out this entire catalog to the 'data' directory
            self.save()

        return

# ----------------------------------------------------------------------------

    def save(self, filename=None):
        '''
        Save the background galaxy catalog in '/data'.
        '''

        if vb == True:
            print 'Pickling this background catalog to disk...'

        # If no filename is given, a default is set using the field numbers (x,y,i,j).
        if filename is None:
            filename = pangloss.DATA_DIR + '/background_catalog_' + str(self.map_x) + '_' + str(
                self.map_y)+'_'+str(self.field_i)+'_'+str(self.field_j)+'.pkl'

        self.filename = filename

        pickle_file = open(self.filename, 'wb')
        pickle.dump(self.__dict__, self.filename, 2)
        pickle_file.close()

        return

# ----------------------------------------------------------------------------

    def load(self, filename=None):
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

# ----------------------------------------------------------------------------

    def bin_to_maps(self, lensed='none', binsize=0.075, center=None, savefile=None, show=False):
        '''
        Bin the background galaxies into WLMaps. We always make both
        a kappa map and a shear map - the 'lensed' kwarg tells us
        which type of lensing, by 'halo's or 'map's, has been applied.
        Binsize is in units of arcmin, and an optional map centroid is
        passed in as the 'center' list [ra,dec] where these world
        coordinates are in degrees.
        ** set binsize=0.586 for same resolution scale as Hilbert maps **
        '''

        maps=['kappa','gamma1','gamma2']

        setup = False
        for map in maps:
            if setup is False:
                # Set up x and y bins - only need to do this the first
                # time through!
                setup = True
                ra,dec = self.galaxies['RA'],self.galaxies['Dec']
                decmin,decmax = np.min(dec),np.max(dec)
                NY = int( np.rad2deg(decmax-decmin)*60.0/binsize )
                ramin,ramax = np.min(ra),np.max(ra)
                NX = int( np.rad2deg(ramax-ramin)*60.0/(binsize*np.cos(np.deg2rad(np.mean(dec)))) )

                rabins = np.linspace(ramin,ramax,NX+1)
                decbins = np.linspace(decmin,decmax,NY+1)

                # Set up some data arrays:
                empty = np.outer(np.zeros(NY),np.zeros(NX))
                kappadata = np.array([empty])
                gammadata = np.array([empty,empty])

            if map == 'kappa':

                if lensed == 'map':
                    values = self.galaxies['kappa']
                elif lensed == 'halo':
                    values = self.galaxies['kappa_halo']
                elif lensed == 'none':
                    values = np.zeros(len(self.galaxies))

            elif map == 'gamma1':

                if lensed == 'map':
                    values = self.galaxies['gamma1']
                elif lensed == 'halo':
                    values = self.galaxies['gamma1_halo']
                elif lensed == 'none':
                    values = np.zeros(len(self.galaxies))

            elif map == 'gamma2':

                if lensed == 'map':
                    values = self.galaxies['gamma2']
                elif lensed == 'halo':
                    values = self.galaxies['gamma2_halo']
                elif lensed == 'none':
                    values = np.zeros(len(self.galaxies))

            # Make 2D histograms, weighted and unweighted. The first
            # gives the sum of the map contributions, the second the
            # number of map contributions. We want the simple average,
            # so we divide one by the other.
            H,x,y = np.histogram2d(ra,dec,weights=values,bins=[rabins,decbins])
            N,x,y = np.histogram2d(ra,dec,bins=[rabins,decbins])

            if map == 'kappa':
                kappadata[0] = 1.0*H/N
            elif map == 'gamma1':
                gammadata[0] = 1.0*H/N
            elif map == 'gamma2':
                gammadata[1] = 1.0*H/N

        # Store these histograms in WLMap objects
        map_xy = [self.map_x, self.map_y]
        kappamap = pangloss.Kappamap(data=[kappadata,self.domain,map_xy])
        shearmap = pangloss.Shearmap(data=[gammadata,self.domain,map_xy])

        # Testing purposes...
        mean_kappa = np.mean(kappamap.values[0])
        print 'mean kappa of type {} is {}'.format(lensed,mean_kappa)

        # Show plot if `show` is True
        if show is True:
            subplot = np.rad2deg([ramax,ramin,decmin,decmax]) # ra flipped b/c left-handed
            kappamap.plot(subplot=subplot,coords='world')

        if savefile is not None:
            plt.savefig(pangloss.DATA_DIR + '/binned_maps/' + savefile, bbox_inches='tight')

        return kappamap, shearmap

# ----------------------------------------------------------------------------

    def calculate_corr(self,corr_type='gg',min_sep=0.1,max_sep=2.0,sep_units='arcmin',binsize=None,N=15.0,lensed='map',foreground=None):
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

            # Return if treecorr is not installed:
            if len(corr_cat.__dict__) == 0:
                print "treecorr is not installed, skipping correlation function calculation."
                return None

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
                config = pangloss.Configuration.example()
                foreground = pangloss.ForegroundCatalog(pangloss.PANGLOSS_MODULE_DIR + '/data/GGL_los_8_' + str(self.map_x) + '_' + str(self.map_y) + '_' + str(self.field_i) + '_' + str(self.field_j) + '_N_4096_ang_4_Guo_galaxies_on_plane_27_to_63.images.txt', config)

            # Create catalog of the foreground galaxy locations
            corr_cat1 = treecorr.Catalog(ra=foreground.galaxies['RA'], dec=foreground.galaxies['Dec'], ra_units='rad', dec_units='rad')

            # Create catalog of the background galaxy ellipticities
            if lensed == 'map':
                corr_cat2 = treecorr.Catalog(ra=galaxies['RA'], dec=galaxies['Dec'], g1=galaxies['e1'], g2=galaxies['e2'], ra_units='rad', dec_units='rad')
            elif lensed == 'halo':
                corr_cat2 = treecorr.Catalog(ra=galaxies['RA'], dec=galaxies['Dec'], g1=galaxies['e1_halo'], g2=galaxies['e2_halo'], ra_units='rad', dec_units='rad')
            elif lensed == 'none':
                corr_cat2 = treecorr.Catalog(ra=galaxies['RA'], dec=galaxies['Dec'], g1=galaxies['e1_int'], g2=galaxies['e2_int'], ra_units='rad', dec_units='rad')

            # Return if treecorr is not installed:
            if len(corr_cat1.__dict__) == 0:
                print "treecorr is not installed, skipping correlation function calculation."
                return None

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

# ----------------------------------------------------------------------------

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
                var1, var2 = corr1.varxi, corr2.varxi

        # For galaxy-galaxy correlation functions:
        elif corr_type == 'ng':
            # Extract the correlation values
            if corr_comp == 'real':
                y1, y2 = corr1.xi, corr2.xi

            elif corr_comp == 'imag':
                y1, y2 = corr1.xi_im, corr2.xi_im

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

        # Propagate the error in the mean percent error calculation
        std_err = 1.0/np.sqrt(np.sum(w)) * 100.0

        return chi2, n_sigma, mean_err, std_err

# ----------------------------------------------------------------------------

    def halo_mass_log_likelihood(self):
        """
            Compare observed ellipticity to ray traced reduced shear.
        """
        std_int = self.std_int
        e1 = self.galaxies['e1'].as_matrix()
        e2 = self.galaxies['e2'].as_matrix()
        std_e1 = e1.std()
        std_e2 = e2.std()
        std_obs = np.mean([std_e1, std_e2])
        sigma = np.sqrt(std_int**2+std_obs**2)

        # Calculate the (log of) normalization constant
        N = 2.*self.galaxy_count
        logZ = (N/2.)*np.log(2.*np.pi*sigma**2)

        # Calculate chi2
        g = self.galaxies['g_halo'].as_matrix()
        g1, g2 = g.real, g.imag
        chi2 = np.sum((e1-g1)**2/sigma**2) + np.sum((e2-g2)**2/sigma**2)

        # Calculate log-likelihood
        log_likelihood = -logZ + (-0.5) * chi2

        return log_likelihood

    def calculate_log_likelihood(self,lensed='halo'):
        '''
        Calculate the log likelihood of the predicted ellipticities given model parameters.
        Taken from Phil's thesis (page 61): http://www.slac.stanford.edu/~pjm/Site/CV_files/Marshall_PhDthesis.pdf
        '''

	if lensed == 'halo':
	    e1_label = 'e1_halo'
	    e2_label = 'e2_halo'
	    g_label = 'g_halo'
	elif lensed == 'map':
	    e1_label = 'e1'
	    e2_label = 'e2'
	    g_label = 'g'
        # Calculate sigma
        std_int = self.std_int
        std_e1 = np.std(self.galaxies[e1_label])
        std_e2 = np.std(self.galaxies[e2_label])
        std_obs = np.mean([std_e1,std_e2])
        sigma = np.sqrt(std_int**2+std_obs**2)

        # Calculate the (log of) normalization constant
        N = 2.*self.galaxy_count
        logZ = (N/2.)*np.log(2.*np.pi*sigma**2)

        # Calculate chi2
        g = self.galaxies[g_label]
        g1, g2 = g.real, g.imag
        e1, e2 = self.galaxies[e1_label], self.galaxies[e2_label]

        chi2 = ( np.sum( (e1 - g1)**2 / sigma**2 ) + np.sum( (e2 - g2)**2 / sigma**2 ) )

        # Calculate log-likelihood
        log_likelihood = -logZ + (-0.5) * chi2

        return log_likelihood

# ----------------------------------------------------------------------------

    def plot(self,subplot=None,mag_lim=[0,24],mass_lim=[0,10**20],z_lim=[0,1.3857],fig_size=10,graph='scatter',lensed='none'):
        '''
        Make scatter plot of generated galaxies.
        '''

        # Get current figure (or make one if it doesn't exist)
        fig = plt.gcf()

        # If there is a Pangloss map open:
        if fig._label == 'Pangloss Map':
            # Adopt axes from the open Kappamap:
            ax = plt.gca()
            ra = ax.coords['ra']
            dec = ax.coords['dec']

            # If the Kappamap subplot was not passed to this Shearmap:
            if subplot == None:
                # Adopt subplot from the open Kappamap:
                subplot = [ax.axis()[0]+.5, ax.axis()[1]+.5, ax.axis()[2]+.5, ax.axis()[3]+.5]
                coords = 'pixel'

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
                # Extract both the lensed-by-map and lensed-by-halo ellipticities
                eMod_halo = galaxies['eMod_halo']
                ePhi_halo = -galaxies['ePhi_halo']
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
        # OLD; before astropy wcs implementation
        #fig.sca(world)
        #world.set_xlim(subplot[0],subplot[1])
        #world.set_ylim(subplot[2],subplot[3])

        if graph == 'scatter':
            # Scale size of point by the galaxy mass
            s = [math.log(mass[i]) for i in range(0,len(mass))]
            ax.scatter(ra,dec,s,alpha=0.5,edgecolor=None,color='blue',transform='wcs')

        elif graph == 'ellipse':
            # Scale galaxy plot size by its mass?
            # scale = ((np.log10(mass)-9.0)/(12.0-9.0))
            scale = 0.1
            floor = 0.01
            size = 0.01*(scale*(scale > 0) + floor)

            # Plot each galaxy as an ellipse
            for i in range(np.shape(galaxies)[0]):
                if lensed == 'none':
                    # Plot intrinsic ellipticities
                    alpha = 0.25
                    pangloss.plot_ellipse(ra[i],dec[i],size,eMod_int[i],ePhi_int[i],world,'blue',alpha)

                elif lensed == 'map':
                    # Plot lensed-by-map ellipticities
                    alpha = 0.3
                    pangloss.plot_ellipse(ra[i],dec[i],size,eMod[i],ePhi[i],world,'green',alpha)

                elif lensed == 'halo':
                    # Plot lensed-by-halo ellipticities
                    alpha = 0.3
                    pangloss.plot_ellipse(ra[i],dec[i],size,eMod_halo[i],ePhi_halo[i],world,'purple',alpha)

                elif lensed == 'both':
                    # Plot both lensed-by-map and intrinsic ellipticities
                    alpha1 = 0.25
                    alpha2 = 0.3
                    pangloss.plot_ellipse(ra[i],dec[i],size,eMod_halo[i],ePhi_halo[i],world,'blue',alpha1)
                    pangloss.plot_ellipse(ra[i],dec[i],size,eMod[i],ePhi[i],world,'green',alpha2)

                elif lensed == 'all':
                    # Plot both types of lensed and intrinsic ellipticities
                    alpha1 = 0.25
                    alpha2 = 0.3
                    alpha3 = 0.3
                    pangloss.plot_ellipse(ra[i],dec[i],size,eMod_int[i],ePhi_int[i],world,'blue',alpha1)
                    pangloss.plot_ellipse(ra[i],dec[i],size,eMod[i],ePhi[i],world,'green',alpha2)
                    pangloss.plot_ellipse(ra[i],dec[i],size,eMod_halo[i],ePhi_halo[i],world,'purple',alpha3)

        elif graph == 'stick':
            if lensed == 'none':
                # Plot intrinsic ellipticity sticks
                pangloss.plot_sticks(ra,dec,eMod_int,ePhi_int,world,'blue')

            elif lensed == 'map':
                # Plot lensed-by-map ellipticity sticks
                pangloss.plot_sticks(ra,dec,eMod,ePhi,world,'green')

            elif lensed == 'halo':
                # Plot lensed-by-halos ellipticity sticks
                pangloss.plot_sticks(ra,dec,eMod_halo,ePhi_halo,world,'purple')

            elif lensed == 'both':
                # Plot both lensed and intrinsic ellipticity sticks
                pangloss.plot_sticks(ra,dec,eMod_halo,ePhi_halo,world,'purple')
                pangloss.plot_sticks(ra,dec,eMod,ePhi,world,'green')

            elif lensed == 'all':
                # Plot both types of lensed and intrinsic ellipticity sticks
                pangloss.plot_sticks(ra,dec,eMod_int,ePhi_int,world,'blue')
                pangloss.plot_sticks(ra,dec,eMod,ePhi,world,'green')
                pangloss.plot_sticks(ra,dec,eMod_halo,ePhi_halo,world,'purple')

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

    def get_sampled_halo_masses_in_lightcones(self):
        row = dict()
        for cone in self.lightcones:
            for index in xrange(len(cone.galaxies)):
                row[cone.galaxies.iloc[index]['GalID']] = [cone.galaxies.iloc[index]['Mh']]
        return DataFrame(row)

    def get_true_halo_masses_in_lightcones(self):
        row = dict()
        for cone in self.lightcones:
            for index in xrange(len(cone.galaxies)):
                row[cone.galaxies.iloc[index]['GalID']] = [cone.galaxies.iloc[index]['Mhalo_obs']]
        return DataFrame(row)