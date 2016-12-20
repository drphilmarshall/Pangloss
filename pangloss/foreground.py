import matplotlib.pyplot as plt
import numpy as np
from utils import int_or_float
import pangloss

# Fast correlation functions:
try:
    import treecorr
except ImportError:
    from pangloss import nocorr as treecorr


# ============================================================================

class ForegroundCatalog(pangloss.Catalog):
    """
    NAME
        ForegroundCatalog

    PURPOSE
        Store, interrogate and plot a collection of foreground galaxy
        data for a given patch of sky.

    COMMENTS
        Inherits from the base class Catalog in catalog.py

    INITIALISATION
        filename:       A string of the catalog filename (likely .txt)
        config:         A config object containing structure of catalog metagalaxies

    METHODS
        read(filename,config): Read in a Hilbert foreground catalog.

        plotForeground: Make a scatterplot of foreground galaxy positions at their
                        respective world coordinates. Only galaxies in the catalog
                        whose attributes are within the optional magnitude, mass,
                        redshift, and coordiante limit arguments are displayed.

    BUGS
        None

    AUTHORS
      This file is part of the Pangloss project, distributed under the
      GPL v2, by Tom Collett (IoA) and  Phil Marshall (Oxford).
      Please cite: Collett et al 2013, http://arxiv.org/abs/1303.6564

    HISTORY
      2015-07-2  Started Everett (SLAC)
    """

    def __init__(self, filename, config, stellar_mass_threshold=None):

        self.filename = filename
        self.type = 'foreground'

        # Structures catalog metadata from configfile and reads in the catalog data
        self.config = config
        self.galaxies = pangloss.read_hilbert_catalog(filename, config)

        if stellar_mass_threshold:
            self.galaxies = self.galaxies[self.galaxies['Mstar_obs'] > stellar_mass_threshold]

        # Add columns to contain free parameters - initialize to true massses:
        self.galaxies['Mh'] = self.galaxies['Mhalo_obs'].copy() * 1.0

        # Parsing the file name
        # 0 <= x,y <= 7, each (i,j) map covers 4x4 square degrees
        input_parse = self.filename.split('_') # Creates list of filename elements separated by '_'
        self.map_x = int_or_float(input_parse[3]) # The x location of the map grid
        self.map_y = int_or_float(input_parse[4]) # The y location of the map grid

        # 0 <= i,j <= 3, each (i,j) field covers 1x1 square degree
        self.field_i = int_or_float(input_parse[5]) # The i location of the field grid in the (x,y) map
        self.field_j = int_or_float(input_parse[6]) # The j location of the field grid in the (x,y) map

        # Calculate mean kappas on the redshift grid
        self.find_mean_kappas()

        # Initialize common catalog attributes:
        pangloss.Catalog.__init__(self)

        self.mass_prior = None

        return

    def __str(self):
        # Add more!
        return 'Foreground catalog with {} galaxies, with redshifts ranging from {} to {}'.format(self.galaxyCount,self.minZ,self.maxZ)

    def snap_to_grid(self,catalog, Grid):
        '''
        Snaps foreground galaxies to redshift grid to aid mean kappa calculation
        NOTE: Different from lightcone snapToGrid
        '''
        z = catalog['z_obs']
        sz, p = Grid.snap(z)
        self.galaxies['z_sz'] = sz
        self.galaxies['sigma_crit'] = Grid.sigma_crit[p] # Units of M_sun/Mpc^2

    def find_mean_kappas(self):
        '''
        Snaps all foreground galaxies to a redshift grid and calculates the mean kappa for each redshfit slice.
        '''

        # Only take the columns from foreground that are needed to calculate densities
        #foreground_galaxies = Table(names=['z_obs','Mhalo_obs'],data=[self.galaxies['z_obs'],self.galaxies['Mhalo_obs']])

        # Set source and strong-lens redshifts
        zl = 0       # There is no strong-lens present
        zs = 1.3857  # All source galaxies are at redshift 1.3857
        planes = 100

        # Setup the grid
        grid = pangloss.Grid(zl,zs,nplanes=planes)
        redshifts = grid.redshifts
        D = pangloss.distances.Distance()

        # Snap all foreground galaxies to grid
        self.snap_to_grid(self.galaxies,grid)

        slice_sigma = np.zeros(len(redshifts))
        slice_sigma_crit = np.zeros(len(redshifts))
        mean_kappa = np.zeros(len(redshifts))

        # Loop over redshift slices and calculate the mass density and critical density
        for i in range(len(redshifts)):

            # Find foreground galaxies in this slice and sum up the masses
            slice_galaxies = self.galaxies[self.galaxies['z_sz']==redshifts[i]]
            slice_mass = 1.0*np.sum(self.galaxies['Mhalo_obs'])

            # Set the scale factor and physical area for this redshift slice
            a = 1.0 / (1 + redshifts[i])
            dtheta = np.deg2rad(1.0) # Foreground catalog is 1'x1'
            area = ( a * D.comoving_transverse_distance(redshifts[i]) * dtheta )**2 # Units of physical Mpc^2

            slice_sigma[i] = slice_mass / area # M_sun/Mpc^2
            slice_sigma_crit[i] = grid.sigma_crit[i] # M_sun/Mpc^2

        # Calculate mean kappa for each redshift slice
        self.mean_kappas = 1.0*slice_sigma / slice_sigma_crit
        self.redshifts = redshifts

        return

    def plot_mean_kappas(self):
        # Set up kappa axis
        fig, ax_kappa = plt.subplots()
        plt.xlabel('Redshift (z)',fontsize='14')
        plt.xticks(fontsize=14)

        # Plot mean kappa
        plt.plot(self.redshifts,self.mean_kappas,'o',zorder=10)
        plt.yticks(fontsize=14)
        plt.yscale('log')
        ax_kappa.set_ylabel('Mean Kappa', color='b',fontsize='14')

        # Set up histogram axis
        ax_hist = ax_kappa.twinx()

        # Plot histogram
        plt.hist(self.galaxies['z_sz'],50,normed=False,facecolor='green', alpha=0.3,zorder=1)
        plt.yscale('log')
        ax_hist.set_ylabel('Galaxy Count', color='g',fontsize='14')
        plt.yticks(fontsize=14)

        # Match axis label colors to plot colors
        for tl in ax_kappa.get_yticklabels(): tl.set_color('b')
        for tl in ax_hist.get_yticklabels():  tl.set_color('g')

        plt.show()

    def plot(self, subplot=None, mag_lim=[0, 24], mass_lim=[0, 10**20],
             z_lim=[0, 1.3857], fig_size=10):
        '''
        Plots the positions of galaxies in the foreground catalog in world coordinates.
        The optional input fig_size is in inches and has a default value of 10.
        The other optional inputs are limits with default values, which limit the number of galaxies that are to be plotted by the respective attribute.
        '''

        # Get current figure (or make one if it doesn't exist)
        fig = plt.gcf()

        # If there is a Pangloss map open:
        if fig._label == 'Pangloss Map':
            # Adopt axes from the open Kappamap:
            # imshow = fig.axes[0]
            # world = fig.axes[1]
            world = fig.gca()
            # This fails though: the halos do not appear...

            # If the Kappamap subplot was not passed to this catalog:
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

            # DEPRECATED: Create new imshow and world axes
            # imshow, world = pangloss.make_axes(fig,subplot)

        ai, af = subplot[0], subplot[1]    # RA limits for subplot
        di, df = subplot[2], subplot[3]    # DEC limits for subplot
        Lx, Ly = abs(ai-af), abs(di-df)    # Length of axes in wcs

        # Find the galaxies that are within the limits, and extract the useful data from them
        ra_lim, dec_lim = [ai, af], [di, df]
        galaxies = self.return_galaxies(mag_lim, mass_lim, z_lim, ra_lim, dec_lim)
        ra = np.rad2deg(galaxies['RA'])
        dec = np.rad2deg(galaxies['Dec'])
        mass = galaxies['Mstar_obs']

        # Set current axis to world coordinates and set the limits
        fig.sca(world)
        world.set_xlim(subplot[0], subplot[1])
        world.set_ylim(subplot[2], subplot[3])

        # Scale galaxy plot size by its mass
        scale = ((np.log10(mass)-9.0)/(12.0-9.0))
        floor = 0.01
        size = 1000.0*(scale*(scale > 0) + floor)

        # Make a scatter plot of the galaxy locations
        plt.scatter(ra,dec,s=size,color='orange',alpha=0.2,edgecolor=None)
        plt.xlabel('Right Ascension / deg')
        plt.ylabel('Declination / deg')

        # Set the correct figure size
        pangloss.set_figure_size(fig,fig_size,Lx,Ly)

        return


    def set_mass_prior(self, distribution):
        self.mass_prior = distribution

    def draw_halo_masses(self):
        self.galaxies['Mh'] = self.mass_prior.draw(self.galaxies)

    @staticmethod
    def guo(stellar_mass_threshold=None):
        return ForegroundCatalog(pangloss.GUO_FILE, pangloss.Configuration.example(), stellar_mass_threshold)