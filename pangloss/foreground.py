import numpy as np
import matplotlib.pyplot as plt
import os, sys
from astropy.table import Table, Column

# Fast correlation functions:
try:
    import treecorr
except ImportError:
    from pangloss import nocorr as treecorr

# Pangloss:
PANGLOSS_DIR = os.path.expandvars("$PANGLOSS_DIR")
sys.path.append(PANGLOSS_DIR)
import pangloss

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

    def __init__(self,filename,config):

        self.filename = filename
        self.type = 'foreground'

        # Structures catalog metadata from configfile and reads in the catalog data
        self.config = config
        self.read(filename,config)

        # Parsing the file name
        # 0 <= x,y <= 7, each (i,j) map covers 4x4 square degrees
        input_parse = self.filename.split('_') # Creates list of filename elements separated by '_'
        self.map_x = eval(input_parse[3]) # The x location of the map grid
        self.map_y = eval(input_parse[4]) # The y location of the map grid

        # 0 <= i,j <= 3, each (i,j) field covers 1x1 square degree
        self.field_i = eval(input_parse[5]) # The i location of the field grid in the (x,y) map
        self.field_j = eval(input_parse[6]) # The j location of the field grid in the (x,y) map

        # Initialize common catalog attributes:
        pangloss.Catalog.__init__(self)

        return

    def __str(self):
        # Add more!
        return 'Foreground catalog with {} galaxies, with redshifts ranging from {} to {}'.format(self.galaxyCount,self.minZ,self.maxZ)

    def read(self,filename,config):
        '''
        Uses astropy.table to read catalog, but with a few specific changes only for
        Hilbert foreground catalogs.
        '''

        self.galaxies = pangloss.io.read_hilbert_catalog(filename,config)

        return

    def plot(self,subplot=None,mag_lim=[0,24],mass_lim=[0,10**20],z_lim=[0,1.3857],fig_size=10):
        '''
        Plots the positions of galaxies in the foreground catalog in world coordinates.
            The optional input fig_size is in inches and has a default value of 10.
        The other optional inputs are limits with default values, which limit
        the number of galaxies that are to be plotted by the respective attribute.
        '''

        # Get current figure (or make one if it doesn't exist)
        fig = plt.gcf()

        # If there is a Pangloss map open:
        if fig._label == 'Pangloss Map':
            # Adopt axes from the open Kappamap:
            imshow = fig.axes[0]
            world = fig.axes[1]

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

            # Create new imshow and world axes
            imshow, world = pangloss.make_axes(fig,subplot)

        ai, af = subplot[0], subplot[1]    # RA limits for subplot
        di, df = subplot[2], subplot[3]    # DEC limits for subplot
        Lx, Ly = abs(ai-af), abs(di-df)    # Length of axes in wcs

        # Find the galaxies that are within the limits, and extract the useful data from them
        ra_lim, dec_lim = [ai, af], [di, df]
        galaxies = self.return_galaxies(mag_lim,mass_lim,z_lim,ra_lim,dec_lim)
        ra = np.rad2deg(galaxies['RA'])
        dec = np.rad2deg(galaxies['Dec'])
        mass = galaxies['Mstar_obs']

        # Set current axis to world coordinates and set the limits
        fig.sca(world)
        world.set_xlim(subplot[0],subplot[1])
        world.set_ylim(subplot[2],subplot[3])

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
