import numpy as np
import matplotlib.pyplot as plt
import os
from pangloss import *
from astropy.table import Table, Column


# ============================================================================

class ForegroundCatalog(Catalog):
    """
    NAME
        ForegroundCatalog

    PURPOSE
        Store, interrogate and plot a collection of foreground galaxy
        data for a given patch of skys.

    COMMENTS
        Inherits from the base class Catalog in catalog.py

    INITIALISATION
        filename:       A string of the catalog filename (likely .txt)
        config:         A config object containing structure of catalog metadata

    METHODS
        plotForeground: Make a scatterplot of foreground galaxy positions at their
                        respective world coordinates. Only galaxies in the catalog
                        whose attributes are within the optional cutoff arguments
                        are displayed.

    BUGS
        ???

    AUTHORS
      This file is part of the Pangloss project, distributed under the
      GPL v2, by Tom Collett (IoA) and  Phil Marshall (Oxford).
      Please cite: Collett et al 2013, http://arxiv.org/abs/1303.6564

    HISTORY
      2015-07-2  Everett (SLAC)
    """

    def __init__(self,filename,config):
        # Calls the superclass init
        Catalog.__init__(self,filename,config)

        # Catalog attributes
        self.galaxyCount = np.shape(self.data)[0]
        self.maxZ = max(self.data['z_obs'])
        self.minZ = min(self.data['z_obs'])
        self.maxM = max(self.data['Mstar_obs'])
        self.minM = min(self.data['Mstar_obs'])
        self.maxMag = max(self.data['mag'])
        self.minMag = min(self.data['mag'])

        # Find world coordinate limits, used for plotting
        self.ra_max = -np.rad2deg(self.data['nRA'].max())
        self.ra_min = -np.rad2deg(self.data['nRA'].min())
        self.dec_max = np.rad2deg(self.data['Dec'].max())
        self.dec_min = np.rad2deg(self.data['Dec'].min())


    def __str(self):
        # Add more!
        return 'Foreground catalog with {} galaxies, '+ \
               'with redshifts ranging from {} to {}'\
                .format(self.galaxyCount,self.minZ,self.maxZ)

# PJM: shouldn't this be just "plot()"? You may want a number of different plot functions
# with different names, but I feel like the default plot should be more interesting
# than just the Catalog.plot() method. Its OK to overload! :-)

    def plotForeground(self,fig_size=10,mag_cutoff=[0,24],mass_cutoff=[0,10**20],z_cutoff=[0,1.3857]):
        '''
        Plots the positions of galaxies in the foreground catalog in world coordinates.
        The optional input fig_size is in inches and has a default value of 10.
        The other optional inputs are cutoffs with default values, which limit
        the number of galaxies that are to be plotted by the respective attribute.
        '''

        # Retrieve list of galaxy world coordinates and magnitudes with values within cutoffs
        ra = -np.rad2deg(self.data['nRA'][(self.data['mag']>mag_cutoff[0]) & (self.data['mag']<mag_cutoff[1]) \
                                        & (self.data['Mstar_obs']>mass_cutoff[0]) & (self.data['Mstar_obs']<mass_cutoff[1]) \
                                        & (self.data['z_obs']>mag_cutoff[0]) & (self.data['z_obs']<mag_cutoff[1])])

        dec = np.rad2deg(self.data['Dec'][(self.data['mag']>mag_cutoff[0]) & (self.data['mag']<mag_cutoff[1]) \
                                        & (self.data['Mstar_obs']>mass_cutoff[0]) & (self.data['Mstar_obs']<mass_cutoff[1]) \
                                        & (self.data['z_obs']>mag_cutoff[0]) & (self.data['z_obs']<mag_cutoff[1])])

        mags = self.data['mag'][(self.data['mag']>mag_cutoff[0]) & (self.data['mag']<mag_cutoff[1]) \
                                        & (self.data['Mstar_obs']>mass_cutoff[0]) & (self.data['Mstar_obs']<mass_cutoff[1]) \
                                        & (self.data['z_obs']>mag_cutoff[0]) & (self.data['z_obs']<mag_cutoff[1])]
        print('ra:',np.shape(ra),'dec:',np.shape(dec),'mags:',np.shape(mags))
        # Scale size of plotted galaxy by the inverse of its magnitude
        size = [500/i for i in mags]

        # Make a scatter plot of the galaxy locations
        plt.scatter(ra,dec,s=size,color='b',alpha=0.5,edgecolor='none')
        plt.xlabel('Right Ascension / deg')
        plt.ylabel('Declination / deg')
        plt.gca().set_xlim((self.ra_min,self.ra_max))
        plt.gca().set_ylim((self.dec_min,self.dec_max))

        # Set figure size to fig_size
        fig = plt.gcf()
        fig.set_size_inches(fig_size,fig_size)

        # The code below can be incorperated with a kappamap to plot an overlay with ForegroundCatalog
        #pix_nra = [K.world2image(a,0)[0] for a in nra]
        #pix_dec = [K.world2image(0,d)[1] for d in dec]
        #plt.scatter(pix_nra,pix_dec)
