import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import os, random, math
from astropy.table import Table, Column

import pangloss

# ============================================================================

class Catalog(object):
    """
    NAME
        Catalog

    PURPOSE
        Read in or generate an astronomical catalog, plot, write out.

    COMMENTS
        This is a base class. Intended subclasses include ForegroundCatalog
        (for deflector halos/galaxies) and BackgroundCatalog (to contain
        weakly lensed source galaxies).

    INITIALISATION
        ??

    METHODS
        findGalaxies: Find all galaxies in the catalog that are within the inputted
                      magnitude, mass, redhsift, and coordinate cutoff ranges,
                      and then return the locations and masses of each galaxy.

        returnGalaxies: Same as the findGalaxies() method, but returns all columns
                        of the catalog for galaxies that satisfy the cutoffs.

    BUGS

    AUTHORS
      This file is part of the Pangloss project, distributed under the
      GPL v2, by Tom Collett (IoA) and  Phil Marshall (Oxford).
      Please cite: Collett et al 2013, http://arxiv.org/abs/1303.6564

    HISTORY
      2015-06-29  Started Everett (SLAC)
    """
    def __init__(self):
        # Common catalog attributes
        self.galaxy_count = np.shape(self.galaxies)[0]
        self.maxZ = max(self.galaxies['z_obs'])
        self.minZ = min(self.galaxies['z_obs'])
        self.maxM = max(self.galaxies['Mstar_obs'])
        self.minM = min(self.galaxies['Mstar_obs'])
        self.maxMag = max(self.galaxies['mag'])
        self.minMag = min(self.galaxies['mag'])

        # Find world coordinate limits, used for plotting
        self.ra_max = np.rad2deg(self.galaxies['RA'].max())
        self.ra_min = np.rad2deg(self.galaxies['RA'].min())
        self.dec_max = np.rad2deg(self.galaxies['Dec'].max())
        self.dec_min = np.rad2deg(self.galaxies['Dec'].min())
        
        return
    
    def __str__(self):
        # Add more!
        return 'General catalog object'

    def write(self,output=os.getcwd()):
        # Writes catalog data to current directory unless otherwise specified
        self.galaxies.write(output,format = 'ascii')
        return
        
    
    def findGalaxies(self,mag_cutoff=[0,24],mass_cutoff=[0,10**20],z_cutoff=[0,1.3857],ra_cutoff=None,dec_cutoff=None):
        '''
        Retrieve list of galaxy world coordinates and their masses with values within inputted cutoffs.
        '''

        # If no ra or dec cutoff are given, use all galaxies
        if ra_cutoff == None: ra_cutoff = [self.ra_max, self.ra_min] # RA flipped because RA is left-handed
        if dec_cutoff == None: dec_cutoff = [self.dec_min, self.dec_max]

        # Convert world coordinate limits to radians
        ra_cutoff, dec_cutoff = np.deg2rad(ra_cutoff), np.deg2rad(dec_cutoff)

        # Select only the ra, dec, and mass values from galaxies that satisfy all cutoffs
        ra = np.rad2deg(self.galaxies['RA'][(self.galaxies['mag']>mag_cutoff[0]) & (self.galaxies['mag']<mag_cutoff[1]) \
                                        & (self.galaxies['Mstar_obs']>mass_cutoff[0]) & (self.galaxies['Mstar_obs']<mass_cutoff[1]) \
                                        & (self.galaxies['z_obs']>mag_cutoff[0]) & (self.galaxies['z_obs']<mag_cutoff[1]) \
                                        & (self.galaxies['RA']>ra_cutoff[1]) & (self.galaxies['RA']<ra_cutoff[0]) \
                                        & (self.galaxies['Dec']>dec_cutoff[0]) & (self.galaxies['Dec']<dec_cutoff[1])])

        dec = np.rad2deg(self.galaxies['Dec'][(self.galaxies['mag']>mag_cutoff[0]) & (self.galaxies['mag']<mag_cutoff[1]) \
                                        & (self.galaxies['Mstar_obs']>mass_cutoff[0]) & (self.galaxies['Mstar_obs']<mass_cutoff[1]) \
                                        & (self.galaxies['z_obs']>mag_cutoff[0]) & (self.galaxies['z_obs']<mag_cutoff[1]) \
                                        & (self.galaxies['RA']>ra_cutoff[1]) & (self.galaxies['RA']<ra_cutoff[0]) \
                                        & (self.galaxies['Dec']>dec_cutoff[0]) & (self.galaxies['Dec']<dec_cutoff[1])])

        mass = self.galaxies['Mstar_obs'][(self.galaxies['mag']>mag_cutoff[0]) & (self.galaxies['mag']<mag_cutoff[1]) \
                                        & (self.galaxies['Mstar_obs']>mass_cutoff[0]) & (self.galaxies['Mstar_obs']<mass_cutoff[1]) \
                                        & (self.galaxies['z_obs']>mag_cutoff[0]) & (self.galaxies['z_obs']<mag_cutoff[1]) \
                                        & (self.galaxies['RA']>ra_cutoff[1]) & (self.galaxies['RA']<ra_cutoff[0]) \
                                        & (self.galaxies['Dec']>dec_cutoff[0]) & (self.galaxies['Dec']<dec_cutoff[1])]

        return ra, dec, mass


    def returnGalaxies(self,mag_cutoff=[0,24],mass_cutoff=[0,10**20],z_cutoff=[0,1.3857],ra_cutoff=None,dec_cutoff=None):
        '''
        Return catalog of galaxies that satisfy the inputted cutoffs.
        '''

        # If no ra or dec cutoff are given, use all galaxies
        if ra_cutoff == None: ra_cutoff = [self.ra_max, self.ra_min] # RA flipped because RA is left-handed
        if dec_cutoff == None: dec_cutoff = [self.dec_min, self.dec_max]

        # Convert world coordinate limits to radians
        ra_cutoff, dec_cutoff = np.deg2rad(ra_cutoff), np.deg2rad(dec_cutoff)

        # Select only galaxies that meet the cutoff criteria
        galaxies = self.galaxies[(self.galaxies['mag']>mag_cutoff[0]) & (self.galaxies['mag']<mag_cutoff[1]) \
                                        & (self.galaxies['Mstar_obs']>mass_cutoff[0]) & (self.galaxies['Mstar_obs']<mass_cutoff[1]) \
                                        & (self.galaxies['z_obs']>mag_cutoff[0]) & (self.galaxies['z_obs']<mag_cutoff[1]) \
                                        & (self.galaxies['RA']>ra_cutoff[1]) & (self.galaxies['RA']<ra_cutoff[0]) \
                                        & (self.galaxies['Dec']>dec_cutoff[0]) & (self.galaxies['Dec']<dec_cutoff[1])]

        return galaxies
