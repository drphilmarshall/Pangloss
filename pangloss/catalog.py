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
        
# ----------------------------------------------------------------------------

    def generate(self,domain=None,N=1000,mag_cut=[24.0,0.0],mass_cut=[10.0**5,10.0**12],z_cut=[0.0,1.3857],plot=False,fig_size=10):
        '''
        Draw N-generated world-coordinate positions of galaxies in the sky inside
        a given domain of the form domain=[ra_init,ra_final,dec_init,dec_final].
        The other optional inputs are value cutoffs; any generated galaxy will
        have attributes within these values. Will make a scatter plot of the
        generated catalog only if plot = True
        '''
        ## *!Note!* Should delete soon - only used in BackgroundCatalog. 
        if domain == None:
            # Make a default domain (shouldn't be used except for testing or demo purposes)
            ra_init = 1    # initial value is larger as ra is left-handed
            dec_init = -1
            ra_final = -1
            dec_final = 1

        else:
            # Set ra and dec limits from domain. domain = [ra_init,ra_final,dec_init,dec_final]
            ra_init = domain[0]
            ra_final = domain[1]
            dec_init = domain[2]
            dec_final = domain[3]

        # Initialize generated variables
        ra = []
        dec = []
        mag = []
        mass = []
        z = []

        # Populate the generated variables
        for i in range(0,N):
            ## NOTE: Not all distributions should be uniform!!!
            ra.append(random.uniform(ra_init,ra_final))
            dec.append(random.uniform(dec_init,dec_final))
            mag.append(random.uniform(mag_cut[0],mag_cut[1]))
            mass.append(random.uniform(mass_cut[0],mass_cut[1]))
            z.append(random.uniform(z_cut[0],z_cut[1]))

        # Save generated catalog as an astropy table
        self.gen_galaxies = Table([ra,dec,mag,mass,z],names=['ra','dec','mag','mass','z'], \
                              meta={'name':'generated catalog','size':N,'mag_cutoff':mag_cut, \
                                    'mass_cutoff':mass_cut,'z_cutoff':z_cut})

        # Make scatter plot of generated galaxies
        if plot == True:
            # Scale size of point by the galaxy mass
            s = [math.log(mass[i]) for i in range(0,len(mass))]
            plt.scatter(ra,dec,s,alpha=0.5,edgecolor=None,color='r')

            plt.xlabel('Right Ascension / deg')
            plt.ylabel('Declination / deg')
            plt.gca().set_xlim(min(ra),max(ra))
            plt.gca().set_ylim(min(dec),max(dec))

            # Set figure size to fig_size
            fig = plt.gcf()
            fig.set_size_inches(fig_size,fig_size)

        return
