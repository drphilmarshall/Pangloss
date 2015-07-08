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
        filename:       A string of the catalog filename (likely .txt)
        config:         A config object containing structure of catalog metadata

    METHODS
        read(filename,config)
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
    def __init__(self,filename,config):
        self.filename = filename
        # Structures catalog metadata from configfile and reads in the catalog data
        self.config = config
        self.read(filename,config)
        return

    def __str__(self):
        # Add more!
        return 'General catalog object'

    def read(self,filename,config):
        # Uses astropy.table to read catalog, but with a few specific changes
        self.data = pangloss.readCatalog(filename,config)
        return

    def write(self,output=os.getcwd()):
        # Writes catalog data to current directory unless otherwise specified
        self.data.write(output,format = 'ascii')
        return

# ----------------------------------------------------------------------------
# Conversions -- This is also in wlmap. We should create a separate file and
# call these functions from both catalog and wlmap.

    # Only approximate WCS transformations - assumes dec=0.0 and small field
    def image2world(self,i,j,mapfile=0):
        a = self.wcs[mapfile]['CRVAL1'] + self.wcs[mapfile]['CD1_1']*(i - self.wcs[mapfile]['CRPIX1'])
        #if a < 0.0: a += 360.0 :We are using nRA instead now
        d = self.wcs[mapfile]['CRVAL2'] + self.wcs[mapfile]['CD2_2']*(j - self.wcs[mapfile]['CRPIX2'])
        return a,d

    def world2image(self,a,d,mapfile=0):
        i = (a - self.wcs[mapfile]['CRVAL1'])/self.wcs[mapfile]['CD1_1'] + self.wcs[mapfile]['CRPIX1']
        # if a negative pixel is returned for i, reinput a as a negative degree
        if i<0:
            a-=360
            i = (a - self.wcs[mapfile]['CRVAL1'])/self.wcs[mapfile]['CD1_1'] + self.wcs[mapfile]['CRPIX1']
        j = (d - self.wcs[mapfile]['CRVAL2'])/self.wcs[mapfile]['CD2_2'] + self.wcs[mapfile]['CRPIX2']
        return i,j

# ----------------------------------------------------------------------------

    def generate(self,domain=None,N=1000,mag_cut=[24.0,0.0],mass_cut=[10.0**5,10.0**12],z_cut=[0.0,1.3857],plot=False,fig_size=10):
        '''
        Draw N-generated world-coordinate positions of galaxies in the sky inside
        a given domain of the form domain=[ra_init,ra_final,dec_init,dec_final].
        The other optional inputs are value cutoffs; any generated galaxy will
        have attributes within these values. Will make a scatter plot of the
        generated catalog only if plot = True
        '''

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
        self.gen_data = Table([ra,dec,mag,mass,z],names=['ra','dec','mag','mass','z'], \
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
