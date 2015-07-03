import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import os
from pangloss import io
from astropy.table import Table, Column


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
      2015-06-29  Started by Everett (SLAC)
    """
    def __init__(self,filename,config):
        self.filename = filename
        # Structures catalog metadata from configfile and reads in the catalog data
        self.config = config
        self.read(filename,config)

    def __str__(self):
        # Add more!
        return 'General catalog object'

    def read(self,filename,config):
        # Uses astropy.table to read catalog, but with a few specific changes
        self.data = io.readCatalog(filename,config)

    def write(self,output=os.getcwd()):
        # Writes catalog data to current directory unless otherwise specified
        self.data.write(output,format = 'ascii')

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

    def generate(self,fig_size=10,mag_cutoff=[24,0],mass_cutoff=[0,10**20],z_cutoff=[0,1.3857]):
        '''
        Draw generated world-coordinate positions of galaxies in the sky in
        world coordinates The optional input fig_size is in inches and has a
        default value of 10. The other optional inputs are value cutoffs; any
        generated galaxy will have attributes within these values.
        '''
