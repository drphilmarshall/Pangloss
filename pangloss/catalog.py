import numpy as np
import scipy as sp
import matplotlib as plt
import os
from pangloss import io, config
from astropy.table import Table, Column
import astropy.wcs


# ============================================================================

class Catalog(object):
    """
    NAME
        Catalog

    PURPOSE
        ???

    COMMENTS
        ???

    INITIALISATION
        filename = file name, config = config object
            
    METHODS
        ???

    BUGS
        ???

    AUTHORS
      This file is part of the Pangloss project, distributed under the
      GPL v2, by Tom Collett (IoA) and  Phil Marshall (Oxford). 
      Please cite: Collett et al 2013, http://arxiv.org/abs/1303.6564

    HISTORY
      2015-06-29  Everett (SLAC)
    """
    def __init__(self,filename,config):        
        self.filename = filename        
        # Structures catalog metadata from configfile and reads in the catalog data
        self.config = config
        self.read(filename,config)
    
    def read(self,filename,config):    
        # Uses astropy.table to read catalog, but with a few specific changes
        self.data = io.readCatalog(filename,config)
    
    def write(self,output=os.getcwd()):
        # Writes catalog data to current directory unless otherwise specified
        self.data.write(output,format = 'ascii')
    
    def generate(self):
        '''
        Draw world-coordinates in the sky in RA and DEC
        '''
        # Will probably need to rewrite - I don't know how the catalog data is actually formatted yet
        self.xmax = self.data['nRA'].max()
        self.xmin = self.data['nRA'].min()
        self.ymax = self.data['Dec'].max()
        self.ymin = self.data['Dec'].min()
        
        # Use previous values to calculate necessary axes and plot in rads
        
        
    