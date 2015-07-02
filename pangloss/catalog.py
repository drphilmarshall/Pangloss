import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import os
from pangloss import io
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
        
        # Find world coordinate limits, used for plotting
        self.nra_max = np.rad2deg(self.data['nRA'].max())
        self.nra_min = np.rad2deg(self.data['nRA'].min())
        self.dec_max = np.rad2deg(self.data['Dec'].max())
        self.dec_min = np.rad2deg(self.data['Dec'].min())
    
    def read(self,filename,config):   
        # Uses astropy.table to read catalog, but with a few specific changes
        self.data = io.readCatalog(filename,config)
    
    def write(self,output=os.getcwd()):
        # Writes catalog data to current directory unless otherwise specified
        self.data.write(output,format = 'ascii')
    
    def generate(self,fig_size=10,mag_cutoff=21.5):
        '''
        Draw world-coordinates in the sky in nRA and DEC. The optional input
        fig_size is in inches and has a default value of 10. The other optional
        input is the magnitude cutoff, which is set to 21.5 by default.
        '''
        
        # Retrieve list of galaxy world coordinates with magnitudes greater than cutoff
        nra = np.rad2deg(self.data['nRA'][self.data['mag']<mag_cutoff])
        dec = np.rad2deg(self.data['Dec'][self.data['mag']<mag_cutoff])
        
        # Scale size of plotted galaxy by the inverse of its magnitude
        mags = self.data['mag'][self.data['mag']<mag_cutoff]
        size = [500/i for i in mags]
        
        # Make a scatter plot of the galaxy locations
        plt.scatter(nra,dec,s=size,color='b',alpha=0.5)
        plt.xlabel('Right Ascension / deg')
        plt.ylabel('Declination / deg')
        plt.gca().set_xlim((self.nra_min,self.nra_max))
        plt.gca().set_ylim((self.dec_min,self.dec_max))
        
        # Set figure size to fig_size
        fig = plt.gcf()
        fig.set_size_inches(fig_size,fig_size)        
    