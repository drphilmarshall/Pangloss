import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import os, random, math
from astropy.table import Table, Column

import pangloss

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

    BUGS

    AUTHORS
      This file is part of the Pangloss project, distributed under the
      GPL v2, by Tom Collett (IoA) and  Phil Marshall (Oxford).
      Please cite: Collett et al 2013, http://arxiv.org/abs/1303.6564

    HISTORY
      2015-06-29  Started Everett (SLAC)
    """
    def __init__(self,domain=None,N=10,mag_cut=[24.0,0.0],mass_cut=[10.0**6,10.0**12],z_cut=[0.0,1.3857],e_mod_cut=[0,0.25]):
        self.type = 'background'
        self.generate(domain,N,mag_cut,mass_cut,z_cut,e_mod_cut)
        
        # Calls the superclass initialization for useful catalog attributes
        pangloss.Catalog.__init__(self)
        
        return

    def __str__(self):
        # *!Need to fix with new attributes!*
        return 'Background catalog with {} galaxies, '+ \
               'with redshifts ranging from {} to {}'\
                .format(self.galaxyCount,self.minZ,self.maxZ)

    def write(self,output=os.getcwd()):
        # Writes catalog data to current directory unless otherwise specified
        self.galaxies.write(output,format = 'ascii')
        return

# ----------------------------------------------------------------------------

    def generate(self,domain=None,N=10,mag_cut=[24.0,0.0],mass_cut=[10.0**6,10.0**12],z_cut=[0.0,1.3857],e_mod_cut=[0,0.25]):
        '''
        Draw N-generated world-coordinate positions of galaxies in the sky per 
        square arcminute inside a given domain of the form 
        domain=[ra_init,ra_final,dec_init,dec_final]. The other optional inputs
        are value cutoffs; any generated galaxy will have attributes within these 
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
        
        # Determine area of domain and the number of generated galaxies contained in it
        # (expecting wcs in degrees)
        self.Lx, self.Ly = abs(np.rad2deg(ra_final)-np.rad2deg(ra_init)), abs(np.rad2deg(dec_final)-np.rad2deg(dec_init))
        area = 3600*self.Lx*self.Ly # square arcminutes
        self.galaxy_count = int(N*area) # N galaxies per square arcminute

        # Initialize generated variables
        ra = []
        dec = []
        mag = []
        mass = []
        z = []
        e_mod = []
        e_phi = []

        # Populate the generated variables
        for i in range(0,self.galaxy_count):
            ## NOTE: Not all distributions should be uniform!!!
            ra.append(random.uniform(ra_init,ra_final))
            dec.append(random.uniform(dec_init,dec_final))
            mag.append(random.uniform(mag_cut[0],mag_cut[1]))
            mass.append(random.uniform(mass_cut[0],mass_cut[1]))
            z.append(random.uniform(z_cut[0],z_cut[1]))
            e_mod.append(random.uniform(e_mod_cut[0],e_mod_cut[1]))
            e_phi.append(random.uniform(0,180))

        # Save generated catalog as an astropy table
        self.galaxies = Table([ra,dec,mag,mass,z,e_mod,e_phi],names=['RA','Dec','mag','Mstar_obs','z_obs','e_mod','e_phi'], \
                              meta={'name':'generated catalog','size':N,'mag_cutoff':mag_cut, \
                                    'mass_cutoff':mass_cut,'z_cutoff':z_cut,'e_mod_cutoff':e_mod_cut})

        return
        
    def lens_by_map(self):
        '''
        Lense background galaxies by the shear and convergence in their respective Kappamaps and Shearmaps. 
        '''
        pass
    
    def add_noise(self):
        '''
        Add shape noise to  
        '''
        pass
    
    def plot(self,subplot=None,mag_cutoff=[0,24],mass_cutoff=[0,10**20],z_cutoff=[0,1.3857],fig_size=10):
        '''
        Make scatter plot of generated galaxies.
        '''
        
        # Get current figure (or make one if it doesn't exist)
        fig = plt.gcf()
        
        # If there is a Pangloss map open:
        if fig._label == 'Pangloss Map':
            # Adopt axes from the open Kappamap:
            imshow = fig.axes[0]
            world = fig.axes[1]
            
            # If the Kappamap subplot was not passed to this Shearmap:
            if subplot == None:
                # Adopt subplot from the open Kappamap:
                fig.sca(world)
                subplot = plt.axis()
                
            # Set RA and Dec cutoffs from subplot
            ra_cutoff = [subplot[0], subplot[1]]
            dec_cutoff = [subplot[2], subplot[3]]
            
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
            
            # Set RA and Dec cutoffs from subplot
            ra_cutoff = [subplot[0], subplot[1]]
            dec_cutoff = [subplot[2], subplot[3]]
                
            # Create new imshow and world axes
            imshow, world = pangloss.make_axes(fig,subplot)

        # Find the galaxies that are within the cutoffs
        ra,dec,mass = pangloss.Catalog.findGalaxies(self,mag_cutoff,mass_cutoff,z_cutoff,ra_cutoff,dec_cutoff)
        
        # Set current axis to world coordinates and set the limits
        fig.sca(world)
        world.set_xlim(subplot[0],subplot[1])
        world.set_ylim(subplot[2],subplot[3])
        
        # Scale size of point by the galaxy mass
        s = [math.log(mass[i]) for i in range(0,len(mass))]
        plt.scatter(ra,dec,s,alpha=0.5,edgecolor=None,color='blue')
        plt.xlabel('Right Ascension / deg')
        plt.ylabel('Declination / deg')

        # Set the correct figure size
        pangloss.set_figure_size(fig,fig_size,self.Lx,self.Ly)
        
        return
