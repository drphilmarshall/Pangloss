import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import os, random, math, cmath
from astropy.table import Table, Column
from matplotlib.patches import Ellipse

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
    def __init__(self,domain=None,N=10,mag_lim=[24.0,0.0],mass_lim=[10.0**6,10.0**12],z_lim=[0.0,1.3857],sigma_e=0.2,M=0.8):
        self.type = 'background'
        self.generate(domain,N,mag_lim,mass_lim,z_lim,sigma_e)
        
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

    def generate(self,domain=None,N=10,mag_lim=[24.0,0.0],mass_lim=[10.0**6,10.0**12],z_lim=[0.0,1.3857],sigma_e=0.2):
        '''
        Draw N-generated world-coordinate positions of galaxies in the sky per 
        square arcminute inside a given domain of the form 
        domain=[ra_init,ra_final,dec_init,dec_final]. The other optional inputs
        are value limits; any generated galaxy will have attributes within these 
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
        ra = np.zeros(self.galaxy_count)
        dec = np.zeros(self.galaxy_count)
        mag = np.zeros(self.galaxy_count)
        mass = np.zeros(self.galaxy_count)
        z = np.zeros(self.galaxy_count)
        e1_int = np.zeros(self.galaxy_count)
        e2_int = np.zeros(self.galaxy_count)
        eMod_int = np.zeros(self.galaxy_count)
        ePhi_int = np.zeros(self.galaxy_count)

        # Populate the generated variables
        ra = np.random.uniform(ra_init,ra_final,self.galaxy_count)
        dec = np.random.uniform(dec_init,dec_final,self.galaxy_count)
        mag = np.random.uniform(mag_lim[0],mag_lim[1],self.galaxy_count)
        mass = np.random.uniform(mass_lim[0],mass_lim[1],self.galaxy_count)
        z = np.random.uniform(z_lim[0],z_lim[1],self.galaxy_count)
        e1_int = np.random.normal(0.0,sigma_e,self.galaxy_count)
        e2_int = np.random.normal(0.0,sigma_e,self.galaxy_count)
        
        # Change any |e|> 1 ellipticity components
        while (e1_int>1.0).any() or (e2_int>1.0).any():
                
            for i in [j for j in range(len(e1_int)) if e1_int[j]>1.0]:
                e1_int[i] = np.random.normal(0.0,sigma_e)
                
            for i in [j for j in range(len(e2_int)) if e2_int[j]>1.0]:
                e2_int[i] = np.random.normal(0.0,sigma_e)

        # Calculate Cartesian components of intrinsic complex ellipticity
        e_int = e1_int+1.0j*e2_int
        eMod_int = abs(e_int)
        ePhi_int = np.rad2deg(np.arctan2(e2_int,e1_int))/2.0

        # Save generated catalog as an astropy table
        self.galaxies = Table([ra,dec,mag,mass,z,eMod_int,ePhi_int,e1_int,e2_int],names=['RA','Dec','mag','Mstar_obs','z_obs','eMod_int','ePhi_int','e1_int','e2_int'], \
                              meta={'name':'generated catalog','size':N,'mag_lim':mag_lim, \
                                    'mass_lim':mass_lim,'z_lim':z_lim,'sigma_e':sigma_e})
                                    
        return
        
    def lens_by_map(self,kappamap,shearmap,subplot=None,mag_lim=[0,24],mass_lim=[0,10**20],z_lim=[0,1.3857]):
        '''
        Lense background galaxies by the shear and convergence in their respective Kappamaps and Shearmaps. 
        '''
        
        # Exctract needed data from catalog galaxies
        #galaxies = pangloss.Catalog.return_galaxies(self,mag_lim,mass_lim,z_lim,ra_lim,dec_lim)
        ra = np.rad2deg(self.galaxies['RA'])
        dec = np.rad2deg(self.galaxies['Dec'])
        e1_int = self.galaxies['e1_int']
        e2_int = self.galaxies['e2_int']
        
        # Initialize new variables
        kappa = np.zeros(self.galaxy_count)
        gamma1 = np.zeros(self.galaxy_count)
        gamma2 = np.zeros(self.galaxy_count)
        g = np.zeros(self.galaxy_count)
        e1 = np.zeros(self.galaxy_count)
        e2 = np.zeros(self.galaxy_count)
        eMod = np.zeros(self.galaxy_count)
        ePhi = np.zeros(self.galaxy_count)
        
        # Extract convergence and shear values at each galaxy location from maps
        for i in range(self.galaxy_count):
            kappa[i] = kappamap.at(ra[i],dec[i],mapfile=0)
            gamma1[i] = shearmap.at(ra[i],dec[i],mapfile=0)
            gamma2[i] = shearmap.at(ra[i],dec[i],mapfile=1)
            
        # Calculate the reduced shear g and its conjugate g_conj
        g = (gamma1 + 1j*gamma2)/(1.0-kappa)
        g_conj = np.array([val.conjugate() for val in g])
        
        # Calculate the observed ellipticity
        e = ((e1_int + 1j*e2_int) + g)/(1.0+g_conj * (e1_int + 1j*e2_int))
        e1, e2 = e.real, e.imag
        eMod = np.abs(e)
        ePhi = np.rad2deg(np.arctan2(e2,e1))/2.0
        #ePhi = np.rad2deg([(cmath.phase(val))/2.0 for val in e])
        
        # Add convergence and shear values to catalog
        self.galaxies['kappa'] = kappa
        self.galaxies['gamma1'] = gamma1
        self.galaxies['gamma2'] = gamma2
        self.galaxies['g'] = g
        self.galaxies['e1'] = e1
        self.galaxies['e2'] = e2
        self.galaxies['eMod'] = eMod
        self.galaxies['ePhi'] = ePhi

        return
        
    def add_noise(self,M=1,sigma_obs=0.1):
        '''
        Add measurement and shape noise to the background galaxy intrinsic shapes.
        '''
        
        # Multiplicatibe shear calibration error:        
        # We tend to systematically underestimate the ellipticity of background galaxies.
        # Multiplying by M < 1 accounts for this error.
        self.galaxies['e1'] = M*self.galaxies['e1']
        self.galaxies['e2'] = M*self.galaxies['e2']
        self.galaxies['eMod'] = M*self.galaxies['eMod']
        
        # Measurement noise:
        self.galaxies['e1'] += np.random.normal(0.0,sigma_obs,self.galaxy_count)
        self.galaxies['e2'] += np.random.normal(0.0,sigma_obs,self.galaxy_count)        
        
        return
    
    def plot(self,subplot=None,mag_lim=[0,24],mass_lim=[0,10**20],z_lim=[0,1.3857],fig_size=10,graph='scatter',lensed=False):
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
        
        # The angles are flipped as we are using a left-handed coordinate reference axis for plotting
        if graph == 'ellipse' or graph == 'stick':
            if lensed == False:
                # Extract intrinsic ellipticity
                eMod_int = galaxies['eMod_int']
                ePhi_int = -galaxies['ePhi_int']               
                
            elif lensed == True:
                # Extract lensed ellipticity
                eMod = galaxies['eMod']
                ePhi = -galaxies['ePhi']
                
            elif lensed == 'both':
                # Extract both the intrinsic and lensed ellipticity
                eMod_int = galaxies['eMod_int']
                ePhi_int = -galaxies['ePhi_int']
                eMod = galaxies['eMod']
                ePhi = -galaxies['ePhi']            
        
        # Set current axis to world coordinates and set the limits
        fig.sca(world)
        world.set_xlim(subplot[0],subplot[1])
        world.set_ylim(subplot[2],subplot[3])
        
        if graph == 'scatter':            
            # Scale size of point by the galaxy mass
            s = [math.log(mass[i]) for i in range(0,len(mass))]
            plt.scatter(ra,dec,s,alpha=0.5,edgecolor=None,color='blue')
        
        elif graph == 'ellipse':             
            # Scale galaxy plot size by its mass?
            # scale = ((np.log10(mass)-9.0)/(12.0-9.0))
            scale = 0.5            
            floor = 0.01
            size = 0.01*(scale*(scale > 0) + floor)
        
            # Plot each galaxy as an ellipse
            for i in range(np.shape(galaxies)[0]):
                if lensed == False:
                    # Plot intrinsic ellipticities
                    alpha = 0.25
                    pangloss.plotting.plot_ellipse(ra[i],dec[i],size,eMod_int[i],ePhi_int[i],world,'blue',alpha)
                
                elif lensed == True:
                    # Plot lensed ellipticities
                    alpha = 0.3
                    pangloss.plotting.plot_ellipse(ra[i],dec[i],size,eMod[i],ePhi[i],world,'green',alpha)
                    
                elif lensed == 'both':
                    # Plot both lensed and intrinsic ellipticities
                    alpha1 = 0.25
                    alpha2 = 0.3
                    pangloss.plotting.plot_ellipse(ra[i],dec[i],size,eMod_int[i],ePhi_int[i],world,'blue',alpha1)
                    pangloss.plotting.plot_ellipse(ra[i],dec[i],size,eMod[i],ePhi[i],world,'green',alpha2)
                    
        elif graph == 'stick':
            if lensed == False:
                # Plot intrinsic ellipticity sticks
                pangloss.plotting.plot_sticks(ra,dec,eMod_int,ePhi_int,world,'blue')
                
            elif lensed == True:
                # Plot lensed ellipticity sticks
                pangloss.plotting.plot_sticks(ra,dec,eMod,ePhi,world,'green')
                    
            elif lensed == 'both':
                # Plot both lensed and intrinsic ellipticity sticks
                pangloss.plotting.plot_sticks(ra,dec,eMod_int,ePhi_int,world,'blue')
                pangloss.plotting.plot_sticks(ra,dec,eMod,ePhi,world,'green')

        # Label axes and set the correct figure size
        plt.xlabel('Right Ascension / deg')
        plt.ylabel('Declination / deg')
        pangloss.set_figure_size(fig,fig_size,Lx,Ly)
        
        return
