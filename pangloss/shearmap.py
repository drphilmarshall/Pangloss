
import numpy as np
import matplotlib.pyplot as plt
from wlmap import WLMap
import cmath

arcmin2rad = (1.0/60.0)*np.pi/180.0
rad2arcmin = 1.0/arcmin2rad
deg2rad = np.pi/180.0
rad2deg = 1.0/deg2rad

vb = False

# ============================================================================

class Shearmap(WLMap):
    """
    NAME
        Shearmap

    PURPOSE
        Read in, store, transform and interrogate a shear map.

    COMMENTS
        A "physical" coordinate system is used, where x = -RA (rad)
        and y = Dec (rad). This is the system favoured by Hilbert et al.

    INITIALISATION
        shearfile      List of files containing a shear map
        FITS           Data file format (def=True)

    METHODS
        plot(self,fig_size=10,subplot=None): Plots either the whole image or a
                                             given sub-image in physical
                                             coordinates

    BUGS
        -Subplots that have the y-axis significantly larger than the x-axis have
        issues with the sticks scaling correctly. Need to look into np.quiver()
        for more options.
        -Alignment of shear sticks is not correct!

    AUTHORS
      This file is part of the Pangloss project, distributed under the
      GPL v2, by Tom Collett (IoA) and  Phil Marshall (Oxford).
      Please cite: Collett et al 2013, http://arxiv.org/abs/1303.6564

    HISTORY
      2015-06-25  Everett (SLAC)
    """

# ----------------------------------------------------------------------------

    def __init__(self,shearfile,FITS=True):

        self.name = 'Shear map kappa from Millenium Simulation, zs = 1.6'
        # Calls the WLMap superclass
        WLMap.__init__(self,shearfile,FITS)

# ----------------------------------------------------------------------------

    def __str__(self):
        ## Add more information!!
        return 'Shear map'

# ----------------------------------------------------------------------------
# Plot the convergence as grayscale:

    def plot(self,fig_size=10,subplot=None,coords='pixel'): # fig_size in inches
        """
        Plot the convergence as a grayscale image.

        Optional arguments:
            fig_size        Figure size in inches
            subplot         List of four plot limits [xmin,xmax,ymin,ymax]
            coords          Type of coordinates inputted for the subplot:
                            'pixel', 'physical', or 'world'
        """
        
        # Use plotting method from WLMap class to calculate values common to both Kappamaps and Shearmaps
        xi,xf,yi,yf,Lx,Ly,xlocs,xlabels,ylocs,ylabels = WLMap.plot(self,fig_size,subplot,coords)

        # Retrieve gamma values in desired subplot
        gamma1 = self.values[0][yi:yf,xi:xf]
        gamma2 = self.values[1][yi:yf,xi:xf]
        
        # Pixel sampling rate for plotting of shear maps

        if Lx >= 40:
            N = np.floor(Lx/40.0)
        else:
            N = 1

        # Set limits and number of points in grid
        X,Y = np.meshgrid(np.arange(0,Lx),np.arange(0,Ly))
        X = X[::N,::N]
        Y = Y[::N,::N]

        # Calculate the modulus and angle of each shear
        mod_gamma = np.sqrt(gamma1*gamma1 + gamma2*gamma2)
        phi_gamma = np.arctan2(gamma2,gamma1)/2.0
        #mod_gamma = cmath.
        #phi_gamma = cmath.phase(gamma1+1j*gamma2)

        # Create the vector components of the shear sticks.
        # We *think* that this transpose is necessary because
        # we are working with maps read in from binary data that
        # have not been transposed. The test is to overlay
        # the foreground galaxy catalogs and see whether the clusters
        # of galaxies line up with the overdensities in kappa/shear...
        dx = mod_gamma * np.sin(phi_gamma)
        dy = mod_gamma * np.cos(phi_gamma)

        # Plot image
        plt.quiver(X,Y,dx[::N,::N],dy[::N,::N],color='r',headwidth=0,pivot='middle')
        plt.title('Shear map of '+self.input[0])

        # Label axes
        plt.xticks(xlocs,xlabels)
        plt.yticks(ylocs,ylabels)
        plt.xlabel('RA (rad)')
        plt.ylabel('DEC (rad)')

        # Set figure size
        fig = plt.gcf()
        if Lx == Ly:
            fig.set_size_inches(fig_size,fig_size)
        elif Lx > Ly:
            fig.set_size_inches(fig_size,fig_size*(1.0*Ly/Lx))
        else:
            fig.set_size_inches(fig_size*(1.0*Lx/Ly),fig_size)

        # Ensures the image is not distorted
        plt.axis('equal')
