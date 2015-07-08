
import numpy as np
import matplotlib.pyplot as plt
import pangloss
import cmath

arcmin2rad = (1.0/60.0)*np.pi/180.0
rad2arcmin = 1.0/arcmin2rad
deg2rad = np.pi/180.0
rad2deg = 1.0/deg2rad

vb = False

# ============================================================================

class Shearmap(pangloss.WLMap):
    """
    NAME
        Shearmap

    PURPOSE
        Read in, store, transform and interrogate a shear map.

    COMMENTS
        A "physical" coordinate system is used, where x = -RA (rad)
        and y = Dec (rad). This is the system favoured by Hilbert et al.

    INITIALISATION
        shearfiles     List of files containing a shear map
        FITS           Data file format (def=True)

    METHODS
        plot(self,fig_size=10,subplot=None): Plots either the whole image or a
                                             given sub-image in physical
                                             coordinates

    BUGS
        -Subplots that have the y-axis significantly larger than the x-axis have
        issues with the sticks scaling correctly. Need to look into np.quiver()
        for more options.

    AUTHORS
      This file is part of the Pangloss project, distributed under the
      GPL v2, by Tom Collett (IoA) and  Phil Marshall (Oxford).
      Please cite: Collett et al 2013, http://arxiv.org/abs/1303.6564

    HISTORY
      2015-06-25  Started Everett (SLAC)
    """

# ----------------------------------------------------------------------------

    def __init__(self,shearfiles,FITS=True):

        self.name = 'Shear map kappa from Millenium Simulation, zs = 1.3857'
        # Calls the WLMap superclass
        pangloss.WLMap.__init__(self,shearfiles,FITS)

# ----------------------------------------------------------------------------

    def __str__(self):
        ## Add more information!!
        return 'Shear map'

# ----------------------------------------------------------------------------
# Plot the convergence as grayscale:

    def plot(self,fig_size=10,subplot=None,coords='world'): # fig_size in inches
        """
        Plot the shear field with shear sticks.

        Optional arguments:
            fig_size        Figure size in inches
            subplot         List of four plot limits [xmin,xmax,ymin,ymax]
            coords          Type of coordinates inputted for the subplot:
                            'pixel', 'physical', or 'world'
        """
        
        # Get current figure and image axes (or make them if they don't exist)
        fig = plt.gcf()
        
# ----------------------------------------------------------------------------
# Note: the following is slightly inelegant. It would be nice to have the following
# be an if-else, but the default subplot is no longer 'None' after calling 
# plot_setupt(), so the check must be done before calling the method. Try to
# fix later.
        
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

        # Use plot_setup method from the base WLMap class:
        pix_xi,pix_xf,pix_yi,pix_yf,Lx,Ly,pix_Lx,pix_Ly,subplot = self.plot_setup(subplot,coords)
        
        # If there is not a Pangloss map open:
        if fig._label != 'Pangloss Map':
            # Create figure and axes from scratch:
            fig._label = "Pangloss Map"
            pangloss.set_figure_size(fig,fig_size,Lx,Ly)            
            imsubplot = [-0.5,pix_Lx-0.5,-0.5,pix_Ly-0.5]
            imshow,world = pangloss.make_axes(fig,subplot,imsubplot)

# ----------------------------------------------------------------------------

        # Set the axes to image, since we'll be computing
        # shear stick positions in world coordinates:
        fig.sca(world)

        # Retrieve gamma values in desired subplot
        gamma1 = self.values[0][pix_yi:pix_yf,pix_xi:pix_xf]
        gamma2 = self.values[1][pix_yi:pix_yf,pix_xi:pix_xf]

        # Create arrays of shear stick positions, one per pixel in world coordinates
        X,Y = np.meshgrid(np.arange(subplot[0],subplot[1],-self.PIXSCALE[0]),np.arange(subplot[2],subplot[3],self.PIXSCALE[0]))

        # Calculate the modulus and angle of each shear
        mod_gamma = np.sqrt(gamma1*gamma1 + gamma2*gamma2)
        phi_gamma = np.arctan2(gamma2,gamma1)/2.0

        # Sticks in world coords need x reversed, to account for left-handed 
        # system:
        dx = mod_gamma * np.cos(phi_gamma)
        dy = mod_gamma * np.sin(phi_gamma)

        # Plot downsampled 2D arrays of shear sticks in current axes.
        # Pixel sampling rate for plotting of shear maps:
        if pix_Lx >= 40:
            N = np.floor(pix_Lx/40.0)
        else:
            N = 1
            
        plt.quiver(X[::N,::N],Y[::N,::N],dx[::N,::N],dy[::N,::N],color='r',headwidth=0,pivot='middle')

        return
