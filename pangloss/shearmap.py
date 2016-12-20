
import numpy as np
import matplotlib.pyplot as plt
import pangloss
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar

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

    def __init__(self, shearfiles, FITS=True):

        self.name = 'Shear map kappa from Millenium Simulation, zs = 1.3857'
        # Calls the WLMap superclass
        pangloss.WLMap.__init__(self, mapfiles=shearfiles, data=None, FITS=FITS)

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
            #ax = fig.axes[0]
            ax = plt.gca()
            ra = ax.coords['ra']
            dec = ax.coords['dec']

            # If the Kappamap subplot was not passed to this Shearmap:
            if subplot == None:
                # Adopt subplot from the open Kappamap:
                subplot = [ax.axis()[0]+.5, ax.axis()[1]+.5, ax.axis()[2]+.5, ax.axis()[3]+.5]
                coords = 'pixel'

            # Adopt figure size from open Kappamap:
            fig_size = plt.gcf().get_size_inches()[0]

        # Use plot_setup method from the base WLMap class:
        pix_xi,pix_xf,pix_yi,pix_yf,Lx,Ly,pix_Lx,pix_Ly,subplot = self.plot_setup(subplot,coords)

        # If there is not a Pangloss map open:
        if fig._label != 'Pangloss Map':
            # Create figure and axes from scratch:
            fig._label = "Pangloss Map"
            pangloss.set_figure_size(fig,fig_size,Lx,Ly)

            # Set the pixel and wcs axes
            imsubplot = [pix_xi, pix_xf, pix_yi, pix_yf]
            ax = pangloss.set_axes(fig,Lx,Ly,self.hdr[0],imsubplot)

# ----------------------------------------------------------------------------

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
        scale = 1.0
        pix_L = np.mean([pix_Lx,pix_Ly])
        L = np.mean([Lx,Ly])
        dx = scale * mod_gamma * np.cos(phi_gamma) * pix_L
        dy = scale * mod_gamma * np.sin(phi_gamma) * pix_L
        # Plot downsampled 2D arrays of shear sticks in current axes.
        # Pixel sampling rate for plotting of shear maps:
        if pix_Lx >= 40:
            N = np.floor(pix_Lx/40.0)
        else:
            N = 1

        ax.quiver(X[::N,::N],Y[::N,::N],dx[::N,::N],dy[::N,::N],color='r',headwidth=0,pivot='middle',transform=ax.get_transform('world'))

        # Add scale bar
        bar = AnchoredSizeBar(ax.transData,L/10.0,'10% Shear',pad=0.5,loc=3,sep=5,borderpad=0.25,frameon=True)
        bar.size_bar._children[0]._linewidth = 2
        bar.size_bar._children[0]._edgecolor = (1,0,0,1)
        ax.add_artist(bar)

        return

    @staticmethod
    def example():
        return Shearmap([pangloss.GAMMA_1_FILE, pangloss.GAMMA_2_FILE], FITS=False)