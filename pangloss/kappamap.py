
import numpy as np
import matplotlib.pyplot as plt

import pangloss

vb = False

# ============================================================================

class Kappamap(pangloss.WLMap):
    """
    NAME
        Kappamap

    PURPOSE
        Read in, store, transform and interrogate a convergence map.

    COMMENTS
        A "physical" coordinate system is used, where x = -RA (rad)
        and y = Dec (rad). This is the system favoured by Hilbert et al.

    INITIALISATION
        kappafile      Name of file containing a convergence map
        FITS           Data file format (def=True)

    METHODS
        plot(self,fig_size=10,subplot=None): Plots either the whole image or a
                                             given sub-image in physical
                                             coordinates

    BUGS

    AUTHORS
      This file is part of the Pangloss project, distributed under the
      GPL v2, by Tom Collett (IoA) and  Phil Marshall (Oxford).
      Please cite: Collett et al 2013, http://arxiv.org/abs/1303.6564

    HISTORY
      2013-03-23  Marshall & Collett (Oxford)
      2015-06-24  Updated to use WLMap class, Everett (SLAC)
    """

# ----------------------------------------------------------------------------

    def __init__(self,kappafile=None,data=None,FITS=True):

        self.name = 'Convergence map kappa from Millenium Simulation, zs = 1.3857'
        # Calls the WLMap superclass
        pangloss.WLMap.__init__(self,mapfiles=kappafile,data=data,FITS=FITS)

# ----------------------------------------------------------------------------

    def __str__(self):
        ## Add more information!!
        return 'Convergence map'

# ----------------------------------------------------------------------------

    def plot(self,fig_size=10,subplot=None,coords='world'):
        """
        Plot the convergence as a grayscale image.

        Optional arguments:
            fig_size        Figure size in inches
            subplot         List of four plot limits [xmin,xmax,ymin,ymax]
            coords          Type of coordinates inputted for the subplot:
                            'pixel', 'physical', or 'world'
        """

        # Always start a figure when plotting kappa maps:
        fig = plt.figure('Pangloss Map')

        # Use plotting method from WLMap class to calculate values common to both Kappamaps and Shearmaps
        pix_xi,pix_xf,pix_yi,pix_yf,Lx,Ly,pix_Lx,pix_Ly,subplot = self.plot_setup(subplot,coords)

        # Set figure size:
        pangloss.set_figure_size(fig,fig_size,Lx,Ly)

        # Set the pixel and wcs axes
        imsubplot = [pix_xi, pix_xf, pix_yi, pix_yf]
        ax = pangloss.set_axes(fig,Lx,Ly,self.hdr[0],imsubplot)

        # Used for colormap scaling
        Kmin = np.min(self.values[0][pix_yi:pix_yf,pix_xi:pix_xf])
        Kmax = np.max(self.values[0][pix_yi:pix_yf,pix_xi:pix_xf])

        # Plot Kappamap image
        ax.imshow(self.values[0],cmap='gray_r',vmin=Kmin,vmax=Kmax,origin='lower')

        return
