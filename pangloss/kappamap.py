
import numpy as np
import matplotlib.pyplot as plt

import pangloss

vb = False

# ============================================================================

class Kappamap(pangloss.WLMap):
    """
    Read in, store, transform and interrogate a convergence map.

    Notes
    -----
    A "physical" coordinate system is used, where x = -RA (rad)
    and y = Dec (rad). This is the system favoured by Hilbert et al.

    Parameters
    ----------
    kappafile : string
        Name of file containing a convergence map
    FITS : boolean
        Data file format (def=True)

    """

# ----------------------------------------------------------------------------

    def __init__(self, kappafile=None, data=None, FITS=True):

        self.name = 'Convergence map kappa from Millenium Simulation, zs = 1.3857'
        # Calls the WLMap superclass
        pangloss.WLMap.__init__(self, mapfiles=kappafile, data=data, FITS=FITS)

# ----------------------------------------------------------------------------

    def __str__(self):
        ## Add more information!!
        return 'Convergence map'

# ----------------------------------------------------------------------------

    def plot(self, fig_size=10, subplot=None, coords='world'):
        """
        Plot the convergence as a grayscale image.

        Parameters
        ----------
        fig_size : float
            Figure size in inches
        subplot : list, float
            Plot limits [xmin,xmax,ymin,ymax]
        coords : string
            Input coordinate system for the subplot: 'pixel', 'physical', or 'world'
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

        # Get the colormap limits
        Kmin = np.min(self.values[0]
                                 [int(pix_yi):int(pix_yf),
                                  int(pix_xi):int(pix_xf)])
        Kmax = np.max(self.values[0]
                                 [int(pix_yi):int(pix_yf),
                                  int(pix_xi):int(pix_xf)])

        # Plot Kappamap image
        ax.imshow(self.values[0],cmap='gray_r',vmin=Kmin,vmax=Kmax,origin='lower')

        return

    @staticmethod
    def example():
        return Kappamap(pangloss.KAPPA_FILE, FITS=False)