
import numpy
import matplotlib.pyplot as plt
from wlmap import WLMap

arcmin2rad = (1.0/60.0)*numpy.pi/180.0
rad2arcmin = 1.0/arcmin2rad
deg2rad = numpy.pi/180.0
rad2deg = 1.0/deg2rad

vb = False

# ============================================================================

class Kappamap(WLMap):
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
      2015-06-24  Everett (SLAC)
    """

# ----------------------------------------------------------------------------

    def __init__(self,kappafile,FITS=True):

        self.name = 'Convergence map kappa from Millenium Simulation, zs = 1.6'
        # Calls the WLMap superclass
        WLMap.__init__(self,kappafile,FITS)

# ----------------------------------------------------------------------------

    def __str__(self):
        ## Add more information!!
        return 'Convergence map'

# ----------------------------------------------------------------------------

    def plot(self,fig_size=10,subplot=None,coords='pixel'):
        """
        Plot the convergence as a grayscale image.

        Optional arguments:
            fig_size        Figure size in inches
            subplot         List of four plot limits [xmin,xmax,ymin,ymax]
            coords          'pixel' or 'physical' ('world' not yet supported)
        """

        # Default subplot is entire image
        if subplot is None:
            subplot = [0,self.NX[0],0,self.NX[0]]

        xi, xf = subplot[0], subplot[1]    # x-limits for subplot
        yi, yf = subplot[2], subplot[3]    # y-limits for subplot
        Lx = xf-xi    # length of x-axis subplot
        Ly = yf-yi    # length of y-axis subplot
        # Number of axis ticks
        if Lx != Ly and Lx/Ly < 0.6:
            N = 5
        else:
            N = 8

        # N-sampled axis values
        xl = numpy.arange(xi,xf+Lx/N*numpy.sign(xf),Lx/N)
        yl = numpy.arange(yi,yf+Ly/N*numpy.sign(yf),Ly/N)

        if coords == 'pixel':
            # Convert axes to physical coordinates, scale correctly with subplot
            xlNew = []; ylNew = [];

            for x in xl:
                xN,yN = self.image2physical(x,0)
                xlNew.append(xN)
            for y in yl:
                xN,yN = self.image2physical(0,y)
                ylNew.append(yN)

            # Format coordinates
            xlabels = ['%.3f' % a for a in xlNew]
            ylabels = ['%.3f' % a for a in ylNew]

        elif coords == 'physical':
            # Label values are already in physical coordinates
            xlabels = ['%.3f' % a for a in xl]
            ylabels = ['%.3f' % a for a in yl]

            # Convert subplot bounds to pixel values
            xi,yi = self.physical2image(xi,yi)
            xf,yf = self.physical2image(xf,yf)
            Lx = xf-xi
            Ly = yf-yi

        else:
            raise IOError('Error: Subplot bounds can only be in pixel or physical coordinates.')

        # Location of tick marks
        xlocs = numpy.arange(0,Lx,Lx/N)
        ylocs = numpy.arange(0,Ly,Ly/N)

        # Plot image
        plt.imshow(self.values[0][yi:yf,xi:xf],cmap = 'gray_r',origin = 'lower')
        plt.title('Convergence map of '+self.input[0])

        # Label axes
        plt.xticks(xlocs,xlabels)
        plt.yticks(ylocs,ylabels)
        plt.xlabel('Physical Coordinate (rad)')
        plt.ylabel('Physical Coordinate (rad)')

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
