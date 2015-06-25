
import numpy, os, string
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
from wlmap import WLMap
#from astropy.io.fits import fits as pyfits

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
        # Calls the map superclass
        WLMap.__init__(self,kappafile,FITS)

# ----------------------------------------------------------------------------

    def __str__(self):
        ## Add more information!!
        return 'Convergence map'

# ----------------------------------------------------------------------------
# Plot the convergence as grayscale:
    
    def plot(self,fig_size=10,subplot=None): # fig_size in inches, default subplot is entire image

        # Maybe add the ability for subplot to be entered in physical coordiantes as well
        if subplot is None:
            subplot = [0,self.NX,0,self.NX]
            
        xi, xf = subplot[0], subplot[1]    # x-limits for subplot
        yi, yf = subplot[2], subplot[3]    # y-limits for subplot
        Lx = xf-xi    # length of x-axis subplot
        Ly = yf-yi    # length of y-axis subplot
        N = 8    #number of ticks
            
        plt.imshow(self.values[yi:yf,xi:xf],cmap = 'gray_r',origin = 'lower')
        plt.title('Convergence map of '+self.input)
        
        # Convert axes to physical coordinates, scale correctly with subplot
        xlNew = []; ylNew = [];
        xl = numpy.arange(xi,xf,Lx/N)
        yl = numpy.arange(yi,yf,Ly/N)
        
        for x in xl:
            xN,yN = self.image2physical(x,0)
            xlNew.append(xN)
        for y in yl:
            xN,yN = self.image2physical(0,y)
            ylNew.append(yN)
        
        # Format coordinates
        xlabels = ['%.3f' % a for a in xlNew]
        ylabels = ['%.3f' % a for a in ylNew]
        xlocs = numpy.arange(0,Lx,Lx/N)
        ylocs = numpy.arange(0,Ly,Ly/N)
        
        # Label axes
        plt.xticks(xlocs,xlabels)
        plt.yticks(ylocs,ylabels)
        plt.xlabel('Physical Coordinate (rad)')
        plt.ylabel('Physical Coordinate (rad)')

        fig = plt.gcf()
        fig.set_size_inches(fig_size,fig_size)
        return None

# ============================================================================
'''
Old Test

if __name__ == '__main__':

    import pylab as plt

    test1=True
    test2=False

# ----------------------------------------------------------------------------

    if test1==True:
    # Self-test: read in map from Stefan, and look up some convergence values.
      vb = True
      FITS = False

      for ext in ( "gamma_1", "fits" ):

          print "Testing ."+ext+" file..."

          # Read in map (and write out as FITS if it doesn't exist):
          kappafile = "/data/tcollett/Pangloss/gammafiles/GGL_los_8_1_1_N_4096_ang_4_rays_to_plane_37_f."+ext
          if ext == "fits": FITS = True
          convergence = Kappamap(kappafile,FITS=FITS)

          # Look up value in some pixel:

          i = 2184 ; j = 2263
          kappa = convergence.at(i,j,coordinate_system='image')
          print "Compare with expected value: 0.0169251"

          # Check WCS / physical coords at same pixel:

          x = 0.00232653752829; y = 0.00367303177544
          kappa = convergence.at(x,y,coordinate_system='physical')
          print "Compare with expected value: 0.0169251"

          print "...done."
          print " "

# ----------------------------------------------------------------------------

    if test2 ==True:

        kappafiles=[]
        for i in range(7):
           for j in range(7):
              kappafiles += ["/data/tcollett/Pangloss/kappafiles/GGL_los_8_%i_%i_N_4096_ang_4_rays_to_plane_37_f.kappa"%(i+1,j+1)]

        l=4096
        U=16

        kappa = numpy.zeros((l/U,l/U,len(kappafiles)))
        print numpy.shape(kappa)

        for k in range(len(kappafiles)):
           convergence = Kappamap(kappafiles[k],FITS=False)
           for i in range(l/U):
              if i % 500 ==0 : print k,",",i
              for j in range(l/U):
                 kappa[i,j,k] = convergence.at(U*i,U*j,coordinate_system='image')

        kappa=kappa.ravel()

        pangloss.writePickle(kappa,'kappalist.dat')
        print numpy.mean(kappa)

        plt.hist(kappa)
        plt.show()

# ============================================================================
'''
