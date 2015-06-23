
import numpy, os, string
import struct
import cPickle
import astropy.io.fits as pyfits
#from astropy.io.fits import fits as pyfits

arcmin2rad = (1.0/60.0)*numpy.pi/180.0
rad2arcmin = 1.0/arcmin2rad
deg2rad = numpy.pi/180.0
rad2deg = 1.0/deg2rad

vb = False

# ============================================================================

class Kappamap:
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
        read_in_fits_data(self):

        read_in_binary_data(self): to cope with hilbert's homegrown format

        setwcs(self): simulated maps often don't have WCS

        get_fits_wcs(self,hdr):

        write_out_to_fits(self):

        image2physical(self,i,j): coord transformation, returns x,y

        physical2image(self,x,y): coord transformation, returns i,j

        image2world(self,i,j): coord transformation, returns a,d

        world2image(self,a,d): coord transformation, returns a,d

        at(self,x,y,coordinate_system='physical'): return pixel values

        lookup(self,i,j): return pixel values given image coords

    BUGS

    AUTHORS
      This file is part of the Pangloss project, distributed under the
      GPL v2, by Tom Collett (IoA) and  Phil Marshall (Oxford). 
      Please cite: Collett et al 2013, http://arxiv.org/abs/1303.6564

    HISTORY
      2013-03-23  Marshall & Collett (Oxford)
    """

# ----------------------------------------------------------------------------

    def __init__(self,kappafile,FITS=True):

        self.name = 'Convergence map kappa from Millenium Simulation, zs = 1.6'
        self.input = kappafile

        # Read in data from file:

        if FITS:
            # Read in FITS image, extract wcs:
            if vb: print "Reading in map from file "+kappafile
            self.read_in_fits_data()

        else:
            # Initialise some necessary WCS parameters:
            self.field = 4.0 # degrees
            self.NX = 4096
            self.PIXSCALE = self.field/(1.0*self.NX) # degrees
            self.setwcs()

            # Read in binary data, to self.values:
            if vb: print "Reading in map from file "+kappafile
            self.read_in_binary_data()

            # If it doesn't already exist, output the map to FITS file:
            # pieces = string.split(self.input,'.')
            # self.output = string.join(pieces[0:len(pieces)-1],'.')+'.fits'
            self.output = self.input+'.fits'
            if os.path.exists(self.output): 
              if vb: print "FITS version already exists: ",self.output
            else:  
              if vb: print "Writing map to "+self.output
              self.write_out_to_fits()
            # This should probably not be in __init__ but hopefully it only gets run once.  

        return None

# ----------------------------------------------------------------------------

    def __str__(self):
        return 'Convergence map'

# ----------------------------------------------------------------------------

    def read_in_fits_data(self):
        hdu = pyfits.open(self.input)[0]
        hdr = hdu.header
        self.get_fits_wcs(hdr)
        self.values = hdu.data
        # This transpose is necessary so that ds9 displays the image correctly.
        self.values = self.values.transpose()
        self.NX = self.values.shape[0]
        self.PIXSCALE = self.wcs['CD1_1']
        self.field = self.NX*self.PIXSCALE

        return None

# ----------------------------------------------------------------------------

    def read_in_binary_data(self):
        file = open(self.input,"rb")
        data = file.read()
        fmt = str(self.NX*self.NX)+'f'
        start = 0
        stop = struct.calcsize(fmt)
        values = struct.unpack(fmt,data[start:stop])
        self.values = numpy.array(values,dtype=numpy.float32).reshape(self.NX,self.NX)
        return None

# ----------------------------------------------------------------------------
# WCS parameters: to allow conversions between
#  image coordinates i,j (pixels)
#  physical coordinates x,y (rad)
#  sky coordinates ra,dec (deg, left-handed system)

    def setwcs(self):
        self.wcs = dict()
        # ra  = CRVAL1 + CD1_1*(i-CRPIX1) 
        # dec = CRVAL2 + CD2_2*(j-CRPIX2) 
        self.wcs['CRPIX1'] = 0.0
        self.wcs['CRPIX2'] = 0.0
        self.wcs['CRVAL1'] =  0.5*self.field + 0.5*self.PIXSCALE
        self.wcs['CRVAL2'] = -0.5*self.field + 0.5*self.PIXSCALE
        self.wcs['CD1_1'] = -self.PIXSCALE
        self.wcs['CD1_2'] = 0.0
        self.wcs['CD2_1'] = 0.0
        self.wcs['CD2_2'] = self.PIXSCALE
        self.wcs['CTYPE1'] = 'RA---TAN'
        self.wcs['CTYPE2'] = 'DEC--TAN'

        # i = LTV1 + LTM1_1*(x/rad) 
        # j = LTV2 + LTM2_2*(y/rad) 
        self.wcs['LTV1'] = 0.5*self.field/self.PIXSCALE - 0.5
        self.wcs['LTV2'] = 0.5*self.field/self.PIXSCALE - 0.5
        self.wcs['LTM1_1'] = 1.0/(self.PIXSCALE*deg2rad)
        self.wcs['LTM2_2'] = 1.0/(self.PIXSCALE*deg2rad)
        return None

# ----------------------------------------------------------------------------

    def get_fits_wcs(self,hdr):
        self.wcs = dict()
        for keyword in hdr.keys():
           self.wcs[keyword] = hdr[keyword]      
        return None

# ----------------------------------------------------------------------------

    def write_out_to_fits(self):

        # Start a FITS header + data unit:
        hdu = pyfits.PrimaryHDU()
        # Add WCS keywords to the FITS header (in apparently random order):
        for keyword in self.wcs.keys():
          hdu.header.update(keyword,self.wcs[keyword])
        # Make image array. The transpose is necessary so that ds9 displays 
        # the image correctly.
        hdu.data = self.values.transpose()
        # Verify and write to file:
        hdu.verify()
        hdu.writeto(self.output)

        return None

# ----------------------------------------------------------------------------
# Interpolating the map to return a single value at a specified point - this 
# is the most important method of this class.

    def at(self,x,y,coordinate_system='physical'):

        if vb: 
            print " "
            print "Looking up kappa value at position",x,",",y," in the "+coordinate_system+" coordinate system"

        # Get pixel indices of desired point, 
        # and also work out other positions for completeness, if verbose:
        if coordinate_system == 'physical':
           i,j = self.physical2image(x,y)
           if vb: print "  - image coordinates:",i,j
        elif coordinate_system == 'image':
           i = x
           j = y
           if vb: 
               x,y = self.image2physical(i,j)
               print "  - physical coordinates:",x,y,"(radians)"

        if vb:
            a,d = self.image2world(i,j)
            print "  - approximate world coordinates:",a,d,"(degrees)"
            print "  - ds9 image coordinates:",i+1,j+1

        # Now look up correct value, doing some bilinear interpolation:

        kappa = self.lookup(i,j)
        if vb: print "  Value of kappa = ",kappa

        return kappa

# ----------------------------------------------------------------------------

    def image2physical(self,i,j):
        x = (i - self.wcs['LTV1'])/self.wcs['LTM1_1'] # x in rad 
        y = (j - self.wcs['LTV2'])/self.wcs['LTM2_2'] # y in rad 
        return x,y

    def physical2image(self,x,y):
        i = self.wcs['LTV1'] + self.wcs['LTM1_1']*x # x in rad 
        j = self.wcs['LTV2'] + self.wcs['LTM2_2']*y # y in rad 
        return i,j

     # Only approximate WCS transformations - assumes dec=0.0 and small field
    def image2world(self,i,j):
        a = self.wcs['CRVAL1'] + self.wcs['CD1_1']*(i - self.wcs['CRPIX1'])
        if a < 0.0: a += 360.0
        d = self.wcs['CRVAL2'] + self.wcs['CD2_2']*(j - self.wcs['CRPIX2'])
        return a,d

    def world2image(self,a,d):
        i = (a - self.wcs['CRVAL1'])/self.wcs['CD1_1'] + self.wcs['CRPIX1']
        j = (d - self.wcs['CRVAL2'])/self.wcs['CD2_2'] + self.wcs['CRPIX2']
        return i,j

# ----------------------------------------------------------------------------

    def lookup(self,i,j):
      
        # Weighted mean of 4 neighbouring pixels, as suggested by Stefan.
        ix = int(i)
        iy = int(j)
        px = i - ix
        py = j - iy

        if ((0 <= ix) and (ix < self.NX-1) and (0 <= iy) and (iy < self.NX-1)):
            mean =   self.values[ix,iy]    *(1.0-px)*(1.0-py) \
                  + self.values[ix+1,iy]  * px     *(1.0-py) \
                  + self.values[ix,iy+1]  *(1.0-px)* py      \
                  + self.values[ix+1,iy+1]* px     * py
        else:
            mean = 0.0        

        return mean
      
# ============================================================================

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

