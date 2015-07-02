
import numpy as np, os, string
import struct
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt

arcmin2rad = (1.0/60.0)*np.pi/180.0
rad2arcmin = 1.0/arcmin2rad
deg2rad = np.pi/180.0
rad2deg = 1.0/deg2rad

vb = False

class WLMap:
    """
    NAME
        WLMap

    PURPOSE
        Superclass for kappamap and gammamap. Used to read in, store,
        transform and interrogate either a convergence or shear map.

    COMMENTS
        A "physical" coordinate system is used, where x = -RA (rad)
        and y = Dec (rad). This is the system favoured by Hilbert et al.

    INITIALISATION
        mapfile      List of file names containing a convergence or shear map.
        FITS           Data file format (def=True)

    METHODS
        read_in_fits_data(self):

        read_in_binary_data(self): to cope with hilbert's homegrown format

        setwcs(self): simulated maps often don't have WCS

        get_fits_wcs(self,hdr,i):

        write_out_to_fits(self,i):

        image2physical(self,i,j,mapfile=0): coord transformation, returns x,y

        physical2image(self,x,y,mapfile=0): coord transformation, returns i,j

        image2world(self,i,j,mapfile=0): coord transformation, returns a,d

        world2image(self,a,d,mapfile=0): coord transformation, returns a,d

        at(self,x,y,mapfile=0,coordinate_system='physical'): return pixel values

        lookup(self,i,j,mapfile): return pixel values given image coords

    BUGS

    AUTHORS
      This file is part of the Pangloss project, distributed under the
      GPL v2, by Tom Collett (IoA) and  Phil Marshall (Oxford).
      Please cite: Collett et al 2013, http://arxiv.org/abs/1303.6564

    HISTORY
      2015-06-24  Everett (SLAC)
    """

    def __init__(self,mapfiles,FITS=True):

        # mapfile should be inputted as list, but is automatically converted to
        # a list if it is a single file string
        if type(mapfiles) != list:
            mapfiles = [mapfiles]
        self.input = mapfiles
        
        # Parsing the file name(s)
        # 0 <= x,y <= 7, each (x,y) map covers 4x4 degrees
        self.field_x = []
        self.field_y = []
        for i in range(0,len(self.input)):
            input_parse = self.input[i].split('_') # Creates list of filename elements separated by '_'
            self.field_x.append(eval(input_parse[3])) # The x location of the map grid
            self.field_y.append(eval(input_parse[4])) # The y location of the map grid

        # Declare needed attributes as lists
        self.values = []
        self.NX = []
        self.PIXSCALE = []
        self.field = []
        self.wcs = []
        self.output = []

        # Read in data from file:
        if FITS:
            # Read in FITS image, extract wcs:
            if vb: print "Reading in map from files "+mapfiles
            self.read_in_fits_data()

        else:
            # Read in binary data, to self.values:
            if vb: print "Reading in map from file "+mapfiles
            self.read_in_binary_data()
        return None

# ----------------------------------------------------------------------------

    def __str__(self):
        return 'abstract map'

# ----------------------------------------------------------------------------

    def read_in_fits_data(self):
        for i in range(0,len(self.input)):
            hdu = pyfits.open(self.input[i])[0]
            hdr = hdu.header
            self.get_fits_wcs(hdr,i)
            self.values.append(hdu.data)
            # This transpose is necessary so that ds9 displays the image correctly.
            self.values[i] = self.values[i].transpose()
            self.NX.append(self.values[i].shape[0])
            self.PIXSCALE.append(self.wcs[i]['CD1_1'])
            self.field.append(self.NX[i]*self.PIXSCALE[i])
            return None

# ----------------------------------------------------------------------------

    def read_in_binary_data(self):
        for i in range(0,len(self.input)):
            # Initialise some necessary WCS parameters:
            self.field.append(4.0) # degrees
            self.NX.append(4096) # pixels
            self.PIXSCALE.append(self.field[i]/(1.0*self.NX[i])) # degrees
            self.setwcs(i)

            mapfile = open(self.input[i],"rb")
            data = mapfile.read()
            fmt = str(self.NX[i]*self.NX[i])+'f'
            start = 0
            stop = struct.calcsize(fmt)
            values = struct.unpack(fmt,data[start:stop])
            self.values.append(np.array(values,dtype=np.float32).reshape(self.NX[i],self.NX[i]))

            # If it doesn't already exist, output the map to FITS file:
            # pieces = string.split(self.input,'.')
            # self.output = string.join(pieces[0:len(pieces)-1],'.')+'.fits'
            self.output.append(self.input[i]+'.fits')
            if os.path.exists(self.output[i]):
              if vb: print "FITS version already exists: ",self.output
            else:
              if vb: print "Writing map to "+self.output[i]
              self.write_out_to_fits(i)
            # This should probably not be in __init__ but hopefully it only gets run once.

# ----------------------------------------------------------------------------
#  WCS parameters: to allow conversions between
#  image coordinates i,j (pixels)
#  physical coordinates x,y (rad)
#  sky coordinates ra,dec (deg, left-handed system)

    def setwcs(self,i):
        self.wcs.append(dict())
        # ra  = CRVAL1 + CD1_1*(i-CRPIX1)
        # dec = CRVAL2 + CD2_2*(j-CRPIX2)
        self.wcs[i]['CRPIX1'] = 0.0
        self.wcs[i]['CRPIX2'] = 0.0
        self.wcs[i]['CRVAL1'] =  0.5*self.field[i] + 0.5*self.PIXSCALE[i] - self.field_x[i]*self.field[i]
        self.wcs[i]['CRVAL2'] = -0.5*self.field[i] + 0.5*self.PIXSCALE[i] + self.field_y[i]*self.field[i]
        self.wcs[i]['CD1_1'] = -self.PIXSCALE[i]
        self.wcs[i]['CD1_2'] = 0.0
        self.wcs[i]['CD2_1'] = 0.0
        self.wcs[i]['CD2_2'] = self.PIXSCALE[i]
        self.wcs[i]['CTYPE1'] = 'RA---TAN'
        self.wcs[i]['CTYPE2'] = 'DEC--TAN'

        # i = LTV1 + LTM1_1*(x/rad)
        # j = LTV2 + LTM2_2*(y/rad)
        self.wcs[i]['LTV1'] = 0.5*self.field[i]/self.PIXSCALE[i] - 0.5
        self.wcs[i]['LTV2'] = 0.5*self.field[i]/self.PIXSCALE[i] - 0.5
        self.wcs[i]['LTM1_1'] = 1.0/(self.PIXSCALE[i]*deg2rad)
        self.wcs[i]['LTM2_2'] = 1.0/(self.PIXSCALE[i]*deg2rad)
        return None

# ----------------------------------------------------------------------------

    def get_fits_wcs(self,hdr,i):
        self.wcs.append(dict())
        for keyword in hdr.keys():
           self.wcs[i][keyword] = hdr[keyword]
        return None

# ----------------------------------------------------------------------------

    def write_out_to_fits(self,i):

        # Start a FITS header + data unit:
        hdu = pyfits.PrimaryHDU()
        # Add WCS keywords to the FITS header (in apparently random order):
        for keyword in self.wcs[i].keys():
          hdu.header.update(keyword,self.wcs[i][keyword])
        # Make image array. The transpose is necessary so that ds9 displays
        # the image correctly.
        hdu.data = self.values[i].transpose()
        # Verify and write to file:
        hdu.verify()
        hdu.writeto(self.output[i])

        return None

# ----------------------------------------------------------------------------
# Interpolating the map to return a single value at a specified point - this
# is the most important method of this class.

    def at(self,x,y,mapfile=None,coordinate_system='physical'):
        # mapfile is the index of the inputted maps that the at() method retrieves values from

        # Shearmaps have two components - which one to look-up *must*
        # be specified! But Kappamaps only have one component, so in
        # this case we can interpret a lack of mapfile kwarg more
        # sympathetically:
        if mapfile is None:
            if len(self.input) > 1:
                return None
            else:
                mapfile = 0

        if vb:
            print " "
            print "Looking up map"+mapfile+" value at position",x,",",y," in the "+coordinate_system+" coordinate system"

        # Get pixel indices of desired point,
        # and also work out other positions for completeness, if verbose:
        if coordinate_system == 'physical':
           i,j = self.physical2image(x,y,mapfile)
           if vb: print "  - image coordinates of map "+mapfile+":",i,j
        elif coordinate_system == 'image':
           i = x
           j = y
           if vb:
               x,y = self.image2physical(i,j,mapfile)
               print "  - physical coordinates of map"+mapfile+":",x,y,"(radians)"

        if vb:
            a,d = self.image2world(i,j,mapfile)
            print "  - approximate world coordinates of map"+mapfile+":",a,d,"(degrees)"
            print "  - ds9 image coordinates:",i+1,j+1

        # Now look up correct value, doing some bilinear interpolation:

        value = self.lookup(i,j,mapfile)
        if vb: print "  Value of map "+mapfile+" at ("+x+", "+y+")= ",value

        return value

# ----------------------------------------------------------------------------

    def image2physical(self,i,j,mapfile=0):
        x = (i - self.wcs[mapfile]['LTV1'])/self.wcs[mapfile]['LTM1_1'] # x in rad
        y = (j - self.wcs[mapfile]['LTV2'])/self.wcs[mapfile]['LTM2_2'] # y in rad
        return x,y

    def physical2image(self,x,y,mapfile=0):
        i = self.wcs[mapfile]['LTV1'] + self.wcs[mapfile]['LTM1_1']*x # x in rad
        j = self.wcs[mapfile]['LTV2'] + self.wcs[mapfile]['LTM2_2']*y # y in rad
        return i,j

     # Only approximate WCS transformations - assumes dec=0.0 and small field
    def image2world(self,i,j,mapfile=0):
        a = self.wcs[mapfile]['CRVAL1'] + self.wcs[mapfile]['CD1_1']*(i - self.wcs[mapfile]['CRPIX1'])
        if a < 0.0: a += 360.0
        d = self.wcs[mapfile]['CRVAL2'] + self.wcs[mapfile]['CD2_2']*(j - self.wcs[mapfile]['CRPIX2'])
        return a,d

    def world2image(self,a,d,mapfile=0):
        i = (a - self.wcs[mapfile]['CRVAL1'])/self.wcs[mapfile]['CD1_1'] + self.wcs[mapfile]['CRPIX1']
        # if a negative pixel is returned for i, reinput a as a negative degree
        if i<0:
            a-=360
            i = (a - self.wcs[mapfile]['CRVAL1'])/self.wcs[mapfile]['CD1_1'] + self.wcs[mapfile]['CRPIX1']
        j = (d - self.wcs[mapfile]['CRVAL2'])/self.wcs[mapfile]['CD2_2'] + self.wcs[mapfile]['CRPIX2']
        return i,j
           
    def physical2world(self,x,y,mapfile=0):
        a = -np.rad2deg(x) - self.field_x[mapfile]*self.field[mapfile]
        if a < 0.0: a += 360.0
        d = np.rad2deg(y) + self.field_y[mapfile]*self.field[mapfile]
        return a,d
    
    def world2physical(self,a,d,mapfile=0):
        x = -np.deg2rad(a + self.field_x[mapfile]*self.field[mapfile])
        y = np.deg2rad(d - self.field_y[mapfile]*self.field[mapfile])
        return x,y

# ----------------------------------------------------------------------------
# In general, we don't know how to plot this map...

    def plot(self,fig_size=10,subplot=None,coords='pixel'):
        '''
        Plot the convergence as a grayscale image.

        Optional arguments:
            fig_size        Figure size in inches
            subplot         List of four plot limits [xmin,xmax,ymin,ymax]
            coords          Type of coordinates inputted for the subplot:
                            'pixel', 'physical', or 'world'
            
        *!NOTE!*: Not a complete plotting method. Ony calculates values common 
        to both Kappamap and Shearmap plots
        '''
        # Default subplot is entire image
        if subplot is None:
            '''
            # Switch these two versions to change default coords from 'pixel' to 'world'
            # coords = 'world':
            ai, di = self.image2world(0,0)
            af, df = self.image2world(self.NX[0],self.NX[0])
            subplot = [ai,af,di,df]
            # coords = 'pixel':
            '''
            subplot = [0,self.NX[0],0,self.NX[0]]
            
        xi, xf = subplot[0], subplot[1]    # x-limits for subplot
        yi, yf = subplot[2], subplot[3]    # y-limits for subplot
        
        Lx = abs(xf-xi)    # length of x-axis subplot
        Ly = abs(yf-yi)    # length of y-axis subplot
        
        #Number of axis ticks
        if (Lx != Ly and Lx/Ly < 0.6) or fig_size < 8:
            tickNum = 5
        else:
            tickNum = 8
            
        # N-sampled axis values
        xl = np.arange(xi,xf+Ly/tickNum*np.sign(yf),Lx/tickNum)     
        yl = np.arange(yi,yf+Ly/tickNum*np.sign(yf),Ly/tickNum) 
        
        '''
        This can be used instead of the above code to set the axis values when 
        the default value for coords is 'world'. However, it is becoming too 
        cumbersome when the previous code was much more clear.
        # N-sampled axis values
        if coords == 'world' and xf > xi:
            # This means that the plot includes the breaking point between 0 and
            # 360, so the axis values must be adjusted
            xf -= 360
            Lx = abs(xf-xi)    # length of x-axis subplot
            xl = np.arange(xi,xf+Lx/tickNum*np.sign(xf),Lx/tickNum)
            xl_above0 = [n for n in xl if n <= 0]
            xl_below0 = [n for n in xl if n > 0]
            xl_below0 = [n+360 for n in xl_below0] 
            print(xl_above0,xl_below0)
            xl = xl_above0+xl_below0
            
        else:
            # For all other cases
            Lx = abs(xf-xi)
            xl = np.arange(xi,xf+Lx/tickNum*np.sign(xf),Lx/tickNum)
        
        Ly = abs(yf-yi)    # length of y-axis subplot
        yl = np.arange(yi,yf+Ly/tickNum*np.sign(yf),Ly/tickNum)        
        '''
        
        if coords == 'pixel':
            # Convert axes to world coordinates, scale correctly with subplot
            xlNew = []; ylNew = [];

            for x in xl:
                xN,yN = self.image2world(x,0)
                xlNew.append(xN)
            for y in yl:
                xN,yN = self.image2world(0,y)
                ylNew.append(yN)

            # Format coordinates
            xlabels = ['%.5f' % a for a in xlNew]
            ylabels = ['%.5f' % a for a in ylNew]

        elif coords == 'physical':
            # Convert axes to world coordinates, scale correctly with subplot
            xlNew = []; ylNew = [];

            for x in xl:
                xN,yN = self.physical2world(x,0)
                xlNew.append(xN)
            for y in yl:
                xN,yN = self.physical2world(0,y)
                ylNew.append(yN)

            # Format coordinates
            xlabels = ['%.5f' % a for a in xlNew]
            ylabels = ['%.5f' % a for a in ylNew]

            # Convert subplot bounds to pixel values
            xi,yi = self.physical2image(xi,yi)
            xf,yf = self.physical2image(xf,yf)
            Lx = xf-xi
            Ly = yf-yi
            
        elif coords == 'world':
            # Label values are already in world coordinates
            xlabels = ['%.5f' % a for a in xl]
            ylabels = ['%.5f' % a for a in yl]

            # Convert subplot bounds to pixel values
            xi,yi = self.world2image(xi,yi)
            xf,yf = self.world2image(xf,yf)
            
            Lx = xf-xi
            Ly = yf-yi

        else:
            raise IOError('Error: Subplot bounds can only be in pixel, physical, or world coordinates.')

        # Location of tick marks
        xlocs = np.arange(0,Lx,Lx/tickNum)
        ylocs = np.arange(0,Ly,Ly/tickNum)
        
        return xi,xf,yi,yf,Lx,Ly,xlocs,xlabels,ylocs,ylabels

# ----------------------------------------------------------------------------

    def lookup(self,i,j,mapfile):

        # Weighted mean of 4 neighbouring pixels, as suggested by Stefan.
        ix = int(i)
        iy = int(j)
        px = i - ix
        py = j - iy

        if ((0 <= ix) and (ix < self.NX[mapfile]-1) and (0 <= iy) and (iy < self.NX[mapfile]-1)):
            mean =  self.values[mapfile][ix,iy]    * (1.0-px) * (1.0-py) \
                  + self.values[mapfile][ix+1,iy]  * px       * (1.0-py) \
                  + self.values[mapfile][ix,iy+1]  * (1.0-px) * py      \
                  + self.values[mapfile][ix+1,iy+1] * px * py
        else:
            mean = None

        return mean
