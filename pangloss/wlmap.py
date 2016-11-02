import os
import struct

import pyfits
import numpy as np
from io import int_or_float

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
      2015-06-24  Started Everett (SLAC)
    """

    def __init__(self, mapfiles=None, data=None, FITS=True):

        assert not (mapfiles is not None and data is not None), \
        "Error: Can only input mapfiles or data; not both!"

        # Make a WLMap by reading in data from a file.
        # Mapfile should be inputted as list, but is automatically
        # converted to a list if it is a single file string
        if mapfiles is not None:
            if type(mapfiles) != list:
                mapfiles = [mapfiles]
            self.input = mapfiles
            self.origin = 'file'

        # Make a WLMap from scratch, given 3 numpy arrays: d,domain,map_xy
        # d must be a 3D array, to allow for shear as well as kappa.
        # x and y must be meshgrids that specify the world coordinates
        # of the pixel centers. domain =[ra_i, ra_f, dec_i, dec_f],
        # map_xy = [map_x, map_y]. (no longer need x and y)
        elif data is not None:
            assert type(data) == list
            assert len(data) == 3
            #d,x,y,domain,map_xy = data[0],data[1],data[2],data[3],data[4] # OLD - don't need x,y if domain is passed
            d,domain,map_xy = data[0],data[1],data[2]
            assert len(d.shape) == 3
            assert len(domain) == 4
            assert len(map_xy) == 2
            #assert d[0].shape[0] == x.shape[0]-1
            #assert d[0].shape[1] == y.shape[0]-1
            #assert x.shape == y.shape # For now, must be square!
            if len(d) == 2: assert d[0].shape == d[1].shape
            self.input = []
            self.origin = 'scratch'

        else:
            print "Empty "+self.__str__()
            return None

        # Declare needed attributes as lists
        self.hdr = []
        self.values = []
        self.NX = []
        self.PIXSCALE = []
        self.field = []
        self.wcs = []
        self.output = []
        self.map_x = []
        self.map_y = []

        if self.origin == 'file':
            # Parsing the file name(s)
            # 0 <= x,y <= 7, each (x,y) map covers 4x4 square degrees
            for i in range(0,len(self.input)):
                input_parse = self.input[i].split('_') # Creates list of filename elements separated by '_'
                self.map_x.append(int_or_float(input_parse[3])) # The x location of the map grid
                self.map_y.append(int_or_float(input_parse[4])) # The y location of the map grid

            # Read in data from file:
            if FITS:
                # Read in FITS image, extract wcs:
                if vb: print "Reading in map from files "+mapfiles
                self.read_in_fits_data()

            else:
                # Read in binary data, to self.values:
                if vb: print "Reading in map from file "+mapfiles
                self.read_in_binary_data()

        elif self.origin == 'scratch':
            # Reformat data into self.values:
            if vb: print "Making map from data arrays"
            for i in range(0,len(d)):
                # The map (x,y) positions corresponding to the Hilbert maps
                self.map_x.append(map_xy[0]) # The x location of the map grid
                self.map_y.append(map_xy[1]) # The y location of the map grid

            self.reformat_data(d,domain)

            # Set up the WCS:
            #self.make_wcs(x,y)

        # Check file consistency
        assert self.PIXSCALE[1:] == self.PIXSCALE[:-1]
        assert self.field[1:] == self.field[:-1]
        assert self.NX[1:] == self.NX[:-1]

        return None

# ----------------------------------------------------------------------------

    def __str__(self):
        return 'Weak Lensing Map object.'

# ----------------------------------------------------------------------------

    def read_in_fits_data(self):
        for i in range(0,len(self.input)):
            hdu = pyfits.open(self.input[i])[0]
            self.hdr.append(hdu.header)
            self.get_fits_wcs(self.hdr[i],i)
            self.values.append(hdu.data)
            hdu.close()
            # This transpose would be necessary so that ds9 displays the image correctly, if we hadn't done it already in read_in_binary_data()
            # self.values[i] = self.values[i].transpose()
            self.NX.append(self.values[i].shape[0])
            self.PIXSCALE.append(-self.wcs[i]['CD1_1'])
            self.field.append(self.NX[i]*self.PIXSCALE[i])
        return None

# ----------------------------------------------------------------------------

    def read_in_binary_data(self):

        for i in range(0,len(self.input)):

            # Initialise some necessary WCS parameters for Stefan
            # Hilbert's binary data files:
            self.field.append(4.0) # degrees
            self.NX.append(4096) # pixels
            self.PIXSCALE.append(self.field[i]/(1.0*self.NX[i])) # degrees
            self.setwcs(i)

            mapfile = open(self.input[i],"rb")
            data = mapfile.read()
            mapfile.close()
            fmt = str(self.NX[i]*self.NX[i])+'f'
            start = 0
            stop = struct.calcsize(fmt)
            values = struct.unpack(fmt,data[start:stop])
            self.values.append(np.array(values,dtype=np.float32).reshape(self.NX[i],self.NX[i]).transpose())

            # If it doesn't already exist, output the map to FITS file:
            # pieces = string.split(self.input,'.')
            # self.output = string.join(pieces[0:len(pieces)-1],'.')+'.fits'
            self.output.append(self.input[i]+'.fits')
            if os.path.exists(self.output[i]):
              if vb: print "FITS version already exists: ",self.output
              self.make_header(i)
            else:
              if vb: print "Writing map to "+self.output[i]
              self.write_out_to_fits(i)
            # This should probably not be in __init__ but hopefully it only gets run once.

        return

# ----------------------------------------------------------------------------

    def reformat_data(self,d,domain):
        #
        for i in range(0,len(d)):
            #
            self.values.append(d[i])

            ra_i, ra_f, dec_i, dec_f = domain[0], domain[1], domain[2], domain[3] # degrees
            dra, ddec = abs(ra_i-ra_f), abs(dec_i-dec_f) # degrees
            assert dra == ddec
            self.field.append(dra) # degrees
            #self.NX.append(int( 4096.0 * (dra / 4.0) ))
            self.NX.append(self.values[i].shape[0])
            self.PIXSCALE.append(self.field[i]/(1.0*self.NX[i])) # degrees
            self.make_wcs(i,domain)
            self.make_header(i)
            self.input.append('')
            self.output.append('')

        return

# ----------------------------------------------------------------------------
    #def make_wcs(self,x,y):
    def make_wcs(self,i,domain):
        '''
        Make wcs for inputted data
        '''
        ra_i, ra_f, dec_i, dec_f = domain[0], domain[1], domain[2], domain[3]
        #ramin = np.min(np.min(x))
        #ramax = np.max(np.max(x))
        #decmin = np.min(np.min(y))
        #decmax = np.max(np.max(y))
        self.wcs.append(dict())
        # ra  = CRVAL1 + CD1_1*(i-CRPIX1)
        # dec = CRVAL2 + CD2_2*(j-CRPIX2)
        self.wcs[i]['CRPIX1'] = 0.0
        self.wcs[i]['CRPIX2'] = 0.0
        self.wcs[i]['CRVAL1'] =  ra_i
        self.wcs[i]['CRVAL2'] = dec_i
        self.wcs[i]['CD1_1'] = -self.PIXSCALE[i]
        self.wcs[i]['CD1_2'] = 0.0
        self.wcs[i]['CD2_1'] = 0.0
        self.wcs[i]['CD2_2'] = self.PIXSCALE[i]
        self.wcs[i]['CTYPE1'] = 'RA---TAN'
        self.wcs[i]['CTYPE2'] = 'DEC--TAN'

        '''
        # i = LTV1 + LTM1_1*(x/rad)
        # j = LTV2 + LTM2_2*(y/rad)
        self.wcs[i]['LTV1'] = 0.5*self.field[i]/self.PIXSCALE[i] - 0.5
        self.wcs[i]['LTV2'] = 0.5*self.field[i]/self.PIXSCALE[i] - 0.5
        self.wcs[i]['LTM1_1'] = 1.0/np.deg2rad(self.PIXSCALE[i])
        self.wcs[i]['LTM2_2'] = 1.0/np.deg2rad(self.PIXSCALE[i])
        '''
        return

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
        self.wcs[i]['CRVAL1'] =  0.5*self.field[i] - self.map_x[i]*self.field[i]
        self.wcs[i]['CRVAL2'] = -0.5*self.field[i] + self.map_y[i]*self.field[i]
        # self.wcs[i]['CRVAL1'] =  0.5*self.field[i] + 0.5*self.PIXSCALE[i] - self.map_x[i]*self.field[i]
        # self.wcs[i]['CRVAL2'] = -0.5*self.field[i] + 0.5*self.PIXSCALE[i] + self.map_y[i]*self.field[i]
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
        self.wcs[i]['LTM1_1'] = 1.0/np.deg2rad(self.PIXSCALE[i])
        self.wcs[i]['LTM2_2'] = 1.0/np.deg2rad(self.PIXSCALE[i])

        return None

# ----------------------------------------------------------------------------

    def get_fits_wcs(self,hdr,i):
        self.wcs.append(dict())
        for keyword in hdr.keys():
            self.wcs[i][keyword] = hdr[keyword]
        return None

# ----------------------------------------------------------------------------

    def make_header(self,i):
        # Start a FITS header + data unit:
        hdu = pyfits.PrimaryHDU()
        # Add WCS keywords to the FITS header (in apparently random order):
        for keyword in self.wcs[i].keys():
            hdu.header.set(keyword,self.wcs[i][keyword])
        # Make image array. The transpose would be necessary so that ds9 displays
        # the image correctly, if we hadn't already done it in read_in_fits_data()
        # hdu.data = self.values[i].transpose()
        hdu.data = self.values[i]
        self.hdr.append(hdu.header)
        return hdu

# ----------------------------------------------------------------------------

    def write_out_to_fits(self,i):
        # Make header
        hdu = self.make_header(i)
        # Verify and write to file:
        hdu.verify()
        hdu.writeto(self.output[i])

        return None

# ----------------------------------------------------------------------------
# Interpolating the map to return a single value at a specified point - this
# is the most important method of this class.

    def at(self,x,y,mapfile=None,coordinate_system='world'):
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
            print "Looking up map",mapfile," value at position",x,",",y," in the ",coordinate_system," coordinate system"

        # Get pixel indices of desired point,
        # and also work out other positions for completeness, if verbose:
        if coordinate_system == 'world':
            i,j = self.world2image(x,y,mapfile)
            if vb: print "  - image coordinates of map ",mapfile,":",i,j
        elif coordinate_system == 'physical':
           i,j = self.physical2image(x,y,mapfile)
           if vb: print "  - image coordinates of map ",mapfile,":",i,j
        elif coordinate_system == 'image':
           i = x
           j = y
           if vb:
               x,y = self.image2physical(i,j,mapfile)
               print "  - physical coordinates of map",mapfile,":",x,y,"(radians)"

        if vb:
            a,d = self.image2world(i,j,mapfile)
            print "  - approximate world coordinates of map",mapfile,":",a,d,"(degrees)"
            print "  - ds9 image coordinates:",i+1,j+1

        # Now look up correct value, doing some bilinear interpolation:

        value = self.lookup(i,j,mapfile)
        if vb: print "  Value of map ",mapfile," at (",x,", ",y,")= ",value

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
        #if a < 0.0: a += 360.0 : We are using WCS with negative RAs in degrees, now
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
        a = -np.rad2deg(x) - self.map_x[mapfile]*self.field[mapfile]
        #if a < 0.0: a += 360.0 :we are using nRA instead now
        d = np.rad2deg(y) + self.map_y[mapfile]*self.field[mapfile]
        return a,d

    def world2physical(self,a,d,mapfile=0):
        x = -np.deg2rad(a + self.map_x[mapfile]*self.field[mapfile])
        y = np.deg2rad(d - self.map_y[mapfile]*self.field[mapfile])
        return x,y

# ----------------------------------------------------------------------------
# In general, we don't know how to plot this map...

    def plot_setup(self,subplot=None,coords='world'):
        '''
        Plot the convergence as a grayscale image.

        Optional arguments:
            subplot         List of four plot limits [xmin,xmax,ymin,ymax]
            coords          Type of coordinates inputted for the subplot:
                            'pixel', 'physical', or 'world'

        *!NOTE!*: Not a complete plotting method. Ony calculates values common
        to both Kappamap and Shearmap plots
        '''

        if subplot is None:
            # Default subplot is entire image
            ai, di = self.image2world(0,0)
            af, df = self.image2world(self.NX[0],self.NX[0])
            subplot = [ai,af,di,df]

        xi, xf = subplot[0], subplot[1]    # x-limits for subplot
        yi, yf = subplot[2], subplot[3]    # y-limits for subplot

        if coords == 'world':
            # Subplot is already in world coordinates
            Lx = 1.0*abs(xf-xi)    # length of x-axis subplot
            Ly = 1.0*abs(yf-yi)    # length of y-axis subplot

        elif coords == 'physical':
            # Convert subplot bounds to world coordinates
            xi,yi = self.physical2world(xi,yi)
            xf,yf = self.physical2world(xf,yf)
            Lx = 1.0*abs(xf-xi)
            Ly = 1.0*abs(yf-yi)
            subplot = [xi,xf,yi,yf]

        elif coords == 'pixel':
            # Convert subplot bounds to world coordinates
            xi,yi = self.image2world(xi,yi)
            xf,yf = self.image2world(xf,yf)
            Lx = 1.0*abs(xf-xi)
            Ly = 1.0*abs(yf-yi)
            subplot = [xi,xf,yi,yf]

        else:
            raise IOError('Error: Subplot bounds can only be in pixel, physical, or world coordinates.')

        # Convert subplot bounds to pixel values
        pix_xi,pix_yi = self.world2image(xi,yi)
        pix_xf,pix_yf = self.world2image(xf,yf)

        # Pixel length of subplot
        pix_Lx = pix_xf-pix_xi
        pix_Ly = pix_yf-pix_yi

        return pix_xi,pix_xf,pix_yi,pix_yf,Lx,Ly,pix_Lx,pix_Ly,subplot

# ----------------------------------------------------------------------------

    def lookup(self,i,j,mapfile):

        # Weighted mean of 4 neighbouring pixels, as suggested by Stefan.
        # Need to transpose x and y, since the maps are transposed...
        # ix = int(i)
        # iy = int(j)
        # px = i - ix
        # py = j - iy
        ix = int(j)
        iy = int(i)
        px = j - ix
        py = i - iy

        if ((0 <= ix) and (ix < self.NX[mapfile]-1) and (0 <= iy) and (iy < self.NX[mapfile]-1)):
            mean =  self.values[mapfile][ix,iy]    * (1.0-px) * (1.0-py) \
                  + self.values[mapfile][ix+1,iy]  * px       * (1.0-py) \
                  + self.values[mapfile][ix,iy+1]  * (1.0-px) * py      \
                  + self.values[mapfile][ix+1,iy+1] * px * py
        else:
            mean = None

        return mean
