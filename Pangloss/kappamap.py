'''
This file is part of the Pangloss project.
Copyright 2012 Tom Collett (IoA) and Phil Marshall (Oxford).

to-do
-----
- check that kappa map is being read in the right way around...
'''

import numpy, os, string
import struct, pyfits

arcmin2rad = (1.0/60.0)*numpy.pi/180.0
rad2arcmin = 1.0/arcmin2rad
deg2rad = numpy.pi/180.0
rad2deg = 1.0/deg2rad

vb = True

# ============================================================================
    
# Read in, write out, look up convergence kappa in map from Stefan Hilbert.

class kappamap:

   def __init__(self,kappafile,FITS=True):

      self.name = 'Convergence map kappa from Millenium Simulation, zs = 1.6'
      self.input = kappafile

      # Read in data from file:
      
      if FITS:
         # Read in FITS image, extract wcs:
         if vb: print "Reading in kappa map from file..."
         self.read_in_fits_data()
         if vb: print "... done"
       
      else:

         # Initialise some necessary WCS parameters:
         self.field = 4.0 # degrees
         self.NX = 4096
         self.PIXSCALE = self.field/(1.0*self.NX) # degrees
         self.setwcs()

         # Read in binary data, to self.values:
         if vb: print "Reading in kappa map from file..."
         self.read_in_binary_data()
         if vb: print "... done"

         # If it doesn't already exist, output the map to FITS file:
         pieces = string.split(self.input,'.')
         self.output = string.join(pieces[0:len(pieces)-1],'.')+'.fits'
         if os.path.exists(self.output): 
           if vb: print "FITS version already exists: ",self.output
         else:  
           self.write_out_to_fits()
           if vb: print "Written kappa map to "+self.output
         # This should probably not be in __init__ but hopefully it only gets run once.  
      
      # Python array storage: need to take transpose of the map!
      self.values = self.values.transpose()
      
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

   # Note: all hdr keywords are copied but only wcs ones are used.
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
      # Make image array:
      hdu.data = self.values
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
      if vb: 
         print "  Value of kappa = ",kappa
      
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

# Self-test: read in map from Stefan, and look up some convergence values.
   
# Read in map (and write out as FITS if it doesn't exist):
#    kappafile = "../../data/GGL_los_8_0_0_N_4096_ang_4_rays_to_plane_37_f.kappa"
#    convergence = kappamap(kappafile,FITS=False)

# Read in map from FITS file:

   kappafile = "../../data/GGL_los_8_0_0_N_4096_ang_4_rays_to_plane_37_f.fits"
   convergence = kappamap(kappafile)
  
# Look up values at some points:
   
   i = 2184 ; j = 2263
   kappa = convergence.at(i,j,coordinate_system='image')
   print "Compare with expected value: 1.17721"

   i = 2263 ; j = 2184
   kappa = convergence.at(i,j,coordinate_system='image')
   print "Compare with expected value: 0.0169251"

   i = 0 ; j = 0
   kappa = convergence.at(i,j,coordinate_system='image')
   print "Compare with expected value: -0.0222247"

# Looking up kappa value at position 2184 , 2263  in the image coordinate system
#   - physical coordinates: 0.00232653752829 0.00367303177544 (radians)
#   - approximate world coordinates: 359.867675781 0.21044921875 (degrees)
#   - ds9 image coordinates: 2185 2264
#   Value of kappa =  1.17720580101
# Compare with expected value: 1.17721
#  
# Looking up kappa value at position 2263 , 2184  in the image coordinate system
#   - physical coordinates: 0.00367303177544 0.00232653752829 (radians)
#   - approximate world coordinates: 359.790527344 0.13330078125 (degrees)
#   - ds9 image coordinates: 2264 2185
#   Value of kappa =  0.0169250555336
# Compare with expected value: 0.0169251
#  
# Looking up kappa value at position 0 , 0  in the image coordinate system
#   - physical coordinates: -0.0348980629244 -0.0348980629244 (radians)
#   - approximate world coordinates: 2.00048828125 -1.99951171875 (degrees)
#   - ds9 image coordinates: 1 1
#   Value of kappa =  -0.0222246814519
# Compare with expected value: -0.0222247

   x = 0.00232653752829; y = 0.00367303177544
   kappa = convergence.at(x,y,coordinate_system='physical')
   print "Compare with expected value: 1.17721"

# Looking up kappa value at position 0.00232653752829 , 0.00367303177544  in the physical coordinate system
#   - image coordinates: 2184.0 2263.0
#   - approximate world coordinates: 359.867675781 0.21044921875 (degrees)
#   - ds9 image coordinates: 2185.0 2264.0
#   Value of kappa =  1.1772058009
# Compare with expected value: 1.17721

# OK, internally all looks good. 

# ============================================================================
