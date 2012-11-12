#!/usr/bin/env python
# ======================================================================

import Pangloss

import sys,getopt,atpy,pyfits
import numpy

# ======================================================================

def hilbert2fits(argv):
   """
   NAME
     hilbert2fits.py

   PURPOSE
     Convert input maps in Hilbert's format to FITS. 

   COMMENTS

   FLAGS
     -h            Print this message [0]
     -v            Verbose operation

   INPUTS
     infiles       Files containing true kappa (or gamma) values.

   OPTIONAL INPUTS

   OUTPUTS
     stdout        Useful information
     outfiles      Output files in FITS format

   EXAMPLES

     hilbert2fits.py  x.*
     
     Outputs are x.kappa.fits, x.gammma_1.fits, x.gamma_2.fits

   BUGS

   AUTHORS
     This file is part of the Pangloss project.
     Copyright 2012 Tom Collett (IoA) and Phil Marshall (Oxford).
     
   HISTORY
     2012-06-15 started Marshall (Oxford)
   """

   # --------------------------------------------------------------------

   try:
      opts, args = getopt.getopt(argv,"hv",["help","verbose"])
   except getopt.GetoptError, err:
      print str(err) # will print something like "option -a not recognized"
      print hilbert2fits.__doc__  # will print the big comment above.
      return

   vb = False
  
   for o,a in opts:
      if o in ("-v", "--verbose"):
         vb = True
      elif o in ("-h", "--help"):
         print hilbert2fits.__doc__
         return
      else:
         assert False, "unhandled option"

   # Files to be converted:
   infiles = args
   if vb: print "Converting",len(infiles),"files to FITS format"
 
   # --------------------------------------------------------------------

   for infile in infiles:
     # Read in map and convert (done automatically):    
     map = Pangloss.kappamap(infile,FITS=False,vb=True)

# ======================================================================

if __name__ == '__main__':
  hilbert2fits(sys.argv[1:])

# ======================================================================
