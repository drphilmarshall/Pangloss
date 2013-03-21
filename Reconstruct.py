#!/usr/bin/env python
# ======================================================================

import pangloss

import getopt,cPickle

# ======================================================================

def Reconstruct(argv):
   """
   NAME
       Reconstruct.py

   PURPOSE
       Read in a lightcone (or list of lightcones) and compute all 
       quantities needed to estimate kappah, the convergence due to
       halos, at the centre of the lightcone. Output is a list of sample
       kappah values drawn from Pr(kappah|D), where D refers to either 
       observed data, or simulated data from a calibration line of
       sight.

   COMMENTS
       The config file contains the list of lightcones to be
       reconstructed, in the form of either a directory or a single
       instance. If a directory is specified, one also gives a number
       ("batchsize") of lightcones to be reconstructed at a time. 
       The number of kappah samples desired must also be given in the
       config file.

   FLAGS
       -h            Print this message [0]

   INPUTS
       configfile    Plain text file containing Pangloss configuration

   OUTPUTS
       stdout        Useful information
       samples       Catalog(s) of samples from Pr(kappah|D)


   EXAMPLE

       Reconstruct.py example.config

   BUGS

   AUTHORS
     This file is part of the Pangloss project, distributed under the
     GPL v2, by Tom Collett (IoA) and  Phil Marshall (Oxford). 
     Please cite: Collett et al 2013, arxiv/###
     
   HISTORY
     2013-03-21 started Collett & Marshall (Oxford)
   """

   # --------------------------------------------------------------------

   try:
      opts, args = getopt.getopt(argv,"h",["help"])
   except getopt.GetoptError, err:
      print str(err) # will print something like "option -a not recognized"
      print Reconstruct.__doc__  # will print the big comment above.
      return

   for o,a in opts:
      if o in ("-h", "--help"):
         print Reconstruct.__doc__
         return
      else:
         assert False, "unhandled option"

   # Check for setup file in array args:
   if len(args) == 1:
      configfile = args[0]
      print "Reconstructing the halo mass in various lightcones, according to instructions in",configfile
   else:
      print Reconstruct.__doc__
      return

   # --------------------------------------------------------------------

   # --------------------------------------------------------------------

   return

# ======================================================================

if __name__ == '__main__':
  Reconstruct(sys.argv[1:])

# ======================================================================
