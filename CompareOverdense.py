#!/usr/bin/env python
# ======================================================================
from matplotlib import rc_file
rc_file('matplotlibrc')

import pangloss

import sys,getopt,cPickle,numpy,glob
import matplotlib.pyplot as plt

from math import pi
from astropy.io import ascii
import os


# ======================================================================

def CompareOverdense(argv):
    """
    NAME
        CompareOverdense.py

    PURPOSE
        Read in a simulation lightcone (or list of lightcones) and compute all 
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
        CompareOverdense.py example.config

    BUGS
        - Code is incomplete.

    AUTHORS
      This file is part of the Pangloss project, distributed under the
      GPL v2, by Tom Collett (IoA) and  Phil Marshall (Oxford). 
      Please cite: Collett et al 2013, http://arxiv.org/abs/1303.6564

    HISTORY
      2013-03-21 started Collett & Marshall (Oxford)
      2014-04-09 modified for BoRG, C Mason (UCSB)
      2014-09-15 C Mason - compare overdensities
    """

    # --------------------------------------------------------------------

    try:
       opts, args = getopt.getopt(argv,"h",["help"])
    except getopt.GetoptError, err:
       print str(err) # will print something like "option -a not recognized"
       print CompareOverdense.__doc__  # will print the big comment above.
       return

    for o,a in opts:
       if o in ("-h", "--help"):
          print CompareOverdense.__doc__
          return
       else:
          assert False, "unhandled option"

    # Check for setup file in array args:
    if len(args) == 1:
        configfile = args[0]
        print pangloss.doubledashedline
        print pangloss.hello
        print pangloss.doubledashedline
        print "CompareOverdense: generating magnification PDFs for lightcones of given density"
        print "CompareOverdense: taking instructions from",configfile
    else:
        print CompareOverdense.__doc__
        return

    # ==============================================================
    # Read in configuration, and extract the ones we need:
    # ==============================================================
    
    experiment = pangloss.Configuration(configfile)

    # Get the experiment name from the configfile name instead?
    EXP_NAME = experiment.parameters['ExperimentName']
    CALIB_DIR = experiment.parameters['CalibrationFolder'][0]
   
    mag = experiment.parameters['MagName']
    
    Ncats = 21
    
    Rc = experiment.parameters['LightconeRadius'] # in arcmin

    # --------------------------------------------------------------------    
    # Load the lightcone pickles
    
    calpickles = []
    Nc = experiment.parameters['NCalibrationLightcones'] * Ncats       ### should be 24

    paths = '%s/*_lightcone.pickle' % (CALIB_DIR)
    found = glob.glob(paths)
    if len(found) > 0: calpickles = found

    print "CompareOverdense: found the lightcones..."
               
    # Reconstruct calibration lines of sight?
    DoCal = experiment.parameters['ReconstructCalibrations']
    
    # ==============================================================    
    # Test how density scales with LoS
    # ==============================================================
    print "CompareOverdense: finding the distribution of lightcones with density..."
   
    area = pi * Rc**2.
  # ndensity_field = 4.2110897902 # galaxies brighter than i 22 in num/arcmin^2   hilbert     
    ndensity_field = 7.87131907218 # galaxies brighter than F125W 22 in num/arcmin^2  Henriques

    lc_dens = []
    
    # --------------------------------------------------------------------
    # Find the overdensity of lightcones cut at m<22 in F125W

    total_galaxies = 0
    lc_galaxies = []
   
    # Sort into lightcones for each field
    for i in xrange(Nc): 
        if i % 1000 == 0 and i !=0:
            print ("Drill: ...on cone %i out of %i..." % (i,Nc))
        
        lc = pangloss.readPickle(calpickles[i])  
        num_galaxies = lc.numberWithin(radius=Rc,cut=[16,22],band=mag,units="arcmin")
        lc_galaxies.append(num_galaxies)   
        # Add to the total number of galaxies
        total_galaxies += num_galaxies
        del lc
    
    del calpickles
    
    lc_dens = [Nc * float(x)/float(total_galaxies) for x in lc_galaxies]
    print 'Mean overdensity in lightcones =',numpy.mean(numpy.array(lc_dens))
    
    # --------------------------------------------------------------------    
    # Find over density in BoRG fields
    borg_dens = numpy.genfromtxt('data/borg_overdensity.txt')[:,2]

    print 'Mean overdensity in BoRG =',numpy.mean(borg_dens)
    
    very_overdense = numpy.where(borg_dens <= 2.5)
    borg_dens = borg_dens[very_overdense]

    # --------------------------------------------------------------------        
    # Plot histogram of overdensities
  
    plt.figure(1)
    
    n_lc, bins_lc, patches_lc = plt.hist(lc_dens, 20, facecolor='DarkOrange', alpha=0.4, normed=True, label='Henriques et al. (2012)')
    n_borg, bins_borg, patches_borg = plt.hist(borg_dens, 10, facecolor='DarkBlue', alpha=0.4, normed=True, label='BoRG fields')
                        
    plt.xlabel(r'$\xi$')
    plt.ylabel(r'$P(\xi)$')
    plt.xlim(0,2.5)
    
    print pangloss.dashedline
    del lc_dens

    plt.legend(loc=1)  
    
    savedfile = "../figs/overdensity_compare.pdf"
    plt.savefig(savedfile,dpi=300)
    print "Plot saved as "+os.getcwd()+"/"+savedfile
    
    return

# ======================================================================

if __name__ == '__main__': 
    CompareOverdense(sys.argv[1:])

# ======================================================================
