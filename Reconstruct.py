#!/usr/bin/env python
# ======================================================================

import pangloss

import sys,getopt,cPickle

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
        - Code is incomplete.

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
        print pangloss.doubledashedline
        print pangloss.hello
        print pangloss.doubledashedline
        print "Reconstruct: assigning halo mass to various lightcones"
        print "Reconstruct: Taking instructions from",configfile
    else:
        print Reconstruct.__doc__
        return

    # --------------------------------------------------------------------
    # Read in configuration, and extract the ones we need:
    
    experiment = pangloss.Configuration(configfile)

    zd = experiment.parameters['StrongLensRedshift']
    zs = experiment.parameters['SourceRedshift']

    calpickles = []
    Nc = experiment.parameters['NCalibrationLightcones']
    for i in range(Nc):
        calpickles.append(experiment.getLightconePickleName('simulated',pointing=i))
    
    obspickle = experiment.getLightconePickleName('real')
    
    # Ray tracing:
    rtscheme = experiment.parameters['RayTracingScheme']
    
    # SHM relation parameters:
    SHMrelation = experiment.parameters['StellarMass2HaloMassRelation']
    
    # Sampling Pr(kappah|D):
    Ns = experiment.parameters['NRealisations']
    
    # --------------------------------------------------------------------
    # Load in stellar mass to halo relation:
    
    if SHMrelation == 'Behroozi':
        # Make the Behroozi pickle if it doesn't exist?
        # ...
        shmr = pangloss.readPickle(CALIB_DIR+'behroozi.pickle')
    else:
        shmr = pangloss.SHMR(SHMrelation)
    
    # --------------------------------------------------------------------
    # Make redshift grid:
    
    grid = pangloss.Grid(zd,zs,nplanes=100)
    
    # --------------------------------------------------------------------
    # Read in lightcones from pickles:

    calcones = []
    for i in range(Nc):
        calcones.append(pangloss.readPickle(calpickles[i]))

    obscone = pangloss.readPickle(obspickle)

    # --------------------------------------------------------------------
    # Make realisations of each lightcone, and store sample kappah vals:

    for lc in calcones+[obscone]:

        p = pangloss.PDF('kappa_halo')
        # Add gamma1, gamma2
 
        lc.define_system(zd,zs)
        lc.load_grid(grid)

#         for j in range(Ns):
#             lc.mimic_photoz_error(sigma=0.1)
#             lc.snap_to_grid(grid)
#             lc.drawMStars(behI)
#             lc.drawMHalos(behT)
#             lc.drawConcentrations(errors=True)
#             lc.Make_kappas(truncationscale=10)
#             lc.Scale_kappas()
#             kappa_halo[j]=lc.kappa_add_total
#             gamma1_halo[j]=lc.gamma1_add_total
#             gamma2_halo[j]=lc.gamma2_add_total
# 
#         # Take Hilbert ray-traced kappa as "truth":
#         p.truth[0] = lc.kappa_hilbert
#         
# 
#         results = kappa_Hilbert,gamma1_Hilbert,gamma2_Hilbert,kappa_halo,gamma1_halo,gamma2_halo
#         #,kappa_T,gamma1_T,gamma2_T
# 
#         return results

    # --------------------------------------------------------------------

    return

# ======================================================================

if __name__ == '__main__': 
    Reconstruct(sys.argv[1:])

# ======================================================================
