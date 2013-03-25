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
        print "Reconstruct: taking instructions from",configfile
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
    RTscheme = experiment.parameters['RayTracingScheme']
    
    # SHM relation parameters:
    SHMrelation = experiment.parameters['StellarMass2HaloMassRelation']
    CALIB_DIR = experiment.parameters['CalibrationFolder'][0]
    SHMfile = CALIB_DIR+'/'+SHMrelation+'.pickle'
    
    # Halo mass function data:
    HMFfile = experiment.parameters['HMFfile'][0]
    
    # Photo-zs:
    zperr = experiment.parameters['PhotozError']
    
    # Stellar mass observations:
    Mserr = experiment.parameters['MstarError']
    Mserr = experiment.parameters['MstarError']
    # Sampling Pr(kappah|D):
    Ns = experiment.parameters['NRealisations']
    
    # --------------------------------------------------------------------
    # Load in stellar mass to halo relation, or make a new one:
    
    try:
        shmr = pangloss.readPickle(SHMfile)
    except IOError:
        print "Reconstruct: generating the stellar mass to halo mass grid."
        print "Reconstruct: this may take a moment..."
        shmr = pangloss.SHMR(method=SHMrelation)
        shmr.makeHaloMassFunction(HMFfile)
        shmr.makeCDFs()
        pangloss.writePickle(shmr,SHMfile)
        print "Reconstruct: SHMR saved to "+SHMfile
    
    # --------------------------------------------------------------------
    # Make redshift grid:
    
    grid = pangloss.Grid(zd,zs,nplanes=100)
    
    # --------------------------------------------------------------------
    # Read in lightcones from pickles:

    calcones = []
    for i in range(Nc):
        calcones.append(pangloss.readPickle(calpickles[i]))
    obscone = pangloss.readPickle(obspickle)
    
    allcones = calcones+[obscone]
    allconefiles = calpickles+[obspickle]

    # --------------------------------------------------------------------
    # Make realisations of each lightcone, and store sample kappah vals:

    for i in range(len(allcones)):

        lc = allcones[i]
        print "Reconstruct: drawing %i samples from Pr(kappah|D)" % (Ns)
        print "Reconstruct: given data in "+allconefiles[i]

        p = pangloss.PDF('kappa_halo')
        # coming soon: gamma1, gamma2...
 
        lc.defineSystem(zd,zs)
        lc.loadGrid(grid)

        # Draw Ns sample realisations of this lightcone, and hence
        # accumulate samples from Pr(kappah|D):
        for j in range(Ns):
            lc.mimicPhotozError(sigma=zperr)
            lc.snapToGrid(grid)
            
            # Simulated lightcones need Mstars drawing from their Mhalos
            if lc.flavor == 'simulated': lc.drawMstars(shmr)
            lc.mimicMstarError(sigmaP=MserrP,sigmaS=MserrS)

            lc.drawMhalos(shmr)
            lc.drawConcentrations(errors=True)

            lc.makeKappas(truncationscale=10)
            lc.combineKappas()
            
            if RTscheme == 'sum':
                p.append(lc.kappa_add_total)
            elif RTscheme == 'keeton':
                p.append(lc.kappa_keeton)
            else:
                raise "Unknown ray-tracing scheme: "+RTscheme
            # also lc.gamma1_add_total, lc.gamma2_add_total

        # Take Hilbert ray-traced kappa as "truth":
        p.truth[0] = lc.kappa_hilbert
        
        # Pickle this lightcone's PDF:
        x = allconefiles[i]
        pfile = x.split('.')[0]+"_PofKappah.pickle"
        pangloss.writePickle(p,pfile)
        print "Reconstruct: Pr(kappah|D) saved to "+pfile

    # --------------------------------------------------------------------

    return

# ======================================================================

if __name__ == '__main__': 
    Reconstruct(sys.argv[1:])

# ======================================================================
