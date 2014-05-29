#!/usr/bin/env python
# ======================================================================

import pangloss

import sys,getopt,cPickle,numpy
import matplotlib.pyplot as plt

from scipy.stats.kde import gaussian_kde


# ======================================================================

def Reconstruct(argv):
    """
    NAME
        Reconstruct.py

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
        Reconstruct.py example.config

    BUGS
        - Code is incomplete.

    AUTHORS
      This file is part of the Pangloss project, distributed under the
      GPL v2, by Tom Collett (IoA) and  Phil Marshall (Oxford). 
      Please cite: Collett et al 2013, http://arxiv.org/abs/1303.6564

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

    # ==============================================================
    # Read in configuration, and extract the ones we need:
    # ==============================================================
    
    experiment = pangloss.Configuration(configfile)

    # Get the experiment name from the configfile name instead?
    EXP_NAME = experiment.parameters['ExperimentName']
    CALIB_DIR = experiment.parameters['CalibrationFolder'][0]
    
    Rc = experiment.parameters['LightconeRadius'] # in arcmin

    zd = experiment.parameters['StrongLensRedshift']
    zs = [1.1, 2.1, 5.7]

    calpickles = []
    Nc = experiment.parameters['NCalibrationLightcones']
    for i in range(Nc):
        calpickles.append(experiment.getLightconePickleName('simulated',pointing=i))
       
    # Ray tracing:
    RTscheme = experiment.parameters['RayTracingScheme']
    
    # Sampling Pr(kappah|D):
    Ns = experiment.parameters['NRealisations']
    
    # Reconstruct calibration lines of sight?
    DoCal = experiment.parameters['ReconstructCalibrations']

    # --------------------------------------------------------------------

    
    # --------------------------------------------------------------------
    # Read in lightcones from pickles:

    calcones = []
    for i in range(Nc):
        calcones.append(pangloss.readPickle(calpickles[i]))

    if DoCal=="False": #must be string type
        calcones=[]
        calpickles=[]

    allcones = calcones
    allconefiles = calpickles 
    
    ndensity_field = 86.0480379792 # in num/arcmin^2    
    

    # ==============================================================
    # Set up kappa_smooth
    # ==============================================================
    #kappa_smooth = 0.188247882068  # BORG
    #kappa_smooth = 0.0716519068853 # Zach

    print pangloss.dashedline
    density = [1.0,0.75,1.25]
    drange = 0.02

    # --------------------------------------------------------------------
    # Select ALL lightcones and find their convergences at each redshift
    
    c = ['r','b','g']       
    
    plt.figure()
    
    for i in range(len(zs)):
        # Make redshift grid:
    
        grid = pangloss.Grid(zd,zs[i],nplanes=100)
        
        pk = []
        pmu =[] 
        lc_density = []   
        for j in range(len(allcones)):        
    
            # Get lightcone, and start PDF for its kappa_halo:
            lc = allcones[j] 
                
            lc_dens = lc.countGalaxies(ndensity = ndensity_field)
            lc_density.append(lc_dens)           
            
            # --------------------------------------------------------------------
            # Calculate mu and kappa for all lightcones
                
            # Redshift scaffolding:
            lc.defineSystem(zd,zs[i])
            lc.loadGrid(grid)
        
            # Figure out data quality etc:
            lc.configureForSurvey(experiment)
        
            if j % 100 == 0 and j !=0:
                print ("Reconstruct: ...on lightcone %i out of %i..." % (j,Nc))
    
                    
            lc.snapToGrid(grid)
                        
            # Draw c from Mhalo:
            lc.drawConcentrations(errors=True)
                        
            # Compute each halo's contribution to the convergence:
            lc.makeKappas(truncationscale=5)
                        
            k_add=lc.combineKappas()
            mu_add=lc.combineMus(weakapprox=False)
                                
            pmu.append([lc.mu_add_total])
            pk.append([lc.kappa_add_total])
                
                        
            x = allconefiles[j]
            
        # Add magnification and convergence to global PDF
        pk = numpy.array(pk)
        pmu = numpy.array(pmu)    
    
            
        # ==============================================================
        # Plot the pdfs
        # ==============================================================
        
        kappa_smooth = numpy.mean(pk)
        mu_smooth = numpy.mean(pmu)
     

        outputfile = "figs/PofKappa_"+EXP_NAME+"_compare_z.png" 
  
        
        par1 = pk[:,0] - kappa_smooth
               
        #par1 = pmu[:,0] - mu_smooth + 1.0
        par1mean = numpy.mean(par1) 
        Nlos = len(par1)
                        
        par1_kde = gaussian_kde(par1)
        x = numpy.linspace(par1.min()-0.2,par1.max()+0.2,3*Nlos)
                                    
        plt.plot(x, par1_kde(x), color=c[i], label=r'$z_s = $'+str(zs[i]))#+r', $\langle \mu \rangle = $%.3f' % par1mean) # distribution function
                        
        
    plt.ticklabel_format(useOffset=False, axis='x')
    
    plt.xlabel(r'$\kappa$')
    plt.ylabel(r"$P(\kappa)$")
                
  #  plt.title(r'PDF for $z_s =$ %.1f' % (zs))
        
    plt.vlines(0.0, 0.0, 35, 'k', linestyle='dashed')
 #   plt.xlim(0.7,1.3)
    plt.xlim(-0.15,0.15)
        
    plt.legend(loc=1)
                        
    plt.savefig(outputfile,dpi=300)
    
    print pangloss.doubledashedline
    print "Reconstruct: saved P(kappa) to",outputfile

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -   
    
    print pangloss.doubledashedline
    return

# ======================================================================

if __name__ == '__main__': 
    Reconstruct(sys.argv[1:])

# ======================================================================
