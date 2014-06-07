#!/usr/bin/env python
# ======================================================================

import pangloss

import sys,getopt,cPickle,numpy
import matplotlib.pyplot as plt

from scipy.stats.kde import gaussian_kde
from math import pi


# ======================================================================

def MakePDF(argv):
    """
    NAME
        MakePDF.py

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
        MakePDF.py example.config

    BUGS
        - Code is incomplete.

    AUTHORS
      This file is part of the Pangloss project, distributed under the
      GPL v2, by Tom Collett (IoA) and  Phil Marshall (Oxford). 
      Please cite: Collett et al 2013, http://arxiv.org/abs/1303.6564

    HISTORY
      2013-03-21 started Collett & Marshall (Oxford)
      2014-04-09 modified for BoRG, C Mason (UCSB)
    """

    # --------------------------------------------------------------------

    try:
       opts, args = getopt.getopt(argv,"h",["help"])
    except getopt.GetoptError, err:
       print str(err) # will print something like "option -a not recognized"
       print MakePDF.__doc__  # will print the big comment above.
       return

    for o,a in opts:
       if o in ("-h", "--help"):
          print MakePDF.__doc__
          return
       else:
          assert False, "unhandled option"

    # Check for setup file in array args:
    if len(args) == 1:
        configfile = args[0]
        print pangloss.doubledashedline
        print pangloss.hello
        print pangloss.doubledashedline
        print "MakePDF: generating magnification PDFs for lightcones of given density"
        print "MakePDF: taking instructions from",configfile
    else:
        print MakePDF.__doc__
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
    zs = experiment.parameters['SourceRedshift']

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
    # Make redshift grid:
    
    grid = pangloss.Grid(zd,zs,nplanes=100)
    
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
    
   # ndensity_field = 4.2110897902 # galaxies brighter than i 22 in num/arcmin^2   hilbert 
    radius_cone = 1.0
    area = pi * radius_cone**2.
    
    ndensity_field = 7.87131907218 # galaxies brighter than F125W 22 in num/arcmin^2  Henriques
    # ==============================================================
    # Test how density scales with LoS
    # ==============================================================
    print "MakePDF: finding the distribution of lightcones with density..."

    lc_dens = []

    for i in range(len(allcones)):
        lc = allcones[i] 
            
        num_galaxies = lc.numberWithin(radius=radius_cone,cut=[16,22],band="F125W",units="arcmin")
        lc_density = num_galaxies/(area * ndensity_field)
            
        lc_dens.append(lc_density)
    
    plt.figure(1)
    
    n, bins, patches = plt.hist(lc_dens, 20, facecolor='k', alpha=0.4)
    plt.setp(patches, 'edgecolor', 'None')
                        
    plt.xlabel(r'Overdensity of halos, $\xi$')
    plt.ylabel('Number of lightcones')
    plt.savefig("figs/"+EXP_NAME+"_lc_density.pdf",dpi=300)
    
    lc_density = numpy.array(lc_dens)
    print 'Mean overdensity =',numpy.mean(lc_density)
    # ==============================================================
    # Set up kappa_smooth
    # ==============================================================
    #kappa_smooth = 0.188247882068  # BORG
    #kappa_smooth = 0.0716519068853 # Zach

    print pangloss.dashedline
#    density = [1.0,0.75,1.25]
    density = [1.0,0.5,1.5]
    drange = 0.05

    # --------------------------------------------------------------------
    # Select ALL lightcones and find their convergences
    
    pk = []
    pmu =[] 

    zmax = zs+0.1
    zbins = numpy.linspace(0.0,zmax,100)
    
    kappa_cont = numpy.zeros((len(allcones), len(zbins)))
    Mh_cont = numpy.zeros((len(allcones), len(zbins))) 
    Mstell_cont = numpy.zeros((len(allcones), len(zbins)))
               
    for j in range(len(allcones)):        
    
        # Get lightcone, and start PDF for its kappa_halo:
        lc = allcones[j] 
               
        # --------------------------------------------------------------------
        # Calculate mu and kappa for all lightcones
            
        # Redshift scaffolding:
        lc.defineSystem(zd,zs)
        lc.loadGrid(grid)
    
        # Figure out data quality etc:
        lc.configureForSurvey(experiment)
    
        if j % 100 == 0 and j !=0:
            print ("MakePDF: ...on lightcone %i out of %i..." % (j,Nc))

                
        lc.snapToGrid(grid)
                    
        # Draw c from Mhalo:
        lc.drawConcentrations(errors=True)
                    
        # Compute each halo's contribution to the convergence:
        lc.makeKappas(truncationscale=5)
                    
        k_add=lc.combineKappas()
        mu_add=lc.combineMus(weakapprox=False)
                               
        pmu.append([lc.mu_add_total])
        pk.append([lc.kappa_add_total])
    
        kappa_cont[j:,] = lc.findContributions('kappa')   
        Mh_cont[j:,] = lc.findContributions('mass') 
        Mstell_cont[j:,] = lc.findContributions('stellarmass') 
            
        # Make a nice visualisation of one of the realisations, in
        # two example cases:
        if j ==0:
            lc.plots('kappa', output=CALIB_DIR+"/example_snapshot_kappa_uncalib.png")
            lc.plots('mu', output=CALIB_DIR+"/example_snapshot_mu_uncalib.png")
        #    lc.plotContributions('kappa', output=CALIB_DIR+"/kappa_contribution_uncalib.png")
        #   lc.plotContributions('mass', output=CALIB_DIR+"/Mh_contribution_uncalib.png")
        #   lc.plotContributions('mu', output=CALIB_DIR+"/mu_contribution_uncalib.png")
        #    lc.plotContributions('stellarmass', output=CALIB_DIR+"/Mstell_contribution_uncalib.png")
    
                    
        x = allconefiles[j]
        
    # Add magnification and convergence to global PDF
    pk = numpy.array(pk)
    pmu = numpy.array(pmu)    

    mean_kappa_cont = numpy.mean(kappa_cont, axis=0)
    mean_Mh_cont = numpy.mean(Mh_cont, axis=0)
    mean_Mstell_cont = numpy.mean(Mstell_cont, axis=0)
    
    plt.clf()
    plt.figure()
    
    ax1 = plt.subplot(2,1,1)
    ax1.plot(zbins, mean_Mh_cont)
    ax1.set_ylabel(r'Cumulative Sum of $M_h$')
    ax2 = plt.subplot(2,1,2)
    ax2.plot(zbins, mean_kappa_cont)
    ax2.set_ylabel(r'Cumulative Sum of $\kappa_h$') 
    plt.tight_layout()
    plt.xlabel('Redshift, z')
    plt.savefig("figs/contribution_z.png",dpi=300)

    # Remove mu outliers

    
    # ==============================================================
    # Plot the pdfs
    # ==============================================================
    
    kappa_smooth = numpy.mean(pk)
    mu_smooth = numpy.mean(pmu)
    print kappa_smooth, mu_smooth

    c = ['r','b','g']
    
    pdf = [{'param':'Kappa', 'name':r'$\kappa$', 'lc':pk[:,0], 'smooth':kappa_smooth, 'mean':0.0, 'height':60, 'min':-0.05, 'max':0.2},
            {'param':'Mu', 'name':r'$\mu$', 'lc':pmu[:,0], 'smooth':mu_smooth, 'mean':1.0, 'height':30, 'min':0.8, 'max':1.5}]
       
    # Plotting convergence and magnification
    for k in range(len(pdf)):

        plt.clf()
        var = pdf[k]
        
        full_pdf = var['lc']
        
        mask = numpy.where((full_pdf > -2.*var['smooth']) & (full_pdf < 2.*var['smooth'])) 
        par = full_pdf[mask]
        smooth_new = numpy.mean(par) 
        full_pdf = par - smooth_new + var['mean']
        par_mean = numpy.mean(par)
        
        outputfile = "figs/Pof"+var['param']+"_"+EXP_NAME+"_test.png" 
        name = var['name']     
        
        plt.figure()
        
        all_LOS = full_pdf - var['smooth'] + var['mean']
        all_kde = gaussian_kde(all_LOS)
        x = numpy.linspace(all_LOS.min()-0.2,all_LOS.max()+0.2,3*Nc)
        plt.plot(x, all_kde(x), color='k', label=r'All LOS, $\langle$'+name+r'$\rangle = $%.3f' % var['smooth']) # distribution function
            
        print pangloss.dashedline
        print "MakePDF: constructing PDF for", var['param'],"..."
            
   
                
        # Select only lightcones within certain number density limits
        for i in range(len(density)):
            sub_pdf = [] 
            
            for j in range(len(full_pdf)):
                if density[i] - drange <= lc_density[mask][j] <= density[i] + drange:
                    sub_pdf.append(full_pdf[j])
            
            Nsub = len(sub_pdf)
            print Nsub
        
            sub_pdf = numpy.array(sub_pdf)
            
            print "MakePDF: sampling %i LoS with number density ~ %.2f the average" % (Nsub, density[i])   
            
            par1 = sub_pdf #- var['smooth'] + var['mean']
            par1mean = numpy.mean(par1) 
            Nlos = len(par1)
            
          #  n, bins, patches = plt.hist(par1, 20, facecolor=c[i], normed=True,alpha=0.4, label=r'$z_s = $'+str(zs))
          #  plt.setp(patches, 'edgecolor', 'None')
                        
            par1_kde = gaussian_kde(par1)
            x = numpy.linspace(par1.min()-0.2,par1.max()+0.2,3*Nlos)
                                    
                        
            plt.plot(x, par1_kde(x), color=c[i], label=r'$\xi = $'+str(density[i])+r', $\langle$'+name+r'$\rangle = $%.3f' % par1mean) # distribution function
                        
        
        plt.ticklabel_format(useOffset=False, axis='x')
    
        plt.xlabel(name)
        plt.ylabel(r"P(%s)"%name)
                
        plt.title(r'PDF for $z_s =$ %.1f' % (zs))
        
     #   plt.vlines(var['mean'], 0.0,var['height'], 'k', linestyle='dashed')
     #   plt.xlim(var['min'],var['max'])
        
        plt.legend(loc=1)
                        
        plt.savefig(outputfile,dpi=300)
    
    print pangloss.doubledashedline
    print "MakePDF: saved PDFs to",outputfile

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -   
    
    print pangloss.doubledashedline
    return

# ======================================================================

if __name__ == '__main__': 
    MakePDF(sys.argv[1:])

# ======================================================================
