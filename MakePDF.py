#!/usr/bin/env python
# ======================================================================

import pangloss

import sys,getopt,cPickle,numpy,glob
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
   
    mag = experiment.parameters['MagName']
    
    Ncats = 21
     
    borg = numpy.genfromtxt('data/borg_overdensity.txt', skip_header=1)

    borg_field = borg[:,0]
    borg_area = borg[:,1]
    borg_overdensity = borg[:,2]
    
    Rc = experiment.parameters['LightconeRadius'] # in arcmin

    zd = experiment.parameters['StrongLensRedshift']
    zs = experiment.parameters['SourceRedshift']

    # --------------------------------------------------------------------    
    # Load the lightcone pickles
    
    calpickles = []
    Nc = experiment.parameters['NCalibrationLightcones'] * Ncats       ### should be 24

    paths = '%s/*_lightcone.pickle' % (CALIB_DIR)
    found = glob.glob(paths)
    if len(found) > 0: calpickles = found
            
    # Ray tracing:
    RTscheme = experiment.parameters['RayTracingScheme']
    
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

    # ==============================================================    
    # Test how density scales with LoS
    # ==============================================================
    print "MakePDF: finding the distribution of lightcones with density..."
   
    area = pi * Rc**2.
  # ndensity_field = 4.2110897902 # galaxies brighter than i 22 in num/arcmin^2   hilbert     
    ndensity_field = 7.87131907218 # galaxies brighter than F125W 22 in num/arcmin^2  Henriques

    lc_dens = []
    
    # --------------------------------------------------------------------
    # Find the overdensity of lightcones cut at m<22 in F125W
    total_galaxies = 0
    lc_galaxies = []
    
    # Sort into lightcones for each field
    for i in range(len(allcones)): 
        lc = allcones[i]   
 #       num_galaxies = lc.numberWithin(radius=Rc,cut=[16,22],band="F125W",units="arcmin")
        num_galaxies = lc.numberWithin(radius=Rc,cut=[16,22],band=mag,units="arcmin")
        lc_galaxies.append(num_galaxies)   
        # Add to the total number of galaxies
        total_galaxies += num_galaxies

    total_density = total_galaxies/Nc
    lc_dens = [Nc * float(x)/float(total_galaxies) for x in lc_galaxies]

    # --------------------------------------------------------------------    
    # Plot histogram of overdensities

    plt.figure(1)
    
    n, bins, patches = plt.hist(lc_dens, 20, facecolor='r', alpha=0.4)
                        
    plt.xlabel(r'Overdensity of halos, $\xi$')
    plt.ylabel('Number of lightcones')
    plt.savefig("figs/"+EXP_NAME+"_lc_density.pdf",dpi=300)
    
    lc_density = numpy.array(lc_dens)
    print 'Mean overdensity in all fields =',numpy.mean(lc_density)
    print pangloss.dashedline

    # ==============================================================
    # Sample all lightcones to make the pdfs
    # ==============================================================

    # --------------------------------------------------------------
    # Set up overdensity range

#    density = borg_overdensity
    density = [1.0,0.75,1.25]
#    density = [1.0,0.5,1.5]

    drange = 0.025   # ~0.02 is what Zach used

    # --------------------------------------------------------------------
    # Find contribution to total kappa and mass at redshift intervals
    zmax = zs+0.1
    zbins = numpy.linspace(0.0,zmax,100)
    
    kappa_cont = numpy.zeros((len(allcones), len(zbins)))
    Mh_cont = numpy.zeros((len(allcones), len(zbins))) 
    Mstell_cont = numpy.zeros((len(allcones), len(zbins)))
               
    # --------------------------------------------------------------------
    # Select ALL lightcones and find their convergences
    
    pk = []
    pmu =[] 
    
    for j in range(len(allcones)):        
    
        # Get lightcone
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
                               
        # Add magnification and convergence to global PDF
        pmu.append([lc.mu_add_total])
        pk.append([lc.kappa_add_total])
    
        kappa_cont[j:,] = lc.findContributions('kappa')   
        Mh_cont[j:,] = lc.findContributions('mass') 
        Mstell_cont[j:,] = lc.findContributions('stellarmass') 
            
        # Make a nice visualisation of one of the lightcones
        if j ==0:
            lc.plots('kappa', output=CALIB_DIR+"/example_snapshot_kappa_uncalib.png")
            lc.plots('mu', output=CALIB_DIR+"/example_snapshot_mu_uncalib.png")
        
        # Name the lightcone       
        x = allconefiles[j]
        
    pk = numpy.array(pk)
    pmu = numpy.array(pmu)    

    # --------------------------------------------------------------------
    # Plot contributions to total kappa and mass at redshifts
    
    mean_kappa_cont = numpy.mean(kappa_cont, axis=0)
    mean_Mh_cont = numpy.mean(Mh_cont, axis=0)
    mean_Mstell_cont = numpy.mean(Mstell_cont, axis=0)
    
    plt.clf()
    plt.figure()
    
    ax1 = plt.subplot(2,1,1)
    ax1.plot(zbins, mean_Mh_cont)
    ax1.set_ylabel(r'Cumulative Sum of $M_h/M_{\odot}$')
    
    ax2 = plt.subplot(2,1,2)
    ax2.plot(zbins, mean_kappa_cont)
    ax2.set_ylabel(r'Cumulative Sum of $\kappa_h$') 
    
    plt.tight_layout()
    plt.xlabel('Redshift, z')
    plt.savefig("figs/contribution_z.png",dpi=300)

    # ==============================================================
    # Plot the pdfs
    # ==============================================================
    
    # --------------------------------------------------------------------
    # Calculate the smooth components
    kappa_smooth = numpy.mean(pk)
    mu_smooth = numpy.mean(pmu)
    print 'Uncalibrated: <kappa> =',kappa_smooth, '<mu> =',mu_smooth

    # --------------------------------------------------------------------
    # Plotting convergence and magnification
    c = ['r','b','g']
    
    pdf = [{'param':'Kappa', 'name':r'$\kappa$', 'lc':pk[:,0], 'smooth':kappa_smooth, 'mean':0.0, 'height':60, 'min':-0.05, 'max':0.2},
            {'param':'Mu', 'name':r'$\mu$', 'lc':pmu[:,0], 'smooth':mu_smooth, 'mean':1.0, 'height':30, 'min':0.7, 'max':1.8}]

    # =====================================================================
    # For only <=3 values of density
    # =====================================================================    
    if len(density) < 4:
        for k in range(len(pdf)):
      
            var = pdf[k]
            full_pdf = var['lc'] - var['smooth'] + var['mean']
            name = var['name']     
           # print 'MakePDF: Old min and max:', full_pdf.min(), full_pdf.max()
            
            # --------------------------------------------------------------------
            # Remove outliers --  only for v high overdensity?
            """
            mask = numpy.where((full_pdf > -1.) & (full_pdf < 2.)) 
            par = full_pdf[mask]
            
            print '         Removing Outliers...'
            print '         New min and max:', par.min(), par.max()
            
            # Recalculate means and the calibrated pdf
            smooth_new = numpy.mean(par) 
            full_pdf = par - smooth_new + var['mean']
            par_mean = numpy.mean(par)
            print 'New mean (this should be',var['mean'],'):', par_mean         
            """
            # --------------------------------------------------------------------
            # Plot PDF for all lines of sight
            
            print pangloss.dashedline
            print "MakePDF: constructing PDF for", var['param'],"..."  
            
            plt.clf()
            plt.figure()
            
            all_LOS = full_pdf

            # Histogram
            n, bins, patches = plt.hist(full_pdf, 20, facecolor='k', normed=True,alpha=0.2)
            plt.setp(patches, 'edgecolor', 'None')            

            # Gaussian KDE            
            all_kde = gaussian_kde(all_LOS)
            x = numpy.linspace(all_LOS.min()-0.2,all_LOS.max()-0.2,3*Nc)
            plt.plot(x, all_kde(x), color='k', label=r'All LOS, $\langle$'+name+r'(uncalib)$\rangle = $%.3f' % var['smooth']) # distribution function
                    
            # --------------------------------------------------------------------
            # Select only lightcones within certain number density limits
            # Need to mask out the outliers
            
            for i in range(len(density)):
                sub_pdf = [] 
                
                for j in range(len(full_pdf)):
                    if density[i] - drange <= lc_density[j] <= density[i] + drange:
    #                if density[i] - drange <= lc_density[mask][j] <= density[i] + drange:
                        sub_pdf.append(full_pdf[j])
                
                sub_pdf = numpy.array(sub_pdf)
                par1 = sub_pdf
                par1mean = numpy.mean(par1) 
                Nlos = len(par1)
                            
                print "MakePDF: sampling %i LoS with number density ~ %.2f the average" % (Nlos, density[i])   
    
                # --------------------------------------------------------------------
                # Plot PDFs
                
                # Histogram
                n, bins, patches = plt.hist(par1, 20, facecolor=c[i], normed=True, alpha=0.2)
                plt.setp(patches, 'edgecolor', 'None')
                            
                # Gaussian KDE            
                par1_kde = gaussian_kde(par1)
                x = numpy.linspace(par1.min()-0.2,par1.max()+0.2,3*Nlos)                        
                plt.plot(x, par1_kde(x), color=c[i], label=r'$\xi = $'+str(density[i])+r', $\langle$'+name+r'$\rangle = $%.3f' % par1mean) # distribution function
                            
            plt.xlim(var['min'], var['max'])
            plt.ticklabel_format(useOffset=False, axis='x')

            plt.axvline(x=var['mean'], color='k', linestyle='dashed')
        
            plt.xlabel(name)
            plt.ylabel(r"P(%s)"%name)
            plt.title(r'PDF for $z_s =$ %.1f' % (zs))
                        
            plt.legend(loc=1)
            
            outputfile = "figs/Pof"+var['param']+"_"+EXP_NAME+"_test.png"                 
            plt.savefig(outputfile,dpi=300)
            print "MakePDF: saved PDFs to",outputfile

    # =====================================================================
    # For many values of overdensity e.g. BORG
    # =====================================================================    

    else:
        for k in range(len(pdf)):
      
            var = pdf[k]
            full_pdf = var['lc']  - var['smooth'] + var['mean']
            name = var['name']     
         #   print 'MakePDF: Old min and max:', full_pdf.min(), full_pdf.max()
            
            # --------------------------------------------------------------------
            # Remove outliers  !!! I think this is only important for v. high overdensity
            """
            mask = numpy.where((full_pdf > -1.) & (full_pdf < 2.)) 
            par = full_pdf[mask]
            
            print '         Removing Outliers...'
            print '         New min and max:', par.min(), par.max()
            
            # Recalculate means and the calibrated pdf
            smooth_new = numpy.mean(par) 
            full_pdf = par - smooth_new + var['mean']
            par_mean = numpy.mean(par)
            print 'New mean (this should be',var['mean'],'):', par_mean        
            """
            # --------------------------------------------------------------------
            # Plot all lines of sight
            
            print pangloss.dashedline
            print "MakePDF: constructing PDF for", var['param'],"..."  
            
            plt.clf()
            plt.figure()
            
            all_LOS = full_pdf

            # Histogram
            n, bins, patches = plt.hist(full_pdf, 20, facecolor='k', normed=True, alpha=0.2)
            plt.setp(patches, 'edgecolor', 'None')            

            # Gaussian KDE            
            all_kde = gaussian_kde(all_LOS)
            x = numpy.linspace(all_LOS.min()-0.2,all_LOS.max()+0.2,3*Nc)
            plt.plot(x, all_kde(x), color='k', label=r'All LOS, $\langle$'+name+r'(uncalib)$\rangle = $%.3f' % var['smooth']) # distribution function
                
            plt.legend(loc=1)
            
            outputfile = "figs/Pof"+var['param']+"_"+EXP_NAME+"_allLOS.png"                 
            plt.savefig(outputfile,dpi=300)
                                                      
            # --------------------------------------------------------------------
            # Select only lightcones within certain number density limits
            # Need to mask out the outliers
            
            for i in range(len(density)):
                sub_pdf = [] 
                
                for j in range(len(full_pdf)):
                    if density[i] - drange <= lc_density[j] <= density[i] + drange:
    #                if density[i] - drange <= lc_density[mask][j] <= density[i] + drange:
                        sub_pdf.append(full_pdf[j])
                
                sub_pdf = numpy.array(sub_pdf)
                par1 = sub_pdf
                par1mean = numpy.mean(par1) 
                Nlos = len(par1)
                            
                print "MakePDF: %s - sampling %i LoS with number density ~ %.2f the average" % (borg_field[j], Nlos, density[i])   
    
                # --------------------------------------------------------------------
                # Write PDFs to pickles
                    
                x = borg_field[j]
                pfile = x+"_Pof"+var['param']+".pickle"
                pangloss.writePickle(par1,pfile)
    
                # --------------------------------------------------------------------
                # Plot PDFs
                
                # Histogram
                n, bins, patches = plt.hist(par1, 20, facecolor='r', normed=True, alpha=0.3)
                plt.setp(patches, 'edgecolor', 'None')
                            
                # Gaussian KDE               
                par1_kde = gaussian_kde(par1)
                x = numpy.linspace(par1.min()-0.2,par1.max()+0.2,3*Nlos)                                        
                plt.plot(x, par1_kde(x), color='r', label=r'$\xi = $'+str(density[i])+r', $\langle$'+name+r'$\rangle = $%.3f' % par1mean) # distribution function
                
                plt.axvline(x=var['mean'], color='k', linestyle='dashed')
                
                plt.xlim(var['min'], var['max'])

                plt.ticklabel_format(useOffset=False, axis='x')
            
                plt.xlabel(name)
                plt.ylabel(r"P(%s)"%name)
                plt.title(r'%s: PDF for $z_s =$ %.1f' % (borg_field[j], zs))
                
                plt.legend(loc=1)
                
                # Make pdf for each field            
                outputfile = "figs/borg/"+str(borg_field[j])+"_Pof"+var['param']+".pdf"                 
                plt.savefig(outputfile,dpi=300)                 

            print "MakePDF: saved PDFs to",outputfile 
            
    print pangloss.doubledashedline

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -   
    
    print pangloss.doubledashedline
    return

# ======================================================================

if __name__ == '__main__': 
    MakePDF(sys.argv[1:])

# ======================================================================
