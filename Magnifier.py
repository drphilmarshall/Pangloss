#!/usr/bin/env python
# ======================================================================
from matplotlib import rc_file
rc_file('matplotlibrc')

import pangloss

import sys,getopt,cPickle,numpy,glob
import matplotlib.pyplot as plt

from scipy.stats.kde import gaussian_kde
from math import pi
from astropy.io import ascii


# ======================================================================

def Magnifier(argv):
    """
    NAME
        Magnifier.py

    PURPOSE
        Read in a simulation lightcone (or list of lightcones) and compute all 
        quantities needed to estimate muh, the convergence due to
        halos, at the centre of the lightcone.
        Output is a list of muh values drawn from Pr(muh|D), 
        where D refers to either observed data, or simulated data from 
        a calibration line of sight.

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
        Magnifier.py example.config --plot_contributions

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
       print Magnifier.__doc__  # will print the big comment above.
       return

    for o,a in opts:
       if o in ("-h", "--help"):
          print Magnifier.__doc__
          return
       else:
          assert False, "unhandled option"

    # Check for setup file in array args:
    if len(args) == 1:
        configfile = args[0]
    elif len(args) == 2:
        configfile = args[0]
        plot_contributions == args[]
        print pangloss.doubledashedline
        print pangloss.hello
        print pangloss.doubledashedline
        print "Magnifier: generating magnification PDFs for lightcones of given density"
        print "Magnifier: taking instructions from",configfile
    else:
        print Magnifier.__doc__
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
    
    field_name = experiment.parameters['FieldName']
    field_overdensity = experiment.parameters['FieldOverdensity']

    field_name = numpy.genfromtxt(str(field_name), comments='#', usecols=0, dtype='S30')
    if len(field_overdensity)==1: field_overdensity = numpy.genfromtxt(str(field_overdensity), comments='#')
    
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

    print "Magnifier: found the lightcones..."
           
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

    for i in xrange(Nc):         
       calcones.append(pangloss.readPickle(calpickles[i]))
       if i==0: print calpickles[i]

    if DoCal=="False": #must be string type
       calcones=[]
       calpickles=[]

    allcones = calcones
    allconefiles = calpickles 

    # ==============================================================    
    # Find the overdensity of lightcones cut at m<22 in F125W
    # ==============================================================
    print "Magnifier: finding the distribution of lightcones with density..."

    lc_dens = []
    lc_galaxies = []
    total_galaxies = 0

    # Sort into lightcones for each field
    for i in xrange(Nc): 
       lc = pangloss.readPickle(calpickles[i])  
       num_galaxies = lc.numberWithin(radius=Rc,cut=[16,22],band=mag,units="arcmin")
       lc_galaxies.append(num_galaxies)   
       # Add to the total number of galaxies
       total_galaxies += num_galaxies
       del lc
    
    lc_dens = [Nc * float(x)/float(total_galaxies) for x in lc_galaxies]

    numpy.savetxt(CALIB_DIR+"/lc_density.txt", lc_dens) 

    print 'Mean overdensity in all fields = %.3f (this should =1)' % numpy.mean(lc_density)
    print 'Lightcone overdensities saved to file'
    print pangloss.dashedline
    del lc_dens
    del lc_density

    # ==============================================================
    # Sample all lightcones to make the pdfs
    # ==============================================================

    # --------------------------------------------------------------
    # Set up overdensity range

    density = field_overdensity
    drange = 0.02   # ~0.02 is what Zach used

    # --------------------------------------------------------------------
    # Find contribution to total kappa and mass at redshift intervals

    zmax = zs+0.1
    zbin = 25
    zbins = numpy.linspace(0.0,zmax,zbin)
    
    kappa_cont = numpy.zeros(Nc, zbin)
    Mh_cont = numpy.zeros(Nc, zbin)
    Mstell_cont = numpy.zeros(Nc, zbin)
          
    # --------------------------------------------------------------------
    # Select ALL lightcones and find their convergences

    pk = []
    pmu =[] 

    for j in xrange(Nc):        

        # Get lightcone
        lc = pangloss.readPickle(calpickles[j])

        # --------------------------------------------------------------------
        # Calculate mu and kappa for all lightcones
           
        # Redshift scaffolding:
        lc.defineSystem(zd,zs)
        lc.loadGrid(grid)

        # Figure out data quality etc:
        lc.configureForSurvey(experiment)

        if j % 1000 == 0 and j !=0:
           print ("Magnifier: ...on lightcone %i out of %i..." % (j,Nc))
               
        lc.snapToGrid(grid)
                   
        # Draw c from Mhalo:
        lc.drawConcentrations(errors=True)
                   
        # Compute each halo's contribution to the convergence:
        lc.makeKappas(truncationscale=5)
                   
        k_add=lc.combineKappas()
        mu_add=lc.combineMus(weakapprox=False)                    
                                                                            
        # Add magnification and convergence to global PDF
        pmu.append(lc.mu_add_total)
        pk.append(lc.kappa_add_total)

        if plot_contributions is True:
            kappa_cont[j:,] = lc.findContributions('kappa')   
            Mh_cont[j:,] = lc.findContributions('mass') 
            Mstell_cont[j:,] = lc.findContributions('stellarmass') 
           
        # Make a nice visualisation of one of the lightcones
        if j ==0:
            lc.plots('kappa', output=CALIB_DIR+"/example_snapshot_kappa_uncalib_z=1.4.png")
            lc.plots('mu', output=CALIB_DIR+"/example_snapshot_mu_uncalibz=1.4.png")
       
       del lc
       
    pk = numpy.array(pk)
    pmu = numpy.array(pmu)    

    # --------------------------------------------------------------------
    # Write PDFs to pickles
                   
    pangloss.writePickle(pk,CALIB_DIR+"/Pofk_z="+str(zs)+".pickle")
    pangloss.writePickle(pmu,CALIB_DIR+"/PofMu_z="+str(zs)+".pickle")

    del pk
    del pmu

    pmu = pangloss.readPickle(CALIB_DIR+"/PofMu_z="+str(zs)+".pickle")
    pk = pangloss.readPickle(CALIB_DIR+"/Pofk_z="+str(zs)+".pickle")

    # --------------------------------------------------------------------
    # Plot contributions to total kappa and mass at redshifts

    if plot_contributions:
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
        plt.savefig("figs/"+EXP_NAME+"contribution_z.pdf",dpi=300)
 
    # ==============================================================
    # Plot the pdfs
    # ==============================================================
    
    # --------------------------------------------------------------------
    # Calculate the smooth components
    kappa_smooth = numpy.mean(pk)
    mu_smooth = numpy.mean(pmu)
    print 'Uncalibrated: <kappa> =',kappa_smooth, '<mu> =',mu_smooth

    pmu = pmu - mu_smooth + 1.
    pk = pk - kappa_smooth + 0.
    print 'mean mu now = ', numpy.mean(pmu)
    # --------------------------------------------------------------------
    # Plotting convergence and magnification
    c = ['purple','g','r','b']
    
#    params = [{'param':'Kappa', 'name':r'$\kappa$', 'lc':pk, 'smooth':kappa_smooth, 'mean':0.0, 'height':60, 'min':-0.2, 'max':0.4},
#            {'param':'Mu', 'name':r'$\mu$', 'lc':pmu, 'smooth':mu_smooth, 'mean':1.0, 'height':30, 'min':0.4, 'max':2.0}]
   
    params = [{'param':'Mu', 'name':r'$\mu$', 'lc':pmu, 'smooth':mu_smooth, 'mean':1.0, 'height':30, 'min':0.4, 'max':2.0}]
    
    lc_density = numpy.genfromtxt(CALIB_DIR+"/lc_density.txt", comments="#")
    
    # =====================================================================
    # For only <=4 values of density
    # =====================================================================    
    if len(density) < 5:
        for k in range(len(params)):
        
            var = params[k]
            full_pdf = var['lc']
            name = var['name']     
            print 'Magnifier: Old min and max:', full_pdf.min(), full_pdf.max()
            
            # --------------------------------------------------------------------
            # Remove outliers --  only for v high overdensity?
            if var['param']=='Kappa':
                mask = numpy.where(full_pdf==full_pdf)
            if var['param']=='Mu':
                mask = numpy.where((full_pdf > var['smooth']-5.) & (full_pdf < 6.+var['smooth']))
                
            par = full_pdf[mask]
            new_density = lc_density[mask]
            
            print '         Removing outliers...'
            print '         New min and max:', par.min(), par.max()
            
            # Recalculate means and the calibrated pdf
            smooth_new = numpy.mean(par) 
            new_pdf = par - smooth_new + var['mean']
            par_mean = numpy.mean(new_pdf)
            print 'New mean (this should be',var['mean'],'):', par_mean         

            # --------------------------------------------------------------------
            # Plot PDF for all lines of sight
            
            print pangloss.dashedline
            print "Magnifier: constructing PDF for", var['param'],"..."  
            
            # --------------------------------------------------------------------
            # Select only lightcones within certain number density limits
            # Need to mask out the outliers
            
            for i in range(len(density)):
                sub_pdf = [] 
                
                for j in xrange(len(new_pdf)):
                    if density[i] - drange <= round(new_density[j],2) <= density[i] + drange:
                        sub_pdf.append(new_pdf[j])
                
                sub_pdf = numpy.array(sub_pdf)
                sub_mean = numpy.mean(sub_pdf) 
                Nlos = len(sub_pdf)
                            
                print "Magnifier: sampling %i LoS with number density ~ %.2f the average" % (Nlos, density[i])   
            
                outputfile = CALIB_DIR+"/"+EXP_NAME+"_Pof"+var['param']+"_"+"_z="+str(zs)+"_overdensity="+str(density[i])+".txt"                 
                numpy.savetxt(outputfile, sub_pdf) 
                print "Magnifier: saved PDFs to",outputfile

    # =====================================================================
    # For many values of overdensity e.g. BORG
    # =====================================================================    

    else:
        for k in range(len(params)):
      
            var = params[k]
            full_pdf = var['lc']  #- var['smooth'] + var['mean']
            name = var['name']     
            print 'Magnifier: Old min and max:', full_pdf.min(), full_pdf.max()
            
            # --------------------------------------------------------------------
            # Remove outliers  !!! I think this is only important for v. high overdensity
            if var['param']=='Kappa':
                mask = numpy.where(full_pdf==full_pdf)
            if var['param']=='Mu':
                mask = numpy.where((full_pdf >= 0.) & (full_pdf < 2.))
#                mask = numpy.where((full_pdf > var['smooth']-5.) & (full_pdf < 6.+var['smooth']))
            
            par = full_pdf[mask]
            new_density = lc_density[mask]

            print '         Removing Outliers...'
            print '         New min and max:', par.min(), par.max()
            
            # Recalculate means and the calibrated pdf
            smooth_new = numpy.mean(par) 
            new_pdf = par - smooth_new + var['mean']
            par_mean = numpy.mean(new_pdf)
            print 'New mean (this should be',var['mean'],'):', par_mean        
            
            # --------------------------------------------------------------------
            # Plot all lines of sight
            
            print pangloss.dashedline
            print "Magnifier: constructing PDF for", var['param'],"..."  
                    
            outputfile = CALIB_DIR+"/"+EXP_NAME+"_Pof"+var['param']+"_"+"_z="+str(zs)+"_allLoS.txt"                 
            numpy.savetxt(outputfile, sub_pdf) 
            print "Magnifier: saved all LoS PDFs to",outputfile
                                                      
            # --------------------------------------------------------------------
            # Select only lightcones within certain number density limits
            # Need to mask out the outliers
            
            means, fieldname = [], []
            
            for i in range(len(density)):
                sub_pdf = [] 
                
                for j in xrange(len(new_pdf)):
                    if density[i] - drange <= round(new_density[j],2)  <= density[i] + drange:
                        sub_pdf.append(new_pdf[j])
                
                sub_pdf = numpy.array(sub_pdf)
                Nlos = len(sub_pdf)
                if Nlos == 0:
                    print "Magnifier: %s - there are NO LoS with number density ~ %.2f the average" % (field_name[i], density[i])    
                    print "Magnifier: %s - no PDF will be made for this field" % (field_name[i])    
                
                else:
 
                    sub_mean = numpy.mean(sub_pdf) 
                                
                    print "Magnifier: %s - sampling %i LoS with number density ~ %.2f the average, mean mu=%.2f" % (field_name[i], Nlos, density[i], sub_mean)
                    
                    if var['param']=='Mu':
                        numpy.savetxt(CALIB_DIR+"/"+EXP_NAME+str(field_name[i])+"_PofMu.txt", sub_pdf)             
                    
                    means.append(sub_mean)
                    fieldname.append(field_name[i])

            meanmu_table = numpy.array([fieldname, means]).T
            ascii.write(meanmu_table, "figs/borg_remove1815/borg_pdf_table_means.txt", names=['#field','mean_mu'])

            print "Mean mu of all the fields = ",numpy.mean(means)        
            print "Magnifier: saved PDFs to",outputfile 
            
    print pangloss.doubledashedline

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -   
    
    print pangloss.doubledashedline
    return

# ======================================================================

if __name__ == '__main__': 
    Magnifier(sys.argv[1:])

# ======================================================================
