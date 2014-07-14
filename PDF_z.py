#!/usr/bin/env python
# ======================================================================

import pangloss

import sys,getopt,cPickle,numpy
import matplotlib.pyplot as plt

import glob

from scipy.stats.kde import gaussian_kde
from math import pi

# ======================================================================

def PDF_z(argv):
    """
    NAME
        PDF_z.py

    PURPOSE
        Read simulation lightcones and produce PDFs for convergence,
        shear, and magnification across all lines of sight used for various
        redshifts
        
    COMMENTS
        The config file contains the list of lightcones to be
        reconstructed, in the form of either a directory or a single
        instance.

    FLAGS
        -h            Print this message [0]

    INPUTS
        configfile    Plain text file containing Pangloss configuration
        parameter     Must be either Kappa, Mu or Gamma

    OUTPUTS

    EXAMPLE
        PDF_z.py example.config Mu

    BUGS

    AUTHORS
      This file is part of the Pangloss project, distributed under the
      GPL v2, by Tom Collett (IoA) and  Phil Marshall (Oxford). 
      Please cite: Collett et al 2013, http://arxiv.org/abs/1303.6564

    HISTORY
      2014-05-29    Started Charlotte Mason (UCSB)
    """

    # --------------------------------------------------------------------

    try:
       opts, args = getopt.getopt(argv,"h",["help"])
    except getopt.GetoptError, err:
       print str(err) # will print something like "option -a not recognized"
       print PDF_z.__doc__  # will print the big comment above.
       return

    for o,a in opts:
       if o in ("-h", "--help"):
          print PDF_z.__doc__
          return
       else:
          assert False, "unhandled option"

    # Check for setup file in array args:
    if len(args) == 2:
        configfile = args[0]
        parameter = args[1]
        print pangloss.doubledashedline
        print pangloss.hello
        print pangloss.doubledashedline
        print "PDF_z: calculating PDF for",parameter,"at various redshifts for all LOS"
        print "PDF_z: taking instructions from",configfile
        print pangloss.dashedline
    else:
        print PDF_z.__doc__
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
  #  zs = numpy.arange(1.0, 8.0, 0.5)
    #zs = [1.1, 2.1, 5.7]
    zs = [8.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0]

    Ncats = 21
    
    calpickles = []
    Nc = experiment.parameters['NCalibrationLightcones'] * Ncats       ### should be 24

    paths = '%s/*_lightcone.pickle' % (CALIB_DIR)
    found = glob.glob(paths)
    if len(found) > 0: calpickles = found
    
#    Nc = experiment.parameters['NCalibrationLightcones']
#    for i in range(Nc):
 #       calpickles.append(experiment.getLightconePickleName('sim_notborg',pointing=i))
       
    # Ray tracing:
    RTscheme = experiment.parameters['RayTracingScheme']
    
    # Sampling Pr(kappah|D):
    Ns = experiment.parameters['NRealisations']
    
    # Reconstruct calibration lines of sight?
    DoCal = experiment.parameters['ReconstructCalibrations']

    
    # --------------------------------------------------------------------
    # Read in lightcones from pickles:

#    calcones = []
#    for i in range(Nc):
#        calcones.append(pangloss.readPickle(calpickles[i]))

    if DoCal=="False": #must be string type
        calcones=[]
        calpickles=[]

#    allcones = calcones
#    allconefiles = calpickles 
    
#    ndensity_field = 4.2110897902 # galaxies brighter than i 22 in num/arcmin^2   hilbert    
    ndensity_field = 7.87131907218 # galaxies brighter than F125W 22 in num/arcmin^2  Henriques  
    
    density = 1.0 # overdensity we are looking at
    drange = 0.03
    # ==============================================================
    # Produce the PDFs
    # ==============================================================

    # --------------------------------------------------------------------
    # Select ALL lightcones and start PDFs for the lensing parameters
    
    c = ['r','b','g','orange','yellow']       
    
    mean_mu=[]
    sd_mu = []
    
    radius_cone = 1.0
    area = pi * radius_cone**2.        
 
    # plt.figure(1)
    
    for i in range(len(zs)):
        print 'Source redshift =',zs[i]
        
        # Make redshift grid
        grid = pangloss.Grid(zd, zs[i], nplanes=100)
        
#        pk = numpy.zeros(len(Nc))
        pmu =numpy.zeros(Nc)
#        pg1=numpy.zeros(len(Nc))
#        pg2=numpy.zeros(len(Nc))
#        pg=numpy.zeros(len(Nc))
        
       # lc_density = []   
                     
        sub_pdf = [] 
                
        for j in range(Nc):        

            # Get lightcone, and fill PDFs
            lc = pangloss.readPickle(calpickles[j])


    #        n_gals = lc.numberWithin(radius=radius_cone,cut=[16,22],band="F125W",units="arcmin")
    #        lc_dens = n_gals/(area * ndensity_field)
                        
            # Redshift scaffolding:
            lc.defineSystem(zd, zs[i])
            lc.loadGrid(grid)
        
            # Figure out data quality etc:
            lc.configureForSurvey(experiment)    
                    
            lc.snapToGrid(grid)
                        
            # Draw c from Mhalo:
            lc.drawConcentrations(errors=True)
                        
            # Compute each halo's contribution to the convergence:
            lc.makeKappas(truncationscale=5)
                        
            k_add=lc.combineKappas()
            mu_add=lc.combineMus(weakapprox=False)
                                
            # Add lensing parameters to global PDFs
            pmu[j] = lc.mu_add_total

#            pk[j]=lc.kappa_add_total
                
#            pg[j]=lc.G1sum
#            pg2[j]=lc.G2sum
#            pg[j]=lc.Gsum
            
            # Select cones by overdensity for the mu vs z plot
#            if density - drange <= lc_dens <= density + drange:
#                sub_pdf.append(mu_add)  
                                                                                    
#            x = allconefiles[j]
        
      #  print lc_density[0]    
#        pk = numpy.array(pk)
        
        pmu = numpy.array(pmu)    
    
    # --------------------------------------------------------------------
    # Write PDFs to pickles
                    
#        pangloss.writePickle(pk,CALIB_DIR+"/Pofk_z="+str(zs)+".pickle")
        pangloss.writePickle(pmu,CALIB_DIR+"/many_z_PofMu_z="+str(zs)+".pickle")
        
#        del pk
        del pmu
    
        #pmu = pangloss.readPickle(CALIB_DIR+"/many_z_PofMu_z="+str(zs)+".pickle")
#        pk = pangloss.readPickle(CALIB_DIR+"/Pofk_z="+str(zs)+".pickle")
    
        
                
#        pg1 = numpy.array(pg1)    
#        pg2 = numpy.array(pg2)    
#        pg = numpy.array(pg)    
#         
#        sub_pdf = numpy.array(sub_pdf) 
#        print 'For fields of overdensity %.2f, we have %i lightcones'%(density,len(sub_pdf))
                            
#        # ==============================================================
#        # Plot the pdfs
#        # ==============================================================
#        
#        # Smooth component corrections
##        kappa_smooth = numpy.mean(pk)
#        mu_smooth = numpy.mean(pmu)
#
#        sub_mu = sub_pdf - mu_smooth + 1.0
##        print numpy.mean(sub_mu)
#        mean_mu.append(numpy.mean(sub_mu))
#        sd_mu.append(numpy.std(sub_mu))
#        
#        if parameter == 'Kappa':
#            pass
##            params = {'param':'Kappa', 'name':r'$\kappa$', 'pdf':pk,
##                        'smooth':kappa_smooth, 'mean':0.0, 'height':35}            
#        
#        elif parameter == 'Gamma':
#            params = {'param':'Gamma', 'name':r'$\gamma$', 'pdf':pg} 
#        
#        else:
#            params = {'param':'Mu', 'name':r'$\mu$', 'pdf':pmu,
#                        'smooth':mu_smooth, 'mean':1.0, 'height':18}
#           
#        param = params['param']
#        name = params['name']
#        pdf = params['pdf']
#        
#        par = pdf                  
#                            
#        if param == 'Kappa':
#            smooth = params['smooth']
#            mean = params['mean']
#            par = par - smooth + mean
#        
#        if param == 'Mu':
#            smooth = params['smooth']
#            mean = params['mean']
#            par1 = par - smooth + mean
#           # print par1.min(), par1.max()
#            mask = numpy.where((par > -1.0) & (par < 2.0)) 
#            par = par[mask]
#            smooth_new = numpy.mean(par) 
#            par = par - smooth_new + mean
#            par_mean = numpy.mean(par) 
##            plt.xlim(0.75,1.25)    
#                
#                        
##        outputfile = "figs/final_"+EXP_NAME+"_compare_z_Pof"+param+"_many.pdf" 
#             
#        par_mean = numpy.mean(par) 
#        
#        Nlos = len(par)            
            
        """    
            par1 = pg1[:,0]
            par2 = pg2[:,0]
            par = pg[:,0]
    
                            
            par1_kde = gaussian_kde(par1)
            par2_kde = gaussian_kde(par2)
            par_kde = gaussian_kde(par)
    
            x = numpy.linspace(par1.min(),par1.max(),3*Nlos)
            y = numpy.linspace(par.min(),par.max(),3*Nlos)
                    
            n1, bins1, patches1 = plt.hist(par1, 8, facecolor='None', histtype='step',linestyle=('dashed'), normed=True, label=r'$\gamma_1, z_s = $'+str(zs[i]))
            plt.setp(patches1, 'edgecolor', c[i])  
            n2, bins2, patches2 = plt.hist(par2, 8, facecolor='None', histtype='step', linestyle=('dotted'), normed=True, label=r'$\gamma_2$')
            plt.setp(patches2, 'edgecolor', c[i])        
            n, bins, patches = plt.hist(par, 8, facecolor=c[i], alpha=0.4, normed=True, label=r'|$\gamma$|')
            plt.setp(patches, 'edgecolor', 'none') 
        
        """                                
#        IQR = numpy.percentile(par, 75) - numpy.percentile(par, 25)
#        Nbins = 6./(2 * IQR * (float(Nlos))**(-1./3.))
#        Nbins = int(round(Nbins,-1))/10
#        print 'Nbins =', Nbins                  
#        
#        n, bins, patches = plt.hist(par, Nbins, facecolor=c[i], normed=True,alpha=0.4, label=r'$z_s = $'+str(zs[i]))
#        plt.setp(patches, 'edgecolor', 'None')
#            
#        par_kde = gaussian_kde(par, bw_method=0.1)#/new_pdf.std(ddof=1))
#        par_kde.covariance_factor = lambda : .05
#        par_kde._compute_covariance()        
#                
#        x = numpy.linspace(par.min()-0.05,par.max()+0.05,1000)
#        plt.plot(x, par_kde(x), color=c[i])
#                                        
#        plt.axvline(x=par_mean, color='k', linestyle='dashed')
#            
#        plt.xlabel(name)
#        plt.ylabel(r'P(%s)' % name)   
#                
#    plt.ticklabel_format(useOffset=False, axis='x')
#    plt.legend(loc=1)                                                        
#
#    outputfile = "figs/final_"+EXP_NAME+"_compare_z_Pof"+param+"_many_test.pdf"
#                           
#    plt.savefig(outputfile,dpi=300)
#
##    plt.figure(2)
#
#    plt.plot(zs, mean_mu, label=r'$\langle \mu \rangle$')
##    plt.plot(zs, sd_mu, label=r'$\sigma_\mu$')
#    plt.xlabel(r'$z_s$')
#    plt.legend(loc=1)
#                            
#    outputfile = "figs/Mu_zs.png" 
#
#    plt.savefig(outputfile,dpi=300)
    
    print pangloss.doubledashedline
    print "PDF_z: saved P("+param+") to",outputfile

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -   
    
    print pangloss.doubledashedline
    return

# ======================================================================

if __name__ == '__main__': 
    PDF_z(sys.argv[1:])

# ======================================================================
