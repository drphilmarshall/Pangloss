#!/usr/bin/env python
# ======================================================================

import pangloss

import sys,getopt,cPickle,numpy

# ======================================================================

def Calibrate(argv,Mode=3):
    """
    NAME
        Calibrate.py

    PURPOSE
        Transform the results of the lightcone reconstruction process,
        Pr(kappah|D), into our target PDF, Pr(kappa|D).

    COMMENTS
        All PDF input is provided as a list of samples. There are two
        modes of operation:

        1) The Pr(kappah|C) for an ensemble of calibration lightcones are
           compressed into a single number (currently the
           median), and then combined with the true kappa values to make
           Pr(kappa,kappah|C). This is written out as a 2D sample list.

        2) The Pr(kappah|D) for a single observed lightcone is compressed
           into a single number (currently the median). This is then used
           to take a slice from Pr(kappa,kappah|C) to make Pr(kappa|D,C).

        Both 1 and 2 can be carried out in series if desired.

    FLAGS
        -h            Print this message [0]

    INPUTS
        configfile    Plain text file containing Pangloss configuration

    OUTPUTS
        stdout        Useful information
        samples       From 1) Pr(kappa,kappah|C) or 2) Pr(kappa|D,C)


    EXAMPLE

        Calibrate.py example.config

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
       print Calibrate.__doc__  # will print the big comment above.
       return

    for o,a in opts:
        if o in ("-h", "--help"):
            print Calibrate.__doc__
            return
        else:
            assert False, "unhandled option"

    # Check for setup file in array args:
    if len(args) == 1:
        configfile = args[0]
        print pangloss.doubledashedline
        print pangloss.hello
        print pangloss.doubledashedline
        print "Calibrate: transforming Pr(kappah|D) to Pr(kappa|D)"
        print "Calibrate: taking instructions from",configfile
    else:
        print Calibrate.__doc__
        return

    # --------------------------------------------------------------------
    # Read in configuration, and extract the ones we need:

    experiment = pangloss.Configuration(configfile)

    Nc = experiment.parameters['NCalibrationLightcones']

    comparator = experiment.parameters['Comparator']
    comparatorType = experiment.parameters['ComparatorType']

    # Figure out which mode is required:
    ModeName = experiment.parameters['CalibrateMode']
    if ModeName=='Joint': Mode = 1
    if ModeName=='Slice': Mode = 2
    if ModeName=='JointAndSlice': Mode = 3

    # --------------------------------------------------------------------
    # Mode 1: generate a joint distribution, eg Pr(kappah,kappa)
    # from the calibration dataset:
    
    if Mode==1 or Mode==3:
    
        #First find the calibration pdfs for kappa_h:
        calpickles = []
        for i in range(Nc):
            calpickles.append(experiment.getLightconePickleName('simulated',pointing=i))

        calresultpickles=[]
        if comparator=="kappa_h" and comparatorType=="median":
            for i in range(Nc):
                x = calpickles[i]
                pfile2 = x.split('.')[0].split("_lightcone")[0]+"_KappaHilbert_Kappah_median.pickle"
                calresultpickles.append(pfile2)

        elif comparator=="kappa_h" and comparatorType!="median": 
            for i in range(Nc):
                x = calpickles[i]
                pfile = x.split('.')[0].split("_lightcone")[0]+"_PofKappah.pickle"

                calresultpickles.append(pfile)
        else:
            print "Calibrate: Unrecognised comparator "+Comparator
            print "Calibrate: If you want to use a comparator other than kappa_h, "
            print "Calibrate: you'll need to code it up!"
            print "Calibrate: (This should be easy, but you can ask tcollett@ast.cam.uk for help)."
            exit()

        # Now read in, or calculate, comparators:
        comparatorlist=numpy.empty(Nc)
        hilbertkappa=numpy.empty(Nc)
        for i in range(Nc):
            C=calresultpickles[i]
            pdf=pangloss.readPickle(C)
            if comparator=="kappa_h":
                if comparatorType=="median":# note we created a special file for this choice of comparator and comparator type. You could also use the comparatortype=="mean" code swapping mean for median.
                    comparatorlist[i]=numpy.median(pdf[1])
                    # BUG: we already took this median in Reconstruct!
                    hilbertkappa[i]=pdf[0]
                elif comparatorType=="mean":
                    hilbertkappa[i]=pdf.truth[0]
                    comparatorlist[i]=numpy.mean(pdf.samples)
                else: 
                    print "I don't know that comparator type. exiting"
                    exit()

        print hilbertkappa,comparatorlist
    
    # --------------------------------------------------------------------
    # Mode 2: calibrate a real line of sight's Pr(kappah|D) using the
    # joint distribution Pr(kappa,<kappah>|D)

    if Mode==2 or Mode==3:
    
        obscone = pangloss.readPickle(obspickle)






    return

# ======================================================================

if __name__ == '__main__':
    Calibrate(sys.argv[1:])

# ======================================================================

# # Old code from MakeCalibrationGuide....
# 
# import cPickle
# import pylab as plt
# import scipy.stats as stats
# from scipy import interpolate,optimize
# import numpy
# import numpy.random as rnd
# import glob
# #------------------------------------------------------------------
# def MakeReconstructionCalibration(CalibrationResultsDirectory,surveyname):
#    results="%s/cone*_%s.result"%(CalibrationResultsDirectory,i,surveyname)
#    resultlist=glob.glob(results)
# 
#    Zl=len(resultlist)
#    
#    Cal_list=numpy.empty((Zl,5))
# 
#    for i in range(len(resultlist)):
#       if i % 1000==0: print "making calibration", i,"of", len(resultlist)
#       F=open(resultlist[i],"rb")
#       res=cPickle.load(F)
#       F.close()
#       
#       kappa_Hilbert,gamma1_Hilbert,gamma2_Hilbert,kappa_halo,gamma1_halo,gamma2_halo=res
# 
#       mu_Hilbert=1/((1-kappa_Hilbert)**2 - gamma_Hilbert**2)
#       mu_halo=1/((1-kappa_halo)**2 - gamma_halo**2)
# 
#       Cal_list[i,0]=kappa_Hilbert
#       Cal_list[i,1]=mu_Hilbert
#       Cal_list[i,2]=numpy.median(kappa_halo)
#       Cal_list[i,3]=numpy.median(mu_halo)
#       Cal_list[i,4]=1
# 
#    F=open("%s/%s.jointdist"%(CalibrationResultsDirectory,surveyname),'wb')
#    cPickle.dump(Cal_list,F,2)
#    F.close()
#    
# def MakeWeightCalibration(CalibrationResultsDirectory,surveyname):
#    results="%s/cone*_%s.result"%(CalibrationResultsDirectory,i,surveyname)
#    resultlist=glob.glob(results)
# 
#    Zl=len(resultlist)
#    
#    Cal_list=numpy.empty((Zl,5))
# 
#    for i in range(len(resultlist)):
#       if i % 1000==0: print "making calibration", i,"of", len(resultlist)
#       F=open(resultlist[i],"rb")
#       res=cPickle.load(F)
#       F.close()
#       
#    #need to work out what goes in here...
#    """
#       kappa_Hilbert,gamma1_Hilbert,gamma2_Hilbert,=res
# 
#       mu_Hilbert=1/((1-kappa_Hilbert)**2 - gamma_Hilbert**2)
#       mu_halo=1/((1-kappa_halo)**2 - gamma_halo**2)
# 
#       Cal_list[i,0]=kappa_Hilbert
#       Cal_list[i,1]=mu_Hilbert
#       Cal_list[i,2]=numpy.median(kappa_halo)
# 
#    F=open("%s/%s.jointdist"%(CalibrationResultsDirectory,surveyname),'wb')
#    cPickle.dump(Cal_list,F,2)
#    F.close()
#    """
#    
# #code to use the above joint distribution...
# 
# def CalibrateResult(CalibrationResultsDirectory,surveyname,Value,width=0.01,Magnification=False):
#    F=open("%s/%s.jointdist"%(CalibrationResultsDirectory,surveyname),'rb')
#    Klist=cPickle.load(F)
#    F.close()
# 
#    if magnification==False:
#       Klist[:,4]=numpy.exp(-((Klist[:,2]-Value)**2)/(2*(width)**2))
#    else:
#       Klist[:,4]=numpy.exp(-((Klist[:,3]-Value)**2)/(2*(width)**2))
# 
#    kappabins=numpy.linspace(-0.1,5.0,5101)
#    Klistdig=numpy.digitize(Klist[:,1],kappabins)
# 
# 
#    Result=numpy.empty((len(Klist[:,0]),3))
#    Result[:,0]=Klist[:,0]
#    Result[:,1]=Klist[:,1]
#    Result[:,2]=Klist[:,4]
#    
#    return Result
# 
