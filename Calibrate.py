#!/usr/bin/env python
# ======================================================================

import pangloss

import sys,getopt,cPickle,numpy

import scipy.stats as stats

# ======================================================================

def Calibrate(argv):
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

        Both 1 and 2 can be carried out in series if desired (Mode=3).

    FLAGS
        -h            Print this message [0]

    INPUTS
        configfile    Plain text file containing Pangloss configuration

    OPTIONAL INPUTS
        --mode        Operating mode 1,2 or 3. See COMMENTS above.

    OUTPUTS
        stdout        Useful information
        samples       From 1) Pr(kappa,kappah|C) or 2) Pr(kappa|D,C)


    EXAMPLE

        Calibrate.py example.config

    BUGS

    AUTHORS
      This file is part of the Pangloss project, distributed under the
      GPL v2, by Tom Collett (IoA) and  Phil Marshall (Oxford).
      Please cite: Collett et al 2013, http://arxiv.org/abs/1303.6564

    HISTORY
      2013-03-21 started Collett & Marshall (Oxford)
    """

    # --------------------------------------------------------------------
    try:
       opts, args = getopt.getopt(argv,"hm:",["help","mode"])
    except getopt.GetoptError, err:
       print str(err) # will print something like "option -a not recognized"
       print Calibrate.__doc__  # will print the big comment above.
       return
    Mode=3
    for o,a in opts:
        if o in ("-h", "--help"):
            print Calibrate.__doc__
            return
        elif o in ("-m", "--mode"):
            Mode = int(a)
            assert Mode < 4 and Mode >0, "unhandled Mode"
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

    EXP_NAME = experiment.parameters['ExperimentName']

    Nc = experiment.parameters['NCalibrationLightcones']

    comparator=experiment.parameters['Comparator'] 
    comparatorType=experiment.parameters['ComparatorType']
    comparatorWidth=experiment.parameters['ComparatorWidth']

    # Figure out which mode is required:
    ModeName = experiment.parameters['CalibrateMode']
    if ModeName=='Joint': Mode = 1
    if ModeName=='Slice': Mode = 2
    if ModeName=='JointAndSlice': Mode = 3

    CALIB_DIR = experiment.parameters['CalibrationFolder'][0]
    jointdistfile= CALIB_DIR+'/'+comparator+'_'+comparatorType+'.pickle'
    jointdistasPDFfile= CALIB_DIR+'/'+comparator+'_'+comparatorType+'_asPDF.pickle'
    

    # Final result is PDF for kappa:
    x = experiment.parameters['ObservedCatalog'][0]
    resultfile = x.split('.')[0]+"_"+EXP_NAME+"_PofKappa.pickle"
    
    # --------------------------------------------------------------------
    # Mode 1: generate a joint distribution, eg Pr(kappah,kappa)
    # from the calibration dataset:
    
    if Mode==1 or Mode==3:
    
        print pangloss.dashedline
    
        # First find the calibration pdfs for kappa_h:
        calpickles = []
        for i in range(Nc):
            calpickles.append(experiment.getLightconePickleName('simulated',pointing=i))

        calresultpickles=[]
        if comparator=="Kappah" and comparatorType=="median":
            for i in range(Nc):
                x = calpickles[i]
                pfile = x.split('.')[0].split("_lightcone")[0]+"_"+EXP_NAME+"_KappaHilbert_Kappah_median.pickle"
                calresultpickles.append(pfile)

        elif comparator=="Kappah" and comparatorType!="median": 
            for i in range(Nc):
                x = calpickles[i]
                pfile = x.split('.')[0].split("_lightcone")[0]+"_"+EXP_NAME+"_KappaHilbert_Kappah_"+comparatorType+".pickle"
                calresultpickles.append(pfile)
        else:
            print "Calibrate: Unrecognised comparator "+Comparator
            print "Calibrate: If you want to use a comparator other than kappa_h, "
            print "Calibrate: you'll need to code it up!"
            print "Calibrate: (This should be easy, but you can ask tcollett@ast.cam.uk for help)."
            exit()

        # Now calculate comparators:
        callist=numpy.empty((Nc,2))
        jd=pangloss.PDF(["kappa_ext",comparator+'_'+comparatorType])

        for i in range(Nc):
            C = calresultpickles[i]
            pdf = pangloss.readPickle(C)
 
            if comparator=="Kappah":

                if comparatorType=="median": 
                    # Recall that we created a special file for this 
                    # choice of comparator and comparator type, in 
                    # Reconstruct. You could also use the 
                    # comparatortype=="mean" code, swapping mean for median.
                    callist[i,0]=pdf[0]
                    callist[i,1]=pdf[1][0]
                
                elif comparatorType=="mean":
                    callist[i,0] = pdf.truth[0]
                    callist[i,1] = numpy.mean(pdf.samples)

                else: 
                    print "Calibrate: Unrecognised comparatorType "+comparatorType
                    print "Calibrate: If you want to use a comparatorType other than median "
                    print "Calibrate: or mean, you'll need to code it up!"
                    print "Calibrate: (This should be easy, but you can ask tcollett@ast.cam.uk for help)."
                    exit()
                jd.append(callist[i])

        pangloss.writePickle(callist,jointdistfile)
        
        # Also store the joint dist as a pangloss pdf:
        pangloss.writePickle(jd,jointdistasPDFfile)
        
        # Plot:
        plotfile = jointdistasPDFfile.split('.')[0]+'.png'
        jd.plot("Kappah_median","kappa_ext",weight=None,output=plotfile,title="The joint distribution of $\kappa_{\mathrm{ext}}$ and calibrator \n\n (more correlated means a better calibrator!)")

        print "Calibrate: calibration joint PDF saved in:"
        print "Calibrate:     "+jointdistfile
        print "Calibrate: and "+jointdistasPDFfile
        print "Calibrate: you can view this PDF in "+plotfile

    # --------------------------------------------------------------------
    # Mode 2: calibrate a real line of sight's Pr(kappah|D) using the
    # joint distribution Pr(kappa,<kappah>|D)

    if Mode==2 or Mode==3:

        print pangloss.dashedline
    
        callibguide = pangloss.readPickle(jointdistfile)

        obspickle = experiment.getLightconePickleName('real')
        pfile = obspickle.split('.')[0].split("_lightcone")[0]+'_'+EXP_NAME+"_PofKappah.pickle"

        pdf=pangloss.readPickle(pfile)

        if comparator=="Kappah":
            if comparatorType=="median":# note we created a special file for this choice of comparator and comparator type. You could also use the comparatortype=="mean" code swapping mean for median.
                RealComparator=numpy.median(pdf.samples)
            elif comparatorType=="mean":
                RealComparator=numpy.mean(pdf.samples)
            else: 
                print "I don't know that comparatorType. exiting"
                exit()

        pdf = pangloss.PDF(["kappa_ext","weight"])

        #print RealComparator
        #print numpy.median(callibguide[:,1]),numpy.std(callibguide[:,1])

        dif=(callibguide[:,1]-RealComparator)
        weights=dif*0.0
        weights[numpy.abs(dif)<comparatorWidth]=1.
        weights/=numpy.sum(weights)
        samples=callibguide[:,0]
        samplesandweights=callibguide.copy()
        samplesandweights[:,1]=weights

        pdf.samples=(samplesandweights)

        plotfile = resultfile.split('.')[0]+".png"
        pdf.plot('kappa_ext',weight='weight',output=plotfile)

        average = numpy.average(samples, weights=weights)
        variance = numpy.dot(weights, (samples-average)**2)/weights.sum()
        average,std=average, variance**.5

        #if step function weights can calculate 68%CL easily:
        included=samples[weights>0]
        onesigconfidence=numpy.abs(\
            stats.scoreatpercentile(included,84)-
            stats.scoreatpercentile(included,16)\
                )/2.
            
        pangloss.writePickle(pdf,resultfile)

        print "Calibrate: your reconstructed lightcone has been calibrated,"
        print "Calibrate: suggesting it has a kappa_ext of",\
            "%.3f +\- %.3f"%(average,onesigconfidence)
        print "Calibrate: the PDF for kappa_ext has been output to "+resultfile
        print "Calibrate: in the form of sample kappa_ext values, and their weights." 
        print "Calibrate: you can view this PDF in "+plotfile
        print
        print "Calibrate: To read and process this file, try:"
        print
        print "   import pangloss"
        print "   pdf = pangloss.readPickle(\"%s\")"%resultfile
        print "   kappa_samples = pdf.getParameter(\"kappa_ext\")"
        print "   kappa_weights = pdf.getParameter(\"weight\")"

    # --------------------------------------------------------------------

    print
    print pangloss.doubledashedline
        
    return resultfile,jointdistasPDFfile

# ======================================================================

if __name__ == '__main__':
    Calibrate(sys.argv[1:])    
       
# ======================================================================
