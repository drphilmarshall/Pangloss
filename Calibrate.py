#!/usr/bin/env python
# ======================================================================

import pangloss

import sys,getopt,cPickle,numpy

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
      Please cite: Collett et al 2013, arxiv/###

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
        print "Pangloss Calibrate: Calibrating lightcones according to instructions in",configfile
    else:
        print Calibrate.__doc__
        return

    # --------------------------------------------------------------------
    # Read in configuration, and extract the ones we need:

    experiment = pangloss.Configuration(configfile)

    Nc = experiment.parameters['NCalibrationLightcones']

    comparator=experiment.parameters['Comparator']
    comparatorType=experiment.parameters['ComparatorType']
    comparatorWidth=experiment.parameters['ComparatorWidth']

    CALIB_DIR = experiment.parameters['CalibrationFolder'][0]
    jointdistfile= CALIB_DIR+'/Calibrationguide/'+comparator+comparatorType+'.pickle'
    
    RES_DIR= experiment.parameters['ResultFolder'][0]
    EXP_NAME=experiment.parameters['ExperimentName']
    resultfile=RES_DIR+'/'+EXP_NAME+'_PofKappaExt.pickle'
    print resultfile
    # --------------------------------------------------------------------
    #Mode1 : generate a jointdistribution from the calibration dataset:
    if Mode==1 or Mode==3:
    #First find the calibration pdfs for kappa_h
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
            print "I don't know that comparator. If you want to use a comparator other than kappa_h, you'll need to code it up (this should be easy but you can ask tcollett@ast.cam.uk for help). Exiting"
            exit()



    #caluclate comparators:
        callist=numpy.empty((Nc,2))
        for i in range(Nc):
            C=calresultpickles[i]
            pdf=pangloss.readPickle(C)
            if comparator=="kappa_h":
                if comparatorType=="median":# note we created a special file for this choice of comparator and comparator type. You could also use the comparatortype=="mean" code swapping mean for median.
                    callist[i,0]=pdf[0]
                    callist[i,1]=numpy.median(pdf[1])
                elif comparatorType=="mean":
                    callist[i,0]=pdf.truth[0]
                    callist[i,1]=numpy.mean(pdf.samples)
                else: 
                    print "I don't know that comparatorType. exiting"
                    exit()

        pangloss.writePickle(callist,jointdistfile)

    # --------------------------------------------------------------------
    
    #Mode2 : calibrate a real line of sight using the joint-distribution
    if Mode==2 or Mode==3:
        callibguide=pangloss.readPickle(jointdistfile)


        obspickle = experiment.getLightconePickleName('real')
        pfile = obspickle.split('.')[0].split("_lightcone")[0]+"_PofKappah.pickle"

        pdf=pangloss.readPickle(pfile)
        if comparator=="kappa_h":
            if comparatorType=="median":# note we created a special file for this choice of comparator and comparator type. You could also use the comparatortype=="mean" code swapping mean for median.
                RealComparator=numpy.median(pdf.samples)
            elif comparatorType=="mean":
                RealComparator=numpy.mean(pdf.samples)
            else: 
                print "I don't know that comparatorType. exiting"
                exit()

        pdf=pangloss.PDF(["kappa_ext","weight"])
        
        dif=(callibguide[:,1]-RealComparator)
        weights=numpy.exp(-(dif**2)/(2*(comparatorWidth)**2))
        weights/=numpy.sum(weights)
        samples=callibguide[:,0]
        samplesandweights=callibguide.copy()
        samplesandweights[:,1]=weights

        pdf.samples=(samplesandweights)


        average = numpy.average(samples, weights=weights)
        variance = numpy.dot(weights, (samples-average)**2)/weights.sum()
        average,std=average, variance**.5

        
        pangloss.writePickle(pdf,resultfile)
        print "Calibrate: Your lightcone has been calibrated."
        print "The reconstruction suggests a kappa_ext of",\
            "%.3f +\- %.3f for this lightcone"%(average,std)
        print "Calibrate: a pdf of kappa_ext has been output to:"
        print resultfile
        print "in the form of sample kappa_ext values and weights." 
        print "To read and process this file, try:"
        print
        print "import pangloss"
        print "pdf=pangloss.readPickle(\"%s\")"%resultfile
        print "kappa_samples=pdf.call(\"kappa_ext\")"
        print "kappa_weights=pdf.call(\"weight\")"
        print
        print "And then do whatever you wanted this for! Goodluck, TEC & PJM"
        
        

    return

# ======================================================================

if __name__ == '__main__':
    Calibrate(sys.argv[1:])    
    test=False
    if test:
        import pangloss
        pdf=pangloss.readPickle("/home/tcollett/Pangloss/results/Example_PofKappaExt.pickle")
        kappa_samples=pdf.call("kappa_ext")
        kappa_weights=pdf.call("weight")
        print kappa_samples,kappa_weights
        print numpy.sum(kappa_weights)

# ======================================================================
