#!/usr/bin/env python
# ======================================================================

import pangloss

import getopt,cPickle

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
         print "Pangloss Calibrate: Calibrating lightcones according to instructions in",configfile
    else:
       print Calibrate.__doc__
       return

    # --------------------------------------------------------------------

    # --------------------------------------------------------------------

    return

# ======================================================================

if __name__ == '__main__':
    Calibrate(sys.argv[1:])

# ======================================================================



