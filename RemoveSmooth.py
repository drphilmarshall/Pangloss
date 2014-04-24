#!/usr/bin/env python
# ======================================================================

import pangloss

import sys,getopt,cPickle,numpy
import matplotlib.pyplot as plt


# ======================================================================

def RemoveSmooth(argv):
    """
    NAME
        RemoveSmooth.py

    PURPOSE


    COMMENTS

    INPUTS
        configfile    Plain text file containing Pangloss configuration

    OUTPUTS


    EXAMPLE

    BUGS

    AUTHORS

    HISTORY
      2014-04-23 started Mason (UCSB)
    """
    # --------------------------------------------------------------------

    print 'hello'
    try:
       opts, args = getopt.getopt(argv,"h",["help"])
    except getopt.GetoptError, err:
       print str(err) # will print something like "option -a not recognized"
       print RemoveSmooth.__doc__  # will print the big comment above.
       return

    for o,a in opts:
       if o in ("-h", "--help"):
          print RemoveSmooth.__doc__
          return
       else:
          assert False, "unhandled option"

    # Check for setup file in array args:
    if len(args) == 1:
        configfile = args[0]
        print pangloss.doubledashedline
        print pangloss.hello
        print pangloss.doubledashedline
        print "RemoveSmooth: removing the smooth component of convergence"
        print "RemoveSmooth: taking instructions from",configfile
    else:
        print RemoveSmooth.__doc__
        return

    # --------------------------------------------------------------------
    # Read in configuration, and extract the ones we need:
    # --------------------------------------------------------------------
    # --------------------------------------------------------------------
    # Calculate kappa_smooth
    experiment = pangloss.Configuration(configfile)

    zd = experiment.parameters['StrongLensRedshift']
    zs = experiment.parameters['SourceRedshift']
    
    grid = pangloss.Grid(zd,zs,nplanes=100)
    
    simcat = experiment.parameters['CalibrationCatalogs']
    simcat = pangloss.readCatalog(simcat,experiment)
    
    simcat.defineSystem(zd,zs)
    simcat.loadGrid(grid) 
    simcat.snapToGrid(grid)

    simcat.drawConcentrations(errors=True)
   
    kappa_slice = simcat.kappa_s
    z_slice = grid.redshifts
    
    plt.figure()
    plt.plot(z_slice, kappa_slice)
    plt.savefig('kappa_smooth.png',dpi=300)

    plt.show()
    
    return
    
    

# ======================================================================

if __name__ == '__main__': 
    RemoveSmooth(sys.argv[1:])

# ======================================================================    