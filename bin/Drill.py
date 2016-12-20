#!/usr/bin/env python
# ======================================================================

import getopt
import numpy
import sys

import pangloss


# ======================================================================

def Drill(argv):
    """
    NAME
        Drill.py

    PURPOSE
        Read in a large catalog of galaxies (either observed or
        simulated), drill out Nc lightcones, at either random or
        user-specified positions, and write them to file as plain text
        catalogs.

    COMMENTS
        The lightcone catalogs propagate some catalog columns (M*, z,
        magnitudes and sky position) and adds some new ones: r, phi
        polar coordinates (in arcmin) relative to the lightcone centre.
        If the catalog is from a simulation, Mhalo is also propagated,
        and the true kappa value for each line of sight is also extracted
        and passed on. Lightcone catalogs are stored as pickles.

    FLAGS
        -h            Print this message [0]

    INPUTS
        configfile    Plain text file containing Pangloss configuration

    OUTPUTS
        stdout        Useful information
        pickle(s)     Lightcone catalog(s)

    EXAMPLE

        Drill.py example.config

    BUGS

    AUTHORS
        This file is part of the Pangloss project, distributed under the
        GPL v2, by Tom Collett (IoA) and  Phil Marshall (Oxford).
        Please cite: Collett et al 2013, http://arxiv.org/abs/1303.6564

    HISTORY
      2013-03-21 started Collett & Marshall (Oxford)
      2015-09-08 updated Everett (SLAC)
    """

    # --------------------------------------------------------------------

    try:
        opts, args = getopt.getopt(argv,"h",["help"])
    except getopt.GetoptError, err:
        print str(err) # will print something like "option -a not recognized"
        print Drill.__doc__  # will print the big comment above.
        return

    for o,a in opts:
        if o in ("-h", "--help"):
            print Drill.__doc__
            return
        else:
            assert False, "unhandled option"

    # Check for setup file in array args:
    if len(args) == 1:
        configfile = args[0]
        print pangloss.doubledashedline
        print pangloss.hello
        print pangloss.doubledashedline
        print "Drill: Drilling out lightcones from input catalogs "
        print "Drill: Taking instructions from",configfile
    else:
        print Drill.__doc__
        return

    # --------------------------------------------------------------------
    # Read in configuration, and extract the ones we need:

    experiment = pangloss.Configuration(configfile)

    Rc = experiment.parameters['LightconeRadius'] # in arcmin
    Nc = experiment.parameters['NCalibrationLightcones']

    # There should only be one calibration folder!
    CALIB_DIR = experiment.parameters['CalibrationFolder'][0]
    EXP_NAME = experiment.parameters['ExperimentName']

    calcats = experiment.parameters['CalibrationCatalogs']
    Ncalcats = len(calcats)
    kappamaps = experiment.parameters['CalibrationKappamaps']


    # There should only be one observed catalog!
    obscat = experiment.parameters['ObservedCatalog'][0]

    # Note nRA - -RA(rad) and Dec is also in rad...
    x0 = experiment.parameters['nRA']
    y0 = experiment.parameters['Dec']

    units = experiment.parameters['Units']

    # Write the observed lightcone to the same directory
    # as its parent catalog:
    obspickle = experiment.getLightconePickleName('real')


    makeNewCalCones=experiment.parameters['MakeNewCalibrations']


    # --------------------------------------------------------------------
    # First, make any calibration lightcones required:

    if (Nc > 0 and Ncalcats > 0 and makeNewCalCones!='False'):
      if makeNewCalCones!= 'false':
        print pangloss.dashedline
        print ("Drill: Making %i calibration lightcones in %i sky patches:" % (Nc,Ncalcats))

        flavor = 'simulated'

        count = 0
        Ncones = Nc/Ncalcats

        for i, catalog in enumerate(calcats):

            print "Drill: Reading in calibration catalog from "+catalog+"..."
            table = pangloss.read_hilbert_catalog(catalog, experiment)

            if units == 'deg':
                table['nRA'] = -table['nRA'] * pangloss.deg2rad
                table['Dec'] = table['Dec'] * pangloss.deg2rad
                table['Mhalo_obs'] = table['Mhalo_obs'] * 1E10
                table['Mstar_obs'] = table['Mstar_obs'] * 1E10
            ###
            #dx = Rc*1.000001*pangloss.arcmin2rad
            #subtable =table.where((table['pos_0[rad]'] > (-dx)) & \
            #                      (table['pos_0[rad]'] < (+dx)) & \
            #                      (table['pos_1[rad]'] > (-dx)) & \
            #                      (table['pos_1[rad]'] < (+dx))   )


            print "Drill: Sampling sky positions in",units,"..."
            x,y = sample_sky(table,Rc,Ncones,method='random')

            if kappamaps is not None:
                print "Drill: Reading in kappa map from "+kappamaps[i]
                MSconvergence = pangloss.Kappamap(kappamaps[i])

            # Coming soon...
            #   gammafile1 = gamma1[i]
            #   MSgamma1 = kappamap.Kappamap(gammafile1)
            #   gammafile2 = gamma2[i]
            #   MSgamma2 = kappamap.Kappamap(gammafile2)

            print "Drill: Pickling lightcones from this patch of sky..."
            for k in range(Ncones):
                if k % 200 == 0 and k !=0:
                    print ("Drill: ...on cone %i out of %i..." % (k,Ncones))

                lc = pangloss.Lightcone(table,'simulated',[x[k],y[k]],Rc)

                if kappamaps is not None:
                    lc.kappa_hilbert = MSconvergence.at(x[k],y[k],coordinate_system='physical')

                # Coming soon...
                #   lc.gamma1_hilbert = MSgamma1.at(x[k],y[k],coordinate_system='physical')
                #   lc.gamma2_hilbert = MSgamma2.at(x[k],y[k],coordinate_system='physical')

                calpickle = experiment.getLightconePickleName('simulated',pointing=count)
                pangloss.writePickle(lc,calpickle)

                count += 1


            # Save memory!
            del table
            del lc
            del catalog

            print "Drill: ...done."

        print ("Drill: All %i calibration lightcones made." % (count))

    # --------------------------------------------------------------------
    # Now, make any observed lightcones required:
    if obscat != 'none':
        if (len(obscat) > 0):

            print pangloss.dashedline
            print "Drill: Reading in observed catalog: "+obscat

            flavor = 'real'

            table = pangloss.read_hilbert_catalog(obscat,experiment)

            xc = [x0,y0]
            lc = pangloss.Lightcone(table,'real',xc,Rc)

            obspickle = experiment.getLightconePickleName('real')
            pangloss.writePickle(lc,obspickle)

            print "Drill: Observed lightcone pickled to "+obspickle

    # --------------------------------------------------------------------

    print pangloss.doubledashedline
    return

# ======================================================================
# Draw calibration sightlines at random:

def sample_sky(table,Rc,Nc,method='random'):

    xmax = table['nRA'].max()
    xmin = table['nRA'].min()
    ymax = table['Dec'].max()
    ymin = table['Dec'].min()

    Rcrad = Rc*pangloss.arcmin2rad

    if method == 'random':
        x = numpy.random.uniform(xmin+Rcrad,xmax-Rcrad,Nc)
        y = numpy.random.uniform(ymin+Rcrad,ymax-Rcrad,Nc)
    else:
        assert False
    return x,y


# ======================================================================

if __name__ == '__main__':
    Drill(sys.argv[1:])

# ======================================================================
