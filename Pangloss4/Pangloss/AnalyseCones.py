import numpy
import numpy.random as rnd
import cPickle
from scipy import interpolate,optimize
import time
import glob
import copy

def AnalyseWeightingScheme(lc,surveyparams):
    print "Code up AnalyseWeightingScheme"
    #note this made need an Nrealisation term depending on if the weight scheme has uncertainties.
    exit()



def AnalyzeReconstruction(lc,zl,zs,Grid,behT,behI,Nrealisation,truncationscale,surveyparams):
    kappa_halo=numpy.empty(Nrealisation)
    gamma1_halo=numpy.empty(Nrealisation)
    gamma2_halo=numpy.empty(Nrealisation)

    kappa_Hilbert=lc.kappa_hilbert
    gamma1_Hilbert=lc.gamma1_hilbert
    gamma2_Hilbert=lc.gamma2_hilbert

    lc.define_system(zl,zs)
    #load grid
    lc.load_grid(Grid)

    #lc.snap_to_grid(Grid)
    #lc.drawConcentrations(errors=False)
    #lc.Make_kappas(truncationscale=10)
    #lc.Scale_kappas()
    #kappa_T=lc.kappa_add_total
    #gamma1_T=lc.gamma1_add_total
    #gamma2_T=lc.gamma2_add_total


    #Set up the survey parameters to give appropriate spectroscopic flag.
    print "set up the survey in Analyse Cones"
    exit()

    for j in range(Nrealisation):
        lc.mimik_photoz_error(sigma=)
        lc.snap_to_grid(Grid)
        lc.drawMStar(behI)
        lc.drawMHalo(behT)
        lc.drawConcentrations(errors=True)
        lc.Make_kappas(truncationscale=10)
        lc.Scale_kappas()
        kappa_halo[j]=lc.kappa_add_total
        gamma1_halo[j]=lc.gamma1_add_total
        gamma2_halo[j]=lc.gamma2_add_total

    results=kappa_Hilbert,gamma1_Hilbert,gamma2_Hilbert,kappa_halo,gamma1_halo,gamma2_halo#,kappa_T,gamma1_T,gamma2_T

    return results
