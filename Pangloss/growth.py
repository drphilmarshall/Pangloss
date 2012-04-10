#A program that outputs the properties of the most important galaxy (largest kappa_keeton) in m lightcones.

import lightcone
import pylab, atpy
import matplotlib.pyplot as plt
import numpy
import numpy.random as rnd
import distances
from mpl_toolkits.axes_grid1 import ImageGrid
from time import clock

D = distances.Distance()
arcmin2rad = (1.0/60.0)*numpy.pi/180.0
rad2arcmin = 1.0/arcmin2rad


datafile = "../../data/GGL_los_8_0_0_1_1_N_4096_ang_4_STARS_SA_galaxies_ANALYTIC_SA_galaxies_on_plane_27_to_63.images.txt"


catalog= atpy.Table(datafile, type='ascii')


for k in range(1):
    plt.clf()
    rmax = 5
    print rmax

    zs=2.0 #should be selecting a source redshift (possibly use forecaster from Collett et al. 2012)

#-0.00633207	-0.0173126
    xpos = -0.008
    ypos = -0.008
    zl =  1.0 
    xc = []
    #xc = [xpos,ypos,zl] #leave as [] to select a lens at random

    print "Initialising lightcone data..."
    lc = lightcone.lightcone(catalog,rmax,zs,lensindex=k+1,position=xc)


    #print "Distributing dark matter in halos..."
    #lc.make_dmhalos()
    #lc.plot(dmhalos=True)
 
    #print "Distributing stars in galaxies..."
    #lc.make_galaxies()
    #lc.plot(galaxies=True)
 
    #print "Computing Keeton (2003) convergence at optical axis due to each halo..."
    lc.make_kappa_contributions()
    i=(numpy.argsort(lc.galaxies.kappa_keeton))[-1]
    if k == 0:
        mostimportant=lc.galaxies.rows([i])
        mostimportant.add_column('kappa_total', numpy.sum(lc.galaxies.kappa_keeton))
    else: 
        dummy=lc.galaxies.rows([i])
        dummy.add_column('kappa_total', numpy.sum(lc.galaxies.kappa_keeton))
        mostimportant.append(dummy)


#mostimportant.write('mostimportant.table',type='ascii')
    

    #plt.scatter(lc.galaxies['z_spec'],(lc.galaxies.kappa_keeton/lc.galaxies.kappa),s=1,edgecolor='none')
    #plt.show()
    
#    print "Total external convergence =",numpy.sum(lc.galaxies.kappa_keeton)#-lc.kappa_expectation
    # Now make an illustrative plot:
    
    if rmax < 8:
       print "Plotting objects in lightcone..."
       lc.plot(starlight=False,dmglow=True,kappa_indiv=False,kappa_keeton=False,observed_light=True)
       

       pngfile = 'view_of_a_lightcone.png'
       plt.savefig(pngfile)
       print "Plot saved in",pngfile

    print "Plotting curve of growth..."
    plt.clf()
    lc.curve_of_growth(ordering="distance",starlight=False,dmglow=False,kappa_indiv=False,kappa_keeton=True,observed_light=True)
    pngfile = 'curve_of_growth.png'
    plt.savefig(pngfile)
    print "Plot saved in",pngfile
    plt.show()
    
