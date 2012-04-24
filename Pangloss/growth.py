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

xmax = catalog['pos_0[rad]'].max()
xmin = catalog['pos_0[rad]'].min()
ymax = catalog['pos_1[rad]'].max()
ymin = catalog['pos_1[rad]'].min()

rmax = 20

iterations = 1
xpos=rnd.uniform(xmin+rmax*arcmin2rad,xmax-rmax*arcmin2rad,iterations)
ypos=rnd.uniform(ymin+rmax*arcmin2rad,ymax-rmax*arcmin2rad,iterations)
xpos=[-0.00645898, -0.01021756, -0.01018463]
ypos=[-0.00366765, -0.00836147, -0.00995278]

print xpos,ypos

for k in range(iterations):
    hardcut=3
    plt.clf()

    zs=1.4 
    zl =  0.6

    xc = [xpos[k],ypos[k],zl]
    print "Initialising lightcone data..."
    lc = lightcone.lightcone(catalog,rmax,zs,position=xc)


    #print "Distributing dark matter in halos..."
    #lc.make_dmhalos()
    #lc.plot(dmhalos=True)
 
    #print "Distributing stars in galaxies..."
    #lc.make_galaxies()
    #lc.plot(galaxies=True)
 
    #print "Computing Keeton (2003) convergence at optical axis due to each halo..."

    lc.make_kappa_contributions(truncation="hard",truncationscale=10)
    i=(numpy.argsort(lc.galaxies.kappa_keeton))[-1]
    if k == 0:
        mostimportant=lc.galaxies.rows([i])
        mostimportant.add_column('kappa_total', numpy.sum(lc.galaxies.kappa_keeton))
    else: 
        dummy=lc.galaxies.rows([i])
        dummy.add_column('kappa_total', numpy.sum(lc.galaxies.kappa_keeton))
        mostimportant.append(dummy)

    """
    if rmax < 8:
       print "Plotting objects in lightcone..."
       lc.plot(starlight=False,dmglow=True,kappa_indiv=False,kappa_keeton=False,observed_light=True)
       

       pngfile = 'view_of_a_lightcone.png'
       plt.savefig(pngfile)
       print "Plot saved in",pngfile

    print "Plotting curve of growth..."
    plt.clf()
    lc.curve_of_growth(ordering="distance",starlight=False,dmglow=False,kappa_indiv=False,kappa_keeton=True,observed_light=True)
    pngfile = 'curve_of_growth_Cone%i_5_Rvir_BM0truncation.png'%(k+1)
    plt.title = '3_Rvir_truncation'
    plt.savefig(pngfile)
    print "Plot saved in",pngfile
    plt.show()
    """
    

    mc = lightcone.lightcone(catalog,rmax,zs,position=xc)
    mc.make_kappa_contributions(truncation="hard",truncationscale=5)

    nc = lightcone.lightcone(catalog,rmax,zs,position=xc)
    nc.make_kappa_contributions(truncation="hard",truncationscale=3)



    zeropoint=0.0
    zero=numpy.ones(1)*zeropoint

    plt.clf()
    plt.subplot(3,1,1)
    args=numpy.argsort((lc.galaxies.r))
    dist=numpy.concatenate((zero,numpy.take(lc.galaxies.r,args)))
    ordered_kappas=numpy.take(lc.galaxies.kappa_keeton,args) 
    cat=numpy.concatenate((zero,(ordered_kappas)))
    cumtot=numpy.cumsum(cat)

    plt.plot(dist,cumtot)
    plt.axhline(y=cumtot[-1], xmin=0, xmax=1000,color='black', ls='dotted')
    plt.axhline(y=0, xmin=0, xmax=1000,color='black', ls='solid')

    plt.subplot(3,1,2)
    args=numpy.argsort((mc.galaxies.r))
    dist=numpy.concatenate((zero,numpy.take(mc.galaxies.r,args)))
    ordered_kappas=numpy.take(mc.galaxies.kappa_keeton,args) 
    cat=numpy.concatenate((zero,(ordered_kappas)))
    cumtot=numpy.cumsum(cat)

    plt.plot(dist,cumtot)
    plt.axhline(y=cumtot[-1], xmin=0, xmax=1000,color='black', ls='dotted')
    plt.axhline(y=0, xmin=0, xmax=1000,color='black', ls='solid')


    plt.subplot(3,1,3)
    args=numpy.argsort((nc.galaxies.r))
    dist=numpy.concatenate((zero,numpy.take(nc.galaxies.r,args)))
    ordered_kappas=numpy.take(nc.galaxies.kappa_keeton,args) 
    cat=numpy.concatenate((zero,(ordered_kappas)))
    cumtot=numpy.cumsum(cat)

    plt.plot(dist,cumtot)
    plt.axhline(y=cumtot[-1], xmin=0, xmax=1000,color='black', ls='dotted')
    plt.axhline(y=0, xmin=0, xmax=1000,color='black', ls='solid')

    plt.xlabel('LOS distance (arcmin)')
    plt.ylabel('$\kappa_{ext}$ (cumulative)')
