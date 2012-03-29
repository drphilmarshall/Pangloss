import pylab, atpy
import matplotlib.pyplot as plt
import numpy
import numpy.random as rnd
import distances
from mpl_toolkits.axes_grid1 import ImageGrid
from time import clock
import lightcone

datafile = "../../data/GGL_los_8_0_0_1_1_N_4096_ang_4_STARS_SA_galaxies_ANALYTIC_SA_galaxies_on_plane_27_to_63.images.txt"
master = atpy.Table(datafile, type='ascii')
print "Read in master table, length",len(master)


plt.clf()
rmax = 2

zs=2.0 #should be selecting a source redshift (possibly use forecaster from Collett et al. 2012)


xpos = -0.00633
ypos = -0.01731
zl =  1.0 
xc = []
#xc = [xpos,ypos,zl] #leave as [] to select a lens at random

iter=30

for j in range(iter):
    lc = lightcone.lightcone(master,rmax,zs,lensindex=j,position=xc)

# calculate a 'blurred' kappa keeton.
    
    lc.make_kappa_contributions(deterministic=True)
    kappatrue=numpy.sum(lc.galaxies.kappa_keeton)
    print kappatrue


# calculate the 'real' kappa_keeton.
    iterations=100
    K=numpy.zeros(iterations)

    for iteration in range(iterations): 
        lc.make_kappa_contributions(deterministic=False,photozerror=True,halomasserror=False)
        K[iteration] = lc.kappa_blurred


    bins=numpy.arange(-0.26,0.29,0.001)
    plt.clf()
    plt.hist(K,bins,normed=True)
    plt.axvline(x=kappatrue, ymin=0,ymax=1, c='k', ls='dashed')
    plt.xlabel("$\kappa_{ext}$")
    plt.ylabel("pdf($\kappa_{ext}$)")
    plt.xlim([kappatrue-0.04,kappatrue+0.02])

#plt.title('Perfect information, but 5% error on z_perturbers, 30% error on M_Halo')
    pngfile = 'photozerror/Kappa_keeton_blurred_%i.png' %j
    plt.savefig(pngfile)
    
    print "Plot saved in",pngfile
