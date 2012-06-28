import pylab, atpy
import matplotlib.pyplot as plt
import numpy
import numpy.random as rnd
import distances
from mpl_toolkits.axes_grid1 import ImageGrid
from time import clock
import lightcone
import cPickle


D = distances.Distance()

arcmin2rad = (1.0/60.0)*numpy.pi/180.0
rad2arcmin = 1.0/arcmin2rad


datafile = "/data/tcollett/Pangloss/catalogs/GGL_los_8_7_7_0_0_N_4096_ang_4_SA_galaxies_on_plane_27_to_63.images.txt"
master = atpy.Table(datafile, type='ascii')

datafile1 = "/data/tcollett/Pangloss/catalogs/GGL_los_8_7_7_3_3_N_4096_ang_4_SA_galaxies_on_plane_27_to_63.images.txt"
datafile2 = "/data/tcollett/Pangloss/catalogs/GGL_los_8_7_7_0_3_N_4096_ang_4_SA_galaxies_on_plane_27_to_63.images.txt"
datafile3 = "/data/tcollett/Pangloss/catalogs/GGL_los_8_7_7_3_0_N_4096_ang_4_SA_galaxies_on_plane_27_to_63.images.txt"

master.append(atpy.Table(datafile1, type='ascii'))
master.append(atpy.Table(datafile2, type='ascii'))
master.append(atpy.Table(datafile3, type='ascii'))

master = master.where(master['M_Subhalo[M_sol/h]']!=0)
Mhalo=master['M_Subhalo[M_sol/h]']

Mstars=master['M_Stellar[M_sol/h]']


#plt.scatter(numpy.log10(Mhalo),numpy.log10(Mstars),s=1,edgecolor='none')
#plt.show()
#
#print numpy.min(Mhalo)
#print numpy.max(numpy.log10(Mhalo))
#plt.hist(numpy.log10(Mhalo))
#plt.show()
#plt.hist(numpy.log10(Mstars))
#plt.show()


breakpoint1=13
breakpoint2=13.5
Mh1 =numpy.linspace(10.,10.4,2,endpoint=False)
Mh2 = numpy.linspace(10.4,breakpoint1,50,endpoint=False)
Mh3 = numpy.linspace(breakpoint1,breakpoint2,5,endpoint=False)
Mh4=(numpy.linspace(breakpoint2,22,11,endpoint=True))

Mh=numpy.append(Mh1,numpy.append(Mh2,numpy.append(Mh3,Mh4)))

Ms = numpy.linspace(4.5,14,30)

redshiftbins=[0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,99]



Mhhist=numpy.empty((len(redshiftbins)-1,len(Mh)-1))
Mshist=numpy.empty((len(redshiftbins)-1,len(Ms)-1))
for i in range(len(redshiftbins)-1):


    cat=master.where((master.z_spec>redshiftbins[i]) & (master.z_spec<redshiftbins[i+1]))
    Mhalo=numpy.log10(cat['M_Subhalo[M_sol/h]'])
    Mstars=numpy.log10(cat['M_Stellar[M_sol/h]'])
    Mhhist[i,:],seph= numpy.histogram(Mhalo,bins=Mh)
    Mshist[i,:],seps= numpy.histogram(Mstars,bins=Ms)
    #correct for density - numpy.histogram is broken, DONT use normed=True.
    for j in range(len(seph)-1):
        Mhhist[i,j]/=Mh[j+1]-Mh[j]
    for j in range(len(seps)-1):
        Mshist[i,j]/=Ms[j+1]-Ms[j]
    

M=Mh,Ms,redshiftbins,Mhhist,Mshist
MHMS=open("MHMS.data","wb")
cPickle.dump(M,MHMS)
