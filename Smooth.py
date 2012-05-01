#!/usr/bin/env python
# ======================================================================

import Pangloss

import sys,getopt,atpy,pyfits
import pylab,matplotlib.pyplot as plt
import numpy
import numpy.random as rnd
import distances
from scipy import optimize
from scipy import stats


arcmin2rad = (1.0/60.0)*numpy.pi/180.0
rad2arcmin = 1.0/arcmin2rad

# ======================================================================

def Smooth(zl,zs,radius,catalogues,truncationscale=3,magnitudecut=99,band='r',nplanes=200):
   lg=Pangloss.lensgrid(zl,zs,nplanes=200)
   lg.populatelensgrid()
   smoothcomponentindiv=numpy.zeros((lg.nplanes,len(catalogues)))
   smoothcomponent=numpy.zeros(lg.nplanes)


   conearea=3.14159*(radius*arcmin2rad)**2

   if band == "u" or band ==  "g" or band == "r" or band ==  "i" or band == "z":
      col = "mag_SDSS_%s" % band
       #elif band == "F814" or band == "F814W" or band == "814" or band == 814:
       #    col = "mag_F814W"
   else:
      col = "mag_%s" % band

   i=0
   for datafile in catalogues:
      cat = atpy.Table(datafile, type='ascii')
      print "Read in catalogue, length",len(cat)


      xmax = cat['pos_0[rad]'].max()
      xmin = cat['pos_0[rad]'].min()
      ymax = cat['pos_1[rad]'].max()
      ymin = cat['pos_1[rad]'].min() 
      catarea=(xmax-xmin)*(ymax-ymin) #in square radians (steradians?)

      normalizedarea=catarea/conearea #how many lightcones could fit in the catalogue?

      cat=cat.where(cat["%s"%col]<magnitudecut)


      M200 = cat['M_Subhalo[M_sol/h]']
      c200 = Pangloss.Relations.MCrelation(M200)
      cat.add_column('NetoC', c200)
      snappedz,planes=lg.snap(cat['z_spec'])
      cat.add_column("snappedplane",planes)
      cat.add_column("snapped_redshift",snappedz)

      rhocrit=lg.rho_crit_p[planes]
      r200 = (3*M200/(800*3.14159*rhocrit))**(1./3)
      cat.add_column('r200TRUE', r200)

      rs = r200/c200
      rhos = Pangloss.Relations.delta_c(c200)*rhocrit 

      R_trunc=truncationscale*r200
    
      mass=4*3.14159*rhos*(rs**3)  *     \
          (  numpy.log(1+(R_trunc)/rs) - \
                R_trunc/(rs+R_trunc)    \
                )
      cat.add_column('Mtrunc', mass)
      #print cat.Mtrunc

      for p in range(lg.nplanes):
         ontheplane=cat.where(cat.snappedplane==p)
         smoothcomponentindiv[p,i]=numpy.sum(ontheplane.Mtrunc)/normalizedarea
      i+=1

   for p in range(lg.nplanes):
      smoothcomponent[p]=numpy.sum(smoothcomponentindiv[p,:])/len(catalogues)
   
   lg.smoothcomponent=smoothcomponent
   lg.kappa=lg.smoothcomponent/lg.Sigma_Crit_p
   lg.kappakeeton=Pangloss.LensingFunc.KappaKeeton_beta(lg.beta_p,lg.kappa,0)

   plt.plot(lg.kappakeeton)
   plt.show()

   print numpy.sum(lg.kappakeeton)




d1= "../data/GGL_los_8_0_0_1_1_N_4096_ang_4_STARS_SA_galaxies_ANALYTIC_SA_galaxies_on_plane_27_to_63.images.txt"
datafile=[d1]
Smooth(0.6,1.4,1,datafile,truncationscale=10)
