#!/usr/bin/env python
# ======================================================================

import LensingFunc as LF
import Relations as Rel
import LensingProfiles as LP
import grid as GRID


import sys,getopt,atpy,pyfits
import pylab,matplotlib.pyplot as plt
import numpy
import numpy.random as rnd
from scipy import optimize
from scipy import stats

#import time

arcmin2rad = (1.0/60.0)*numpy.pi/180.0
rad2arcmin = 1.0/arcmin2rad

# ======================================================================

def smooth(zl,zs,catalogues,truncationscale=10,magnitudecut=99,band='r',nplanes=100,hardcut="RVir",cosmo=[0.25,0.75,0.73],scaling="add",errors=True,grid=None,centralsonly=False, Mh2Mh=False, perfectsatellites=False):
   #print time.clock()
   if grid==None:
      lg=GRID.lensgrid(zl,zs,nplanes=nplanes,cosmo=cosmo)
      lg.populatelensgrid()
   else: lg=grid


   smoothcomponentindiv=numpy.zeros((lg.nplanes,len(catalogues)))
   smoothcomponent=numpy.zeros(lg.nplanes)
   
   #print "grid",time.clock() 


   if band == "u" or band ==  "g" or band == "r" or band ==  "i" or band == "z":
      col = "mag_SDSS_%s" % band
       #elif band == "F814" or band == "F814W" or band == "814" or band == 814:
       #    col = "mag_F814W"
   else:
      col = "mag_%s" % band


   i=0
   for datafile in catalogues:
      cat = datafile#atpy.Table(datafile, type='ascii')
      #print "cat read in", time.clock()
      cat=cat.where(cat["%s"%col] < magnitudecut)
      if centralsonly==True:
         cat=cat.where(cat.Type==0)
      #cat=cat.where(cat.Type!=2)

      xmax = cat['pos_0[rad]'].max()
      xmin = cat['pos_0[rad]'].min()
      ymax = cat['pos_1[rad]'].max()
      ymin = cat['pos_1[rad]'].min() 
      catarea=(xmax-xmin)*(ymax-ymin) #in square radians (steradians?)
      physicalarea=catarea*(lg.Da_p+1e-99)**2

      #Msub=cat['M_Subhalo[M_sol/h]']
      #Mhal=cat['M_Halo[M_sol/h]']
      #Msub[Msub==0]=Mhal[Msub==0]
      #cat.remove_columns("M_Subhalo[M_sol/h]")
      #cat.add_column("M_Subhalo[M_sol/h]",Msub)
      #print numpy.sum(cat["M_Subhalo[M_sol/h]"][cat.Type==0])
      #print numpy.sum(cat["M_Subhalo[M_sol/h]"][cat.Type==1])
      #print numpy.sum(cat["M_Subhalo[M_sol/h]"][cat.Type==2])




      cat['M_Subhalo[M_sol/h]'][cat['M_Subhalo[M_sol/h]']==0.0]=1.0

      if Mh2Mh:
         Mstar = self.grid.drawMstar(self.galaxies['M_Subhalo[M_sol/h]'],self.galaxies['z_spec'])

         Mstore=Mstar*1.0

         #Now give the Mstars a scatter for observational errors:
         Mstar = 10**(numpy.log10(Mstar)+rnd.normal(0,0.15,len(Mstar)))

         Mstar[numpy.log10(Mstar)>12.45]=10**12.46

         self.galaxies.add_column("Mstar",Mstar)
         #trim out very low stellar mass halos:
         cat=cat.where(numpy.log10(Mstar)>7)

         M200  = self.grid.drawMhalo(self.galaxies.Mstar,self.galaxies['z_spec'])
         

         # This should never print anything.
         if len(M200[numpy.isnan(M200)==True]!=0):
            M200[numpy.isnan(M200)==True]=self.grid.drawMhalo(10**(numpy.log10(self.galaxies.Mstar[numpy.isnan(M200)==True])),self.galaxies['z_spec'][numpy.isnan(M200)==True]+0.001)
            print "boo"
            if len(M200[numpy.isnan(M200)==True]!=0):
               print numpy.log10(Mstar[numpy.isnan(M200)==True])
               print numpy.log10(Mstore[numpy.isnan(M200)==True])


      else: M200 = cat['M_Subhalo[M_sol/h]']

      if perfectsatellites==True:
                    M200[cat.Type==1]=cat['M_Subhalo[M_sol/h]'][cat.Type==1]
                    M200[cat.Type==2]=cat['M_Subhalo[M_sol/h]'][cat.Type==2]

      #weed out stefans 'massless' halos.
      cat=cat.where(M200>1e1)
      M200=M200[M200>1e1]





      c200 = Rel.MCrelation(M200,MCerror=False)


      #cat.add_column('NetoC', c200)
      snappedz,planes=lg.snap(cat['z_spec'])
      #cat.add_column("snappedplane",planes)
      #cat.add_column("snapped_redshift",snappedz)


      rhocrit=lg.rho_crit_p[planes]
      r200 = (3*M200/(800*3.14159*rhocrit))**(1./3)
      #cat.add_column('r200TRUE', r200)

      rs = r200/c200
      rhos = LP.delta_c(c200)*rhocrit 

      if hardcut == "Rvir" or hardcut=="RVir" or hardcut == "r_vir" or hardcut == "rvir":
         R_trunc=truncationscale*r200
      elif hardcut == "rs" or hardcut=="Rs" or hardcut == "R_s" or hardcut == "r_s":
         R_trunc=truncationscale*rs
      elif hardcut== "subhalo":
         R_trunc=3*r200
         R_trunc[cat.Type == 1]=(r200[cat.Type == 1])*truncationscale
      else: print "what hardcut did you mean?"  

      #spherical mass - not coded in yet.
      xtrunc=R_trunc/rs

      sphmass=4*3.14159*rhos*(rs**3)  *     \
          (  numpy.log(1+xtrunc) - \
                xtrunc/(1 + xtrunc)    \
                )


      #cylindrical mass
      cylmass=4*rhos*3.14159*(rs**3)*(
         (2/(xtrunc**2-1)**.5)
         *
         numpy.arctan(((xtrunc-1.)/(xtrunc+1))**.5)
         +
         numpy.log(xtrunc/2)
         )


      mass=cylmass

      cat.add_column('Mtrunc', mass)
 

      for p in range(lg.nplanes):
         ontheplane=cat.where(planes==p)
         smoothcomponentindiv[p,i]=numpy.sum(ontheplane.Mtrunc)/physicalarea[p]
      i+=1

   for p in range(lg.nplanes):
      smoothcomponent[p]=numpy.sum(smoothcomponentindiv[p,:])/len(catalogues)
   
   lg.smoothcomponent=smoothcomponent
   lg.kappa=lg.smoothcomponent/lg.sigma_crit_p


   lg.kappaScaled=LF.KappaScale_beta(lg.beta_p,lg.kappa,0,scaling=scaling)

   cat.remove_columns('Mtrunc')
   
   return numpy.sum(lg.kappaScaled)
   

#-------------------------------------------
if __name__ == '__main__':
  test1=True
  test2=False


  if test1 ==True:
   #print "running"
# Kappa Smooth as a function of magnitude cut:
   d1= "/data/tcollett/Pangloss/catalogs/GGL_los_8_7_7_3_0_N_4096_ang_4_SA_galaxies_on_plane_27_to_63.images.txt"
   datafile=[atpy.Table(d1, type='ascii')]
   """
   maglist=numpy.linspace(19,25,15,endpoint=True)
   smoothlist=numpy.zeros(len(maglist))
   for i in range(len(maglist)):
      magcut=maglist[i]
      smoothlist[i]=smooth(0.6,1.4,datafile,truncationscale=5,magnitudecut=magcut,band="i",nplanes=50)
   
   plt.plot(maglist,smoothlist)
   plt.xlabel("i band depth")
   plt.ylabel("$\kappa_{\mathrm{smooth}}$")
   plt.title("Truncation at 5 R$\mathrm{_{Vir}}$, $z_l = 0.6$, $z_s = 1.4$")
   plt.savefig("../figure3.png")
   plt.show()
   """
   n=10

   smoothlist0=smooth(0.6,1.4,datafile,truncationscale=3,nplanes=200,errors=False)
   print smoothlist0

   """
   smoothlist1=numpy.zeros(n)
   smoothlist2=numpy.zeros(n)
   smoothlist3=numpy.zeros(n)
   smoothlist4=numpy.zeros(n)
   for i in range(n):

      smoothlist1[i]=smooth(0.6,1.4,datafile,truncationscale=3,nplanes=20,scaling='add')
      #print smoothlist1[i]
      if i%(n/10)==0:print i, "of", n

   a= numpy.std(smoothlist1)
   b= numpy.std(smoothlist2)
   c= numpy.std(smoothlist3)
   d= numpy.std(smoothlist4)

   bins=numpy.linspace(0,1.0,60)
   plt.hist(smoothlist1,bins,label="$\sigma\mathrm{[M*]}=0.15$")
   plt.legend(title="$\sigma(\kappa)$=%.3f"%a)
   plt.xlabel("$\kappa_{\mathrm{smooth}}$")
   plt.savefig("smoothvsmasserror")
   plt.show()
   """

  if test2 ==True:
   print "running"
