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

def smooth(zl,zs,catalogues,truncationscale=10,magnitudecut=99,band='r',nplanes=100,hardcut="RVir",cosmo=[0.25,0.75,0.73],scaling="add",errors=True,grid=None,BehrooziSpline=None,eBer=1e-99,centralsonly=False, Mh2Mh=False, Mstar2Mh=False, perfectsatellites=False):
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
            if BehrooziSpline==None:
               HALOSTARlowz,STARHALOlowz,HALOSTARhighz,STARHALOhighz=Rel.Behroozi_Spline()
            else:
               HALOSTARlowz,STARHALOlowz,HALOSTARhighz,STARHALOhighz=BehrooziSpline[0],BehrooziSpline[1],BehrooziSpline[2],BehrooziSpline[3]

            M200 = Rel.Mhalo_to_Mhalo(cat['M_Subhalo[M_sol/h]'],cat['z_spec'],HALOSTARlowz,STARHALOlowz,HALOSTARhighz,STARHALOhighz,eBer=eBer)
      elif Mstar2Mh:
            Mstar=10**(numpy.log10(cat['M_Stellar[M_sol/h]'])+rnd.normal(0,eBer,len(cat.z_spec)))
            M200=Rel.Mstar_to_M200(Mstar,cat['z_spec'],scatter=False)

            #Mt=Rel.Mstar_to_M200(cat['M_Stellar[M_sol/h]'],cat['z_spec'],scatter=False)
            #print numpy.max(numpy.log10(Mt/M200)[Mt>10**13])


      else: M200 = cat['M_Subhalo[M_sol/h]']

      if perfectsatellites==True:
                    M200[cat.Type==1]=cat['M_Subhalo[M_sol/h]'][cat.Type==1]
                    M200[cat.Type==2]=cat['M_Subhalo[M_sol/h]'][cat.Type==2]

 



      #print numpy.std((numpy.log10(M200[M200>1e15]/cat['M_Subhalo[M_sol/h]'][M200>1e15])))



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

      #circular mass - not coded in yet.
      #mass=4*3.14159*rhos*(rs**3)  *     \
      #    (  numpy.log(1+(R_trunc)/rs) - \
      #          R_trunc/(rs+R_trunc)    \
      #          )

      xtrunc=R_trunc/rs
      print c200.max()
      print c200[M200>1e10].min()

      print xtrunc.max()
      print xtrunc[M200>1e11].min()

      mass=4*rs*rhos*(
         (2/(xtrunc**2-1)**.5)
         *
         numpy.arctan(((xtrunc-1.)/(xtrunc+1))**.5)
         +
         numpy.log(xtrunc/2)
         )

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
   

#d1= "../../data/GGL_los_8_0_0_1_1_N_4096_ang_4_STARS_SA_galaxies_ANALYTIC_SA_galaxies_on_plane_27_to_63.images.txt"
#datafile=[d1]
#Smooth(0.6,1.4,datafile,truncationscale=5,nplanes=20)
   #print "finished",time.clock()

#-------------------------------------------
if __name__ == '__main__':
  test1=True
  test2=False


  if test1 ==True:
   #print "running"
# Kappa Smooth as a function of magnitude cut:
   d1= "/data/tcollett/Pangloss/catalogs/GGL_los_8_7_7_0_0_N_4096_ang_4_SA_galaxies_on_plane_27_to_63.images.txt"
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

   #smoothlist0=smooth(0.6,1.4,datafile,truncationscale=3,nplanes=50,eBer=1e-99)

   smoothlist1=numpy.zeros(n)
   smoothlist2=numpy.zeros(n)
   smoothlist3=numpy.zeros(n)
   smoothlist4=numpy.zeros(n)
   for i in range(n):

      smoothlist1[i]=smooth(0.6,1.4,datafile,truncationscale=3,nplanes=20,eBer=0.18,scaling='add')
      print smoothlist1[i]
      #smoothlist2[i]=smooth(0.6,1.4,datafile,truncationscale=3,nplanes=50,eBer=0.19)
      #smoothlist3[i]=smooth(0.6,1.4,datafile,truncationscale=3,nplanes=50,eBer=0.20)
      #smoothlist4[i]=smooth(0.6,1.4,datafile,truncationscale=3,nplanes=50,eBer=0.21)
      if i%(n/10)==0:print i, "of", n

   a= numpy.std(smoothlist1)
   b= numpy.std(smoothlist2)
   c= numpy.std(smoothlist3)
   d= numpy.std(smoothlist4)

   bins=numpy.linspace(0,1.0,60)
   #plt.subplot(411)
   plt.hist(smoothlist1,bins,label="$\sigma\mathrm{[M*]}=0.15$")
   plt.legend(title="$\sigma(\kappa)$=%.3f"%a)
   #plt.subplot(412)
   #plt.hist(smoothlist2,bins,label="$\sigma\mathrm{[M*]}=0.19$")
   #plt.legend(title="$\sigma(\kappa)$=%.3f"%b)
   #plt.subplot(413)
   #plt.hist(smoothlist3,bins,label="$\sigma\mathrm{[M*]}=0.20$")
   #plt.legend(title="$\sigma(\kappa)$=%.3f"%c)
   #plt.subplot(414)
   #plt.hist(smoothlist4,bins,label="$\sigma\mathrm{[M*]}=0.21$")
   #plt.legend(title="$\sigma(\kappa)$=%.3f"%d)


   plt.xlabel("$\kappa_{\mathrm{smooth}}$")
   plt.savefig("smoothvsmasserror")
   plt.show()


  if test2 ==True:
   print "running"
