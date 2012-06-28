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
import cPickle, copy

arcmin2rad = (1.0/60.0)*numpy.pi/180.0
rad2arcmin = 1.0/arcmin2rad

N_45mean = 41.8
N_60mean = 74.475
N_90mean = 166.864
N_30mean = 18.831
N_15mean = 4.812

# ======================================================================

def MScompare(argv):
   """
   NAME
     MScompare.py

   PURPOSE
     Compare reconstructed convergence with "true" convergence for a set of 
     lightcones. Truth comes from Stefan Hilbert's ray tracing work in 
     the Millenium Simulation, in the form of a pixellated convergence
     map. Reconstructions are made from simulated catalogs in the same
     field.
   COMMENTS

   FLAGS
     -u            Print this message [0]

   INPUTS
     catalog       Catalog containing halos, stellar mass distributions etc
     kappafile     FITS file containing true kappa values.

   OPTIONAL INPUTS
     -n Ncones     No. of cones to compare
     -r Rcone      Cone radius used in reconstruction
     -zd zd        Deflector redshift (to match kappafile)
     -zs zs        Source redshift (to match kappafile)
     -t truncationscale number of virial radii to truncate NFW halos at.

   OUTPUTS
     stdout        Useful information
     pngfiles      Output plots in png format

   EXAMPLES

     MScompare.py -n 1000 catalog.txt kappa.fits

   BUGS

   AUTHORS
     This file is part of the Pangloss project.
     Copyright 2012 Tom Collett (IoA) and Phil Marshall (Oxford).
     
   HISTORY
     2012-04-06 started Marshall (Oxford)
   """

   # --------------------------------------------------------------------

   try:
      opts, args = getopt.getopt(argv,"hvn:r:",["help","verbose","zd","zs"])
   except getopt.GetoptError, err:
      print str(err) # will print something like "option -a not recognized"
      print MScompare.__doc__  # will print the big comment above.
      return

   vb = False
   Ncones = 1000
   Rcone = 3 # arcmin
   truncationscale=3   # *R_200 halo truncation
  
   scaling = "add" 
   #scaling = "tom"
   # Defaults are for B1608 (CHECK):
   zl = 0.62
   zs = 1.39
  
   for o,a in opts:
      if o in ("-v", "--verbose"):
         vb = True
      elif o in ("-n"):
         Ncones = int(a)
      elif o in ("-r"):
         Rcone = a
      elif o in ("-t"):
         truncationscale = a 

      elif o in ("--zl"):
         zl = a
      elif o in ("--zs"):
         zs = a
      elif o in ("-h", "--help"):
         print MScompare.__doc__
         return
      else:
         assert False, "unhandled option"

   # Check for datafiles in array args:
   if len(args) == 2:
      catalog = args[0]
      kappafile = args[1]
      if vb:
         print "Reconstructing convergence in lightcones from:",catalog
         print "Comparing to convergence in:",kappafile
         print "Number of lightcones to be reconstructed:",Ncones
   else:
      print MScompare.__doc__
      return

   # --------------------------------------------------------------------

   # Read in master catalog, and kappa map:
    
   master = atpy.Table(catalog, type='ascii')
   if vb: print "Read in master catalog, length",len(master)
   print "Read in master catalog, length",len(master)
   xmax = master['pos_0[rad]'].max()
   xmin = master['pos_0[rad]'].min()
   ymax = master['pos_1[rad]'].max()
   ymin = master['pos_1[rad]'].min()
   if vb: print "Catalog extent:",xmin,xmax,ymin,ymax
    
   MSconvergence = Pangloss.kappamap(kappafile)
   if vb: print "Read in true kappa map, dimension",MSconvergence.NX

   #build grid, populate it and invert the behroozi relation on each plane.
   grid=Pangloss.lensgrid(zl,zs,nplanes=100,cosmo=[0.25,0.75,0.73])
   grid.populatelensgrid()
   FILE=open("/home/tcollett/Pangloss/Pangloss/MHMS.data")
   MFs=cPickle.load(FILE)
   grid.Behroozigrid(MFs)


   # --------------------------------------------------------------------

   # Generate Ncones random positions
   
   x = rnd.uniform(xmin+Rcone*arcmin2rad,xmax-Rcone*arcmin2rad,Ncones)
   y = rnd.uniform(ymin+Rcone*arcmin2rad,ymax-Rcone*arcmin2rad,Ncones)
   

   #------------------------------------------------------------------
   kappa_Scaled = numpy.zeros(Ncones)
   kappa_Scaled_truth = numpy.zeros(Ncones)
   kappa_tom = numpy.zeros(Ncones)
   kappa_hilbert = numpy.zeros(Ncones)
   N_45 = numpy.zeros(Ncones)
   N_60 = numpy.zeros(Ncones)
   N_90 = numpy.zeros(Ncones)
   N_30 = numpy.zeros(Ncones)
   N_15 = numpy.zeros(Ncones)
   other = numpy.zeros(Ncones)
   other2 = numpy.zeros(Ncones)
   K_0 = numpy.zeros(Ncones)
   K_1 = numpy.zeros(Ncones)
   K_2 = numpy.zeros(Ncones)
   N_0 = numpy.zeros(Ncones)
   N_1 = numpy.zeros(Ncones)
   N_2 = numpy.zeros(Ncones)
   M_0 = numpy.zeros(Ncones)
   M_1 = numpy.zeros(Ncones)
   M_2 = numpy.zeros(Ncones)
   large005 = numpy.zeros(Ncones)
   large01 = numpy.zeros(Ncones)
   large02 = numpy.zeros(Ncones)
   kap=[]
   dist=[]
   Rbins=numpy.linspace(0,Rcone,30,endpoint=True)
   kappa_Scaled_R=numpy.zeros((len(Rbins),Ncones))
   delta_kappa_R=numpy.zeros((len(Rbins),Ncones))


   errors=True

   Mstar2Mh=False
   Mh2Mh=True
   perfectsatellites=False
   centralsonly=False


#  reconstruct kappa_scaled and look up true kappa_hilbert at that position:
   smooth=False
   mean0=True

   if mean0==True and smooth == True: print "make your mind up! smooth or mean=0?"

   kappa_empty = 0
   kappa_empty_truth = 0

   if smooth==True:
      kappa_empty_truth = Pangloss.smooth(zl,zs,[master],truncationscale=truncationscale,hardcut="RVir",scaling=scaling,grid=grid,errors=False,Mstar2Mh=False,perfectsatellites=True,Mh2Mh=False,centralsonly=False)

   if errors==False and smooth == True:
         kappa_empty = Pangloss.smooth(zl,zs,[master],truncationscale=truncationscale,hardcut="RVir",scaling=scaling,errors=errors,grid=grid,centralsonly=centralsonly,Mstar2Mh=Mstar2Mh,perfectsatellites=perfectsatellites,Mh2Mh=Mh2Mh)

   for k in range(Ncones):
      if k % 50 == 0: print ("evaluating cone %i of %i" %(k,Ncones))
      xc = [x[k],y[k]]

      # Truth:
      kappa_hilbert[k] = MSconvergence.at(x[k],y[k],coordinate_system='physical')
      # THE CATALOGUES NEED TRANSPOSING!!!! (don't make that mistake again!)

      # Reconstruction:
      lc = Pangloss.lens_lightcone(master,Rcone,xc,zl,zs,grid=grid)
      col="mag_SDSS_i"
      magcutcat=lc.galaxies.where((lc.galaxies["%s"%col] < 22))
      other[k]=numpy.min(magcutcat.rphys)
      other2[k]=numpy.min(lc.galaxies.rphys)

      lc.make_kappa_contributions(hardcut="RVir",truncationscale=truncationscale,scaling=scaling,errors=errors,centralsonly=centralsonly,perfectsatellites=perfectsatellites,Mh2Mh=Mh2Mh,Mstar2Mh=Mstar2Mh)

      tc= Pangloss.lens_lightcone(master,Rcone,xc,zl,zs,grid=grid)
      tc.make_kappa_contributions(hardcut="RVir",truncationscale=truncationscale,scaling=scaling,errors=False,centralsonly=False,perfectsatellites=True,Mh2Mh=False,Mstar2Mh=False)

      if errors==True and smooth == True:
         kappa_empty = Pangloss.smooth(zl,zs,[master],truncationscale=truncationscale,hardcut="RVir",scaling=scaling,errors=errors,grid=grid,centralsonly=centralsonly,Mstar2Mh=Mstar2Mh,perfectsatellites=perfectsatellites,Mh2Mh=Mh2Mh)


      kappa_Scaled[k] = (lc.kappa_Scaled_total-kappa_empty)
      kappa_Scaled_truth[k]=tc.kappa_Scaled_total-kappa_empty_truth

      #plt.scatter(lc.galaxies.M200,tc.galaxies.M200)
      #plt.show()

      #print  lc.kappa_Scaled_total,1
      #if  numpy.isnan(lc.kappa_Scaled_total)==True:
      #   jg =lc.galaxies.where(numpy.isnan(lc.galaxies.M200)==True)
      #   print jg.Mstar,jg.z_spec
      #print  tc.kappa_Scaled_total,2


      K_0[k] = numpy.sum(lc.galaxies.kappa_Scaled[lc.galaxies.Type==0])/numpy.sum(lc.galaxies.kappa_Scaled)
      K_1[k] = numpy.sum(lc.galaxies.kappa_Scaled[lc.galaxies.Type==1])/numpy.sum(lc.galaxies.kappa_Scaled)
      K_2[k] = numpy.sum(lc.galaxies.kappa_Scaled[lc.galaxies.Type==2])/numpy.sum(lc.galaxies.kappa_Scaled)


      N_0[k] = numpy.size(lc.galaxies.Type[lc.galaxies.Type==0])*1.0/numpy.size(lc.galaxies.Type)
      N_1[k] = numpy.size(lc.galaxies.Type[lc.galaxies.Type==1])*1.0/numpy.size(lc.galaxies.Type)
      N_2[k] = numpy.size(lc.galaxies.Type[lc.galaxies.Type==2])*1.0/numpy.size(lc.galaxies.Type)

      N_1[k]+=N_2[k]


      M_0[k] = numpy.sum(lc.galaxies['M_Subhalo[M_sol/h]'][lc.galaxies.Type==0])*1.0/numpy.sum(lc.galaxies['M_Subhalo[M_sol/h]'])
      M_1[k] = numpy.sum(lc.galaxies['M_Subhalo[M_sol/h]'][lc.galaxies.Type==1])*1.0/numpy.sum(lc.galaxies['M_Subhalo[M_sol/h]'])
      M_2[k] = numpy.sum(lc.galaxies['M_Subhalo[M_sol/h]'][lc.galaxies.Type==2])*1.0/numpy.sum(lc.galaxies['M_Subhalo[M_sol/h]'])


      large005[k]= numpy.size(lc.galaxies.kappa[lc.galaxies.kappa>0.002])*1.0
      large01[k] = numpy.size(lc.galaxies.kappa[lc.galaxies.kappa>0.005])*1.0
      large02[k] = numpy.size(lc.galaxies.kappa[lc.galaxies.kappa>0.01])*1.0


      for zaz in range(len(lc.galaxies.r)):
         dist.append(lc.galaxies.r[zaz])
         kap.append(lc.galaxies.kappa[zaz])

      # calculate for smaller Rcone.

 
      for j in range(len(Rbins)):
         mc=lc.galaxies.where(lc.galaxies.r<Rbins[j])
         kappa_Scaled_R[j,k]=numpy.sum(mc.kappa_Scaled)-kappa_empty
         delta_kappa_R[j,k]=numpy.sum(mc.kappa_Scaled)-kappa_empty-kappa_hilbert[k]

  

   if mean0==True:
      aleph=numpy.mean(kappa_Scaled)
      print numpy.mean(kappa_Scaled)
      kappa_Scaled-=numpy.mean(kappa_Scaled)
      print numpy.mean(kappa_Scaled)
      kappa_Scaled_truth-=numpy.mean(kappa_Scaled_truth)
      for j in range(len(Rbins)):
          kappa_Scaled_R[j]-=aleph
          delta_kappa_R[j]-=aleph

 
   print numpy.mean(kappa_hilbert)
   print numpy.mean(kappa_Scaled)

   #for j in range(len(Rbins)):
   #      plt.hist(kappa_hilbert)
   #      plt.hist(kappa_Scaled_R[j,:])
   #      plt.show()



   # --------------------------------------------------------------------
   # Basic statistics of two arrays:
   
   difference = kappa_Scaled - kappa_hilbert
   bias = numpy.average(difference)
   scatter = numpy.std(difference)
   #difference2 = kappa_tom - kappa_hilbert
   #scatter2 = numpy.std(difference2)
   #bias2 = numpy.average(difference2)

   y=difference
   x=kappa_Scaled

   #scatterprime = numpy.std(kappa_hilbert)
   #print scatterprime

   print "$\kappa_{\mathrm{Scaled}}-\kappa_{\mathrm{Hilbert}}$ = ",bias,"+/-",scatter

   bias_R=numpy.zeros(len(Rbins))
   scatter_R=numpy.zeros(len(Rbins))
   centre_R=numpy.zeros(len(Rbins))
   upper_R=numpy.zeros(len(Rbins))
   lower_R=numpy.zeros(len(Rbins))


   for j in range(len(Rbins)):
      scatter_R[j]=numpy.std(delta_kappa_R[j,:])
      bias_R[j]=numpy.mean(delta_kappa_R[j,:])

      centre_R[j]=numpy.median(delta_kappa_R[j,:])
      upper_R[j]=stats.scoreatpercentile(delta_kappa_R[j,:],16)
      lower_R[j]=stats.scoreatpercentile(delta_kappa_R[j,:],84)


   FILE = open("radcutbias.txt", 'w')
   FILE2 = open("radcutscatter.txt", 'w')
   cPickle.dump(bias_R , FILE )
   cPickle.dump(scatter_R , FILE2)

#========================================================================
   #Now lots of plotting routines
#========================================================================  
   figprops = dict(figsize=(10.0, 10.0), dpi=128)                                          # Figure properties
   fig = pylab.figure(**figprops)

  # Need small space between subplots to avoid deletion due to overlap...
   adjustprops = dict(\
         left=0.1,\
         bottom=0.1,\
         right=0.95,\
         top=0.95,\
         wspace=0.1,\
         hspace=0.2)
   fig.subplots_adjust(**adjustprops)

  # Font sizes:
   params = { 'axes.labelsize': 10,
              'text.fontsize': 10,
            'legend.fontsize': 8,
            'xtick.labelsize': 8,
            'ytick.labelsize': 8}
   pylab.rcParams.update(params)
  
   

   bins = numpy.linspace(0,1,20)
   #plt.subplot(331)
   plt.subplot(311)
   plt.hist(K_0,bins,alpha=0.5,fc='r',hatch=' ',label="Central galaxies")
   #plt.hist(K_2,bins,alpha=0.5,fc='g',hatch=' ',label="Non-Central main galaxies")
   plt.hist(K_1,bins,alpha=0.5,fc='b',label="Satellite galaxies")
   plt.xlabel('$\kappa$ fraction')
   plt.legend(loc=9)


   #plt.subplot(332)
   #plt.hist(K_0[kappa_Scaled>0.03],bins,fc='b',alpha=0.5)
   #plt.hist(K_2[kappa_Scaled>0.03],bins,fc='g',alpha=0.5)
   #plt.hist(K_1[kappa_Scaled>0.03],bins,fc='r',alpha=0.5)
   #plt.xlabel('$\kappa > 0.03$')

   #plt.subplot(333)
   #plt.hist(K_0[kappa_Scaled<0.0],bins,fc='b',alpha=0.5)
   #plt.hist(K_2[kappa_Scaled<0.0],bins,fc='g',alpha=0.5)
   #plt.hist(K_1[kappa_Scaled<0.0],bins,fc='r',alpha=0.5)
   #plt.xlabel('$\kappa < 0.0$')


   #plt.subplot(334)
   plt.subplot(312)
   plt.hist(N_0,bins,fc='r',alpha=0.5,hatch='/ ')
   #plt.hist(N_2,bins,fc='g',alpha=0.5,hatch='/ ')
   plt.hist(N_1,bins,fc='b',alpha=0.5,hatch='/ ')
   plt.xlabel('Number fraction')

   #plt.subplot(335)
   #plt.hist(N_0[kappa_Scaled>0.03],bins,fc='b',alpha=0.5,hatch='/ ')
   #plt.hist(N_2[kappa_Scaled>0.03],bins,fc='g',alpha=0.5,hatch='/ ')
   #plt.hist(N_1[kappa_Scaled>0.03],bins,fc='r',alpha=0.5,hatch='/ ')
   #plt.xlabel('$\kappa > 0.03$')

   #plt.subplot(336)
   #plt.hist(N_0[kappa_Scaled<0.0],bins,fc='b',alpha=0.5,hatch='/ ')
   #plt.hist(N_2[kappa_Scaled<0.0],bins,fc='g',alpha=0.5,hatch='/ ')
   #plt.hist(N_1[kappa_Scaled<0.0],bins,fc='r',alpha=0.5,hatch='/ ')
   #plt.xlabel('$\kappa < 0.0$')

   #plt.subplot(337)
   plt.subplot(313)
   plt.hist(M_0,bins,fc='r',alpha=0.5,hatch='\ ')
   #plt.hist(M_2,bins,fc='g',alpha=0.5,hatch='\ ')
   plt.hist(M_1,bins,fc='b',alpha=0.5,hatch='\ ')
   plt.xlabel('Mass fraction')


   #plt.subplot(338)
   #plt.hist(M_0[kappa_Scaled>0.03],bins,fc='b',alpha=0.5,hatch='\ ')
   #plt.hist(M_2[kappa_Scaled>0.03],bins,fc='g',alpha=0.5,hatch='\ ')
   #plt.hist(M_1[kappa_Scaled>0.03],bins,fc='r',alpha=0.5,hatch='\ ')
   #plt.xlabel('$\kappa > 0.03$')

   #plt.subplot(339)
   #plt.hist(M_0[kappa_Scaled<0.0],bins,fc='b',alpha=0.5,hatch='\ ')
   #plt.hist(M_2[kappa_Scaled<0.0],bins,fc='g',alpha=0.5,hatch='\ ')
   #plt.hist(M_1[kappa_Scaled<0.0],bins,fc='r',alpha=0.5,hatch='\ ')
   #plt.xlabel('$\kappa < 0.0$')
   
   plt.savefig("contributions.png")
   plt.show()
   plt.clf()

   plt.subplot(111)
   plt.scatter(K_0,difference,c='b',s=2,edgecolor='none',label="Central galaxies in main halo")
   plt.scatter(K_2,difference,c='g',s=2,edgecolor='none',label="Non-Central main galaxies")
   #plt.scatter(K_1,difference,c='r',s=2,edgecolor='none',label="Satellite galaxies")
   plt.xlim([0,1])
   plt.xlabel('$\kappa$ fraction')
   plt.legend(loc=9)

   plt.savefig("contributionsvsdelta.png")
   plt.show()
   plt.clf()




   """
   plt.plot(Rbins,scatter_R,c='k',label="scatter")
   plt.xlabel("Radius out to which LoS is reconstructed(arcmin)")
   plt.ylabel("$\kappa_{\mathrm{TC}}-\kappa_{\mathrm{Hilbert}}$")
   plt.plot(Rbins,bias_R,c='k', ls = 'dashed',label="bias")
   plt.text(3,0.35,"Truncation at 5 R_Vir")
   plt.axhline(y=0,xmin=0,xmax=1,ls="dotted")
   plt.legend(loc=4,title ="%i LoS"%Ncones)
   plt.savefig("DeltaKappa_vs_RconeOLD.png")
   plt.show()
   """
   
   plt.xlabel("Radius out to which LoS is reconstructed(arcmin)")
   plt.ylabel("$\kappa_{\mathrm{%s}}-\kappa_{\mathrm{Hilbert}}$"%scaling)
   plt.plot(Rbins,centre_R,c='k',label="Median")
   plt.plot(Rbins,upper_R,c='k', ls = 'dashed',label="84th percentile")
   plt.plot(Rbins,lower_R,c='k', ls = 'dashed',label="16th percentile")
   plt.plot(Rbins,bias_R,c='k', ls = '-.',label="Mean")
   plt.text(3,0.35,"Truncation at 3 R_Vir")
   plt.axhline(y=0,xmin=0,xmax=1,ls="dotted")
   plt.legend(loc=4,title ="%i LoS"%Ncones)
   plt.savefig("DeltaKappa_vs_Rcone.png")
   plt.show()
   
   plt.clf()
   bins = numpy.linspace(0,15,16)
   plt.hist(large005,bins,alpha=0.5,fc='b',label="N_gal with $\kappa>0.002$")
   plt.hist(large01 ,bins,hatch='///',fc='none',label="N_gal with $\kappa>0.005$")
   plt.hist(large02 ,bins,alpha=0.5,fc='r',label="N_gal with $\kappa>0.01$")
   plt.xlabel('Number of important galaxies')
   plt.legend(loc=9)
   plt.savefig("howmanyhalosareimportant.png")
   plt.show()
   plt.clf()
   plt.scatter(large005,kappa_Scaled,edgecolor='none',c='k',s=1)
   plt.show()

   #print len(dist)
   #plt.clf()
   #plt.scatter(dist,kap,edgecolor='none',c='k',s=0.5)
   #plt.show()
   bins=numpy.linspace(-5,1,30)
   plt.hist(numpy.log10(kap),bins,alpha=0.5)
   plt.title("Individual log($\kappa$) contributions")
   plt.ylabel("Number of halos")
   plt.xlabel("log($\kappa$)")
   nnn=0
   for i in range(len(kap)):
      if kap[i] > 0: nnn+=1
   plt.figtext(0.5,0.8," %.2f percent of halos within 3\' \n contribute no convergence" \
                  %  (100-(nnn*100./len(kap)))\
                  )
   plt.savefig("kaphist.png")
   plt.show()

   """
   list=numpy.linspace(-0.05,0.25,30)
   plt.subplot(311)
   plt.title("$\kappa_{\mathrm{TC}}-\kappa_{\mathrm{Hilbert}}$ = %.3f +/- %.3f" % (bias2,scatter2))
   plt.hist(kappa_tom, bins=list,normed=True, label="$\kappa_{\mathrm{TC}}$")
   plt.legend(title='Cut at %.0f R_vir.'%truncationscale, loc=1)
   plt.subplot(312)
   plt.hist(kappa_hilbert, bins=list,normed=True,label="$\kappa_{\mathrm{Hilbert}}$")
   plt.legend(loc=1)
   plt.subplot(313)
   plt.hist(difference2, bins=20,normed=True,label="$\kappa_{\mathrm{TC}}-\kappa_{\mathrm{Hilbert}}$")
   plt.legend(loc=2)
   plt.savefig("kappatom.png")
   plt.show()
   """

   plt.clf()
   list=numpy.linspace(-0.25,0.25,30)
   plt.subplot(311)
   #plt.title("$\kappa_{\mathrm{%s}}-\kappa_{\mathrm{Hilbert}}$ = %.3f +/- %.3f" % (scaling,bias,scatter))
   plt.hist(kappa_Scaled, bins=list,alpha=0.8,normed=True, label="re-inferred halo mass")
   plt.hist(kappa_Scaled_truth, bins=list,alpha=0.3,normed=True, label="perfect halo mass")

   plt.legend(title='$\kappa_{\mathrm{reconstruction}}$', loc=1)
   plt.subplot(312)
   plt.hist(kappa_hilbert, bins=list,normed=True,label="$\kappa_{\mathrm{Hilbert}}$",color='r',alpha=0.5)
   plt.legend(loc=1)
   plt.subplot(313)
   plt.hist(difference, bins=list,alpha=0.8,normed=True,label="re-inferred halo mass")
   plt.hist(kappa_Scaled_truth-kappa_hilbert, bins=list,alpha=0.3,normed=True,label="perfect halo mass")

   plt.legend(loc=1,title="difference")
   plt.savefig("comparison.png")
   plt.show()
   


   """
   plt.subplot(121)
   plt.scatter(difference,other2,s=1,c='k',edgecolor='none', label="Any halo")
   plt.xlabel("$\kappa_{\mathrm{Keeton}}-\kappa_{\mathrm{Hilbert}}$")
   #plt.ylabel("LOS distance to most important object (arcmin)")
   plt.ylabel("Closest halo to LOS (arcmin)")
   plt.xlim([-0.4,0.2])
   plt.subplot(121).xaxis.set_ticklabels(["",-0.3,"",-0.1,"",0.1])
   plt.ylim([0,0.8])
   plt.legend()

   plt.subplot(122)
   plt.scatter(difference,other,s=1,c='g',edgecolor='none', label="Halo with i<22")
   plt.legend()
   plt.xlabel("$\kappa_{\mathrm{Keeton}}-\kappa_{\mathrm{Hilbert}}$")
   plt.xlim([-0.4,0.2])
   plt.ylim([0,0.8])
   plt.subplot(122).xaxis.set_ticklabels(["",-0.3,"",-0.1,"",0.1])

   plt.savefig("other.png")
   plt.show()  
   
   plt.subplot(121)
   plt.scatter(difference,other2,s=1,c='k',edgecolor='none', label="Any halo")
   plt.xlabel("$\kappa_{\mathrm{TC}}-\kappa_{\mathrm{Hilbert}}$")
   plt.ylabel("Closest halo to LOS (Mpc)")
   #plt.xlim([-0.4,0.2])
   #plt.subplot(121).xaxis.set_ticklabels(["",-0.3,"",-0.1,"",0.1])
   #plt.ylim([0,0.8])
   plt.legend()

   plt.subplot(122)
   plt.scatter(difference,other,s=1,c='g',edgecolor='none', label="Halo with i<22")
   plt.legend()
   plt.xlabel("$\kappa_{\mathrm{TC}}-\kappa_{\mathrm{Hilbert}}$")
   #plt.xlim([-0.4,0.2])
   #plt.ylim([0,0.8])
   #plt.subplot(122).xaxis.set_ticklabels(["",-0.3,"",-0.1,"",0.1])

   plt.savefig("other.png")
   plt.show()    
   
   plt.clf()
   plt.subplot(111)
   plt.scatter(kappa_Scaled,difference,s=1,c='k',edgecolor='none')
   plt.ylabel("$\kappa_{\mathrm{%s}}-\kappa_{\mathrm{Hilbert}}$"%scaling)
   plt.xlabel("$\kappa_{\mathrm{%s}}$"%scaling)
   plt.xlim([-0.08,0.25])
   plt.ylim([-0.2,0.1])
   plt.savefig("Fig3.png")

   
   """
   plt.clf()
   plt.subplot(111)
   plt.scatter(kappa_hilbert,kappa_Scaled,s=1,c='k',edgecolor='none')
   plt.ylabel("$\kappa_{\mathrm{%s}}$"%scaling)
   plt.xlabel("$\kappa_{\mathrm{Hilbert}}$")
   plt.savefig("Fig1.png")
   plt.show()  
   """

   #print numpy.std(scatter)
   plt.clf()
   plt.subplot(111)
   z=numpy.linspace(-0.05,0.2,30)
   plt.plot(z, fitfunc(pfinal,z))
   plt.scatter(x,y,s=1,c='k',edgecolor='none')
   plt.ylabel("$\kappa_{\mathrm{Keeton}}-\kappa_{\mathrm{Hilbert}}$")
   plt.xlabel("$\kappa_{\mathrm{Keeton}}$")
   plt.xlim([-0.08,0.25])
   plt.ylim([-0.2,0.1])
   plt.savefig("Fig4.png")
   plt.show() 
   


   plt.clf()
   plt.subplot(111)
   z=numpy.linspace(-0.05,0.2,30)
   plt.scatter(N_45,y,s=1,c='k',edgecolor='none')
   plt.ylabel("$\kappa_{\mathrm{Keeton}}-\kappa_{\mathrm{Hilbert}}$")
   plt.xlabel("N$_{45}$")
   plt.savefig("Fig5.png")
   plt.show() 

   plt.clf()
   plt.subplot(111)
   z=numpy.linspace(-0.05,0.2,30)
   plt.scatter(N_45,kappa_hilbert,s=1,c='k',edgecolor='none')
   plt.ylabel("$\kappa_{\mathrm{Hilbert}}$")
   plt.xlabel("N$_{45}$")
   plt.savefig("Fig6.png")
   plt.show()
   


   plt.clf()
   plt.subplot(111)
   z=numpy.linspace(-0.05,0.2,30)
   plt.scatter(N_60,kappa_hilbert,s=1,c='k',edgecolor='none')
   plt.ylabel("$\kappa_{\mathrm{Hilbert}}$")
   plt.xlabel("N$_{60}$")
   plt.savefig("Fig6b.png")
   plt.show()


   plt.clf()
   plt.subplot(111)
   z=numpy.linspace(-0.05,0.2,30)
   plt.scatter(N_90,kappa_hilbert,s=1,c='k',edgecolor='none')
   plt.ylabel("$\kappa_{\mathrm{Hilbert}}$")
   plt.xlabel("N$_{90}$")
   plt.savefig("Fig6c.png")
   plt.show()

   plt.clf()
   plt.subplot(111)
   z=numpy.linspace(-0.05,0.2,30)
   plt.scatter(N_30,kappa_hilbert,s=1,c='k',edgecolor='none')
   plt.ylabel("$\kappa_{\mathrm{Hilbert}}$")
   plt.xlabel("N$_{30}$")
   plt.savefig("Fig6d.png")
   plt.show()

   plt.clf()
   plt.subplot(111)
   z=numpy.linspace(-0.05,0.2,30)
   plt.scatter(N_15,kappa_hilbert,s=1,c='k',edgecolor='none')
   plt.ylabel("$\kappa_{\mathrm{Hilbert}}$")
   plt.xlabel("N$_{15}$")
   plt.savefig("Fig6e.png")
   plt.show()
   """
# ======================================================================



if __name__ == '__main__':
  MScompare(sys.argv[1:])
  #print "check your fiducial cosmology?"


# ======================================================================
