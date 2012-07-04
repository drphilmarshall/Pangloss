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
import time

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
     Reconstruct a pdf for the convergence along a LoS and compare with
     "true" convergence for a set of lightcones. Truth comes from 
     Stefan Hilbert's ray tracing work in the Millenium Simulation, in
     the form of a pixellated convergence map. Reconstructions are made
     from simulated catalogs in the same field.
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
   Rcone = 5 # arcmin
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


   # --------------------------------------------------------------------

   #build grid, populate it and invert the behroozi relation on each plane.
   grid=Pangloss.lensgrid(zl,zs,nplanes=5,cosmo=[0.25,0.75,0.73])
   grid.populatelensgrid()
   FILE=open("/home/tcollett/Pangloss/Pangloss/MHMS.data")
   MFs=cPickle.load(FILE)
   grid.Behroozigrid(MFs)


   # --------------------------------------------------------------------

   # Generate Ncones random positions
   
   x = rnd.uniform(xmin+Rcone*arcmin2rad,xmax-Rcone*arcmin2rad,Ncones)
   y = rnd.uniform(ymin+Rcone*arcmin2rad,ymax-Rcone*arcmin2rad,Ncones)
   

   #------------------------------------------------------------------
   # Make a distribution for 1000 random cones - that way we can estimate <kappa>, to make a void corrrection.

   kappa_Scaled = numpy.zeros(Ncones)
   kappa_Scaled_truth = numpy.zeros(Ncones)
   kappa_hilbert = numpy.zeros(Ncones)
   errors=True


   #!#!#!# Don't change from these defaults, without good reason!
   #{Mstar2Mh=False,Mh2Mh=True,perfectsatellites=False,centralsonly=False}
   Mstar2Mh=False
   Mh2Mh=True             
   perfectsatellites=False
   centralsonly=False
   #!#!#!#

   mean0=True

   for k in range(Ncones):
      if k % 200 == 0: print ("evaluating cone %i of %i, for the smooth correction" %(k,Ncones))
      xc = [x[k],y[k]]

      # Hilbert Truth:
      kappa_hilbert[k] = MSconvergence.at(x[k],y[k],coordinate_system='physical')
      # THE CATALOGUES NEED TRANSPOSING!!!! (don't make that mistake again!)


      # Reconstruction Truth, with perfect knowledge of halo mass: 
      tc= Pangloss.lens_lightcone(master,Rcone,xc,zl,zs,grid=grid)
      tc.make_kappa_contributions(hardcut="RVir",truncationscale=truncationscale,scaling=scaling,errors=False,centralsonly=False,perfectsatellites=True,Mh2Mh=False,Mstar2Mh=False)

      # Reconstruction, with Behroozi blurred halo masses:
      lc = Pangloss.lens_lightcone(master,Rcone,xc,zl,zs,grid=grid)
      lc.make_kappa_contributions(hardcut="RVir",truncationscale=truncationscale,scaling=scaling,errors=errors,centralsonly=centralsonly,perfectsatellites=perfectsatellites,Mh2Mh=Mh2Mh,Mstar2Mh=Mstar2Mh)

 
      kappa_Scaled[k] = lc.kappa_Scaled_total
      kappa_Scaled_truth[k]=tc.kappa_Scaled_total


   kappa_empty = numpy.mean(kappa_Scaled)
   kappa_empty_truth = numpy.mean(kappa_Scaled_truth)

   kappa_Scaled       -=  kappa_empty 
   kappa_Scaled_truth -=  kappa_empty_truth


   #sanity check - do hilbert and reconstruction look correlated?
   plt.scatter(kappa_hilbert,kappa_Scaled_truth)
   plt.show()

   blurfactor=numpy.std(kappa_hilbert-kappa_Scaled_truth)
   plt.scatter(kappa_hilbert-kappa_Scaled_truth,kappa_Scaled)
   plt.show()

   #------------------------------------------------------------------

   #Now we have the void corrections, for both kappa_truth, and kappa_reconstruction

   #so now lets draw the cones that we want to creat pdfs for.

   Npdf=1000 #how many pdf's do we want?
   #Ncones=1000 #How many realisations? 

   xn = rnd.uniform(xmin+Rcone*arcmin2rad,xmax-Rcone*arcmin2rad,Npdf)
   yn = rnd.uniform(ymin+Rcone*arcmin2rad,ymax-Rcone*arcmin2rad,Npdf)


   kappa_hilbert_n = numpy.zeros(Npdf)
   kappa_Truth_n = numpy.zeros(Npdf)
   kappa_Reconstruct_n = numpy.zeros((Npdf,Ncones))

   #True values on each cone:
   for n in range(Npdf):
      xc = [xn[n],yn[n]]

      # Hilbert Truth:
      kappa_hilbert_n[n] = MSconvergence.at(xn[n],yn[n],coordinate_system='physical')
      #Reconstruction Truth:
      tc= Pangloss.lens_lightcone(master,Rcone,xc,zl,zs,grid=grid)
      tc.make_kappa_contributions(hardcut="RVir",truncationscale=truncationscale,scaling=scaling,errors=False,centralsonly=False,perfectsatellites=True,Mh2Mh=False,Mstar2Mh=False)
      kappa_Truth_n[n]=tc.kappa_Scaled_total-kappa_empty_truth

      #pdf:
      for k in range(Ncones):
         if k % 200 == 0:  print ("realization %i of %i for pdf number %i of %i" %(k,Ncones,n,Npdf))
         lc = Pangloss.lens_lightcone(master,Rcone,xc,zl,zs,grid=grid)
         lc.make_kappa_contributions(hardcut="RVir",truncationscale=truncationscale,scaling=scaling,errors=errors,centralsonly=centralsonly,perfectsatellites=perfectsatellites,Mh2Mh=Mh2Mh,Mstar2Mh=Mstar2Mh)
         kappa_Reconstruct_n[n,k] = lc.kappa_Scaled_total-kappa_empty+rnd.normal(0,blurfactor)
      #last term comes from the intrinsic difference in kappa_hilbert vs kappa_reconstruct_truth

   #------------------------------------------------------------------
         """
      #Now lets plot the pdfs!
      colors=["b","g","r","y","purple"]
      if n < 5:
         #print numpy.std(kappa_Reconstruct_n[n,:])
         #print blurfactor
         #lc.plot()
         #plt.savefig("plot.png")
         ##plt.show()
         #plt.clf()
         #AX=plt.subplot(1,1,1, aspect ='equal')
         #lc.massplot(AX)
         #plt.savefig("massplot.png")
         #plt.show()
         #plt.clf()

         bins=numpy.linspace(-0.2,0.6,41)
         plt.hist(kappa_Reconstruct_n[n,:],bins=bins,color=colors[n],alpha=0.3)
         plt.axvline(x=kappa_hilbert_n[n],ymin=0,ymax=1,c=colors[n],ls='dashed',label="$\kappa$ Hilbert")
         plt.axvline(x=kappa_Truth_n[n],ymin=0,ymax=1,c=colors[n],ls='dotted',label="$\kappa$ Perfect Mass")
         plt.legend(loc=1)
         plt.savefig("singlepdf.png")
   plt.show()
   """


   #output results
   F=open("pdfs1000e.dat","wb")
   results = kappa_hilbert_n,kappa_Truth_n,kappa_Reconstruct_n
   cPickle.dump(results,F,protocol=2)


# ======================================================================



if __name__ == '__main__':
  MScompare(sys.argv[1:])
  #print "check your fiducial cosmology?"


# ======================================================================
