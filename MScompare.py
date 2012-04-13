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

N_45mean = 41.8

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
   Rcone = 2 # arcmin
  
   # Defaults are for B1608 (CHECK):
   zd = 0.62
   zs = 1.39
  
   for o,a in opts:
      if o in ("-v", "--verbose"):
         vb = True
      elif o in ("-n"):
         Ncones = int(a)
      elif o in ("-r"):
         Rcone = a
      elif o in ("--zd"):
         zd = a
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

   xmax = master['pos_0[rad]'].max()
   xmin = master['pos_0[rad]'].min()
   ymax = master['pos_1[rad]'].max()
   ymin = master['pos_1[rad]'].min()
   if vb: print "Catalog extent:",xmin,xmax,ymin,ymax
    
   MSconvergence = Pangloss.kappamap(kappafile)
   if vb: print "Read in true kappa map, dimension",MSconvergence.NX



   # --------------------------------------------------------------------

   # Generate Ncones random positions, and reconstruct kappa_keeton in
   # each one. Also look up true kappa_hilbert at that position:
   
   x = rnd.uniform(xmin,xmax,Ncones)
   y = rnd.uniform(ymin,ymax,Ncones)
      
   kappa_keeton = numpy.zeros(Ncones)
   kappa_hilbert = numpy.zeros(Ncones)
   N_45 = numpy.zeros(Ncones)
   N_60 = numpy.zeros(Ncones)
   N_90 = numpy.zeros(Ncones)
   N_30 = numpy.zeros(Ncones)
   N_15 = numpy.zeros(Ncones)
   other = numpy.zeros(Ncones)  
   for k in range(Ncones):
      if k % 100 == 0 and k !=0: print ("evaluating cone %i of %i" %(k,Ncones))
      xc = [x[k],y[k],zd]

      # Truth:
      kappa_hilbert[k] = MSconvergence.at(x[k],y[k],coordinate_system='physical')

      # Reconstruction:
      lc = Pangloss.lightcone(master,Rcone,zs,position=xc)
      #lc.N_radius()
      lc.make_kappa_contributions()
      ind=(numpy.argsort(lc.galaxies.kappa_keeton))[-1]
      #other[k] = lc.galaxies.r[ind]
      #N_45[k]=lc.N_radius(45,cut=[18.5,22.5])
      #N_60[k]=lc.N_radius(45,cut=[18.5,23.5])
      #N_90[k]=lc.N_radius(45,cut=[18.5,24.5])
      #N_30[k]=lc.N_radius(45,cut=[18.5,25.5])
      #N_15[k]=lc.N_radius(45,cut=[18.5,26.5])
 
      N_45[k]=lc.N_radius(45,cut=[18.5,24.5])
      N_60[k]=lc.N_radius(60,cut=[18.5,24.5])
      N_90[k]=lc.N_radius(90,cut=[18.5,24.5])
      N_30[k]=lc.N_radius(30,cut=[18.5,24.5])
      N_15[k]=lc.N_radius(15,cut=[18.5,24.5])
      
      

      kappa_keeton[k] = numpy.sum(lc.galaxies.kappa_keeton)

   print numpy.average(N_45)
   print numpy.average(N_60)
   print numpy.average(N_90)
   print numpy.average(N_30)
   print numpy.average(N_15)


   # --------------------------------------------------------------------

   # Basic statistics of two arrays:
   
   difference = kappa_keeton - kappa_hilbert
   bias = numpy.average(difference)
   scatter = numpy.std(difference)

   print "$\kappa_{\mathrm{Keeton}}-\kappa_{\mathrm{Hilbert}}$ = ",bias,"+/-",scatter

   h=kappa_hilbert
   print stats.linregress(N_45,h)[2]
   print stats.linregress(N_60,h)[2]
   print stats.linregress(N_90,h)[2]
   print stats.linregress(N_30,h)[2]
   print stats.linregress(N_15,h)[2]

#========================================================================
   #Now lots of plotting routines!!!!
#========================================================================  
   """
   list=numpy.linspace(-0.05,0.25,30)
   plt.subplot(311)
   plt.title("$\kappa_{\mathrm{Keeton}}-\kappa_{\mathrm{Hilbert}}$ = %.3f +/- %.3f" % (bias,scatter))
   plt.hist(kappa_keeton, bins=list,normed=True, label="$\kappa_{\mathrm{Keeton}}$")
   plt.legend(loc=1)
   plt.subplot(312)
   plt.hist(kappa_hilbert, bins=list,normed=True,label="$\kappa_{\mathrm{Hilbert}}$")
   plt.legend(loc=1)
   plt.subplot(313)
   plt.hist(difference, bins=20,normed=True,label="$\kappa_{\mathrm{Keeton}}-\kappa_{\mathrm{Hilbert}}$")
   plt.legend(loc=2)
   plt.savefig("Fig2.png")
   plt.show()
   plt.clf()
   plt.subplot(111)
   plt.scatter(kappa_keeton,difference,s=1,c='k',edgecolor='none')
   plt.ylabel("$\kappa_{\mathrm{Keeton}}-\kappa_{\mathrm{Hilbert}}$")
   plt.xlabel("$\kappa_{\mathrm{Keeton}}$")
   plt.xlim([-0.08,0.25])
   plt.ylim([-0.2,0.1])
   plt.savefig("Fig3.png")
   plt.show()
   plt.scatter(difference,other,s=1,c='k',edgecolor='none')
   plt.xlabel("$\kappa_{\mathrm{Keeton}}-\kappa_{\mathrm{Hilbert}}$")
   #plt.ylabel("LOS distance to most important object (arcmin)")
   plt.ylabel("Largest individual $\kappa_{\mathrm{Keeton}}$ in LOS")
   plt.savefig("other.png")
   plt.show()   
   plt.clf()
   plt.subplot(111)
   plt.scatter(kappa_hilbert,kappa_keeton,s=1,c='k',edgecolor='none')
   plt.ylabel("$\kappa_{\mathrm{Keeton}}$")
   plt.xlabel("$\kappa_{\mathrm{Hilbert}}$")
   plt.savefig("Fig1.png")
   plt.show()  

   y=difference
   x=kappa_keeton
   yerr=0.01 #assume an error on y
   fitfunc= lambda p,x: p[0]+p[1]*x  
   errfunc= lambda p,x,y, err: (y-fitfunc(p,x))/err
   pinit=[0.03,-1.]
   out =optimize.leastsq(errfunc,pinit,args=(x,y,yerr),full_output=1)
   pfinal=out[0]
   covar=out[1]
   #print covar
   #print pfinal

   scatter=y-pfinal[0]-pfinal[1]*x
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

# ======================================================================
