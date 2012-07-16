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
   try:
      opts, args = getopt.getopt(argv,"hvn:r:",["help","verbose","zd","zs"])
   except getopt.GetoptError, err:
      print str(err) # will print something like "option -a not recognized"
      print MScompare.__doc__  # will print the big comment above.
      return

   vb = False
   Rcone = 5 # arcmin
   truncationscale=numpy.linspace(0.5,6,20) # *R_200 halo truncation
   truncationscale=numpy.append(truncationscale,[7,8,9,10,12,14,16,18,20])
   adoptedtruncationscale=10
   Ntruncs=len(truncationscale)
   Ncones=1000
   magcuts=numpy.linspace(14,26,50)


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
   if len(args) == 4:
      catalog = args[0]
      kappafile = args[1]
      gammafile1 = args[2]
      gammafile2 = args[3]

      if vb:
         print "Reconstructing convergence in lightcones from:",catalog
         print "Comparing to convergence in:",kappafile
         print "Number of lightcones to be reconstructed:",Ncones
   else:
      print MScompare.__doc__
      return

   # --------------------------------------------------------------------
   
   # Read in master catalog, and kappa map:
   if everything == True:
      for i in range(4):
         for j in range(4):
            catij="/data/tcollett/Pangloss/catalogs/GGL_los_8_%i_%i_%i_%i_N_4096_ang_4_SA_galaxies_on_plane_27_to_63.images.txt"%(oneseven,oneseven,i,j)
            if i == 0 and j == 0:
               master = atpy.Table(catij, type='ascii')
            else:
               additional=atpy.Table(catij, type='ascii')
               master.append(additional)
   else:
      master = atpy.Table(catalog, type='ascii')
   if vb: print "Read in master catalog, length",len(master)
   print "Read in master catalog, length",len(master)
   xmax = master['pos_0[rad]'].max()
   xmin = master['pos_0[rad]'].min()
   ymax = master['pos_1[rad]'].max()
   ymin = master['pos_1[rad]'].min()
   if vb: print "Catalog extent:",xmin,xmax,ymin,ymax
   


   kappafile = "/data/tcollett/Pangloss/kappafiles/GGL_los_8_%i_%i_N_4096_ang_4_rays_to_plane_37_f.fits"%(oneseven,oneseven)
   MSconvergence = Pangloss.kappamap(kappafile)
   aa=numpy.arange(0,4096,4)
   bb=numpy.arange(0,4096,4)
   c=0
   for a in aa:
      for b in bb:
         c+= MSconvergence.at(a,b,coordinate_system='image')
   print c/(len(aa)**2)

   if vb: print "Read in true kappa map, dimension",MSconvergence.NX
   MSshear1 = Pangloss.kappamap(gammafile1)
   MSshear2 = Pangloss.kappamap(gammafile2)
   if vb: print "Read in true gamma map, dimension",MSshear.NX

   # --------------------------------------------------------------------

   #build grid, populate it and invert the behroozi relation on each plane.
   grid=Pangloss.lensgrid(zl,zs,nplanes=200,cosmo=[0.25,0.75,0.73])
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

   Rcuts=numpy.linspace(0,Rcone,50)

   kappa_Scaled = numpy.zeros(Ncones)
   kappa_Scaled_truth = numpy.zeros((Ncones,Ntruncs))
   kappa_Scaled_truth_M = numpy.zeros((Ncones,len(magcuts)))
   kappa_Scaled_truth_R = numpy.zeros((Ncones,len(Rcuts)))
   kappa_difference_truth_M= numpy.zeros((Ncones,len(magcuts)))
   kappa_difference_truth_R = numpy.zeros((Ncones,len(Rcuts)))

   kappa_hilbert = numpy.zeros(Ncones)
   gamma_hilbert_2 = numpy.zeros(Ncones)
   gamma_hilbert_1 = numpy.zeros(Ncones)

   errors=True

   #!#!#!# Don't change from these defaults, without good reason!
   #{Mstar2Mh=False,Mh2Mh=True,perfectsatellites=False,centralsonly=False}
   hardcut="BMO"
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
      gamma_hilbert_1[k] =MSshear1.at(x[k],y[k],coordinate_system='physical')
      gamma_hilbert_2[k] =MSshear2.at(x[k],y[k],coordinate_system='physical')


      # Reconstruction Truth, with perfect knowledge of halo mass, but different truncation scales:
      tc= Pangloss.lens_lightcone(master,Rcone,xc,zl,zs,grid=grid)
      if investigateT:
         for T in range(Ntruncs):
            tc.make_kappa_contributions(hardcut=hardcut,truncationscale=truncationscale[T],scaling=scaling,errors=False,centralsonly=False,perfectsatellites=True,Mh2Mh=False,Mstar2Mh=False)
            kappa_Scaled_truth[k,T]=tc.kappa_Scaled_total
         #keep only the columns that we don't need to overwrite
            tc.galaxies.keep_columns(['GalID','HaloID','SubhaloID','Type','PlaneNumber','z_spec','pos_0[rad]','pos_1[rad]','Dc_los[Mpc/h]','M_Halo[M_sol/h]','M_Subhalo[M_sol/h]','M_Stellar[M_sol/h]','mag_SDSS_u','mag_SDSS_g','mag_SDSS_r','mag_SDSS_i','mag_SDSS_z','mag_J','mag_H','mag_K','x','y','r','zsnapped','psnapped','Da_d','Da_ds','Da_dl','beta','rho_crit','sigma_crit','rphys',"mag_F814W"])

      #now do magnitude cut catalogues
      if investigateM or investigateR:
         tc.make_kappa_contributions(hardcut=hardcut,truncationscale=adoptedtruncationscale,scaling=scaling,errors=False,centralsonly=False,perfectsatellites=True,Mh2Mh=False,Mstar2Mh=False)
      if investigateM:
         for M in range(len(magcuts)):
            mc=tc.galaxies.where(tc.galaxies["mag_SDSS_i"]<magcuts[M])
            kappa_Scaled_truth_M[k,M]=numpy.sum(mc.kappa_Scaled)
            kappa_difference_truth_M[k,M]= kappa_Scaled_truth_M[k,M] - kappa_hilbert[k]


      #now do radius cut catlogues:
      if investigateR:
         for R in range(len(Rcuts)):
            rc=tc.galaxies.where(tc.galaxies.r<Rcuts[R])
            kappa_Scaled_truth_R[k,R]=numpy.sum(rc.kappa_Scaled)
            kappa_difference_truth_R[k,R]= kappa_Scaled_truth_R[k,R] - kappa_hilbert[k]



      #keep only the columns that we don't need to overwrite
      tc.galaxies.keep_columns(['GalID','HaloID','SubhaloID','Type','PlaneNumber','z_spec','pos_0[rad]','pos_1[rad]','Dc_los[Mpc/h]','M_Halo[M_sol/h]','M_Subhalo[M_sol/h]','M_Stellar[M_sol/h]','mag_SDSS_u','mag_SDSS_g','mag_SDSS_r','mag_SDSS_i','mag_SDSS_z','mag_J','mag_H','mag_K','x','y','r','zsnapped','psnapped','Da_d','Da_ds','Da_dl','beta','rho_crit','sigma_crit','rphys',"mag_F814W"])

   #------------------------------------------------------------------
   #zero correct, and do some statistics

   print numpy.mean(kappa_hilbert)

   #T quantities
   if investigateT:
      kappa_empty_truth= numpy.zeros(Ntruncs)
      kappa_scatter_truth = numpy.zeros(Ntruncs)
      for T in range(Ntruncs):
         kappa_empty_truth[T] = numpy.mean(kappa_Scaled_truth[:,T])
         kappa_scatter_truth[T] = numpy.std(kappa_Scaled_truth[:,T]-kappa_hilbert[:])
         kappa_Scaled_truth[:,T] -=  kappa_empty_truth[T]

   #M quantities
   if investigateM:
      kappa_empty_truth_M =numpy.zeros(len(magcuts))
      kappa_scatter_truth_M = numpy.zeros(len(magcuts))
      kappa_16_M=numpy.zeros((len(magcuts)))
      kappa_50_M=numpy.zeros((len(magcuts)))
      kappa_84_M=numpy.zeros((len(magcuts)))
      kappa_dif_mean_M=numpy.zeros((len(magcuts)))
      for M in range(len(magcuts)):
         kappa_empty_truth_M[M] = numpy.mean(kappa_Scaled_truth_M[:,M])
         kappa_scatter_truth_M[M] = numpy.std(kappa_Scaled_truth_M[:,M]-kappa_hilbert[:])
         kappa_Scaled_truth_M[:,M] -=  kappa_empty_truth_M[M]
         kappa_difference_truth_M[:,M] -= kappa_empty_truth_M[M]
         kappa_dif_mean_M[M]=numpy.mean(kappa_difference_truth_M[:,M])
         kappa_16_M[M]=stats.scoreatpercentile(kappa_difference_truth_M[:,M],16)
         kappa_50_M[M]=stats.scoreatpercentile(kappa_difference_truth_M[:,M],50)
         kappa_84_M[M]=stats.scoreatpercentile(kappa_difference_truth_M[:,M],84)
   #R quantities
   if investigateR:
      kappa_empty_truth_R =numpy.zeros(len(Rcuts))
      kappa_scatter_truth_R = numpy.zeros(len(Rcuts))
      kappa_16_R=numpy.zeros((len(Rcuts)))
      kappa_50_R=numpy.zeros((len(Rcuts)))
      kappa_84_R=numpy.zeros((len(Rcuts)))
      kappa_dif_mean_R=numpy.zeros((len(Rcuts)))
      for R in range(len(Rcuts)):
         kappa_empty_truth_R[R] = numpy.mean(kappa_Scaled_truth_R[:,R])
         kappa_scatter_truth_R[R] = numpy.std(kappa_Scaled_truth_R[:,R]-kappa_hilbert[:])
         kappa_Scaled_truth_R[:,R] -=  kappa_empty_truth_R[R]
         kappa_difference_truth_R[:,R] -= kappa_empty_truth_R[R]
         kappa_dif_mean_R[R]=numpy.mean(kappa_difference_truth_R[:,R])
         kappa_16_R[R]=stats.scoreatpercentile(kappa_difference_truth_R[:,R],16)
         kappa_50_R[R]=stats.scoreatpercentile(kappa_difference_truth_R[:,R],50)
         kappa_84_R[R]=stats.scoreatpercentile(kappa_difference_truth_R[:,R],84)

   #------------------------------------------------------------------
   #set plotting sizes
   figprops = dict(figsize=(4.0, 3.0), dpi=128)          # Figure properties
   fig = plt.figure(**figprops)

   # Need small space between subplots to avoid deletion due to overlap...
   adjustprops = dict(\
      left=0.1,\
         bottom=0.12,\
         right=0.9,\
         top=0.95,\
         wspace=0.08,\
         hspace=0.1)
   fig.subplots_adjust(**adjustprops)

   fig_width_pt = 246.*2  # Get this from LaTeX using \showthe\columnwidth
   inches_per_pt = 1.0/72.27               # Convert pt to inch
   golden_mean = ((5**0.5)-1.0)/2.0         # Aesthetic ratio
   fig_size=[10, 10*golden_mean]
   params = {'backend': 'ps',
          'axes.labelsize': 16,
          'text.fontsize': 14,
          'legend.fontsize': 14,
          'xtick.labelsize': 16,
          'ytick.labelsize': 16,
          'text.usetex': True,
          'figure.figsize': fig_size}
   plt.rcParams.update(params)
   
# Plot data
   plt.figure(1)
   plt.clf()
   alp=0.1
   bet=0.15
   plt.axes([alp,bet,0.95-alp,0.95-bet])
#pylab.legend()

   #------------------------------------------------------------------

   #sanity check - do hilbert and reconstruction look correlated?
   if investigateT:
      plt.scatter(kappa_hilbert,kappa_Scaled_truth[:,-1])
      plt.show()
   elif investigateM:
      plt.scatter(kappa_hilbert,kappa_Scaled_truth_M[:,-1])
      plt.show() 
   #------------------------------------------------------------------
   #truncation plot
   if investigateT:
      plt.plot(truncationscale,kappa_scatter_truth,c='k')
      plt.xlabel(r"Truncation scale ($R_{200}$)")
      plt.ylabel(r"$\sigma(\kappa_{\mathrm{reconstruction}} - \kappa_{\mathrm{raytrace}})$",size='x-large')
      plt.savefig("truncation_scatter.png")
      plt.savefig("doc/truncation_scatter.eps")
      plt.show()

   #------------------------------------------------------------------
   #magnitude plot
   if investigateM:
      plt.plot(magcuts,kappa_16_M,c='k', ls = 'dashed',label="16th percentile")
      plt.plot(magcuts,kappa_50_M,c='k',label="Median")
      plt.plot(magcuts,kappa_84_M,c='k', ls = 'dashed',label="84th percentile")
      #plt.plot(magcuts,kappa_dif_mean_M,c='k', ls = '-.',label="Mean")                             
      plt.axhline(y=0,xmin=0,xmax=1,ls="dotted",c='k')
      plt.xlabel(r"Limiting i magnitude of reconstruction")
      plt.ylabel(r"$\kappa_{\mathrm{reconstruction}} - \kappa_{\mathrm{raytrace}}$",size='x-large')
      plt.ylim([-0.04,0.03])
      plt.legend(loc=4,title ="%i LoS"%Ncones)
      plt.savefig("mag_scatter.png")
      plt.savefig("doc/mag_scatter.eps")
      plt.show()

   #------------------------------------------------------------------
   #radius plot
   if investigateR:
      plt.plot(Rcuts,kappa_16_R,c='k', ls = 'dashed',label="16th percentile")
      plt.plot(Rcuts,kappa_50_R,c='k',label="Median")
      plt.plot(Rcuts,kappa_84_R,c='k', ls = 'dashed',label="84th percentile")
      #plt.plot(Rcuts,kappa_dif_mean_R,c='k', ls = '-.',label="Mean")                             
      plt.axhline(y=0,xmin=0,xmax=1,ls="dotted",c='k')
      plt.xlabel("Radius out to which halos are reconstructed")
      plt.ylabel("$\kappa_{\mathrm{reconstruction}} - \kappa_{\mathrm{raytrace}}$",size='x-large')
      plt.legend(loc=4,title ="%i LoS"%Ncones)
      plt.savefig("radius_scatter.png")
      plt.savefig("doc/radius_scatter.eps")
      plt.show()

   #------------------------------------------------------------------


   #output results
   F=open("truncresults%i.dat"%oneseven,"wb")
   G=open("magresults%i.dat"%oneseven,"wb")
   H=open("radresults%i.dat"%oneseven,"wb")
   if investigateT:
      results =  x,y,kappa_hilbert,truncationscale,kappa_scatter_truth,kappa_Scaled_truth
      cPickle.dump(results,F,protocol=2)
      F.close()
   if investigateM:
      results =  x,y,kappa_hilbert,magcuts,kappa_16_M,kappa_50_M,kappa_84_M,kappa_Scaled_truth_M
      cPickle.dump(results,G,protocol=2)
      G.close()
   if investigateR:
      results =  x,y,kappa_hilbert,Rcuts,kappa_16_R,kappa_50_R,kappa_84_R,kappa_Scaled_truth_R
      cPickle.dump(results,H,protocol=2)
      H.close()
# ======================================================================

if __name__ == '__main__':
   everything = True
   investigateT=True
   investigateR=True
   investigateM=True
   oneseven=1 #calculate on patch 1_1 or 7_7?
   MScompare(sys.argv[1:])
   #print "check your fiducial cosmology?"

# ======================================================================
