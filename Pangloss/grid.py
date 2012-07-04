'''
This file is part of the Pangloss project.
Copyright 2012 Tom Collett (IoA) and Phil Marshall (Oxford).

description:
------------
Code to make a lensgrid object, which lightcone objects and smooth comoponent objects can be 'snapped' to.

to-do:
------
It should work now. TC #noguarantees

issues:
-------



'''

# ==========================================================================

import pylab
import matplotlib.pyplot as plt
import numpy
import numpy.random as rnd
import distances
from mpl_toolkits.axes_grid1 import ImageGrid
from time import clock
import LensingProfiles as LP
import LensingFunc as LF
import cPickle
import Relations as Rel
from scipy import interpolate,optimize
import copy

#t0=time.clock()    

arcmin2rad = (1.0/60.0)*numpy.pi/180.0
rad2arcmin = 1.0/arcmin2rad

vb = False

# ============================================================================
    

class grid(distances.Distance):
   def __init__(self,zmax=3,nplanes=100,cosmo=[0.25,0.75,0.73]): 
       distances.Distance.__init__(self,cosmo=cosmo)
       self.name = '1D Redshift grid, with precalculated quantities '
       self.zmax=zmax
       self.nplanes=nplanes    
       #These are the planes:
       self.zplane,self.dz=self.redshiftbins=numpy.linspace(0.0,self.zmax,self.nplanes,endpoint=True,retstep=True)

       self.Da_p=numpy.zeros(len(self.zplane))
       for i in range(len(self.zplane)):
           self.Da_p[i]= self.Da(self.zplane[i])+1e-100
       self.rho_crit_p=self.rho_crit_univ(self.zplane)

       self.cosmo=cosmo

   def snap(self,z):
      snapped_p=numpy.digitize(z,self.zplane-self.dz)-1
      snapped_z=self.zplane[snapped_p]
      return snapped_z,snapped_p

# ----------------------------------------------------------------------------

   def __str__(self):
        return '1-D Grid of %i planes seperated in redshift by dz= %f' % (self.nplanes,self.dz)

# ============================================================================
# ============================================================================

class lensgrid(grid):

   def __init__(self,zl,zs,zmax=[],nplanes=100,cosmo=[0.25,0.75,0.73]):
      if zmax==[]:zmax=zs
      grid.__init__(self,zmax,nplanes,cosmo=cosmo)
      self.name = '1D Redshift grid, with precalculated quantities for a fixed lens and source'


      # SNAP LENS AND SOURCE TO GRID:
      self.zl=self.snap([zl])[0][0]
      self.zs=self.snap([zs])[0][0]
      return None
# ----------------------------------------------------------------------------

   def __str__(self):
        return '1-D Grid of %i planes seperated in redshift by dz= %f, for a lens with zl=%.2f and zs=%.2f' % (self.nplanes,self.dz, self.zl,self.zs)

# ----------------------------------------------------------------------------

   def populatelensgrid(self):
       self.Da_l=self.Da(self.zl)
       self.Da_s=self.Da(self.zs)
       self.Da_ls=self.Da(self.zl,self.zs)
       #Calculate angular diameter distances for each plane.
       self.Da_ps=numpy.zeros(len(self.zplane))
       self.Da_pl=numpy.zeros(len(self.zplane))
       for i in range(len(self.zplane)):
           self.Da_ps[i] = self.Da(self.zplane[i],self.zs)+1e-100
           self.Da_pl[i] = self.Da(self.zplane[i],self.zl)+1e-100
      
       #Calculate sigma crit on each plane. #will raise a div0 warning, but this should be ignored
       self.sigma_crit_p=LF.SigmaCrit_Da(self.Da_p,self.Da_s,self.Da_ps)

       #Set angular diameter distances for beta calculation:
       D1s=numpy.zeros(len(self.zplane))
       D2=numpy.zeros(len(self.zplane))
       for i in range(len(self.zplane)):
          if self.zplane[i]>self.zl: #1 is lens, 2 is perturber
             D1s[i] = self.Da_ls
             D2[i]  = self.Da_p[i]
          else: #1 is perturber, 2 is lens
             D1s[i] = self.Da_ps[i]
             D2[i]  = self.Da_l
       D12= self.Da_pl
       Ds=self.Da_s*numpy.ones(len(self.zplane))

       self.beta_p=LF.beta_Da(D12,Ds,D2,D1s)

# ----------------------------------------------------------------------------


   def Behroozigrid(self,MFs):
       #Calculate the Behroozi relation and scatter on each plane.
       
       Mh = numpy.linspace(9.,25.,1601)
       Ms = numpy.linspace(5.5,12.5,801)
       sigma = 0.15
       self.Mh=Mh
       self.Ms=Ms
       


       #MFs contains the binned up stellar and halo mass functions from the catalogue (currently only one square)
       MhMF,MsMF,redshiftbinsMF,MhhistMF,MshistMF=MFs
       Mhdigitize=numpy.digitize(Mh,MhMF)-1
       Msdigitize=numpy.digitize(Ms,MsMF)-1


       ZPD=-1
       MhaloFromMstar_given_z = {}
       MstarFromMhalo_given_z = {}
       for p in range(len(self.zplane)):
          # Behroozi confusingly gives M_halo(M_*) but scatter for M_*(M_halo)
          TentotheMhMean = Rel.Mstar_to_M200(10**Ms,numpy.ones(len(Ms))*self.zplane[p])
          if vb: print MhMean
          MhMean=numpy.log10(TentotheMhMean)


          # Invert the relationship for a reasonable scatter
          invModel = interpolate.splrep(MhMean,Ms,s=0)

          if vb: print invModel

          # Calculate the mean M_* at fixed M_halo
          MsMean = interpolate.splev(Mh,invModel)
          if vb: print MsMean

          # Evaluate the distribution on the grid determined by Mh,Ms
          norm = sigma*(2*numpy.pi)**0.5
          pdf = numpy.empty((Ms.size,Mh.size))



          #Evaluate Halo + Stellar Mass functions
          #calculate correct redshift slice:
          zpd=numpy.digitize([self.zplane[p]],redshiftbinsMF)[0]
          #print zpd

          """
          SMfunc=MshistMF[zpd,:]
          SMF=SMfunc[Msdigitize]
          plt.plot(Ms,numpy.log10(SMF))
          plt.show(
          SMF/=numpy.sum(SMF)
          """
          #Fit a powerlaw to the binned HMF, if this hasn't already been done.
          TCHM=MhhistMF[zpd]
          TCM = numpy.array([(MhMF[i]+MhMF[i+1])/2. for i in range(MhMF.size-1)])
          aa=numpy.where(TCHM==TCHM.max())[0]
          bb=3
          C = TCHM[aa:-bb]>0
          def getPL(p,getM=False):
             N = 10**(p[0]+TCM*p[1])
             if getM:
                return N
             return (N-TCHM)[aa:-bb][C]/TCHM[aa:-bb][C]**0.5
          if ZPD != zpd:
             coeff,ier = optimize.leastsq(getPL,[14.56,-1.])


          HMF = 10**(coeff[0]+Mh*coeff[1])

          HMfunc=MhhistMF[zpd,:]  
          
          for i in range(Mh.size):

             pdf[:,i] = numpy.exp(-0.5*(Ms-MsMean[i])**2/sigma**2)/norm

          #multiply pdf by Halo Mass Function


          plotting=False
          if plotting:
             if zpd==3 and ZPD != 3:
                pdfm=pdf[50:751,100:601]
                plt.imshow(pdfm,origin="lower",extent = [Mh[100:601].min(), Mh[100:601].max(), Ms[50:751].min(), Ms[50:751].max()])
                plt.xlabel("log$_{10}$(M$_{\mathrm{halo}}$/M$_{\odot}$)")
                plt.ylabel("log$_{10}$(M$_{\mathrm{stellar}}$/M$_{\odot}$)")
                plt.title("P(M$_{\mathrm{stellar}}$|M$_{\mathrm{halo}}$) - Behroozi Relation")
                plt.savefig("Behroozi.png")
                plt.show()

          for i in range(Mh.size):
              pdf[:,i]*=HMF[i]
          if zpd==3 and ZPD != 3:
          #   #plt.subplot(212)
          #   plt.imshow(pdf,origin="lower",extent = [Mh.min(), Mh.max(), Ms.min(), Ms.max()])
          #   plt.xlabel("log$_{10}$(M$_{\mathrm{halo}}$/M$_{\odot}$)")
          #   plt.ylabel("log$_{10}$(M$_{\mathrm{stellar}}$/M$_{\odot}$)")
          #   plt.show()

             if plotting:
                pdfl=pdf*1.0
                for i in range(len(numpy.sum(pdf,1))):
                   pdfl[i,:]/=numpy.sum(pdf,1)[i]
                   
                pdfl=pdfl[50:751,100:601]
                plt.title("P(M$_{\mathrm{halo}}$|M$_{\mathrm{stellar}}$)")
                plt.imshow(pdfl,origin="lower",extent = [Mh[100:601].min(), Mh[100:601].max(), Ms[50:751].min(), Ms[50:751].max()])
                plt.xlabel("log$_{10}$(M$_{\mathrm{halo}}$/M$_{\odot}$)")
                plt.ylabel("log$_{10}$(M$_{\mathrm{stellar}}$/M$_{\odot}$)")
                plt.savefig("BehrooziInverse.png")
                plt.show()



          cdf = numpy.cumsum(pdf,1).astype(numpy.float32)
          cdf/= numpy.max(cdf)
          cdf = (cdf.T-cdf[:,0]).T
          cdf = (cdf.T/cdf[:,-1]).T

          if vb: print cdf
          MhaloFromMstar = {}

          
          for i in range(Ms.size):
             tmp = numpy.round(cdf[i]*1e5).astype(numpy.int64)/1e5
             lo = tmp[tmp==0].size-1
             hi = tmp[tmp<1].size+1
             tmp = int(numpy.round(Ms[i]*100))/100.
             MhaloFromMstar[tmp] = interpolate.splrep(cdf[i][lo:hi],Mh[lo:hi],s=0,k=1)



          ZPD=zpd*1.0
          MstarFromMhalo_given_z[self.zplane[p]] = invModel
          MhaloFromMstar_given_z[self.zplane[p]] = MhaloFromMstar
       self.MhaloDist = MhaloFromMstar_given_z
       self.Mstar_given_halo= MstarFromMhalo_given_z
       # Use like: MhaloFromMstar_given_z[z][Ms]



   def drawMhalo(self,Mstar,z,r=None):
      if len(Mstar)!=len(z): print "Draw MHalo takes a list of Mstars and redshifts. The redshifts must be snapped."
      Mst=numpy.log10(Mstar)
      Mh=numpy.zeros(len(Mstar))
      z=self.snap(z)[0]
      for i in range(len(Mstar)):
         if r==None: r = numpy.random.random()
         Mstr = numpy.round(Mst[i]*100.0)/100.
         #print Mstr
         zed=z[i]
         
         #print self.MhaloDist.keys()


         Mh[i]=10**(interpolate.splev(r,self.MhaloDist[zed][int(numpy.round(Mst[i]*100))/100.]))

         #hack since the cdf interpolator has a small bug
         k=0
         if numpy.isnan(Mh[i])==True:
            for j in range(20):
               if numpy.isnan(Mh[i])==True:
                  if j % 2 == 0:
                     Mh[i]=10**(interpolate.splev(r,self.MhaloDist[zed][int(numpy.round(Mst[i]*100+(j+2)/2))/100.]))
                  if j % 2 == 1:
                     Mh[i]=10**(interpolate.splev(r,self.MhaloDist[zed][int(numpy.round(Mst[i]*100-(j+1)/2))/100.]))
                  if k > j: k=j*1
            if k > 0: print "I had to cheat with the stellar mass by",(k+1.)/1.,"cex"

         if numpy.isnan(Mh[i])==True:   
            print numpy.log10(Mstar[i]),z[i],Mh[i]

            #Mast=10**numpy.linspace(9,12,20000)
            #Mah=self.drawMhalo(Mast,(numpy.ones(len(Mast))*zed))
            #plt.scatter(numpy.log10(Mah),numpy.log10(Mast),s=1,edgecolor='none')
            #plt.show()

      return Mh

   def drawMstar(self,Mhalo,z,r=None,scatter=False):
      if len(Mhalo)!=len(z): print "DrawMstar takes a list of Mhalos and redshifts. The redshifts must be snapped."
      Mh=numpy.log10(Mhalo)
      Mst=numpy.zeros(len(Mhalo))
      z=self.snap(z)[0]
      for i in range(len(Mhalo)):
         zed=z[i]
         if scatter==False:
            Mst[i]=10**(interpolate.splev(Mh[i],self.Mstar_given_halo[zed]))
         else:
            Mst[i]=10**(interpolate.splev(Mh[i],self.Mstar_given_halo[zed])+rnd.normal(0,0.15))
            

      #print numpy.log10(Mst)
      #print numpy.log10(Mhalo)
      #plt.hist(numpy.log10(Mst))
      #plt.show()
      #plt.scatter(Mst)

      return Mst
# ============================================================================

if __name__ == '__main__':
    zl,zs=0.6,1.4    
    lg=lensgrid(zl,zs,nplanes=20)
    lg.populatelensgrid()
    #filename='lensgrid_zl%.2f_zs%.2f.pcl'%(zl,zs)
    #filename='test1.test'

    #f=open(filename,'wb')
    #cPickle.dump(lg,f,2)
    #f.close()


    FILE=open("/home/tcollett/Pangloss/Pangloss/MHMS.data")
    MFs=cPickle.load(FILE)
    lg.Behroozigrid(MFs)
    
    MhTrue=10**(numpy.linspace(10,13,2000))
    print MhTrue
    """
    for i in range(len(lg.zplane)):
       Ms=lg.drawMstar(MhTrue,(numpy.ones(len(MhTrue))*lg.zplane[i]))
       #plt.scatter(numpy.log10(MhTrue),numpy.log10(Ms))
       #plt.show()


       Mh=lg.drawMhalo(Ms,(numpy.ones(len(MhTrue))*lg.zplane[i]))


       print len(numpy.isnan(Mh)[numpy.isnan(Mh)==True])

       plt.scatter(numpy.log10(Mh),numpy.log10(Ms),s=1,edgecolor='none')
       plt.show()
       #plt.scatter(numpy.log10(Mh),numpy.log10(MhTrue),s=1,edgecolor='none')
       #plt.show()
    
    #print lg.drawMhalo([  5.75568442e+10,   1.06243864e+11], [ 1.01593,   0.963173])
    """

    i=4
    Ms=lg.drawMstar(MhTrue,(numpy.ones(len(MhTrue))*lg.zplane[i]))
    plt.plot(numpy.log10(MhTrue),numpy.log10(Ms),c='k')
    plt.plot(numpy.log10(MhTrue),numpy.log10(Ms)+0.15,c='k',ls="dashed")
    plt.plot(numpy.log10(MhTrue),numpy.log10(Ms)-0.15,c='k',ls="dashed")
    plt.plot(numpy.log10(MhTrue),numpy.log10(Ms)+0.30,c='k',ls="dotted")
    plt.plot(numpy.log10(MhTrue),numpy.log10(Ms)-0.30,c='k',ls="dotted")
    plt.xlabel("log$_{10}$(M$_{\mathrm{halo}}$/M$_{\odot}$)")
    plt.ylabel("log$_{10}$(M$_{\mathrm{stellar}}$/M$_{\odot}$)")
    plt.title("P(M$_{\mathrm{stellar}}$|M$_{\mathrm{halo}}$) - Behroozi Relation")
    plt.savefig("BehrooziConf.png")
    plt.show()

    MsTrue=10**(numpy.linspace(6,12,200))
    Mh=lg.drawMhalo(MsTrue,(numpy.ones(len(MsTrue))*lg.zplane[i]),r=0.5)
    Mh1=lg.drawMhalo(MsTrue,(numpy.ones(len(MsTrue))*lg.zplane[i]),r=0.16)
    Mh2=lg.drawMhalo(MsTrue,(numpy.ones(len(MsTrue))*lg.zplane[i]),r=0.84)
    Mh3=lg.drawMhalo(MsTrue,(numpy.ones(len(MsTrue))*lg.zplane[i]),r=0.025)
    Mh4=lg.drawMhalo(MsTrue,(numpy.ones(len(MsTrue))*lg.zplane[i]),r=0.975)

    plt.plot(numpy.log10(MsTrue),numpy.log10(Mh),c='k')
    plt.plot(numpy.log10(MsTrue),numpy.log10(Mh1),c='k',ls="dashed")
    plt.plot(numpy.log10(MsTrue),numpy.log10(Mh2),c='k',ls="dashed")
    plt.plot(numpy.log10(MsTrue),numpy.log10(Mh3),c='k',ls="dotted")
    plt.plot(numpy.log10(MsTrue),numpy.log10(Mh4),c='k',ls="dotted")
    plt.ylabel("log$_{10}$(M$_{\mathrm{halo}}$/M$_{\odot}$)")
    plt.xlabel("log$_{10}$(M$_{\mathrm{stellar}}$/M$_{\odot}$)")
    plt.title("P(M$_{\mathrm{halo}}$|M$_{\mathrm{stellar}}$)")
    plt.savefig("InvBehrooziConf.png")
    plt.show()





# ============================================================================
