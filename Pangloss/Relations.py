import numpy
from numpy import exp,pi
from scipy.special import gamma
import numpy.random as rnd
import atpy,copy
import matplotlib,pylab as plt
from scipy import interpolate

#--------------------------------------------------------------
def logerr(l,m,s):
      c=10**rnd.normal(m,s)
      return c

def MCrelation(M200,MCerror=False):
      if MCerror==False:
         c_200 = 4.67*(M200/(10**14))**0.11 #Neto et al. equation 5
      if MCerror==True:
       c_200=4.67*(M200/(10**14))**0.11
       logc_200=numpy.log10(c_200)
       lM200 = numpy.log10(M200)
       c_new = c_200*10**(rnd.normal(1,0.2,len(M200)))/10**(1+(0.2**2)/2)
       #print numpy.max(c_new),numpy.min(c_new)

       return c_new


       for i in range(len(M200)):    #best fit scatter parameters of neto et al (double log normal)
         if lM200[i]<11.875: ### FILLER
            f=0.205 ###FILLER
            if rnd.rand()<f:
               c_200[i]=logerr(logc_200[i],0.683,0.147) ###FILLER
            else:
               c_200[i]=logerr(logc_200[i],0.920,0.106) ###FILLER
         if lM200[i]<12.125:
            f=0.205
            if rnd.rand()<f:
               c_200[i]=logerr(logc_200[i],0.683,0.147)
            else:
               c_200[i]=logerr(logc_200[i],0.920,0.106)
         elif lM200[i]<12.375:
            f=0.171
            if rnd.rand()<f:
               c_200[i]=logerr(logc_200[i],0.658,0.150)
            else:
               c_200[i]=logerr(logc_200[i],0.903,0.108)
         elif lM200[i]<12.625:
            f=0.199
            if rnd.rand()<f:
               c_200[i]=logerr(logc_200[i],0.646,0.139)
            else:
               c_200[i]=logerr(logc_200[i],0.881,0.099)        
         elif lM200[i]<12.875:
            f=0.229
            if rnd.rand()<f:
               c_200[i]=logerr(logc_200[i],0.605,0.158)
            else:
               c_200[i]=logerr(logc_200[i],0.838,0.101)          
         elif lM200[i]<13.125:
            f=0.263
            if rnd.rand()<f:
               c_200[i]=logerr(logc_200[i],0.603,0.136)
            else:
               c_200[i]=logerr(logc_200[i],0.810,0.100)          
         elif lM200[i]<13.375:
            f=0.253
            if rnd.rand()<f:
               c_200[i]=logerr(logc_200[i],0.586,0.140)
            else:
               c_200[i]=logerr(logc_200[i],0.793,0.099)         
         elif lM200[i]<13.625:
            f=0.275
            if rnd.rand()<f:
               c_200[i]=logerr(logc_200[i],0.566,0.142)
            else:
               c_200[i]=logerr(logc_200[i],0.763,0.095)            
         elif lM200[i]<13.875:
            f=0.318
            if rnd.rand()<f:
               c_200[i]=logerr(logc_200[i],0.543,0.140)
            else:
               c_200[i]=logerr(logc_200[i],0.744,0.094)
         elif lM200[i]<14.125:
            f=0.361
            if rnd.rand()<f:
               c_200[i]=logerr(logc_200[i],0.531,0.131)
            else:
               c_200[i]=logerr(logc_200[i],0.716,0.088)            
         elif lM200[i]<14.375:
            f=0.383
            if rnd.rand()<f:
               c_200[i]=logerr(logc_200[i],0.510,0.121)
            else:
               c_200[i]=logerr(logc_200[i],0.689,0.095)       
         elif lM200[i]<14.625:
            f=0.370
            if rnd.rand()<f:
               c_200[i]=logerr(logc_200[i],0.490,0.133)
            else:
               c_200[i]= logerr(logc_200[i],0.670,0.094)           
         elif lM200[i]<14.875:
            f=0.484
            if rnd.rand()<f:
               c_200[i]=logerr(logc_200[i],0.519,0.121)
            else:
               c_200[i]=logerr(logc_200[i],0.635,0.091)    
         elif lM200[i]<15.125:
            f=0.578
            if rnd.rand()<f:
               c_200[i]=logerr(logc_200[i],0.493,0.094)
            else:
               c_200[i]=logerr(logc_200[i],0.661,0.061)
      return c_200

#----------------------------------------------------------

def binMS(cat=None):
   #binning up the MS catalogue and using that relation:
      if cat == None:
            d1= "../../data/GGL_los_8_0_0_1_1_N_4096_ang_4_STARS_SA_galaxies_ANALYTIC_SA_galaxies_on_plane_27_to_63.images.txt"
            cat=[atpy.Table(d1, type='ascii')]
      for i in range(len(cat)):
            cat_j=copy.copy(cat[i])
            cat_j.keep_columns(['M_Subhalo[M_sol/h]','M_Stellar[M_sol/h]','z_spec'])
            if i == 0:
                  c=copy.copy(cat_j)
            else:
                  c.append(cat_j)
            
            Mhalo=c['M_Subhalo[M_sol/h]']
            Mstars=c['M_Stellar[M_sol/h]']

            LH=numpy.log10(Mhalo)
            LS=numpy.log10(Mstars)
            z=c.z_spec

            massbin,delta=numpy.linspace(7,12,20,retstep=True)
            zbin=[0,99]
          
            MB=numpy.digitize(LH,massbin)
            zB=numpy.digitize(z,zbin)

            plt.scatter(LH,LS-LH,s=0.2,c='k',edgecolor='')
            plt.xlabel('log(M_Halo/M$_\odot$)')
            plt.ylabel('log(M_Stellar/M_Halo)')
            plt.xlim([10,14])
            plt.ylim([-3,0.5])
            plt.savefig("starformationefficiency.png")
            plt.show()

            mean=numpy.zeros((len(massbin),len(zbin)))
            for i in range(len(massbin)):
              for j in range(len(zbin)):
                M=c.where(LH>massbin[i]-delta/2. &\
                          LH<massbin[i]+delta/2.)
                mean[i,j]=numpy.mean(numpy.log10(M['M_Subhalo[M_sol/h]']))
      return None

#--------------------------------------------------------------
def Mstar_to_M200(M_Star,redshift,scatter=True,Behroozi=True,cat=None):
   if Behroozi==True:
      #Following Behroozi et al. 2010.
      M_200=numpy.zeros(len(M_Star))
      #parameters:
      for i in range(len(M_Star)):
         z=redshift[i]
         if z<0.9:
            Mstar00 = 10.72
            Mstar0a = 0.55
            Mstar0aa=0.0
            M_10 = 12.35
            M_1a = 0.28
            beta0 = 0.44
            betaa = 0.18
            delta0 = 0.57
            deltaa = 0.17
            gamma0 = 1.56
            gammaa = 2.51
         else:
            Mstar00 = 11.09
            Mstar0a = 0.56
            Mstar0aa= 6.99
            M_10 = 12.27
            M_1a = -0.84
            beta0 = 0.65
            betaa = 0.31
            delta0 = 0.56
            deltaa = -0.12
            gamma0 = 1.12
            gammaa = -0.53


      #scaled parameters:
         a=1./(1.+z)
         M_1=10**(M_10+M_1a*(a-1))
         beta=beta0+betaa*(a-1)
         Mstar0=10**(Mstar00+Mstar0a*(a-1)+Mstar0aa*(a-0.5)**2)
         delta=delta0+deltaa*(a-1)
         gamma=gamma0+gammaa*(a-1)

      #reltationship ****NO SCATTER****

         M_200[i] =10.0**(numpy.log10(M_1)+beta*numpy.log10(M_Star[i]/Mstar0)+((M_Star[i]/Mstar0)**delta)/(1.+(M_Star[i]/Mstar0)**-gamma)-0.5)

         if scatter==True:
             M_200[i]=10.0**(numpy.log10(M_200[i])*rnd.normal(0,0.2))
      return M_200 

def Behroozi_Spline():
      points = numpy.linspace(6,12,500)
      redshift=0*points+0.1
      masses= 10**points

      MH=Mstar_to_M200(masses,redshift,scatter=False)
      BS1=interpolate.splrep(MH,masses)
      BS2=interpolate.splrep(masses,MH)

      redshift=0*points+1.5

      MH2=Mstar_to_M200(masses,redshift,scatter=False)

      BS3=interpolate.splrep(MH2,masses)
      BS4=interpolate.splrep(masses,MH2)
 
      return BS1,BS2,BS3,BS4


#----------------------------------------------------------
def Mhalo_to_Mhalo(OLD,z,HALOSTARlowz,STARHALOlowz,HALOSTARhighz,STARHALOhighz,eBer=0.212):
      Mstar=OLD*0.0
      NEW=OLD*0.0
      Mstar[numpy.where(z<0.9)]=(interpolate.splev(OLD[numpy.where(z<0.9)],HALOSTARlowz))
      Mstar[numpy.where(z>0.9)]=(interpolate.splev(OLD[numpy.where(z>0.9)],HALOSTARhighz))

      Mstar=10**(numpy.log10(Mstar)+rnd.normal(0,eBer,size=len(OLD)))

      NEW[numpy.where(z<0.9)]=(interpolate.splev(Mstar[numpy.where(z<0.9)],STARHALOlowz))
      NEW[numpy.where(z>0.9)]=(interpolate.splev(Mstar[numpy.where(z>0.9)],STARHALOhighz))

      return NEW

#Mstar_to_M200(1e10,0.2,scatter=True,Behroozi=False,cat=None)#--------------------------------------------------------------
def reffFromMass(mstar):
     # See http://adsabs.harvard.edu/abs/2009MNRAS.394.1978H
     logm = numpy.log10(mstar)
     return 10**(7.55-1.84*logm+0.11*logm**2) #units are kiloparsecs
