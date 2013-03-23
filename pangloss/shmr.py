# ===========================================================================

import numpy

# ============================================================================

class SHMR(object):
    """
    TEST DOCSTRING
    """

# ----------------------------------------------------------------------------

    def __init__(self,method='Behroozi',directory):
        
        self.name = self.__str__()
        self.method = method
        self.nMs,self.nMh,self.nz = 
        
        
        
        
        
        return None

# ----------------------------------------------------------------------------

    def __str__(self):
        return 'Stellar Mass to Halo Mass relation'

# ----------------------------------------------------------------------------
def Mstar_to_M200(M_Star,redshift,relationship='Behroozi'):
    #Takes an array of stellar mass and an array of redshifts and gives the best fit halo mass of {behroozi}.
    M_Star=10**(M_Star)
    if relationship == 'Behroozi':
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
 
          M_200[i] =(numpy.log10(M_1)+beta*numpy.log10(M_Star[i]/Mstar0)+((M_Star[i]/Mstar0)**delta)/(1.+(M_Star[i]/Mstar0)**-gamma)-0.5)
       return M_200 

# ----------------------------------------------------------------------------

    def drawMstar(self):
        return Mstar

# ----------------------------------------------------------------------------

    def drawMhalo(self):
        return Mhalo

# ----------------------------------------------------------------------------

    def makeCDFs(self):
        # Make model and model2 and name them informatively:
        self.CDF = numpy.zeros(self.nMs,self.nMh,self.nz)
        return

#=============================================================================

if __name__ == '__main__':
    shmr = SHMR('Behroozi')
    print shmr

#=============================================================================
# PofMgivenMcommaz.py:
# 
# import numpy,cPickle
# from scipy import interpolate,optimize
# import ndinterp
# import pylab as plt
# 
# 

# 
# # Data from lightcones
# inhalomass,inhaloZ = numpy.load('MassRedshift.cat')
# #inhaloZ[inhaloZ<0]=0
# # Define the grid over which we'll work
# Mh = numpy.linspace(10.,20.,1001)
# Ms = numpy.linspace(8.,13.,501)
# zeds,dz  = numpy.linspace(-0.2,1.6,10,retstep=True)
# dz/=2.
# print dz
# 
# model = numpy.empty((Ms.size,Mh.size,zeds.size))
# model2 = numpy.empty((Mh.size,zeds.size))
# 
# 
# for k in range(len(zeds)): 
#     z=zeds[k]
# # Behroozi confusingly gives M_halo(M_*) but scatter for M_*(M_halo)
#     MhMean = Mstar_to_M200(Ms,numpy.ones(len(Ms))*z)
# 
# # Invert the relationship for a reasonable scatter
#     invModel = interpolate.splrep(MhMean,Ms,s=0)
# 
# # Calculate the mean M_* at fixed M_halo
#     MsMean = interpolate.splev(Mh,invModel)
#     model2[:,k]=MsMean
# 
# # Evaluate the distribution on the grid determined by Mh,Ms
#     sigma=0.15
#     norm = sigma*(2*numpy.pi)**0.5
#     pdf = numpy.empty((Ms.size,Mh.size))
#     for i in range(Mh.size):
#         pdf[:,i] = numpy.exp(-0.5*(Ms-MsMean[i])**2/sigma**2)/norm
# 
#     
# # Deal with the halo mass function; here we fit a powerlaw
#     #bin by redshift:
#     Mhalos=inhalomass[inhaloZ<z+0.1]
#     Z=inhaloZ[inhaloZ<(z+dz)]
#     #This is all the halos in the relevant redshift bin:
#     Mhalos=Mhalos[Z>(z-dz)]
# 
#     Massbins=numpy.linspace(10,20,101)
#     hist,bins=numpy.histogram(Mhalos,Massbins)
#     MOD = interpolate.splrep(Massbins[:-1],hist,s=0,k=1)
#     HMF = interpolate.splev(Mh,MOD)     # This is an emperical HMF
#     #plt.plot(Mh,HMF)
#     #plt.show()
# 
#     TCM = Mh[HMF.argmax()+1:]
#     TCHM = HMF[HMF.argmax()+1:]
#     ""
#     def getPL(p,getM=False):
#         N = 10**(p[0]+TCM*p[1])
#         if getM:
#             return N
#         return (N-TCHM)/TCHM**0.5
#     coeff,ier = optimize.leastsq(getPL,[14.56,-1.])
#     HMF1 = 10**(coeff[0]+Mh*coeff[1])   # This is the powerlaw fit
#     ""
# 
# # Perform P(M*|Mh)*P(Mh)
#     pdf *= HMF1
# 
# # Calculate the CDF for P(Mh|M*)
#     pdf /= pdf.sum()
#     cdf = numpy.cumsum(pdf,1).astype(numpy.float32)
#     cdf = (cdf.T-cdf[:,0]).T
#     cdf = (cdf.T/cdf[:,-1]).T
# 
#     CDF = numpy.empty((cdf.shape[0],Mh.size))
#     X = numpy.linspace(0.,1.,Mh.size)
#     for i in range(Ms.size):
#     # Some hacks for numerical stability...
#         tmp = numpy.round(cdf[i]*1e5).astype(numpy.int64)/1e5
#         lo = tmp[tmp==0].size-1
#         hi = tmp[tmp<1].size+1
#     # Re-evaulate the CDF on a regular grid X
#         mod = interpolate.splrep(cdf[i][lo:hi],Mh[lo:hi],s=0,k=1)
#         q = interpolate.splev(X,mod)
#         CDF[i] = interpolate.splev(X,mod)
#     model[:,:,k]=CDF
# 
# 
# 
# # Form Mh(M*,X)
# axes = {}
# axes[0] = interpolate.splrep(Ms,numpy.arange(Ms.size),k=1)
# axes[1] = interpolate.splrep(X,numpy.arange(X.size),k=1)
# axes[2] = interpolate.splrep(zeds,numpy.arange(zeds.size),k=1)
# 
# model = ndinterp.ndInterp(axes,model)
# 
# 
# axes2 = {}
# axes2[0] = interpolate.splrep(Mh,numpy.arange(Mh.size),k=1)
# axes2[1] = interpolate.splrep(zeds,numpy.arange(zeds.size),k=1)
# model2 = ndinterp.ndInterp(axes2,model2)
# 
# 
# MODI=open("/data/tcollett/Pangloss/inverse.behroozi","wb")
# MODB=open("/data/tcollett/Pangloss/truth.behroozi","wb")
# cPickle.dump(model2,MODI,2)
# cPickle.dump(model,MODB,2)
# 
# 
# def drawMHalo(model,Mslist,redshiftList):
#    R = numpy.random.random(Mslist.size)
#    return model.eval(numpy.array([Mslist,R,redshiftList]).T)
# 
# def drawMStar(model,Mhlist,redshiftList):
#    return model.eval(numpy.array([Mhlist,redshiftList]).T)
# 
# 
# # Tom gave lots of data!!
# #inhalomass = inhalomass[::20]
# #inhaloZ=inhaloZ[::20]
# inhalomass=inhalomass[inhaloZ<1.5]
# inhaloZ=inhaloZ[inhaloZ<1.5]
# 
# 
# instarmass = drawMStar(model2,inhalomass,inhaloZ)
# 
# 
# plt.scatter(inhalomass,instarmass,edgecolor='',c=inhaloZ)
# plt.colorbar()
# plt.show()
# 
# outhalomass= drawMHalo(model,instarmass,inhaloZ)
# plt.scatter(inhalomass,outhalomass,edgecolor='',c=inhaloZ)
# x=numpy.linspace(9,16,101)
# plt.plot(x,x,c='r')
# plt.show()
# 
# """
# # Does the distribution Pr(Mh|M*obs) look Gaussian? Meh....
# O = drawMhalo(11.+numpy.random.randn(instarmass.size)*0.45)
# pylab.hist(O[O>1])
# pylab.show()
# 
# masses = []
# for i in range(100):
#     outhalomass = drawMhalo(instarmass+numpy.random.randn(instarmass.size)*0.45)
#     masses.append(outhalomass)
# masses = numpy.array(masses)
# 
# # Sometimes the stellar masses scatter out of the pre-defined stellar mass
# #   grid; in practice we strictly require 8 < M* < 13.
# masses[masses==0] = numpy.nan
# from scipy import stats
# M = stats.stats.nanmean(masses,0)
# e = stats.stats.nanstd(masses,0)
# pylab.errorbar(inhalomass,M,e,fmt='ko')
# pylab.plot([0.,20.],[0.,20.],'b')
# pylab.xlim(11.,17.)
# pylab.ylim(11.,17.)
# pylab.show()
# 
# 
# #MODI=open("/data/tcollett/Pangloss/inverse.behroozi","wb")
# MODB=open("/data/tcollett/Pangloss/mattruth.behroozi","wb")
# #cPickle.dump(invModel,MODI,2)
# cPickle.dump(model,MODB,2)
# 
# """
