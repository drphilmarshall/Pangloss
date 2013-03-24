# ===========================================================================

import numpy
from scipy import interpolate,optimize
import ndinterp

# ============================================================================

class SHMR(object):
    """
    TEST DOCSTRING
    """

# ----------------------------------------------------------------------------

    def __init__(self,method='Behroozi'):
        
        self.name = self.__str__()
        self.method = method
        self.nMh,self.nMs,self.nz = 501,251,10
        
        # Define the grid over which we'll work
        self.Mh_axis = numpy.linspace(10.,20.,self.nMh)
        self.Ms_axis = numpy.linspace(8.,13.,self.nMs)
        self.zed_axis,self.dz  = numpy.linspace(0.,1.6,self.nz,retstep=True)
        
        return None

# ----------------------------------------------------------------------------

    def __str__(self):
        return 'Stellar Mass to Halo Mass relation'

# ----------------------------------------------------------------------------

    def drawMstars(self,Mh,z):
        assert len(Mh)==len(z)
        MstarBest=self.S2H_model.eval(([Mh,z]).T)
        Mstar=MstarBest+numpy.random.randn(Mh.size)*0.15#intrinsic Mstar scatter of the behroozi relation.
        return Mstar

# ----------------------------------------------------------------------------

    def drawMhalos(self,Ms,z,X=None):
        assert len(Ms) == len(z)
        if X != None: assert len(X) == len(Ms)
        else: X=numpy.random.random(Ms.size)
        return self.H2S_model.eval(numpy.array([Ms,R,z]).T)

      
# ----------------------------------------------------------------------------

    def makeHaloMassFunction(self,catalog):
        #Infer halo mass function from Millenium Mh,z catalogue ; we use a power-law for this.
        zeds,dz  = numpy.linspace(0,1.8,10,retstep=True)#coarse redshift bin. these thingschage slowly.
        self.HMF={}
        self.HMFzkeys,self.HMFdz=zeds-dz,dz
        
        infer_from_data=True
        if infer_from_data:
            #load in the catalog's list of halo masses and redshift. Note this is not included with the 
            #gitrelease of pangloss, but it's pretty easy to get hold of. (e.g. ask tcollett@ast.cam.ac.uk!)


            import cPickle
            F=open('/data/tcollett/Pangloss/HaloMassRedshift.catalog','rb')
            inhalomass,inhaloZ=cPickle.load(F)
            F.close()
            inhaloZ[inhaloZ<0]=0

            for i in range(len(zeds)):
                z=zeds[i]+dz/2.
                Masses=inhalomass[inhaloZ>z-dz/2]
                newinZ=inhaloZ[inhaloZ>z-dz/2]
                Mhalos=Masses[newinZ<z+dz/2]
                Massbins=numpy.linspace(10,20,101)  
                hist,bins=numpy.histogram(Mhalos,Massbins)
                MOD = interpolate.splrep(Massbins[:-1],hist,s=0,k=1)
                HMF = interpolate.splev(self.Mh_axis,MOD)
                       #This is an emperical HMF
                self.TCM = self.Mh_axis[HMF.argmax()+1:]
                self.TCHM = HMF[HMF.argmax()+1:]


                #fit a powerlaw to the HMF
                PLcoeff,ier = optimize.leastsq(self.getPL,[14.56,-1.])
                print z, PLcoeff
                                #bug here...

                exit() 

                self.HMF[i]=PLcoeff
        # We've already fit a powerlaw to millenium: it's parameters as a function of z are 
        # included here.  
#        elif catalog='Millennium':
#            for Z in zeds:
#                z=Z+dz/2.
#                if z>0 and z<:self.HMF[z]=
#                if z> and z<:self.HMF[z]=
#                if z> and z<:self.HMF[z]=
#                if z> and z<:self.HMF[z]=
#                if z> and z<:self.HMF[z]=
#                if z> and z<:self.HMF[z]=
#                if z> and z<:self.HMF[z]=
#                if z> and z<:self.HMF[z]=
#                if z> and z<:self.HMF[z]=
#                if z> and z<:self.HMF[z]=
                
        return
# ----------------------------------------------------------------------------

    def getPL(self,p,getM=False):
        N = 10**(p[0]+self.TCM*p[1])    
        print N
        if getM:
            return N
        return (N-self.TCHM)/self.TCHM**0.5
        
# ----------------------------------------------------------------------------

    def getHaloMassFunction(self,z,HMFcatalog='Millennium'):
        try: self.HMF[HMFcatalog]
        except AttributeError:
            self.makeHaloMassFunction(HMFcatalog)
            
        for i in range(len(self.HMFzkeys)):
            key=self.HMFzkeys[i]
            if key<=z+self.HMFdz/2. and key>=z-self.HMFdz/2.:
                zkey=i
            
        return 10**(self.HMF[zkey][0]+self.Mh_axis*self.HMF[zkey][1])

# ----------------------------------------------------------------------------

    def makeCDFs(self):
        # make the models of the SHMR.
        
        #create the empty models that we will populate:
        S2H_grid = numpy.empty((self.Ms_axis.size,self.Mh_axis.size,self.zed_axis.size))
        H2S_grid = numpy.empty((self.Mh_axis.size,self.zed_axis.size))
        
        
        #1) invert the analytic behroozi MS->Mh relation
        Mh,Ms,zeds,dz=self.Mh_axis, self.Ms_axis, self.zed_axis,self.dz 
        
        for k in range(self.nz):
            z=zeds[k]
        
            MhMean = self.Mstar_to_M200(Ms,numpy.ones(len(Ms))*z)
            
            #fit a spline to the inverse of the behroozi relation
            invModel_z = interpolate.splrep(MhMean,Ms,s=0)
        
            # Calculate the mean M_* at fixed M_halo
            MsMean = interpolate.splev(Mh,invModel_z)
            H2S_grid[:,k]=MsMean
        
        
            #Now we can make P(Ms|Mh)
            sigma=0.15
            norm = sigma*(2*numpy.pi)**0.5
            pdflist = numpy.empty((Ms.size,Mh.size))
            for i in range(Mh.size):
                pdf = numpy.exp(-0.5*(Ms-MsMean[i])**2/sigma**2)/norm
                pdflist[:,i] = pdf

            #now we can convert this into a joint distribution, P(Ms,Mh) by multiplying by the halo
            #massfunction at this redshift (Bayes...)
 
            # Perform P(M*|Mh)*P(Mh)
            pdf *= self.getHaloMassFunction(z,HMFcatalog='Millennium')
            
            
            # # Calculate the CDF for P(Mh|M*)
            pdf /= pdf.sum()
            cdf = numpy.cumsum(pdf,1).astype(numpy.float32)
            cdf = (cdf.T-cdf[:,0]).T
            cdf = (cdf.T/cdf[:,-1]).T

            CDF = numpy.empty((cdf.shape[0],Mh.size))
            X = numpy.linspace(0.,1.,Mh.size)
            for i in range(Ms.size):
            # Some hacks for numerical stability...
                tmp = numpy.round(cdf[i]*1e5).astype(numpy.int64)/1e5
                lo = tmp[tmp==0].size-1
                hi = tmp[tmp<1].size+1
            # Re-evaulate the CDF on a regular grid
            mod = interpolate.splrep(cdf[i][lo:hi],Mh[lo:hi],s=0,k=1)
            q = interpolate.splev(X,mod)
            CDF[i] = interpolate.splev(X,mod)
            S2H_grid[:,:,k]=CDF
            
        # Form Mh(M*,X)
        axes = {}
        axes[0] = interpolate.splrep(Ms,numpy.arange(Ms.size),k=1)
        axes[1] = interpolate.splrep(X,numpy.arange(X.size),k=1)
        axes[2] = interpolate.splrep(zeds,numpy.arange(zeds.size),k=1)

        self.S2H_model = ndinterp.ndInterp(axes,S2H_grid)
            
        
        #Make the zero-scatter halo to stellar mass relation.
        # but first we have to do some more stuff...
        axes2 = {}
        axes2[0] = interpolate.splrep(Mh,numpy.arange(Mh.size),k=1)
        axes2[1] = interpolate.splrep(zeds,numpy.arange(zeds.size),k=1)
        self.H2S_model = ndinterp.ndInterp(axes2,H2S_grid)        
        
        return
# ----------------------------------------------------------------------------
    def Mstar_to_M200(self,M_Star,redshift):
    #Takes an array of stellar mass and an array of redshifts and gives the best fit halo mass of {behroozi}.
        M_Star=10**(M_Star)
        if self.method == 'Behroozi':
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

#=============================================================================

if __name__ == '__main__':
    shmr = SHMR('Behroozi')
    shmr.makeCDFs()

    print shmr.drawMstars([12],[0.1])

#=============================================================================
# PofMgivenMcommaz.py:
# 
# import numpy,cPickle

# import pylab as plt
# 
# 

# 
# # Data from lightcones


# 
# 
# MODI=open("/data/tcollett/Pangloss/inverse.behroozi","wb")
# MODB=open("/data/tcollett/Pangloss/truth.behroozi","wb")
# cPickle.dump(model2,MODI,2)
# cPickle.dump(model,MODB,2)
# 
# 

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
