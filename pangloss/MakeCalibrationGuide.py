import cPickle
import pylab as plt
import scipy.stats as stats
from scipy import interpolate,optimize
import numpy
import numpy.random as rnd
import glob
#------------------------------------------------------------------
def MakeReconstructionCalibration(CalibrationResultsDirectory,surveyname):
   results="%s/cone*_%s.result"%(CalibrationResultsDirectory,i,surveyname)
   resultlist=glob.glob(results)

   Zl=len(resultlist)
   
   Cal_list=numpy.empty((Zl,5))

   for i in range(len(resultlist)):
      if i % 1000==0: print "making calibration", i,"of", len(resultlist)
      F=open(resultlist[i],"rb")
      res=cPickle.load(F)
      F.close()
      
      kappa_Hilbert,gamma1_Hilbert,gamma2_Hilbert,kappa_halo,gamma1_halo,gamma2_halo=res

      mu_Hilbert=1/((1-kappa_Hilbert)**2 - gamma_Hilbert**2)
      mu_halo=1/((1-kappa_halo)**2 - gamma_halo**2)

      Cal_list[i,0]=kappa_Hilbert
      Cal_list[i,1]=mu_Hilbert
      Cal_list[i,2]=numpy.median(kappa_halo)
      Cal_list[i,3]=numpy.median(mu_halo)
      Cal_list[i,4]=1

   F=open("%s/%s.jointdist"%(CalibrationResultsDirectory,surveyname),'wb')
   cPickle.dump(Cal_list,F,2)
   F.close()
   
def MakeWeightCalibration(CalibrationResultsDirectory,surveyname):
   results="%s/cone*_%s.result"%(CalibrationResultsDirectory,i,surveyname)
   resultlist=glob.glob(results)

   Zl=len(resultlist)
   
   Cal_list=numpy.empty((Zl,5))

   for i in range(len(resultlist)):
      if i % 1000==0: print "making calibration", i,"of", len(resultlist)
      F=open(resultlist[i],"rb")
      res=cPickle.load(F)
      F.close()
      
   #need to work out what goes in here...
   """
      kappa_Hilbert,gamma1_Hilbert,gamma2_Hilbert,=res

      mu_Hilbert=1/((1-kappa_Hilbert)**2 - gamma_Hilbert**2)
      mu_halo=1/((1-kappa_halo)**2 - gamma_halo**2)

      Cal_list[i,0]=kappa_Hilbert
      Cal_list[i,1]=mu_Hilbert
      Cal_list[i,2]=numpy.median(kappa_halo)

   F=open("%s/%s.jointdist"%(CalibrationResultsDirectory,surveyname),'wb')
   cPickle.dump(Cal_list,F,2)
   F.close()
   """
   
#code to use the above joint distribution...

def CalibrateResult(CalibrationResultsDirectory,surveyname,Value,width=0.01,Magnification=False):
   F=open("%s/%s.jointdist"%(CalibrationResultsDirectory,surveyname),'rb')
   Klist=cPickle.load(F)
   F.close()

   if magnification==False:
      Klist[:,4]=numpy.exp(-((Klist[:,2]-Value)**2)/(2*(width)**2))
   else:
      Klist[:,4]=numpy.exp(-((Klist[:,3]-Value)**2)/(2*(width)**2))

   kappabins=numpy.linspace(-0.1,5.0,5101)
   Klistdig=numpy.digitize(Klist[:,1],kappabins)


   Result=numpy.empty((len(Klist[:,0]),3))
   Result[:,0]=Klist[:,0]
   Result[:,1]=Klist[:,1]
   Result[:,2]=Klist[:,4]
   
   return Result
