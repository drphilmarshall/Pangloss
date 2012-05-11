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

#t0=time.clock()    

D = distances.Distance()

arcmin2rad = (1.0/60.0)*numpy.pi/180.0
rad2arcmin = 1.0/arcmin2rad

vb = False

# ============================================================================
    

class grid(object):
   def __init__(self,zmax=3,nplanes=200,cosmo=[0.25,0.75,0.73]): 
       D=distances.Distance(cosmo=cosmo)
       self.name = '1D Redshift grid, with precalculated quantities '

       self.zmax=zmax
       self.nplanes=nplanes    
       #These are the planes:
       self.zplane,self.dz=self.redshiftbins=numpy.linspace(0.,self.zmax,self.nplanes,endpoint=True,retstep=True)

       self.Da_p=numpy.zeros(len(self.zplane))
       for i in range(len(self.zplane)):
           self.Da_p[i]= D.Da(self.zplane[i])+1e-100
       self.rho_crit_p=D.rho_crit_univ(self.zplane)

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

   def __init__(self,zl,zs,zmax=[],nplanes=200,cosmo=[0.25,0.75,0.73]):
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
       self.Da_l=D.Da(self.zl)
       self.Da_s=D.Da(self.zs)
       self.Da_ls=D.Da(self.zl,self.zs)
       #Calculate angular diameter distances for each plane.
       self.Da_ps=numpy.zeros(len(self.zplane))
       self.Da_pl=numpy.zeros(len(self.zplane))
       for i in range(len(self.zplane)):
           self.Da_ps[i] = D.Da(self.zplane[i],self.zs)+1e-100
           self.Da_pl[i] = D.Da(self.zplane[i],self.zl)+1e-100
      
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
       



# ============================================================================

if __name__ == '__main__':
    zl,zs=0.6,1.4
    lg=lensgrid(zl,zs)
    lg.populatelensgrid()
    filename='lensgrid_zl%.2f_zs%.2f.pcl'%(zl,zs)
    #filename='test1.test'

    f=open(filename,'wb')
    cPickle.dump(lg,f,2)
    f.close()
# ============================================================================
