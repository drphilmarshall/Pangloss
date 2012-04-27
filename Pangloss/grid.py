'''
This file is part of the Pangloss project.
Copyright 2012 Tom Collett (IoA) and Phil Marshall (Oxford).

description:
------------
Code to make a lensgrid object, which lightcone objects and smooth comoponent objects can be 'snapped' to.

to-do:
------
Make it work

issues:
-------



'''

# ======================================================================

import pylab
import matplotlib.pyplot as plt
import numpy
import numpy.random as rnd
import distances
from mpl_toolkits.axes_grid1 import ImageGrid
from time import clock
import LensingProfiles as LP
import LensingFunc as LF

#import time
#t0=time.clock()    

D = distances.Distance()
D.h = 0.7


arcmin2rad = (1.0/60.0)*numpy.pi/180.0
rad2arcmin = 1.0/arcmin2rad

vb = False

# ============================================================================
    

class grid:

   def __init__(self,zmax=3,nplanes=300): 
       self.name = '1D Redshift grid, with precalculated quantities '

       self.zmax=zmax
       self.nplanes=nplanes
       
       #These are the planes:
       self.zplane,self.dz=self.redshiftbins=numpy.linspace(0,self.zmax,self.nplanes,endpoint=True,retstep=True)

       self.Da_p=numpy.zeros(len(self.zplane))
       for i in range(len(self.zplane)):
           self.Da_p[i]= D.Da(self.zplane[i])
       self.rho_crit_p=D.rho_crit_univ(self.zplane)

       return None
# ----------------------------------------------------------------------------

   def __str__(self):
        return '1-D Grid of %i planes seperated in redshift by dz= %f' % (self.nplanes,self.dz)

# ============================================================================
# ============================================================================

class lensgrid(grid):

   def __init__(self,zl,zs): 
      self.name = '1D Redshift grid, with precalculated quantities for a fixed lens and source'
      self.zl=zl
      self.zs=zs

      return None
# ----------------------------------------------------------------------------

   def __str__(self):
        return '1-D Grid of %i planes seperated in redshift by dz= %f, for a lens with zl=%.2f and zs=%.2f' % (self.nplanes,self.dz, self.zl,self.zs)

# ----------------------------------------------------------------------------

   def populatelensgrid(self):
       self.Da_l=D.Da(self.zl)
       self.Da_s=D.Da(self.zs)
       
       #Calculate angular diameter distances for each plane.
       #self.Da_d=numpy.zeros(len(self.zplane))
       self.Da_ps=numpy.zeros(len(self.zplane))
       self.Da_pl=numpy.zeros(len(self.zplane))
       for i in range(len(self.zplane)):
           #self.Da_d[i]= D.Da(self.zplane[i])
           self.Da_ps[i]= D.Da(self.zplane[i],self.zs)
           self.Da_pl[i]= D.Da(self.zplane[i],self.zl)


# ============================================================================

# TESTS:

def test(zl,zs):
    g=lensgrid(zl,zs)
    print g.Da_p
    return

#-------------------------------------------------------------------------

def test2(catalog): #plots a kappa distribution:
###not finished###
#we could consider feeding in our lens selection parameters here)


    return

#-------------------------------------------------------------------------

def test3(catalog):
    return

# ============================================================================

if __name__ == '__main__':
    test(0.6,1.4)

# ============================================================================
