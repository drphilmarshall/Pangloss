import pylab
import matplotlib.pyplot as plt
import numpy, numpy.random as rnd, atpy
import distances
from mpl_toolkits.axes_grid1 import ImageGrid
from time import clock

#import time
#t0=time.clock()    

D = distances.Distance()
D.h = 0.7

arcmin2rad = (1.0/60.0)*numpy.pi/180.0
rad2arcmin = 1.0/arcmin2rad
# ------------------------------------------------------------------------
def F(x):
    return numpy.arccos(1/x))/((x**2-1)**.5
def L(x,t):
    return numpy.log(x/((t**2+x**2)**.5+t))
#=========================================================================
# NFW profile functions.
# ------------------------------------------------------------------------
    # Function needed to calculate kappa for an NFW halo. 
def Ffunc(x):
       z=numpy.zeros(len(x))
       for i in range(len(x)):
          if x[i]==-1:
              z[i]=0.0
          elif x[i]>1:
             z[i]= (1-(2./(x[i]**2-1)**.5)*numpy.arctan(((x[i]-1.)/(x[i]+1))**.5))/(x[i]**2-1.)
          else: 
             #print "WARNING You are very close to a halo"
             if x[i]==1:
                z[i] =1./3
             else:
                y=(-x[i]**2+1)**.5
                z[i] = (1.-(2./(1-x[i]**2)**.5)*numpy.arctanh(((1.-x[i])/(x[i]+1))**.5))/(x[i]**2-1)


          if z[i] < 0: print 'warning Ffunc'
       return 2*z
# ----------------------------------------------------------------------------

   # Function needed to calculate gamma for an NFW halo. 
   # Form is  long, but follows http://arxiv.org/pdf/astro-ph/9908213v1.pdf
def Gfunc(self,x):
       z=numpy.zeros(len(x))
       for i in range(len(x)):
          X=x[i]
          if x[i]==-1:
              z[i]=0.0
          elif x[i]>1:
             y=(((X-1)/(X+1))**.5)
             z[i]= (8* numpy.arctan(y) / (X**2*(X**2-1)**0.5)) +\
                 (4/X**2)*numpy.log(X/2) - \
                 2/(X**2-1) +\
                 4*numpy.arctan(y)/(((X**2)-1)**(3./2))
          else: 
             #print "WARNING You are very close to a halo"
             if x[i]==1:
                z[i] =(10./3+4*numpy.log(0.5))
             else:
                y=(((1-X)/(X+1))**.5)
                z[i]= (8* numpy.arctanh(y) / (X**2*(1-X**2)**0.5)) +\
                    (4/X**2)*numpy.log(X/2) - \
                    2/(X**2-1) +\
                    4*numpy.arctanh(y)/((X**2-1)*(1-X**2)**(1./2))
          if z[i]<0: print 'warning Gfunc'; print x[i]; print z[i]

       return z

# ========================================================================
# BMO profile functions.
# ------------------------------------------------------------------------
def BMO1Ffunc(x,t):
    z=numpy.zeros(len(x))
    for i in range(len(x)):
        if x[i]>1:
            z[i]= (t**2/(t**2+1)**2)              * ( \
                (t**2+1)/(x**2-1)*(1-F(x))+2*F(x) -   \
                3.14159/(t**2+x**2)               +   \ 
                ((t**2-1)/(t*(t**2+x**2)**.5))*L(x,t) )
        else: 
            z[i]=Ffunc([x[i]])
        if z[i] < 0: print 'warning BMO2Ffunc'
       return z

def BMO1GFunc(x,t)
