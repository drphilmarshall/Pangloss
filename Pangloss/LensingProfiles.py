import pylab
import matplotlib.pyplot as plt
import numpy, numpy.random as rnd, atpy
import distances
from mpl_toolkits.axes_grid1 import ImageGrid
from time import clock
from numpy import exp,pi
from scipy.special import gamma,gammainc

#import time
#t0=time.clock()    

D = distances.Distance()
D.h = 0.7

arcmin2rad = (1.0/60.0)*numpy.pi/180.0
rad2arcmin = 1.0/arcmin2rad
# ------------------------------------------------------------------------
def F(x):
    z=numpy.ones(len(x))
    z[x>1]=numpy.arccos(1/x[x>1])/((x[x>1]**2-1)**.5)
    lowx=x[x<1]
    z[x<1]=numpy.arccosh(1/x[x<1])/((1-x[x<1]**2)**.5)
    A=1.+1e-5
    z[x==1]=numpy.arccos(1/A)/((A**2-1)**.5)
    return z

def L(x,t):
    return numpy.log(x/(((t**2+x**2)**.5)+t))
#=========================================================================
# NFW profile functions.

#--------------------------------------------------------------
def delta_c(c):
    return (200./3)*(c**3)/(numpy.log(1+c)-c/(1+c))
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
def Gfunc(x): 
       z=numpy.zeros(len(x))
       for i in range(len(x)):
          X=x[i]
          if x[i]==-1:
              z[i]=0.0 ###THIS IS NOT TRUE, SHOULD BE GFUNC OF A POINT MASS HERE.####
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
    #t=float(t)
    z=numpy.zeros(len(x))
    z[x!=1]=t**2/(2*(t**2+1)**2)*(
        ((t**2+1)/(x[x!=1]**2-1))*(1-F(x[x!=1]))
        +
        2*F(x[x!=1])
        -
        3.14159/(t**2+x[x!=1]**2)**.5
        +
        (t**2-1)*L(x[x!=1],t)
        /
        (t*(t**2+x[x!=1]**2)**.5)
        )
    delta=numpy.array([1.+1e-5])
    z[x==1]=t**2/(2*(t**2+1)**2)*(
        ((t**2+1)/(delta**2-1))*(1-F(delta))
        +
        2*F(delta)
        -
        3.14159/(t**2+delta**2)**.5
        +
        (t**2-1)*L(delta,t)
        /
        (t*(t**2+delta**2)**.5)
        )

    return 4*z
# ------------------------------------------------------------------------

# ========================================================================
# BMO2 profile functions.
# ------------------------------------------------------------------------




# ========================================================================
# de Vaucelour profile functions.
# ------------------------------------------------------------------------

#def sersic(r,re,amp=1.,n=4.):
#     k = 2.*n-1./3+4./(405.*n)+46/(25515.*n**2)
#     amp = amp/((re**2)*exp(k)*n*(k**(-2*n))*gamma(2*n)*2*pi)
#     R = r/re
#     return amp*exp(-k*(R**(1./n)-1.))

def sersic(r,re,amp=1.,n=4.):
     k = 2.*n-1./3+4./(405.*n)+46/(25515.*n**2)
     amp = amp/((re**2)*exp(k)*n*(k**(-2*n))*gamma(2*n)*2*pi)
     R = r/re
     kappa = amp*exp(-k*(R**(1./n)-1.))
     kbar = 2*n*amp*exp(k)*gamma(2*n)*gammainc(2*n,k*R**(1./n))/(R**2*k**(2*n))
     return kappa,kbar-kappa




# ========================================================================
# 
# ------------------------------------------------------------------------
if __name__ == '__main__':
    l=numpy.linspace(0.001,48,1001)
    t=numpy.ones(len(l))*10
    #print F(l)
    #print L(l,t)
    t0=clock()
    f1=BMO1Ffunc(l,t)
    t1=clock()
    print t1-t0
    f2=Ffunc(l)
    t2=clock()
    print t2-t1
    plt.plot(l,numpy.log10(f2))
    plt.plot(l,numpy.log10(f1))
    plt.show()
    plt.plot(l,f1/f2)
    plt.show()

