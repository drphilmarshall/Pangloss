
import pangloss

import numpy
from scipy.special import gamma,gammainc # for Sersic profile.

#=========================================================================

"""
    NAME
        lensing

    PURPOSE
        Compute gravitational lensing quantities.

    COMMENTS
            
    FUNCTIONS
        NFW profile:
            delta_c(c):
            F(x):
            L(x,t):
            F2(x):
            Ffunc(x):
            Gfunc(x): 
        Baltz, Marshall & Oguri truncated NFW profile:
            BMO1Ffunc(x,t):
            BMO1Gfunc(x,t):
            BMO2Ffunc(x,t):
            BMO2Gfunc(x,t):
        Sersic profile:
            sersic(r,re,amp=1.,n=4.):

    BUGS
        - Some redundancy between F() and Ffunc() etc?

    AUTHORS
      This file is part of the Pangloss project, distributed under the
      GPL v2, by Tom Collett (IoA) and  Phil Marshall (Oxford). 
      Please cite: Collett et al 2013, http://arxiv.org/abs/1303.6564

    HISTORY
      2013-03-25  Collett & Marshall (Oxford)
"""

#=========================================================================
# NFW profile functions:

def delta_c(c):
    return (200./3)*(c**3)/(numpy.log(1+c)-c/(1+c))

def F(x):
    z=numpy.ones(len(x))
    z[x>1]=numpy.arccos(1/x[x>1])/((x[x>1]**2-1)**.5)
    z[x<1]=numpy.arccosh(1/x[x<1])/((1-x[x<1]**2)**.5)
    z[x==1]=numpy.log(2)
    return z

def L(x,t):
    return numpy.log(x/(((t**2+x**2)**.5)+t))

def F2(x):
    z=numpy.ones(len(x))
    z[x>1]=numpy.arctan((x[x>1])**2-1)/((x[x>1]**2-1)**.5)
    z[x<1]=numpy.arctanh(1-(x[x<1])**2)/((1-x[x<1]**2)**.5)
    z[x==1]=numpy.log(2)
    return z

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
          if x[i]==1:
             z[i] =1./3
          else:
             y=(-x[i]**2+1)**.5
             z[i] = (1.-(2./(1-x[i]**2)**.5)*numpy.arctanh(((1.-x[i])/(x[i]+1))**.5))/(x[i]**2-1)
       if z[i] < 0: print 'warning Ffunc' # BUG - non-informative alert
    return 2*z
    
# ------------------------------------------------------------------------
# Function needed to calculate gamma for an NFW halo. 
# Form is  long, but follows http://arxiv.org/pdf/astro-ph/9908213v1.pdf

def Gfunc(x): 
    z=numpy.zeros(len(x))
    for i in range(len(x)):
        X=x[i]
        if x[i]>1:
            y=(((X-1)/(X+1))**.5)
            z[i]= (8* numpy.arctan(y) / (X**2*(X**2-1)**0.5)) +\
                (4/X**2)*numpy.log(X/2) - \
                2/(X**2-1) +\
                4*numpy.arctan(y)/(((X**2)-1)**(3./2))
        else: 
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
# Baltz, Marshall & Oguri truncated profile functions.

def BMO1Ffunc(x,t):
    x[x==1]=1.+1e-5
    z=numpy.zeros(len(x))
    z[x!=1]=t**2/(2*(t**2+1)**2)*(
        ((t**2+1)/((x[x!=1])**2-1))*(1-F(x[x!=1]))
        +
        2*F(x[x!=1])
        -
        3.14159/(t**2+x[x!=1]**2)**.5
        +
        (t**2-1)*L(x[x!=1],t)
        /
        (t*(t**2+x[x!=1]**2)**.5)
        )
    return 4*z

# ------------------------------------------------------------------------

def BMO1Gfunc(x,t):
    z=numpy.zeros(len(x))
    x[x==1]=1.+1e-5
    z[x!=1]=t**2/((t**2+1)**2)*(
        ((t**2+1)+2*(x[x!=1]**2-1))*(F(x[x!=1])) #possibly need -1 here!!
        +
        t*3.14159
        +
        (t**2-1)*numpy.log(t)
        +
        ((t**2+x[x!=1]**2)**.5)*(-3.14159+(t**2-1)*L(x[x!=1],t)/t)
        )
    return 4*z/(x**2)

# ------------------------------------------------------------------------

def BMO2Ffunc(x,t):
    #t=float(t)
    x[x==1]=1.+1e-5
    z=numpy.zeros(len(x))
    z=t**4/(4*(t**2+1)**3)*(
        (2*(t**2+1)/((x)**2-1))*(1-F(x))
        +
        8*F(x)
        +
        (t**4-1)/((t**2)*(t**2+x**2))
        -
        3.14159*(4*(t**2+x**2)+t**2+1)/(t**2+x**2)**1.5
        +
        ((t**2)*(t**4-1)+(t**2+x**2)*(3*(t**4)-6*(t**2)-1))*L(x,t)
        /
        ((t**3)*(t**2+x**2)**1.5)
        )
    return 4*z

# ------------------------------------------------------------------------

def BMO2Gfunc(x,t):
    #t=float(t)
    z=numpy.zeros(len(x))
    x[x==1]=1.+1e-5
    z=(t**4/(2*(t**2+1)**3))*(
        (t**2+1+4*(x**2-1))*(2*F(x)) 
        +
        (1/t)*(3.14159*(3*t**2-1)+2*t*(t**2-3)*numpy.log(t))
        +
        ((-3.14159*t**3)*(4*(t**2+x**2)-t**2-1))
        /
        (t**3*(t**2+x**2)**.5)
        +
        (-t**2*(t**4-1)+(t**2+x**2)*(3*t**4-6*t**2-1))*L(x,t)
        )
    return 4*z/(x**2)

# ========================================================================
# de Vaucelour profile functions.

#def sersic(r,re,amp=1.,n=4.):
#     k = 2.*n-1./3+4./(405.*n)+46/(25515.*n**2)
#     amp = amp/((re**2)*numpy.exp(k)*n*(k**(-2*n))*gamma(2*n)*2*numpy.pi)
#     R = r/re
#     return amp*numpy.exp(-k*(R**(1./n)-1.))

def sersic(r,re,amp=1.,n=4.):
     k = 2.*n-1./3+4./(405.*n)+46/(25515.*n**2)
     amp = amp/((re**2)*numpy.exp(k)*n*(k**(-2*n))*gamma(2*n)*2*numpy.pi)
     R = r/re
     kappa = amp*numpy.exp(-k*(R**(1./n)-1.))
     kbar = 2*n*amp*numpy.exp(k)*gamma(2*n)*gammainc(2*n,k*R**(1./n))/(R**2*k**(2*n))
     return kappa,kbar-kappa

# ========================================================================

