import pangloss
import numpy as np, gc
from scipy.interpolate import RectBivariateSpline
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
        Unused Spencer functions:
            FSpencer(x,xMask):
            BMO1FSpencerFunc(x,t):
            BMO1GSpencerFunc(x,t):

    BUGS
        - Some redundancy between F() and Ffunc() etc?

    AUTHORS
      This file is part of the Pangloss project, distributed under the
      GPL v2, by Tom Collett (IoA) and  Phil Marshall (Oxford).
      Please cite: Collett et al 2013, http://arxiv.org/abs/1303.6564

    HISTORY
      2013-03-25  Collett & Marshall (Oxford)
      2015-08-07  Everett (SLAC)
"""

class LensingTable():
    '''
    This small class creates a lookup table for the lensing values to save CPU time
    for lensing caluclations.
    '''

    def __init__(self,x_lim=[1e-5,500.0],t_lim=[50.0,200.0],N=1000):

        # Create lookup table
        self.generate_lookup_table(x_lim,t_lim,N)
        return

    def generate_lookup_table(self,x_lim=[0.0,500.0],t_lim=[50.0,200.0],N=1000):

        # Set x and t limits
        self.x_min = x_lim[0]
        self.x_max = x_lim[1]
        self.t_min = t_lim[0]
        self.t_max = t_lim[1]

        # Create stepsize based on N
        self.x_stepsize = abs(self.x_max-self.x_min)/(1.0*N)
        self.t_stepsize = abs(self.t_max-self.t_min)/(1.0*N)

        # Set the x and t axis
        x = np.arange(self.x_min,self.x_max,self.x_stepsize)
        t = np.arange(self.t_min,self.t_max,self.t_stepsize)

        # Create the lookup grid domain
        self.x_grid, self.t_grid = np.meshgrid(x,t)

        # Create the lookup grids
        self.BMO1F_grid = BMO1Ffunc(self.x_grid,self.t_grid)
        self.BMO1G_grid = BMO1Gfunc(self.x_grid,self.t_grid)

        # Create the spline interpolater
        self.BMO1F_spline = RectBivariateSpline(x,t,self.BMO1F_grid)
        self.BMO1G_spline = RectBivariateSpline(x,t,self.BMO1G_grid)

        return

#=========================================================================
# Lookup functions for BMO1 tables

    def lookup_BMO1F(self,x,t):
        assert (x > self.x_min).all() and (x < self.x_max).all()
        assert (t > self.t_min).all() and (t < self.t_max).all()

        # Found problem! Can use numpy arrays, but must be of form
        # input = ([x1,x2,x3,...],[y1,y2,y3,...])
        F = self.BMO1F_spline.ev(x,t)
        assert not np.isnan(F).any()

        return F

    def lookup_BMO1G(self,x,t):
        assert (x > self.x_min).all() and (x < self.x_max).all()
        assert (t > self.t_min).all() and (t < self.t_max).all()

        G = self.BMO1G_spline.ev(x,t)
        assert not np.isnan(G).any()

        return G


#=========================================================================
# NFW profile functions:

def delta_c(c):
    return (200./3)*(c**3)/(np.log(1+c)-c/(1+c))

def F(x):
    z=np.ones(np.shape(x))
    z[x>1]=np.arccos(1/x[x>1])/((x[x>1]**2-1)**.5)
    z[x<1]=np.arccosh(1/x[x<1])/((1-x[x<1]**2)**.5)
    z[x==1]=0.69314718 #np.log(2)

    assert not (np.isnan(z).any())
    return z

def FSpencer(x,xMask):

    # Calculate z in one step using masks
    z1 = np.arccos(1.0 / x*xMask) / (((x*xMask)**2-1.0)**.5 )
    z2 = np.arccosh(1.0 / x*~xMask) / ((1.0 - (x*~xMask)**2)**.5)

    # Convert any NaN's to zeros
    #z1 = np.nan_to_num(z1)
    z1[np.isnan(z1)] = 0.0
    #assert not (np.isnan(z1).any())
    #z2 = np.nan_to_num(z2)
    z2[np.isnan(z2)] = 0.0
    #assert not (np.isnan(z2).any())

    z = z1 + z2

    assert not (np.isnan(z).any())

    return z

def L(x,t):
    return np.log(x/(((t**2+x**2)**.5)+t))

def F2(x):
    z=np.ones(len(x))
    z[x>1]=np.arctan((x[x>1])**2-1)/((x[x>1]**2-1)**.5)
    z[x<1]=np.arctanh(1-(x[x<1])**2)/((1-x[x<1]**2)**.5)
    z[x==1]=0.69314718 #np.log(2)
    return z

# ------------------------------------------------------------------------
# Function needed to calculate kappa for an NFW halo.

def Ffunc(x):
    z=np.zeros(len(x))
    for i in range(len(x)):
       if x[i]==-1:
           z[i]=0.0
       elif x[i]>1:
          z[i]= (1-(2./(x[i]**2-1)**.5)*np.arctan(((x[i]-1.)/(x[i]+1))**.5))/(x[i]**2-1.)
       else:
          if x[i]==1:
             z[i] =1./3
          else:
             y=(-x[i]**2+1)**.5
             z[i] = (1.-(2./(1-x[i]**2)**.5)*np.arctanh(((1.-x[i])/(x[i]+1))**.5))/(x[i]**2-1)
       if z[i] < 0: print 'warning Ffunc' # BUG - non-informative alert
    return 2*z

# ------------------------------------------------------------------------
# Function needed to calculate gamma for an NFW halo.
# Form is  long, but follows http://arxiv.org/pdf/astro-ph/9908213v1.pdf

def Gfunc(x):
    z=np.zeros(len(x))
    for i in range(len(x)):
        X=x[i]
        if x[i]>1:
            y=(((X-1)/(X+1))**.5)
            z[i]= (8* np.arctan(y) / (X**2*(X**2-1)**0.5)) +\
                (4/X**2)*np.log(X/2) - \
                2/(X**2-1) +\
                4*np.arctan(y)/(((X**2)-1)**(3./2))
        else:
            if x[i]==1:
                z[i] =(10./3+4*np.log(0.5))
            else:
                y=(((1-X)/(X+1))**.5)
                z[i]= (8* np.arctanh(y) / (X**2*(1-X**2)**0.5)) +\
                    (4/X**2)*np.log(X/2) - \
                    2/(X**2-1) +\
                    4*np.arctanh(y)/((X**2-1)*(1-X**2)**(1./2))
        if z[i]<0: print 'warning Gfunc'; print x[i]; print z[i]
    return z

# ========================================================================
# Baltz, Marshall & Oguri truncated profile functions.

def BMO1FSpencerFunc(x,t):
    '''
    New BMO1Ffunc method using masks. NOTE: Didn't end up being more efficient.
    '''

    x[x==1]=1.+1e-5
    xMask = x > 1
    Fs = FSpencer(x,xMask)

    z = t**2 / (2*(t**2+1)**2) * (
        ((t**2+1) / ((x)**2-1)) * (1-Fs)
        + 2 * Fs
        - 3.14159 / (t**2+x*2)**.5
        + (t**2-1) * L(x,t)
        / (t * (t**2 + x**2)**.5)
        )

    assert not (np.isnan(z).any())

    return 4.0*z

def BMO1GSpencerFunc(x,t):
    '''
    New BMO1Gfunc method using masks. NOTE: Didn't end up being more efficient.
    '''

    z = np.zeros(len(x))
    x[x==1] = 1.+1e-5
    xMask = x > 1

    z = t**2 / ((t**2+1)**2) * (
        ((t**2+1)+2 * (x**2-1))*(FSpencer(x,xMask)) #possibly need -1 here!!
        + t * 3.14159
        + (t**2-1) * np.log(t)
        + ((t**2 + x**2)**.5) * (-3.14159 + (t**2-1) * L(x,t) / t)
        )

    assert not (np.isnan(z).any())

    return 4*z/(x**2)


def BMO1Ffunc(x,t):
    x[x==1]=1.+1e-5
    f = F(x)
    z = t**2/(2*(t**2+1)**2)*(
        ((t**2+1)/((x)**2-1))*(1-f)
        +
        2*f
        -
        3.14159/(t**2+x**2)**.5
        +
        (t**2-1)*L(x,t)
        /
        (t*(t**2+x**2)**.5)
        )
    if np.isnan(z).any():
        print 'x,t:',x,t
        print 'z:',z
    assert not (np.isnan(z).any())

    return 4*z

# ------------------------------------------------------------------------

def BMO1Gfunc(x,t):
    x[x==1]=1.+1e-5
    z = t**2/((t**2+1)**2)*(
        ((t**2+1)+2*(x**2-1))*F(x) #possibly need -1 here!!
        +
        t*3.14159
        +
        (t**2-1)*np.log(t)
        +
        ((t**2+x**2)**.5)*(-3.14159+(t**2-1)*L(x,t)/t)
        )
    return 4*z/(x**2)

# ------------------------------------------------------------------------

def BMO2Ffunc(x,t):
    #t=float(t)
    x[x==1]=1.+1e-5
    z=np.zeros(len(x))
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
    z=np.zeros(len(x))
    x[x==1]=1.+1e-5
    z=(t**4/(2*(t**2+1)**3))*(
        (t**2+1+4*(x**2-1))*(2*F(x))
        +
        (1/t)*(3.14159*(3*t**2-1)+2*t*(t**2-3)*np.log(t))
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
#     amp = amp/((re**2)*np.exp(k)*n*(k**(-2*n))*gamma(2*n)*2*np.pi)
#     R = r/re
#     return amp*np.exp(-k*(R**(1./n)-1.))

def sersic(r,re,amp=1.,n=4.):
     k = 2.*n-1./3+4./(405.*n)+46/(25515.*n**2)
     amp = amp/((re**2)*np.exp(k)*n*(k**(-2*n))*gamma(2*n)*2*np.pi)
     R = r/re
     kappa = amp*np.exp(-k*(R**(1./n)-1.))
     kbar = 2*n*amp*np.exp(k)*gamma(2*n)*gammainc(2*n,k*R**(1./n))/(R**2*k**(2*n))
     return kappa,kbar-kappa

# ========================================================================
