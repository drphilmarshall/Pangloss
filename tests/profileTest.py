# profileTest script

import os, sys, math, numpy as np
import line_profiler
import cProfile

# Pangloss:
#PANGLOSS_DIR = os.path.expandvars("$PANGLOSS_DIR")
#sys.path.append(PANGLOSS_DIR)
#import pangloss

def nVidia_acos(x):
    #negate = float(x<0)
    negate = np.zeros(np.size(x))
    negate[x<0] = 1.0
    x=np.abs(x)
    ret = -0.0187293
    ret = ret * x
    ret = ret + 0.0742610
    ret = ret * x
    ret = ret - 0.2121144
    ret = ret * x
    ret = ret + 1.5707288
    ret = ret * np.sqrt(1.0-x)
    ret = ret - 2 * negate * ret
    return negate * 3.14159265358979 + ret

def fast_acos(x):
    '''
    From http://stackoverflow.com/questions/3380628/fast-arc-cos-algorithm
    '''
    return (-0.69813170079773212 * x * x - 0.87266462599716477) * x + 1.5707963267948966

#@profile
'''
def F(x):
    z=np.ones(len(x))
    s1=x[x>1]
    m1=1/s1
    #f1=np.arccos(m1)
    f1=fast_acos(m1)
    s2=x[x>1]
    m2=((s2**2-1)**.5)
    m3=f1/m2
    s3=z[x>1]

    s4=x[x<1]
    m4=1/s4
    f2=np.arccosh(m4)
    s5=x[x<1]
    m5=((1-s5**2)**.5)
    m6=f2/m5
    s6=z[x<1]


    m7 = np.log(2)
    s7=z[x==1]

    print np.size(s1),np.size(s4)

    return z
'''

def BMO1FSpencerFunc(x,t):
    '''
    New BMO1Ffunc method that is more efficient.
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

def F(x):
    z=np.ones(len(x))
    z[x>1]=np.arccos(1/x[x>1])/((x[x>1]**2-1)**.5)
    z[x<1]=np.arccosh(1/x[x<1])/((1-x[x<1]**2)**.5)
    z[x==1]=0.69314718 #np.log(2)
    return z

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
    return 4*z

x = np.arange(0.001,2,.001)
t = np.ones(np.size(x),dtype='double')

profiler = line_profiler.LineProfiler()
profiler.add_function(BMO1FSpencerFunc)
profiler.add_function(FSpencer)
profiler.add_function(L)
profiler.add_function(F)
profiler.add_function(BMO1Ffunc)
profiler.enable_by_count()
BMO1FSpencerFunc(x,t)
BMO1Ffunc(x,t)
profiler.print_stats()

cProfile.run('BMO1Ffunc(x,t); print')
cProfile.run('BMO1FSpencerFunc(x,t); print')
