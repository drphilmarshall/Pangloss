# -*- coding: utf-8 -*-
# A quadratic fit
#———————————————————–
import numpy as np
import pymc
import pylab as plt

from math import pi
from scipy.stats import lognorm


# create some test data
s = 0.95
mean, var, skew, kurt = lognorm.stats(s, moments='mvsk')
x = np.linspace(lognorm.ppf(0.001, s), lognorm.ppf(0.99, s), 100)
sample = lognorm.rvs(s, size=500)
f = lognorm.pdf(x, s)
plt.figure()

plt.plot(x,f, label='pure data')
np.random.seed(76523654)
noise = np.random.normal(size=100) * .01     # create some Gaussian noise
f = f + noise                                # add noise to the data
#plt.plot(x,f, label='add noise')

shape, loc, scale = lognorm.fit(sample) # Fit a curve to the variates
mean = np.log(scale) # Mean of log(X)
sigma = shape # Standard deviation of log(X)
M = np.exp(mean) # Geometric mean == median
s = np.exp(sigma) # Geometric standard deviation

# Plot figure of results
plt.plot(x, lognorm.pdf(x, shape, loc, scale=scale), 'k--',label=r'$\Chi^2$ fit')
plt.legend(loc=2)
plt.show()


#priors
sig = pymc.Uniform('sig', 0.0, 100.0, value=1.)

a = pymc.Uniform('a', -10.0, 10.0, value= 0.0)
b = pymc.Uniform('b', -10.0, 10.0, value= 0.0)
x0 = pymc.Uniform('x0', -10.0, 10.0, value= 0.0)

#model
@pymc.deterministic(plot=False)
def model_lognormal(x=x, x0=x0, a=a, b=b):
      return np.exp((-(np.log(x-x0)-a)**2.)/(2*b**2.))/(np.sqrt(2*pi)*(x-x0)*b)


#likelihood
y = pymc.Normal('y', mu=model_lognormal, tau=1.0/sig**2, value=f, observed=True)
#———————————————————–

#Now, go to command line and run the following (or alternatively put them in a file):

import pymc, test              # load the model file
R = pymc.MCMC(test)    #  build the model
R.sample(10000)              # populate and run it
print 'a   ', R.a.value        # print outputs
print 'b    ', R.b.value
print 'x0    ', R.x0.value