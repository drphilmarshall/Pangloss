import pymc as mc
import numpy as np
import pylab as pl
def GaussFunc(x, amplitude, centroid, sigma):
    return amplitude * np.exp(-0.5 * ((x - centroid) / sigma)**2)

wavelength = np.arange(5000, 5050, 0.02)

# Profile 1
centroid_one = 5025.0
sigma_one = 2.2
height_one = 0.8
profile1 = GaussFunc(wavelength, height_one, centroid_one, sigma_one, )

# Profile 2
centroid_two = 5027.0
sigma_two = 1.2
height_two = 0.5
profile2 = GaussFunc(wavelength, height_two, centroid_two, sigma_two, )

# Measured values
noise = np.random.normal(0.0, 0.02, len(wavelength))
combined = profile1 + profile2 + noise

# If you want to plot what this looks like
pl.plot(wavelength, combined, label="Measured")
pl.plot(wavelength, profile1, color='red', linestyle='dashed', label="1")
pl.plot(wavelength, profile2, color='green', linestyle='dashed', label="2")
pl.title("Feature One and Two")
pl.legend()

### Suggested one runs the above code first.
### Unknowns we are interested in

est_centroid_one = mc.Uniform("est_centroid_one", 5000, 5050 )
est_centroid_two = mc.Uniform("est_centroid_two", 5000, 5050 )

est_sigma_one = mc.Uniform( "est_sigma_one", 0, 5 )
est_sigma_two = mc.Uniform( "est_sigma_two", 0, 5 )

est_height_one = mc.Uniform( "est_height_one", 0, 5 ) 
est_height_two = mc.Uniform( "est_height_two", 0, 5 ) 

#std deviation of the noise, converted to precision by tau = 1/sigma**2
precision= 1./mc.Uniform("std", 0, 1)**2

#Set up the model's relationships.

@mc.deterministic( trace = False) 
def est_profile_1(x = wavelength, centroid = est_centroid_one, sigma = est_sigma_one, height= est_height_one):
    return GaussFunc( x, height, centroid, sigma )


@mc.deterministic( trace = False) 
def est_profile_2(x = wavelength, centroid = est_centroid_two, sigma = est_sigma_two, height= est_height_two):
    return GaussFunc( x, height, centroid, sigma )


@mc.deterministic( trace = False )
def mean( profile_1 = est_profile_1, profile_2 = est_profile_2 ):
    return profile_1 + profile_2


observations = mc.Normal("obs", mean, precision, value = combined, observed = True)


model = mc.Model([est_centroid_one, 
              est_centroid_two, 
                est_height_one,
                est_height_two,
                est_sigma_one,
                est_sigma_two,
                precision])

#always a good idea to MAP it prior to MCMC, so as to start with good initial values
map_ = mc.MAP( model )
map_.fit()

mcmc = mc.MCMC( model )
mcmc.sample( 50000,40000 ) 