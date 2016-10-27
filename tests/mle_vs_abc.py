# Script to test feasibility of maximum likelihood estimation (MLE) vs approximate
# Bayesian computation (ABC) for pangloss inference

import numpy as np
import scipy as sp
import os,sys,timeit
import matplotlib.pyplot as plt
import cmath, cProfile
import copy
import cPickle as pickle
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar

# Turn on for verbose
vb = True

# Turn on for recording times
time = True

# Turn on for plotting correlation function plots
corr_plots = True

# Pangloss:
PANGLOSS_DIR = os.path.expandvars("$PANGLOSS_DIR")
sys.path.append(PANGLOSS_DIR)
import pangloss

# Import Shear and Convergence maps
if vb is True: print('Loading Kappa and Shear maps...')
K = pangloss.Kappamap(PANGLOSS_DIR+'/data/GGL_los_8_0_0_N_4096_ang_4_rays_to_plane_37_f.kappa',FITS=False)
S = pangloss.Shearmap([PANGLOSS_DIR+'/data/GGL_los_8_0_0_N_4096_ang_4_rays_to_plane_37_f.gamma_1',PANGLOSS_DIR+'/data/GGL_los_8_0_0_N_4096_ang_4_rays_to_plane_37_f.gamma_2'],FITS=False)

# Import Foreground Catalog
if vb is True: print('Loading foreground catalog...')
config = pangloss.Configuration(PANGLOSS_DIR+'/example/example.config')
F = pangloss.ForegroundCatalog(PANGLOSS_DIR+'/data/GGL_los_8_0_0_0_0_N_4096_ang_4_Guo_galaxies_on_plane_27_to_63.images.txt',config)
#if smooth_plots is True: F.plot_mean_kappas()

# Generate Background Catalog in the middle of the (0,0,0,0) field
if vb is True: print('Generating background catalog...')

# Can pick one of the domains below
#d = [1.85,1.15,-1.85,-1.15]
#d = [1.75,1.25,-1.75,-1.25]
#d = [1.65,1.35,-1.65,-1.35]
d = [1.6,1.4,-1.6,-1.4] # 1440 galaxies
#d = [1.55,1.45,-1.55,-1.45]
#d = [1.55,1.48,-1.55,-1.48] # 176 galaxies
#d = [1.55,1.52,-1.61,-1.59] # only galaxies in subplot
#d = [1.55,1.54,-1.61,-1.6] # ~3 galaxies

# time data is stored in these arrays...
T = 20 # number of time data points
likelihood_time = []
likelihood_total_time = []
abc_time = []
abc_total_time = []

#
area = 0.2*0.2*(3600) #arcmin^2
n = 3*10**2 # number of sources
N = n / area # number of soruces per arcmin^2

for i in range(T):
    print 'Run {}'.format(i)
    # Start all clocks
    if time is True: start_time = timeit.default_timer()

    B = pangloss.BackgroundCatalog(N=N,sigma_e=0.01,domain=d,field=[0,0,0,0])

    # Lens the background catalog by map
    if vb is True: print('Lensing background by map...')
    B.lens_by_map(K,S)
    print 'Background catalog has',B.galaxy_count,'galaxies'

    ''' OLD: can use lens-by-map!
    # Drill the lightcones
    if vb is True: print('Drilling lightcones...')
    lc_radius = 0.5
    B.drill_lightcones(radius=lc_radius,foreground=F,save=False)

    # Calculate mean/std galaxies per lightcone
    galaxy_counts = [lightcone.galaxy_count for lightcone in B.lightcones]
    mean_galaxies = np.mean(galaxy_counts)
    std_galaxies = np.std(galaxy_counts)
    print 'Lightcones have {0:.2f} +/- {1:.2f} galaxies'.format(mean_galaxies,std_galaxies)

    # Lens the background catalog by foreground halos
    if vb is True: print('Lensing background by halos..')
    #relevance_lim = 0.0
    #relevance_lim = 0.0001
    relevance_lim = 0.000005
    smooth_corr = False
    #cProfile.run('B.lens_by_halos(relevance_lim=relevance_lim,lookup_table=True); print')
    B.lens_by_halos(relevance_lim=relevance_lim,lookup_table=True,smooth_corr=smooth_corr)
    print 'Lightcones have {0:.2f} +/- {1:.2f} relevant galaxies'.format(B.mean_relevant_halos,B.std_relevant_halos)
    '''

    # Calculate likelihood
    start_likelihood = timeit.default_timer()
    likelihood = B.calculate_log_likelihood(lensed='map')
    #cProfile.run('B.calculate_log_likelihood(lensed="map"); print')
    B.calculate_log_likelihood(lensed='map')
    likelihood_time.append(timeit.default_timer() - start_likelihood)
    likelihood_total_time.append(timeit.default_timer() - start_time)
    print 'lens_by_halos likelihood = {}; likelihood time was {} s, Total time taken was {} s'.format(likelihood,likelihood_time[i],likelihood_total_time[i])

    # Calculate e-e corr for ABC
    start_abc = timeit.default_timer()
    gg_halo = B.calculate_corr(corr_type='gg',lensed='map',foreground=F)
    # Calculate distance metric and make selection decision
    abc_time.append(timeit.default_timer() - start_abc)
    abc_total_time.append(timeit.default_timer() - start_time - likelihood_time[i])
    abc_decision = 'Rejected' # TEMPORARY!!
    print 'abc decision = {}; abc decision time was {} s, Total time taken (excluding likelihood estimation) was {} s'.format(abc_decision,abc_time[i],abc_total_time[i])

    del B

mean_likelihood = np.mean(likelihood_time)
std_likelihood = np.std(likelihood_time)

mean_likelihood_t = np.mean(likelihood_total_time)
std_likelihood_t = np.std(likelihood_total_time)

mean_abc = np.mean(abc_time)
std_abc = np.std(abc_time)

mean_abc_t = np.mean(abc_total_time)
std_abc_t = np.std(abc_total_time)
print 'likelihood_time = {}+/-{},\nlikelihood_total_time = {}+/-{},\nabc_time = {}+/-{},\nabc_total_time = {}+/-{}'.format(mean_likelihood,std_likelihood,mean_likelihood_t,std_likelihood_t,mean_abc,std_abc,mean_abc_t,std_abc_t)
print '{}'


if vb is True: print('Closing...')
