# Script to test the 'bin_to_maps()' method
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import cPickle as pickle
import pyfits
import os, sys, cmath, cProfile, copy

# Turn on for verbose
vb = True

# Pangloss:
PANGLOSS_DIR = os.path.expandvars("$PANGLOSS_DIR")
sys.path.append(PANGLOSS_DIR)
import pangloss

# Import Shear and Convergence maps
if vb is True: print('Loading Kappa and Shear maps...')
K = pangloss.Kappamap(kappafile=[PANGLOSS_DIR+'/data/GGL_los_8_0_0_N_4096_ang_4_rays_to_plane_37_f.kappa.fits'],FITS=True)
S = pangloss.Shearmap(shearfiles=[PANGLOSS_DIR+'/data/GGL_los_8_0_0_N_4096_ang_4_rays_to_plane_37_f.gamma_1.fits',PANGLOSS_DIR+'/data/GGL_los_8_0_0_N_4096_ang_4_rays_to_plane_37_f.gamma_2.fits'],FITS=True)


# Import Foreground Catalog
if vb is True: print('Loading foreground catalog...')
config = pangloss.Configuration(PANGLOSS_DIR+'/example/example.config')
F = pangloss.ForegroundCatalog(PANGLOSS_DIR+'/data/GGL_los_8_0_0_0_0_N_4096_ang_4_Guo_galaxies_on_plane_27_to_63.images.txt',config)

# Generate Background Catalog in the middle of the (0,0,0,0) field
if vb is True: print('Generating background catalog...')
# Can pick one of the domains below
#d = [2.,1.,-2.,-1.]
#d = [1.,0.,-1.,0.]
#d = [1.85,1.15,-1.85,-1.15]
#d = [1.75,1.25,-1.75,-1.25]
#d = [1.65,1.35,-1.65,-1.35]
#d = [1.6,1.4,-1.6,-1.4] # 1440 galaxies
#d = [1.55,1.45,-1.55,-1.45]
#d = [1.55,1.48,-1.55,-1.48] # 176 galaxies
d = [1.20,1.15,-1.50,-1.45]
#d = [1.55,1.52,-1.61,-1.59] # only galaxies in subplot
#d = [1.55,1.54,-1.61,-1.6] # ~3 galaxies

K.plot(subplot=d,coords='world')
plt.savefig(PANGLOSS_DIR+'/data/binned_maps/'+'kappa_hilbert', bbox_inches='tight')
plt.close(plt.gcf())

B = pangloss.BackgroundCatalog(N=60.0,sigma_e=0.01,domain=d,field=[0,0,0,0],spacing=np.deg2rad(K.PIXSCALE[0]))

# Lens the background catalog by map
if vb is True: print('Lensing background by map...')
B.lens_by_map(K,S)
print 'Background catalog has',B.galaxy_count,'galaxies'

# Binsize
binsize = 0.075
#binsize = 0.05

# Bin kappa and shear maps
Kmap, Smap = B.bin_to_maps(lensed='map',binsize=binsize,savefile='kappa_map')

# Drill the lightcones
if vb is True: print('Drilling lightcones...')
lc_radius = 6.0
smooth_corr = True
B.drill_lightcones(radius=lc_radius,foreground=F,save=False,smooth_corr=smooth_corr)

# Calculate mean/std galaxies per lightcone
galaxy_counts = [lightcone.galaxy_count for lightcone in B.lightcones]
mean_galaxies = np.mean(galaxy_counts)
std_galaxies = np.std(galaxy_counts)
print 'Lightcones have {0:.2f} +/- {1:.2f} galaxies'.format(mean_galaxies,std_galaxies)

# Lens the background catalog by foreground halos
if vb is True: print('Lensing background by halos..')
relevance_lim = 0.0
#relevance_lim = 10**-5
B.lens_by_halos(relevance_lim=relevance_lim,lookup_table=True,smooth_corr=smooth_corr)
print 'Lightcones have {0:.2f} +/- {1:.2f} relevant galaxies'.format(B.mean_relevant_halos,B.std_relevant_halos)

# Bin kappa and shear maps
Khalo, Shalo = B.bin_to_maps(lensed='halo',binsize=binsize,savefile='kappa_halo')

# Make difference kappamap
kappadata = np.array([Kmap.values[0] - Khalo.values[0]])
assert Kmap.map_x == Khalo.map_x
assert Kmap.map_y == Khalo.map_y
map_xy = [Kmap.map_x, Kmap.map_y]
ra_i, dec_i = Kmap.image2world(0,0)
ra_f, dec_f = Kmap.image2world(Kmap.NX[0],Kmap.NX[0])
domain = [ra_i, ra_f, dec_i, dec_f]
Kdiff = pangloss.Kappamap(data=[kappadata,domain,map_xy])
Kdiff.plot(subplot=d,coords='world')
plt.savefig(PANGLOSS_DIR+'/data/binned_maps/'+'kappa_difference', bbox_inches='tight')

### Testing
print 'mean Khalo = {}'.format(np.mean(Khalo.values[0]))
Khalo.values[0] = Khalo.values[0] - np.mean(Khalo.values[0])

# Make 3 subplots of all maps
pangloss.plotting.compare_binned_maps(Kmap,Khalo,fig_size=20,savefile='kappamap_plots')
