# Script to run the `LensByHalo` and `CalculateCorr` demos outside of the notebook

import numpy as np
import scipy as sp
import os,sys
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import cmath, cProfile
import cPickle as pickle
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar

# Turn on for verbose
vb = True

# Turn on for pickling correlation data
pickle = False

# Turn on for plotting maps
maps = True

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

# Generate Background Catalog in the middle of the (0,0,0,0) field
if vb is True: print('Generating background catalog...')
# Can pick one of the domains below
#d = [1.85,1.15,-1.85,-1.15]
#d = [1.75,1.25,-1.75,-1.25]
#d = [1.65,1.35,-1.65,-1.35]
#d = [1.6,1.4,-1.6,-1.4]
#d = [1.55,1.45,-1.55,-1.45]
d = [1.55,1.48,-1.55,-1.48]
#d = [1.55,1.52,-1.61,-1.59] # only galaxies in subplot
#d = [1.55,1.54,-1.61,-1.6] # ~3 galaxies
B = pangloss.BackgroundCatalog(N=10.0,sigma_e=0.01,domain=d,field=[0,0,0,0])

# Lens the background catalog by map
if vb is True: print('Lensing background by map...')
B.lens_by_map(K,S)
print 'Background catalog has',B.galaxy_count,'galaxies'

# Drill the lightcones
if vb is True: print('Drilling lightcones...')
B.drill_lightcones(radius=2.0,foreground=F,save=False)

# Calculate mean/std galaxies per lightcone
galaxy_counts = [B.lightcones[i].galaxy_count for i in range(len(B.lightcones))]
mean_galaxies = np.mean(galaxy_counts)
std_galaxies = np.std(galaxy_counts)
#important_counts = [B.lightcones[i] for i in range(len(B.lightcones))]
#mean_important = np.mean(galaxy_counts)
#std_important = np.std(galaxy_counts)
print 'Lightcones have',mean_galaxies,'+/-',std_galaxies,'galaxies'
#print 'Lightcones have',mean_important,'+/-',std_important,'important galaxies'

# Lens the background catalog by foreground halos
if vb is True: print('Lensing background by halos..')
cProfile.run('B.lens_by_halos(importance_lim=0.0); print')
#B.lens_by_halos(importance_lim=0.0)
'''
# Calculate the correlation function for each lensing type
if vb is True: print('Calculating correlation...')
gg_map = B.calculate_corr(corr_type='gg',lensed='map')
gg_halo = B.calculate_corr(corr_type='gg',lensed='halo')
ng_map = B.calculate_corr(corr_type='ng',lensed='map')
ng_halo = B.calculate_corr(corr_type='ng',lensed='halo')

# Plot the correlation functions
pangloss.plotting.plot_corr(gg_map,corr_type='gg',corr_comp='plus',lensed='map',color='green',galaxy_count=B.galaxy_count)
pangloss.plotting.plot_corr(gg_halo,corr_type='gg',corr_comp='plus',lensed='halo',color='purple')
pangloss.plotting.plot_corr(gg_map,corr_type='gg',corr_comp='cross',lensed='map',color='green')
pangloss.plotting.plot_corr(gg_halo,corr_type='gg',corr_comp='cross',lensed='halo',color='purple')
plt.gcf().set_size_inches(10,10)
plt.show()

pangloss.plotting.plot_corr(ng_map,corr_type='ng',corr_comp='real',lensed='map',color='green',galaxy_count=B.galaxy_count)
pangloss.plotting.plot_corr(ng_halo,corr_type='ng',corr_comp='real',lensed='halo',color='purple')
plt.gcf().set_size_inches(10,10)
plt.show()

# Compare the correlation functions
chi2,n_sigma,percent_err,std_err = B.compare_corr(gg_halo,gg_map,corr_type='gg',corr_comp='plus')
print 'Ellipticity-Ellipticity correlation difference intrinsic to mapped:','chi^2: ',chi2,'n_sigma: ',n_sigma,'percent_err: ',percent_err,'+\-',std_err

chi2,n_sigma,percent_err,std_err = B.compare_corr(ng_halo,ng_map,corr_type='ng',corr_comp='real')
print 'Galaxy-Galaxy correlation difference intrinsic to mapped:','chi^2: ',chi2,'n_sigma: ',n_sigma,'percent_err: ',percent_err,'+\-',std_err

# Plot a map near a lens
if maps is True:
    K.plot(fig_size=15,subplot=[1.55,1.52,-1.61,-1.59])
    S.plot()
    B.plot(lensed='all',graph='stick')
    plt.show()

# Save gg_map and gg_halo
if pickle is True:

    if vb is True: print('Pickling gg_map...')
    filename1 = PANGLOSS_DIR+'/data/tests/gg_map_'+str(B.galaxy_count)+'.pickle'
    gg1 = [gg_map.logr,gg_map.xip,gg_map.xim,gg_map.xip_im,gg_map.xim_im,gg_map.varxi]
    pickle_file1 = open(filename1, 'wb')
    pickle.dump(gg1,pickle_file1)
    pickle_file1.close()

    if vb is True: print('Pickling gg_halo...')
    filename2 = PANGLOSS_DIR+'/data/tests/gg_halo_'+str(B.galaxy_count)+'.pickle'
    gg2 = [gg_halo.logr,gg_halo.xip,gg_halo.xim,gg_halo.xip_im,gg_halo.xim_im,gg_halo.varxi]
    pickle_file2 = open(filename2, 'wb')
    pickle.dump(gg2,pickle_file2)
    pickle_file2.close()
'''
if vb is True: print('Closing...')
