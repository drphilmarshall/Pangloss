# Script used to create correlation function progression plots

import numpy as np
import scipy as sp
import os,sys,timeit
import matplotlib.pyplot as plt
import cmath, cProfile
import copy
import cPickle as pickle
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar

# Fast correlation functions:
try:
    import treecorr
except ImportError:
    import pangloss.nocorr as treecorr

# Turn on for verbose
vb = True

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
B = pangloss.BackgroundCatalog(N=10.0,sigma_e=0.01,domain=d,field=[0,0,0,0])

# Lens the background catalog by map
if vb is True: print('Lensing background by map...')
B.lens_by_map(K,S)
print 'Background catalog has',B.galaxy_count,'galaxies'

# The lightcone radii used for comparison
radius = [1, 2, 4, 8]
line_style = [':','-.','--','-']

# Loop over process for different lightcone radius
for lc_radius in radius:
    print 'Radius = {}'.format(lc_radius)
    # Drill the lightcones
    if vb is True: print('Drilling lightcones...')
    B.drill_lightcones(radius=lc_radius,foreground=F,save=False)
    # Calculate mean/std galaxies per lightcone
    galaxy_counts = [lightcone.galaxy_count for lightcone in B.lightcones]
    mean_galaxies = np.mean(galaxy_counts)
    std_galaxies = np.std(galaxy_counts)
    print 'Lightcones have {0:.2f} +/- {1:.2f} galaxies'.format(mean_galaxies,std_galaxies)

    # Lens the background catalog by foreground halos
    if vb is True: print('Lensing background by halos..')
    relevance_lim = 0.0
    #relevance_lim = 0.00001
    smooth_corr = True
    #cProfile.run('B.lens_by_halos(relevance_lim=relevance_lim,lookup_table=True); print')
    B.lens_by_halos(relevance_lim=relevance_lim,lookup_table=True,smooth_corr=smooth_corr)
    print 'Lightcones have {0:.2f} +/- {1:.2f} relevant galaxies'.format(B.mean_relevant_halos,B.std_relevant_halos)

    # Calculate the correlation function for each lensing type
    if vb is True: print('Calculating correlation...')
    if lc_radius == radius[-1]:
        gg_map = B.calculate_corr(corr_type='gg',lensed='map',foreground=F)
        ng_map = B.calculate_corr(corr_type='ng',lensed='map',foreground=F)
    gg_halo = B.calculate_corr(corr_type='gg',lensed='halo',foreground=F)
    ng_halo = B.calculate_corr(corr_type='ng',lensed='halo',foreground=F)

    # Save correlation data...
    if lc_radius == radius[-1]:
        if vb is True: print('Pickling gg_map...')
        filename1 = PANGLOSS_DIR+'/data/tests/gg_map_'+str(B.galaxy_count)+'_'+str(lc_radius)+'.pickle'
        gg1 = [gg_map.logr,gg_map.xip,gg_map.xim,gg_map.xip_im,gg_map.xim_im,gg_map.varxi]
        pickle_file1 = open(filename1, 'wb')
        pickle.dump(gg1,pickle_file1)
        pickle_file1.close()

    if vb is True: print('Pickling gg_halo...')
    filename2 = PANGLOSS_DIR+'/data/tests/gg_halo_'+str(B.galaxy_count)+'_'+str(lc_radius)+'.pickle'
    gg2 = [gg_halo.logr,gg_halo.xip,gg_halo.xim,gg_halo.xip_im,gg_halo.xim_im,gg_halo.varxi]
    pickle_file2 = open(filename2, 'wb')
    pickle.dump(gg2,pickle_file2)
    pickle_file2.close()

    if lc_radius == radius[-1]:
        if vb is True: print('Pickling ng_map...')
        filename1 = PANGLOSS_DIR+'/data/tests/ng_map_'+str(B.galaxy_count)+'_'+str(lc_radius)+'.pickle'
        ng1 = [ng_map.logr,ng_map.xi,ng_map.xi_im,ng_map.varxi]
        pickle_file1 = open(filename1, 'wb')
        pickle.dump(ng1,pickle_file1)
        pickle_file1.close()

    if vb is True: print('Pickling ng_halo...')
    filename2 = PANGLOSS_DIR+'/data/tests/ng_halo_'+str(B.galaxy_count)+'_'+str(lc_radius)+'.pickle'
    ng2 = [ng_halo.logr,ng_halo.xi,ng_halo.xi_im,ng_halo.varxi]
    pickle_file2 = open(filename2, 'wb')
    pickle.dump(ng2,pickle_file2)
    pickle_file2.close()

# Plot parameters
N = 15.0
min_sep = 0.1
max_sep = 2.0
sep_units = 'arcmin'
binsize = binsize = np.log10(1.0*max_sep/min_sep)/(1.0*N)

# Make the gg correlation progression plot
filename = PANGLOSS_DIR+'/data/tests/gg_map_'+str(B.galaxy_count)+'_'+str(radius[-1])+'.pickle'
gg1 = pickle.load(open(filename, 'rb'))
gg_map = treecorr.GGCorrelation(bin_size=binsize, min_sep=min_sep, max_sep=max_sep, sep_units=sep_units, bin_slop=0.05/binsize)
gg_map.logr,gg_map.xip,gg_map.xim,gg_map.xip_im,gg_map.xim_im,gg_map.varxi = gg1[0],gg1[1],gg1[2],gg1[3],gg1[4],gg1[5]
pangloss.plotting.plot_corr(gg_map,corr_type='gg',corr_comp='plus',lensed='map',color='green',galaxy_count=B.galaxy_count)

for i in range(len(radius)-1,-1,-1):
    filename = PANGLOSS_DIR+'/data/tests/gg_halo_'+str(B.galaxy_count)+'_'+str(radius[i])+'.pickle'
    gg2 = pickle.load(open(filename, 'rb'))
    # Plot the gg correlation function
    gg_halo = treecorr.GGCorrelation(bin_size=binsize, min_sep=min_sep, max_sep=max_sep, sep_units=sep_units, bin_slop=0.05/binsize)
    gg_halo.logr,gg_halo.xip,gg_halo.xim,gg_halo.xip_im,gg_halo.xim_im,gg_halo.varxi = gg2[0],gg2[1],gg2[2],gg2[3],gg2[4],gg2[5]
    ls = line_style[i]
    pangloss.plotting.plot_corr(gg_halo,corr_type='gg',corr_comp='plus',lensed='halo',color='purple',line_style=ls,radius=radius[i])

# Append legend labels
handles, labels = plt.gca().get_legend_handles_labels()

# remove the errorbars
handles = [h[0] for h in handles]

rad = [radius]
for i in range(len(labels)):
    if i == 0: rad = radius[-1]
    else: rad = radius[-i]
    labels[i] += r' ($R=${})'.format(rad)

plt.legend(handles, labels)

# Save, show, and close plot
plt.gcf().set_size_inches(10,10)
plt.savefig(PANGLOSS_DIR+'/data/tests/gg_progression', bbox_inches='tight')
if corr_plots is True: plt.show()
plt.close(plt.gcf())

# Make the ng correlation progression plot
filename = PANGLOSS_DIR+'/data/tests/ng_map_'+str(B.galaxy_count)+'_'+str(radius[-1])+'.pickle'
ng1 = pickle.load(open(filename, 'rb'))
ng_map = treecorr.GGCorrelation(bin_size=binsize, min_sep=min_sep, max_sep=max_sep, sep_units=sep_units, bin_slop=0.05/binsize)
ng_map.logr,ng_map.xi,ng_map.xi_im,ng_map.varxi = ng1[0],ng1[1],ng1[2],ng1[3]
pangloss.plotting.plot_corr(ng_map,corr_type='ng',corr_comp='real',lensed='map',color='green',galaxy_count=B.galaxy_count)

for i in range(len(radius)-1,-1,-1):
    filename = PANGLOSS_DIR+'/data/tests/ng_halo_'+str(B.galaxy_count)+'_'+str(radius[i])+'.pickle'
    ng2 = pickle.load(open(filename, 'rb'))
    # Plot the gg correlation function
    ng_halo = treecorr.NGCorrelation(bin_size=binsize, min_sep=min_sep, max_sep=max_sep, sep_units=sep_units, bin_slop=0.05/binsize)
    ng_halo.logr,ng_halo.xi,ng_halo.xi_im,ng_halo.varxi = ng2[0],ng2[1],ng2[2],ng2[3]
    ls = line_style[i]
    pangloss.plotting.plot_corr(ng_halo,corr_type='ng',corr_comp='real',lensed='halo',color='purple',line_style=ls,radius=radius[i])

# Append legend labels
handles, labels = plt.gca().get_legend_handles_labels()

# remove the errorbars
handles = [h[0] for h in handles]

rad = [radius]
for i in range(len(labels)):
    if i == 0: rad = radius[-1]
    else: rad = radius[-i]
    labels[i] += r' ($R=${})'.format(rad)

plt.legend(handles, labels)

# Save, show, and close plot
plt.gcf().set_size_inches(10,10)
plt.savefig(PANGLOSS_DIR+'/data/tests/ng_progression', bbox_inches='tight')
if corr_plots is True: plt.show()
plt.close(plt.gcf())
