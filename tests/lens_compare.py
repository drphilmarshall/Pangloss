# Script to run the `LensByHalo` and `CalculateCorr` demos outside of the notebook

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

# Turn on for likelihood estimation time
time = False

# Turn on for pickling correlation data
pickle = False

# Turn on for relevance plots
rel_plots = False

# Turn on for smooth-component correction plots
smooth_plots = False

# Turn on for plotting correlation function plots
corr_plots = True

# Turn on for plotting correlation function comparrison plots for relevant halos
rel_compare = False

# Turn on for plotting correlation function comparrison plots for smooth-component correctoin
smooth_compare = False

# Turn on for plotting maps
maps = False

# time should only be used if `rel_compare` and `smooth_compare` are False!
if time is True: assert rel_compare is False and smooth_compare is False

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

# Start likelihood clock
if time is True: start_time = timeit.default_timer()

# Lens the background catalog by map
if vb is True: print('Lensing background by map...')
B.lens_by_map(K,S)
print 'Background catalog has',B.galaxy_count,'galaxies'

# Drill the lightcones
if vb is True: print('Drilling ligh7cones...')
lc_radius = 4.0
smooth_corr = True
B.drill_lightcones(radius=lc_radius,foreground=F,save=False,smooth_corr=smooth_corr)

# Save copy of background lightcones for relevant halos:
if rel_compare is True: B_rel = copy.deepcopy(B)
# Save copy of background lightcones for smooth-comp correction halos:
if smooth_compare is True: B_smooth = copy.deepcopy(B)

# Calculate mean/std galaxies per lightcone
galaxy_counts = [lightcone.galaxy_count for lightcone in B.lightcones]
mean_galaxies = np.mean(galaxy_counts)
std_galaxies = np.std(galaxy_counts)
print 'Lightcones have {0:.2f} +/- {1:.2f} galaxies'.format(mean_galaxies,std_galaxies)

# Lens the background catalog by foreground halos
if vb is True: print('Lensing background by halos..')
#relevance_lim = 0.0
#relevance_lim = 0.00001
relevance_lim = 10**-8
#cProfile.run('B.lens_by_halos(relevance_lim=relevance_lim,lookup_table=True,smooth_corr=smooth_corr); print')
B.lens_by_halos(relevance_lim=relevance_lim,lookup_table=True,smooth_corr=smooth_corr)
print 'Lightcones have {0:.2f} +/- {1:.2f} relevant galaxies'.format(B.mean_relevant_halos,B.std_relevant_halos)

# Calculate likelihood
likelihood = B.calculate_log_likelihood()
if time is True:
    elapsed_time = timeit.default_timer() - start_time
    print 'lens_by_halos likelihood = {}; Total time taken was {} s'.format(likelihood,elapsed_time)

if rel_compare is True:
    # Lens the background catalog using only relevant foreground halos
    if vb is True: print('Lensing background by relevant halos..')
    relevance_lim2 = 0.00001 # ~60 galaxies/lightcone
    B_rel.lens_by_halos(relevance_lim=relevance_lim2,lookup_table=True,smooth_corr=smooth_corr)
    mean_rel =  B_rel.mean_relevant_halos
    std_rel = B_rel.std_relevant_halos
    print 'Lightcones have {0:.2f} +/- {1:.2f} relevant galaxies'.format(mean_rel,std_rel)

if smooth_compare is True:
    # Lens the background catalog using smooth-component correction
    if vb is True: print('Lensing background with smooth-component correction..')
    B_smooth.lens_by_halos(relevance_lim=relevance_lim,lookup_table=True,smooth_corr=True)

# Plot 'relevance' distribution
if relevance_lim == 0.0 and rel_plots is True:
    mean_relevance = [np.mean(lightcone.galaxies['relevance']) for lightcone in B.lightcones]
    plt.subplot(2, 1, 1)
    plt.hist(mean_relevance,100,alpha=0.75,log=True)
    plt.xlabel('Mean Relevance per Lightcone (M=10^12 Sol Mass, R=10 kpc)',fontsize=16)
    plt.ylabel('Lightcone Count ({} total )'.format(B.galaxy_count),fontsize=16)
    #plt.show()

    plt.subplot(2, 1, 2)
    max_relevance = [np.max(lightcone.galaxies['relevance']) for lightcone in B.lightcones]
    plt.hist(max_relevance,100,alpha=0.75,log=True)
    plt.xlabel('Max Relevance per Lightcone (M=10^12 Sol Mass, R=10 kpc)',fontsize=16)
    plt.ylabel('Lightcone Count ({} total )'.format(B.galaxy_count),fontsize=16)
    plt.show()

if smooth_plots is True:
    lc = B.lightcones[np.random.randint(0,B.galaxy_count)]
    lc.plot_kappas()
    pangloss.plotting.plot_densities(F,lc,density_type='surface')
    pangloss.plotting.plot_densities(F,lc)

if corr_plots is True:
    # Calculate the correlation function for each lensing type
    if vb is True: print('Calculating correlation...')
    gg_map = B.calculate_corr(corr_type='gg',lensed='map',foreground=F)
    gg_halo = B.calculate_corr(corr_type='gg',lensed='halo',foreground=F)
    ng_map = B.calculate_corr(corr_type='ng',lensed='map',foreground=F)
    ng_halo = B.calculate_corr(corr_type='ng',lensed='halo',foreground=F)

    if rel_compare is True:
        # Calculate the correlation function using only most relevant halos
        gg_halo_r = B_rel.calculate_corr(corr_type='gg',lensed='halo',foreground=F)
        ng_halo_r = B_rel.calculate_corr(corr_type='ng',lensed='halo',foreground=F)

    if smooth_compare is True:
        # Calculate the correlation function using the smooth-component correction
        gg_halo_s = B_smooth.calculate_corr(corr_type='gg',lensed='halo',foreground=F)
        ng_halo_s = B_smooth.calculate_corr(corr_type='ng',lensed='halo',foreground=F)

    # Plot the correlation functions
    pangloss.plotting.plot_corr(gg_map,corr_type='gg',corr_comp='plus',lensed='map',color='green',galaxy_count=B.galaxy_count)
    pangloss.plotting.plot_corr(gg_halo,corr_type='gg',corr_comp='plus',lensed='halo',color='purple')
    pangloss.plotting.plot_corr(gg_map,corr_type='gg',corr_comp='cross',lensed='map',color='green')
    pangloss.plotting.plot_corr(gg_halo,corr_type='gg',corr_comp='cross',lensed='halo',color='purple')
    plt.gcf().set_size_inches(10,10)
    plt.show()

    # Calculate RMSE
    gg_nrmse = np.sqrt( np.sum( (gg_halo.xip - gg_map.xip)**2 ) / np.size(gg_map.xip) ) / np.sqrt( np.sum( gg_halo.xip**2 ) / np.size(gg_map.xip) )
    print 'The shear-shear correlation function NRMSE = {}'.format(gg_nrmse)

    pangloss.plotting.plot_corr(ng_map,corr_type='ng',corr_comp='real',lensed='map',color='green',galaxy_count=B.galaxy_count)
    pangloss.plotting.plot_corr(ng_halo,corr_type='ng',corr_comp='real',lensed='halo',color='purple')
    plt.gcf().set_size_inches(10,10)
    plt.show()

    # Calculate RMSE
    ng_nrmse = np.sqrt( np.sum( (ng_halo.xi - ng_map.xi)**2 ) / np.size(ng_map.xi) ) / np.sqrt( np.sum( ng_halo.xi**2 ) / np.size(ng_map.xi) )
    print 'The galaxy-mass correlation function NRMSE = {}'.format(ng_nrmse)

    # Compare the correlation functions
    chi2,n_sigma,percent_err,std_err = B.compare_corr(gg_halo,gg_map,corr_type='gg',corr_comp='plus')
    print 'Ellipticity-Ellipticity correlation difference intrinsic to mapped:','chi^2: ',chi2,'n_sigma: ',n_sigma,'percent_err: ',percent_err,'+\-',std_err

    chi2,n_sigma,percent_err,std_err = B.compare_corr(ng_halo,ng_map,corr_type='ng',corr_comp='real')
    print 'Galaxy-Galaxy correlation difference intrinsic to mapped:','chi^2: ',chi2,'n_sigma: ',n_sigma,'percent_err: ',percent_err,'+\-',std_err

    if rel_compare is True:
        # Compare the correlations for intrinsic, halos, and relevant halos
        if vb is True: print('Comparing correlations (rel)...')
        pangloss.plotting.compare_relevant_halos(gg_map,gg_halo,gg_halo_r,corr_type='gg',galaxy_count=B.galaxy_count,radius=lc_radius,rel_halos=[mean_rel,std_rel])
        pangloss.plotting.compare_relevant_halos(ng_map,ng_halo,ng_halo_r,corr_type='ng',galaxy_count=B.galaxy_count,radius=lc_radius,rel_halos=[mean_rel,std_rel])

    if smooth_compare is True:
        # Compare the correlations for intrinsic, halos, and halos with smooth-component correction
        if vb is True: print('Comparing correlations (smooth)...')
        if rel_compare is False:
            mean_rel =  None
            std_rel = None
        pangloss.plotting.compare_smooth_component(gg_map,gg_halo,gg_halo_s,corr_type='gg',galaxy_count=B.galaxy_count,radius=lc_radius,rel_halos=[mean_rel,std_rel])
        pangloss.plotting.compare_smooth_component(ng_map,ng_halo,ng_halo_s,corr_type='ng',galaxy_count=B.galaxy_count,radius=lc_radius,rel_halos=[mean_rel,std_rel])

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

if vb is True: print('Closing...')
