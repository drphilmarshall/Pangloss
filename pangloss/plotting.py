import numpy as np
import matplotlib.pyplot as plt
import os, sys
from astropy.table import Table, Column
from matplotlib.patches import Ellipse
from matplotlib.collections import LineCollection
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar

# Import Pangloss:
PANGLOSS_DIR = os.path.expandvars("$PANGLOSS_DIR")
sys.path.append(PANGLOSS_DIR)
import pangloss

viewport = [0.1,0.1,0.8,0.8]

# ============================================================================

"""
NAME
    plotting

PURPOSE
    A collection of functions related to coordinate conversions and plotting,
    useful for maps, catalogs, and lightcones.

COMMENTS

FUNCTIONS
    make_axes(fig,subplot,imsubplot):    Creates axes for plotting maps and catalogs.

    set_figure_size(fig,fig_size,Lx,Ly):    Changes figure size.

    plot_ellipse(ra,dec,size,mod,phi,axis,color,alpha):    Plot an ellipse at (ra,dec) with
                                                           given modulus and angle phi.

    plot_sticks(ra,dec,mod,phi,axis,color):    Plot ellipticity stick at (ra,dec) with given
                                               modulus, angle, and color.

    plot_corr(corr,corr_type,corr_comp,sep_units,lensed,color,
                            line_style,fig_size,M,galaxy_count):    Plot the correlation function component of
                                                                    the given correlation type.

    The following are used to create a demonstrative plot of how correlation functions describe
    the lensing of galaxies:
    calc_corr_demo(units):
    plot_corr_demo(corr,corr_type,units):
    plot_lensed_colors(B,subplot,center,fig_size):
    calc_corr_components(points):
    plot_corr_component(r,corr,corr_type,c):
    plot_corr_color_demo(N):

BUGS
    None

AUTHORS
  This file is part of the Pangloss project, distributed under the
  GPL v2, by Tom Collett (IoA) and  Phil Marshall (Oxford).
  Please cite: Collett et al 2013, http://arxiv.org/abs/1303.6564

HISTORY
  2015-07-5  Collected and extended by Everett (SLAC)
"""


def make_axes(fig,subplot,imsubplot=[0,1,0,1]):
    '''
    Creates axes for plotting maps and catalogs. Both sets of axes are
    contained in the current figure instance, and also returned for
    ease of use.
    '''

    '''
    # Works for shearmaps...
    # Make image axes for plotting maps in:
    imshow = fig.add_axes(viewport,label='imshow')
    imshow.set_xlim(imsubplot[0],imsubplot[1])
    imshow.set_ylim(imsubplot[2],imsubplot[3])
    imshow.get_xaxis().set_visible(False)
    imshow.get_yaxis().set_visible(False)

    # Make wcs axes for plotting catalogs in:
    world = fig.add_axes(viewport,label='world')
    world.set_xlabel('Right Ascension / deg')
    world.set_ylabel('Declination / deg')
    world.set_title('')
    world.set_xlim(subplot[0],subplot[1])
    world.set_ylim(subplot[2],subplot[3])
    '''


    # works for kappamaps...
    # Make wcs axes for plotting catalogs in:
    world = fig.add_axes(viewport,label='world')
    world.set_xlabel('Right Ascension / deg')
    world.set_ylabel('Declination / deg')
    world.set_title('')
    world.set_xlim(subplot[0],subplot[1])
    world.set_ylim(subplot[2],subplot[3])

    # Make image axes for plotting maps in:
    imshow = fig.add_axes(viewport,label='imshow')
    imshow.set_xlim(imsubplot[0],imsubplot[1])
    imshow.set_ylim(imsubplot[2],imsubplot[3])
    imshow.get_xaxis().set_visible(False)
    imshow.get_yaxis().set_visible(False)

    return imshow,world


def set_figure_size(fig,fig_size,Lx,Ly):
    '''
    Takes inputted figure and changes its size according to Lx and Ly given in
    inches.
    '''

    if Lx == Ly:
        fig.set_size_inches(fig_size,fig_size)
    elif Lx > Ly:
        fig.set_size_inches(fig_size,fig_size*(1.0*Ly/Lx))
    else:
        fig.set_size_inches(fig_size*(1.0*Lx/Ly),fig_size)
    return

# ----------------------------------------------------------------------------

def plot_ellipse(ra,dec,size,mod,phi,axis,color,alpha):
    '''
    Plot an ellipse centered at (ra,dec) with width 'size' and height<width dependent
    on 'q'. The ellipse is rotated by angle 'phi' which must be in deg. The ellipse
    has color set to the string 'color', and plots on the inputted axis.
    '''
    q = (1-mod)/(1+mod)
    ellipse = Ellipse(xy=[ra,dec],width=size,height=np.sqrt(q)*size,angle=phi)
    axis.add_artist(ellipse)
    ellipse.set_clip_box(axis.bbox)
    ellipse.set_alpha(alpha)
    ellipse.set_facecolor(color)

    return


def plot_sticks(ra,dec,mod,phi,axis,color):
    '''
    Plot an ellipticity stick at (ra,dec) with modulus 'mod' and orientation angle 'phi'.
    The stick is plotted on the inputted axis and with given 'color'.
    '''

    # Make sure that all of the columns are the same length
    assert np.size(ra) == np.size(dec) and np.size(dec) == np.size(phi)
    N = np.size(ra)

    # If N is 1, convert inputs into lists
    if N == 1:
        ra = np.array([ra])
        dec = np.array([dec])
        mod = np.array([mod])
        phi = np.array([phi])

    # Convert the angles phi to rad
    phi = np.deg2rad(phi)

    # Preallocation
    pt1 = np.zeros(N,dtype=tuple)
    pt2 = np.zeros(N,dtype=tuple)
    lines = np.zeros(N,dtype=tuple)

    # Set scale size of sticks
    xi, xf = axis.get_xlim()
    yi, yf = axis.get_ylim()
    Lx, Ly = abs(xf-xi), abs(yf-yi)
    L = np.mean([Lx,Ly])

    # Need this to see sticks (weak lensing)
    scale = 1.0
    size = scale*mod*L

    # For every object, create a line centered at (ra[i],dec[i]) with appropriate size
    # and orientation angle phi
    for i in range(N):
        pt1[i] = (ra[i]-size[i]*0.5*np.cos(phi[i]), dec[i]-size[i]*0.5*np.sin(phi[i]))
        pt2[i] = (ra[i]+size[i]*0.5*np.cos(phi[i]), dec[i]+size[i]*0.5*np.sin(phi[i]))
        #pt1[i] = (ra[i]-s*np.cos(phi[i]), dec[i]-s*np.sin(phi[i]))
        #pt2[i] = (ra[i]+s*np.cos(phi[i]), dec[i]+s*np.sin(phi[i]))
        lines[i] = tuple([pt1[i],pt2[i]])

    # Turn the array of lines into tuples
    lines = tuple(lines)

    # Create the collection of sticks from these lines, and add them to the inputted axis
    sticks = LineCollection(lines,linestyles='solid',color=color,lw=2)
    axis.add_collection(sticks)

    return

# ----------------------------------------------------------------------------

def plot_corr(corr,corr_type='gg',corr_comp='plus',sep_units='arcmin',lensed='map',color=None,line_style=None,fig_size=10,M=None,galaxy_count=None,radius=None,rel_halos=None):
    '''
    Plot the correlation component 'corr_comp' of type 'corr_type' in separation units of 'sep_units'.
    By default the method will plot the lensed-by-map correlation, but can also plot without lensing
    (lensed='none'), lensed-by-halo correlation (lensed='halo'), or with the lensed-by-halo correlation
    using only the most relevant halos (lensed='halo_rel').
    '''

    # If no color is inputted, plot in black
    if color is None:
        color = 'black'

    # Get current figure (or make one if it doesn't exist)
    fig = plt.gcf()

    if fig._label != 'Correlation':

        # Set the figure name
        fig._label = 'Correlation'

        # Set max and min dtheta
        min_sep = corr.min_sep
        max_sep = corr.max_sep
        #min_sep = np.deg2rad(0.1/60.0)
        #max_sep = min_sep * 100.0 # 10 arcmin
        if sep_units == 'arcmin':
            min_sep = np.rad2deg(min_sep)*60
            max_sep = np.rad2deg(max_sep)*60

        elif sep_units == 'deg':
            min_sep = np.rad2deg(min_sep)
            max_sep = np.rad2deg(max_sep)

        elif sep_units == 'rad':
            min_sep = min_sep
            max_sep = max_sep

        # Set figure size
        plt.gcf().set_size_inches(fig_size,fig_size)

        # Set plot settings
        plt.xscale('log')
        plt.gca().set_xlim(min_sep,max_sep)

        # Set axes labels
        plt.xlabel(r'$\Delta\theta$ (arcmin)',fontsize=20)

        if corr_type == 'gg': plt.ylabel(r'Ellipticity-Ellipticity Correlation $\xi(\Delta\theta)$',fontsize=20)
        elif corr_type == 'ng': plt.ylabel(r'Galaxy-Mass Correlation $\left<\gamma_t\right>$',fontsize=20)
        else: pass # Can incorporate other correlation functions later if needed

        # Plot xi=0 for reference
        plt.plot([min(np.exp(corr.logr)),max(np.exp(corr.logr))],[0,0],c='k',linestyle='dashed',linewidth=2)

        # Make axis ticks larger
        plt.gca().tick_params('both', length=5, width=1, which='major')
        plt.gca().tick_params('both', length=5, width=1, which='minor')

    if corr_type == 'gg':
        # For shear-shear (or ellipticity-ellipticity) correlation

        # Mark first label as intrinsic, ray-traced, or halo model
        if lensed == 'none': label1 = 'Intrinsic '
        elif lensed == 'map': label1 = 'Ray-Traced '
        elif lensed == 'halo': label1 = 'Halo Model '
        elif lensed == 'halo_rel': label1 = 'Most Relevant Halos '
        elif lensed == 'halo_smooth': label1 = 'Smooth-Component Correction '

        # Create the correct component label and linestyle to be used in legend
        if corr_comp == 'plus':
            # Plot xi_+ with a solid line
            correlation = corr.xip
            err = np.sqrt(corr.varxi)
            ls = '-'
            label2 = '+'

        elif corr_comp == 'minus':
            # Plot xi_- with a dotted line
            correlation = corr.xim
            err = np.sqrt(corr.varxi)
            ls = ':'
            label2 = '-'

        elif corr_comp == 'cross':
            # Plot xi_x with a dashed line
            correlation = 0.5*(corr.xim_im-corr.xip_im)
            err = 0.5*np.sqrt(2.0*corr.varxi)
            ls = '--'
            label2 = r'\times'

        elif corr_comp == 'cross_prime':
            # Plot xi_x prime with a dot-dash line
            correlation = corr.xip_im
            err = np.sqrt(corr.varxi)
            ls = '-.'
            label2 = r'\times^\prime'

        # If a linestyle is passed, overwrite default
        if not (line_style is None):
            ls = line_style

        # If the multiplicative error M is not None, add it to the legend label
        if M is not None:
            label2 += r',\,M={}'.format(M)

        # Plot the inputted shear-shear (or ellipticity-ellipticity) correlation function component
        plt.errorbar(np.exp(corr.logr),correlation,err,c=color,linestyle=ls,label=label1+r'$\xi_'+label2+'$',linewidth=2)

    elif corr_type == 'ng':
        # Plot the galaxy-mass correlation function (xi_gm)

        # Mark first label as intrinsic, observed, or predicted
        if lensed == 'none': label1 = 'Intrinsic '
        elif lensed == 'map': label1 = 'Ray-Traced '
        elif lensed == 'halo': label1 = 'Halo Model '
        elif lensed == 'halo_rel': label1 = 'Most Relevant Halos '
        elif lensed == 'halo_smooth': label1 = 'Smooth-Component Correction '

        if corr_comp == 'real':
            # Plot real(xi) with a solid line
            correlation = corr.xi
            ls = '-'
            label2 = r'$\operatorname{Re}$'

        elif corr_comp == 'imag':
            # Plot the imag(xi) with a dashed line
            correlation = corr.xi_im
            ls = '--'
            label2 = r'$\operatorname{Im}$'

        # Error is the same for either case
        err = np.sqrt(corr.varxi)

        # If a linestyle is passed, overwrite default
        if not (line_style is None):
            ls = line_style

        # If the multiplicative error M is not None, add it to the legend label
        if M is not None:
            label2 += r',\,M={}'.format(M)

        # Plot the inputted galaxy-galaxy correlation function component
        plt.errorbar(np.exp(corr.logr),correlation,err,c=color,linestyle=ls,label=label1+label2+r'$\left(\xi\right)$',linewidth=2)

    else:
        # Can incorporate other correlation functions here if needed
        pass

    # get legend handles
    handles, labels = plt.gca().get_legend_handles_labels()
    # remove the errorbars
    handles = [h[0] for h in handles]
    plt.legend(handles,labels,fontsize=18)

    # Print the number of lightcones and their radii
    if galaxy_count and radius is not None:
        # place a text box in upper left in axes coords
        plt.gca().text(0.3,0.95,'{} lightcones with radius {} arcmin'.format(galaxy_count,radius),transform=plt.gca().transAxes,fontsize=18,verticalalignment='top')
    elif galaxy_count is not None:
        # place a text box in upper left in axes coords
        plt.gca().text(0.3,0.95,'{} lightcones'.format(galaxy_count),transform=plt.gca().transAxes,fontsize=18,verticalalignment='top')

    # Print the number of relevant halos in each lightcone if applicable
    if lensed == 'halo_rel' and rel_halos is not None:
        # place a text box in upper left in axes coords
        mean_rel, std_rel = rel_halos[0], rel_halos[1]
        plt.gca().text(0.3,0.85,'{:.3}+/-{:.3} relevant halos per lightcone'.format(mean_rel,std_rel),transform=plt.gca().transAxes,fontsize=18,verticalalignment='top')

    return

def compare_relevant_halos(corr_map,corr_halo,corr_rel,corr_type='gg',sep_units='arcmin',galaxy_count=None,radius=None,rel_halos=None):
    '''
    Plots a comparrison of the ray-traced/map correlation, halo model correlation, and relevant-halos correlation.
    '''

    # Set max and min dtheta
    min_sep = corr_map.min_sep
    max_sep = corr_map.max_sep

    if sep_units == 'arcmin':
        min_sep = np.rad2deg(min_sep)*60
        max_sep = np.rad2deg(max_sep)*60

    elif sep_units == 'deg':
        min_sep = np.rad2deg(min_sep)
        max_sep = np.rad2deg(max_sep)

    elif sep_units == 'rad':
        min_sep = min_sep
        max_sep = max_sep

    dtheta = np.exp(corr_map.logr)

    if corr_type == 'gg':
        # Calculate percent error of relevant halos vs all halos
        gg_percent_err = np.abs((corr_rel.xip-corr_halo.xip)/corr_halo.xip)*100
        # Reject outliers (error can get very large when corr is near 0)
        gg_percent_err, indices = reject_outliers(gg_percent_err)
        dtheta = dtheta[indices]
        # Calculate mean/std percent error of relevant halos vs all halos
        mean_gg_err = np.mean(gg_percent_err)
        std_gg_err = np.std(gg_percent_err)
        # Calculate mean Fraction for gg
        gg_frac = corr_rel.xip/corr_halo.xip
        mean_gg_frac = np.mean(gg_frac)
        std_gg_frac = np.std(gg_frac)
        # Calculate RMSE
        nrmse = np.sqrt( np.sum( (corr_rel.xip - corr_halo.xip)**2 ) / np.size(corr_halo.xip) ) / np.sqrt( np.sum( corr_halo.xip**2 ) / np.size(corr_halo.xip) )
        # Set correlation component for plotting
        corr_comp = 'plus'

    elif corr_type == 'ng':
        # Calculate percent error of relevant halos vs all halos
        ng_percent_err = np.abs((corr_rel.xi-corr_halo.xi)/corr_halo.xi)*100
        # Reject outliers (error can get very large when corr is near 0)
        ng_percent_err, indices = reject_outliers(ng_percent_err)
        dtheta = dtheta[indices]
        # Calculate mean/std percent error of relevant halos vs all halos
        mean_ng_err = np.mean(ng_percent_err)
        std_ng_err = np.std(ng_percent_err)
        # Calculate mean Fraction for ng
        ng_frac = corr_rel.xi/corr_halo.xi
        mean_ng_frac = np.mean(ng_frac)
        std_ng_frac = np.std(ng_frac)
        # Calculate RMSE
        nrmse = np.sqrt( np.sum( (corr_rel.xi - corr_halo.xi)**2 ) / np.size(corr_halo.xi) ) / np.sqrt( np.sum( corr_halo.xi**2 ) / np.size(corr_halo.xi) )
        # Set correlation component for plotting
        corr_comp = 'real'

    # Plot the correlation functions
    plt.subplot(2,1,1)
    plot_corr(corr_map,corr_type=corr_type,corr_comp=corr_comp,lensed='map',color='green',galaxy_count=galaxy_count,radius=radius)
    plot_corr(corr_halo,corr_type=corr_type,corr_comp=corr_comp,lensed='halo',color='purple')
    plot_corr(corr_rel,corr_type=corr_type,corr_comp=corr_comp,lensed='halo_rel',color='blue',rel_halos=rel_halos)

    # Plot the percent error of all halos vs relevant halos at each separation distance
    plt.subplot(2,1,2)
    if corr_type == 'gg':
        plt.plot(dtheta,gg_percent_err,'ob',label='Individual Error')
        plt.plot([min_sep, max_sep],[mean_gg_err, mean_gg_err],'--b',label='Mean Error: {:.3}%'.format(mean_gg_err),linewidth=2)
    if corr_type == 'ng':
        plt.plot(dtheta,ng_percent_err,'ob',label='Individual Error')
        plt.plot([min_sep, max_sep],[mean_ng_err, mean_ng_err],'--b',label='Mean Error: {:.3}%'.format(mean_ng_err),linewidth=2)
    # Set plot settings
    plt.xscale('log')
    plt.gca().set_xlim(min_sep,max_sep)

    # Make axis ticks larger
    plt.gca().tick_params('both', length=5, width=1, which='major')
    plt.gca().tick_params('both', length=5, width=1, which='minor')

    # Set axes labels
    plt.xlabel(r'$\Delta\theta$ (arcmin)',fontsize=20)
    plt.ylabel('Percent Error',fontsize=20)
    plt.legend(fontsize=18)

    # State the NRMSE
    plt.gca().text(0.025,0.95,'Normalized RMSE: {:.3}'.format(nrmse),transform=plt.gca().transAxes,fontsize=18,verticalalignment='top')

    plt.gcf().set_size_inches(10,10)
    plt.show()

    return

def compare_smooth_component(corr_map,corr_halo,corr_smooth,corr_type='gg',sep_units='arcmin',galaxy_count=None,radius=None,rel_halos=None):
    '''
    Plots a comparrison of the ray-traced/map correlation, halo model correlation, and halo model with smooth-component correction correlation.
    '''

    # Set max and min dtheta
    min_sep = corr_map.min_sep
    max_sep = corr_map.max_sep

    if sep_units == 'arcmin':
        min_sep = np.rad2deg(min_sep)*60
        max_sep = np.rad2deg(max_sep)*60

    elif sep_units == 'deg':
        min_sep = np.rad2deg(min_sep)
        max_sep = np.rad2deg(max_sep)

    elif sep_units == 'rad':
        min_sep = min_sep
        max_sep = max_sep

    dtheta = np.exp(corr_map.logr)

    if corr_type == 'gg':
        # Calculate percent error of relevant halos vs all halos
        gg_percent_err = np.abs((corr_smooth.xip-corr_halo.xip)/corr_halo.xip)*100
        # Reject outliers (error can get very large when corr is near 0)
        gg_percent_err, indices = reject_outliers(gg_percent_err)
        dtheta = dtheta[indices]
        # Calculate mean/std percent error of relevant halos vs all halos
        mean_gg_err = np.mean(gg_percent_err)
        std_gg_err = np.std(gg_percent_err)
        # Calculate mean Fraction for gg
        gg_frac = corr_smooth.xip/corr_halo.xip
        mean_gg_frac = np.mean(gg_frac)
        std_gg_frac = np.std(gg_frac)
        # Calculate RMSE
        nrmse = np.sqrt( np.sum( (corr_smooth.xip - corr_halo.xip)**2 ) / np.size(corr_halo.xip) ) / np.sqrt( np.sum( corr_halo.xip**2 ) / np.size(corr_halo.xip) )
        # Set correlation component for plotting
        corr_comp = 'plus'

    elif corr_type == 'ng':
        # Calculate percent error of relevant halos vs all halos
        ng_percent_err = np.abs((corr_smooth.xi-corr_halo.xi)/corr_halo.xi)*100
        # Reject outliers (error can get very large when corr is near 0)
        ng_percent_err, indices = reject_outliers(ng_percent_err)
        dtheta = dtheta[indices]
        # Calculate mean/std percent error of relevant halos vs all halos
        mean_ng_err = np.mean(ng_percent_err)
        std_ng_err = np.std(ng_percent_err)
        # Calculate mean Fraction for ng
        ng_frac = corr_smooth.xi/corr_halo.xi
        mean_ng_frac = np.mean(ng_frac)
        std_ng_frac = np.std(ng_frac)
        # Calculate RMSE
        nrmse = np.sqrt( np.sum( (corr_smooth.xi - corr_halo.xi)**2 ) / np.size(corr_halo.xi) ) / np.sqrt( np.sum( corr_halo.xi**2 ) / np.size(corr_halo.xi) )
        # Set correlation component for plotting
        corr_comp = 'real'

    # Plot the correlation functions
    #plt.subplot(2,1,1)
    plot_corr(corr_map,corr_type=corr_type,corr_comp=corr_comp,lensed='map',color='green',galaxy_count=galaxy_count,radius=radius)
    plot_corr(corr_halo,corr_type=corr_type,corr_comp=corr_comp,lensed='halo',color='purple')
    plot_corr(corr_smooth,corr_type=corr_type,corr_comp=corr_comp,lensed='halo_smooth',color='blue',rel_halos=rel_halos)

    '''
    **Eventually change this to plot the percent change in correlation in respect to existing difference between maps
    and halos without correction**
    # Plot the percent error of all halos vs relevant halos at each separation distance
    plt.subplot(2,1,2)
    if corr_type == 'gg':
        plt.plot(dtheta,gg_percent_err,'ob',label='Individual Error')
        plt.plot([min_sep, max_sep],[mean_gg_err, mean_gg_err],'--b',label='Mean Error: {:.3}%'.format(mean_gg_err),linewidth=2)
    if corr_type == 'ng':
        plt.plot(dtheta,ng_percent_err,'ob',label='Individual Error')
        plt.plot([min_sep, max_sep],[mean_ng_err, mean_ng_err],'--b',label='Mean Error: {:.3}%'.format(mean_ng_err),linewidth=2)
    # Set plot settings
    plt.xscale('log')
    plt.gca().set_xlim(min_sep,max_sep)

    # Make axis ticks larger
    plt.gca().tick_params('both', length=5, width=1, which='major')
    plt.gca().tick_params('both', length=5, width=1, which='minor')

    # Set axes labels
    plt.xlabel(r'$\Delta\theta$ (arcmin)',fontsize=20)
    plt.ylabel('Percent Error',fontsize=20)
    plt.legend(fontsize=18)

    # State the NRMSE
    plt.gca().text(0.025,0.95,'Normalized RMSE: {:.3}'.format(nrmse),transform=plt.gca().transAxes,fontsize=18,verticalalignment='top')
    '''

    plt.gcf().set_size_inches(10,10)
    plt.show()

    return


def reject_outliers(data, m=5):
    '''
    Used in compare_relevant_halos() to reject percent error outliers that are
    very high due to the correlation values being very near zero. Returns non-outliers
    and their indices
    '''
    return data[abs(data - np.mean(data)) < m * np.std(data)], abs(data - np.mean(data)) < m * np.std(data)

def plot_densities(foreground,lightcone,density_type='volume'):
    '''
    Plots the universe's mass density, halo mass density, and smooth-component mass density
    as a function of redshift.
    '''

    # Use the cosmology defined by the lightcone
    cosmo = lightcone.cosmo
    redshifts = lightcone.redshifts
    N = np.size(redshifts)

    D = pangloss.Distance(cosmo)

    rho_m = np.zeros(N)
    rho_h = np.zeros(N)
    rho_s = np.zeros(N)
    sigma_m = np.zeros(N)
    sigma_h = np.zeros(N)
    sigma_s = np.zeros(N)

    # Calculate the universe's mean mass density at each redshift slice
    for i in range(N):
        z = redshifts[i]
        rho_m[i] = D.mean_mass_density(z) # M_sol/Mpc^3
        rho_h[i] = lightcone.halo_density(z) # M_sol/Mpc^3
        rho_s[i] = rho_m[i] - rho_h[i] # M_sol/Mpc^3

        if density_type == 'surface':
            sigma_m[i] = lightcone.slice_sigma(rho_m[i],z) # M_sol/Mpc^2
            sigma_h[i] = lightcone.slice_sigma(rho_h[i],z) # M_sol/Mpc^2
            sigma_s[i] = lightcone.slice_sigma(rho_s[i],z) # M_sol/Mpc^2

    # Set up density axis
    fig, ax_rho = plt.subplots()
    plt.xlabel('Redshift (z)',fontsize='14')
    plt.xticks(fontsize=14)

    # Plot densities
    if density_type == 'volume':
        plt.plot(redshifts,rho_m,'--k',linewidth=2,label='Universe Mean Density',zorder=10)
        plt.plot(redshifts,rho_h,'--og',label='Foreground Halos',zorder=10)
        plt.plot(redshifts,rho_s,'--ob',label='Smooth-Component',zorder=10)
        degree = 3
    elif density_type == 'surface':
        # Critical surface density
        sigma_crit = lightcone.sigma_crit
        plt.plot(redshifts,sigma_m,'--k',linewidth=2,label='Universe Mean Density',zorder=10)
        plt.plot(redshifts,sigma_crit,'--r',linewidth=2,label='Sigma_crit',zorder=10)
        plt.plot(redshifts,sigma_s,'--ob',label='Sigma_s',zorder=10)
        plt.plot(redshifts,sigma_h,'--og',label='Sigma_h',zorder=10)
        degree = 2

    plt.yticks(fontsize=14)
    plt.yscale('log')
    ax_rho.set_ylabel(r'Mass Density ($M_\odot$/Mpc$^{}$)'.format(degree),fontsize='16')
    plt.legend(loc=2)

    # Set up histogram axis
    ax_hist = ax_rho.twinx()

    # Plot histogram
    plt.hist(lightcone.galaxies['z_sz'],50,normed=False,facecolor='green', alpha=0.1,zorder=1)
    plt.yscale('log')
    ax_hist.set_ylabel('Galaxy Count', color='g',fontsize='16')
    plt.yticks(fontsize=14)

    # Match axis label colors to plot colors
    #for tl in ax_rho.get_yticklabels(): tl.set_color('b')
    for tl in ax_hist.get_yticklabels():  tl.set_color('g')

    plt.gca().text(0.25,0.95,'Lightcone with radius {} arcmin and {} galaxies'.format(lightcone.rmax,lightcone.galaxy_count),transform=plt.gca().transAxes,fontsize=18,verticalalignment='top')

    plt.show()

    return


#-------------------------------------------------------------------------------------------------------------------
# This section is for code only used to create demo plots. None of these are currently used for actual pangloss use.

# ----------------------------------------------------------------------------
# Create phase plots for correlation coefficients. Used in VisualizingCorrelationFunction.ipyn.
def calc_corr_demo(units='deg'):
    '''
    Calculates the plus, minus, cross, and cross-prime components of the ellipticity-ellipticity
    correlation function for ellipticity sticks with a range of orientations from 0 to 180 degrees.
    This is only used for demonstrations to show the connection between the difference in galaxy
    orientation and the different correlation function values.
    '''

    if units == 'deg':
        # Colormap locations
        i = np.arange(-20,135+25,.1)
        j = np.arange(-20,135+25,.1)

        # Create grid
        I,J = np.meshgrid(i,j)

        # Calculate shear components
        gt_i = -np.cos(np.deg2rad(2*I))
        gt_j = -np.cos(np.deg2rad(2*(I+J)))
        gx_i = -np.sin(np.deg2rad(2*I))
        gx_j = -np.sin(np.deg2rad(2*(I+J)))

    elif units == 'rad':
        # Colormap locations
        i = np.arange(-0.35,np.pi+0.35,.01)
        j = np.arange(-0.35,np.pi+0.35,.01)

        # Create grid
        I,J = np.meshgrid(i,j)

        # Calculate shear components
        gt_i = -np.cos(2*I)
        gt_j = -np.cos(2*(I+J))
        gx_i = -np.sin(2*I)
        gx_j = -np.sin(2*(I+J))

    # Calculate correlations
    gtt = gt_i*gt_j
    gxx = gx_i*gx_j
    gtx = gt_i*gx_j
    gxt = gx_i*gt_j

    # Calculate xi_plus, xi_minus, xi_cross, and xi_cross'
    xi_p = gtt+gxx
    xi_m = gtt-gxx
    xi_x = gtx
    xi_xp = gtx-gxt

    return xi_p,xi_m,xi_x,xi_xp


def plot_corr_demo(corr,corr_type='plus',units='deg'):
    '''
    Plots the inputted correlation function type as a color map with the ellipticity sticks at their
    corresponding orientations overlaid. This is only used for demonstrations to show the connection
    between the difference in galaxy orientation and the different correlation function values.
    '''

    if units == 'deg':
        # Colormap locations
        i = np.arange(-20,135+25,.1)
        j = np.arange(-20,135+25,.1)
        # Stick locations
        x = np.arange(0,180,45)
        y = np.arange(0,180,45)

        # Create grid
        I,J = np.meshgrid(i,j)
        X,Y = np.meshgrid(x,y)

        # Stick offset
        dx = 7

    elif units == 'rad':
        # Colormap locations
        i = np.arange(-0.35,np.pi+0.35,.01)
        j = np.arange(-0.35,np.pi+0.35,.01)
        # Stick locations
        x = np.arange(0,np.pi,np.pi/4)
        y = np.arange(0,np.pi,np.pi/4)

        # Create grid
        I,J = np.meshgrid(i,j)
        X,Y = np.meshgrid(x,y)

        # Stick offset
        dx = 0.125

    # Create plot
    plt.imshow(corr,origin='lower',cmap='bwr',extent=(np.min(I), np.max(I), np.min(J), np.max(J)))
    plt.xlabel(r'$\theta_i\,(degree)$',fontsize=20)
    plt.ylabel(r'$\Delta\theta\,(|\theta_i-\theta_j|)\,(degree)$',fontsize=20)
    plt.gcf().set_size_inches(10,10)
    plt.colorbar()
    plt.grid()

    # Set title
    if corr_type == 'plus': plt.title(r'$\xi_+$',fontsize=25)
    elif corr_type == 'minus': plt.title(r'$\xi_-$',fontsize=25)
    elif corr_type == 'cross': plt.title(r'$\xi_\times$',fontsize=25)
    elif corr_type == 'crossp': plt.title(r'$\xi^\prime_\times$',fontsize=25)

    # Get axes and set stick modulus
    ax = plt.gca()
    mod = np.ones(len(I[0,:]))/20.0

    # Create sticks at dtheta = 0 deg:
    pangloss.plotting.plot_sticks(X[0,:]-dx,Y[0,:],mod,X[0,:],ax,'k')
    pangloss.plotting.plot_sticks(X[0,:]+dx,Y[0,:],mod,X[0,:]+Y[0,:],ax,'k')
    #plt.scatter(X[0,:]-dx,Y[0,:],c='k')
    #plt.scatter(X[0,:]+dx,Y[0,:],c='k')
    plt.scatter(X[0,:],Y[0,:],marker='+',s=30,c='k')

    # Create sticks at dtheta = 45 deg:
    pangloss.plotting.plot_sticks(X[1,:]-dx,Y[1,:],mod,X[1,:],ax,'k')
    pangloss.plotting.plot_sticks(X[1,:]+dx,Y[1,:],mod,X[1,:]+Y[1,:],ax,'k')
    #plt.scatter(X[1,:]-dx,Y[1,:],c='k')
    #plt.scatter(X[1,:]+dx,Y[1,:],c='k')
    plt.scatter(X[1,:],Y[1,:],marker='+',s=30,c='k')

    # Create sticks at dtheta = 90 deg:
    pangloss.plotting.plot_sticks(X[2,:]-dx,Y[2,:],mod,X[2,:],ax,'k')
    pangloss.plotting.plot_sticks(X[2,:]+dx,Y[2,:],mod,X[2,:]+Y[2,:],ax,'k')
    #plt.scatter(X[2,:]-dx,Y[2,:],c='k')
    #plt.scatter(X[2,:]+dx,Y[2,:],c='k')
    plt.scatter(X[2,:],Y[2,:],marker='+',s=30,c='k')

    # Create sticks at dtheta = 135 deg:
    pangloss.plotting.plot_sticks(X[3,:]-dx,Y[3,:],mod,X[3,:],ax,'k')
    pangloss.plotting.plot_sticks(X[3,:]+dx,Y[3,:],mod,X[3,:]+Y[3,:],ax,'k')
    #plt.scatter(X[3,:]-dx,Y[3,:],c='k')
    #plt.scatter(X[3,:]+dx,Y[3,:],c='k')
    plt.scatter(X[3,:],Y[3,:],marker='+',s=30,c='k')

    # Set axis limits
    plt.xlim([np.min(I),np.max(I)])
    plt.ylim([np.min(J),np.max(J)])

    return

# ----------------------------------------------------------------------------
# Plot color-coded galaxies around a lens and a scatter plot of the correlation coefficients using
# the same color scheme. Used in VisualizingCorrelationFunction.ipyn.

def plot_lensed_colors(B,subplot=[1.1175,1.0925,-1.54,-1.5175],center=[1.10425,-1.52875],fig_size=10):
    '''
    Plot the background galaxies in the catalog B contained in the inputted subplot in different colors
    as a function of their distance from the inputted center.
    '''

    # Create new figure with correct axes
    fig = plt.gcf()
    pangloss.plotting.make_axes(fig,subplot=[1.1175,1.0925,-1.54,-1.5175])
    ax = plt.gca()

    # Extract data from galaxies contained in subplot
    galaxies = B.return_galaxies(ra_lim=[subplot[0],subplot[1]],dec_lim=[subplot[2],subplot[3]])
    ra = np.rad2deg(galaxies['RA'])
    dec = np.rad2deg(galaxies['Dec'])
    mod = galaxies['eMod']
    phi = galaxies['ePhi']

    # Set center of lens by hand
    ra_c = center[0]
    dec_c = center[1]

    # Distance bin cutoffs
    c1 = 0.0025
    c2 = 0.005
    c3 = 0.0075
    c4 = 0.01
    c5 = 0.0125

    for i in range(np.size(ra)):
        # Don't plot if galaxy is strongly lensed
        if galaxies['strong_flag'][i] == 1: continue

        # Calculate galaxy distance to center
        dra_c = abs(ra[i]-ra_c)
        ddec_c = abs(dec[i]-dec_c)
        r_c = np.sqrt(dra_c**2+ddec_c**2)

        # Set galaxy color based upon distance to center
        if (r_c < c1):  color = 'red'
        elif (r_c >= c1) and (r_c < c2): color = 'orange'
        elif (r_c >= c2) and (r_c < c3): color = 'yellow'
        elif (r_c >= c3) and (r_c < c4): color = 'green'
        elif (r_c >= c4) and (r_c < c5): color = 'blue'
        else: color = 'purple'

        # Plot galaxy ellipticity
        #plt.subplot(1,5,1)
        pangloss.plot_sticks(ra[i],dec[i],mod[i],-phi[i],ax,color=color)

    # Plot lens center
    #plt.subplot(1,5,1)
    plt.scatter(1.10425,-1.52875,s=200)

    # Set figure size
    fig.set_size_inches(fig_size,fig_size)

    # Add scale bar
    Lx, Ly = abs(subplot[0]-subplot[1]), abs(subplot[2]-subplot[3])
    L = np.mean([Lx,Ly])
    bar = AnchoredSizeBar(ax.transData,L/10.0,'10% Ellipticity',pad=0.5,loc=3,sep=5,borderpad=0.25,frameon=True)
    bar.size_bar._children[0]._linewidth = 2
    #bar.size_bar._children[0]._edgecolor = (1,0,0,1)
    ax.add_artist(bar)

    plt.show()

    return


def calc_corr_components(points):
    '''
    Calculate the plus, minus, cross, and cross^prime components of the correlation between each galaxy pair.
    '''

    del_xi_p = []
    del_xi_m = []
    del_xi_x = []
    del_xi_xp = []

    r = []
    c = []

    # Center of lens
    ra_c = 1.10425
    dec_c = -1.52875

    # Distance bin cutoffs
    c1 = 0.0025
    c2 = 0.005
    c3 = 0.0075
    c4 = 0.01
    c5 = 0.0125

    for point1 in points:
        for point2 in points:
            if point1 == point2:
                # Don't calculate correlation for a point with itself
                continue
            else:
                # Extract galaxy locations and ellipticities
                ra1 = point1[0]
                ra2 = point2[0]
                dec1 = point1[1]
                dec2 = point2[1]
                el1 = point1[2]
                el2 = point2[2]

                # Calculate separation distance and angle
                dra = ra2-ra1
                ddec = dec2-dec1
                r.append(np.sqrt(dra**2+ddec**2))
                phi = np.arctan2(ddec,dra)

                # Calculate separation distance of each object to center of lens
                dra1_c = abs(ra1-ra_c)
                dra2_c = abs(ra2-ra_c)
                ddec1_c = abs(dec1-dec_c)
                ddec2_c = abs(dec2-dec_c)

                r1_c = np.sqrt(dra1_c**2+ddec1_c**2)
                r2_c = np.sqrt(dra2_c**2+ddec2_c**2)

                # Determine color of correlation scatter point based upon distance from lens
                if (r1_c < c1) and (r2_c < c1): c.append((1,0,0,1))
                elif (r1_c >= c1) and (r1_c < c2) and (r2_c >= c1) and (r2_c < c2): c.append((1,0.5,0,0.8))
                elif (r1_c >= c2) and (r1_c < c3) and (r2_c >= c2) and (r2_c < c3): c.append((1,1,0,0.6))
                elif (r1_c >= c3) and (r1_c < c4) and (r2_c >= c3) and (r2_c < c4): c.append((0,1,0,0.4))
                elif (r1_c >= c4) and (r1_c < c5) and (r2_c >= c4) and (r2_c < c5): c.append((0,0,1,0.2))
                else: c.append((0.5,0,0.5,0.01))

                # Calculate shear components
                g1t = -(el1*np.e**(-2j*phi)).real
                g1x = -(el1*np.e**(-2j*phi)).imag

                g2t = -(el2*np.e**(-2j*phi)).real
                g2x = -(el2*np.e**(-2j*phi)).imag

                # Calculate correlation components
                del_tt = g1t*g2t
                del_xx = g1x*g2x
                del_tx = g1t*g2x
                del_xt = g1x*g2t

                # Calculate correlation
                del_xi_p.append(del_tt+del_xx)
                del_xi_m.append(del_tt-del_xx)
                del_xi_x.append(del_tx)
                del_xi_xp.append(del_xt)

    return r,del_xi_p,del_xi_m,del_xi_x,del_xi_xp,c


def plot_corr_component(r,corr,corr_type,c):
    '''
    Plot a single correlation component.
    '''

    # Create another subplot
    #plt.subplot(1,5,i)

    plt.scatter(np.array(r)*60.0,corr,color=c)
    plt.plot([min(60.0*np.array(r)),max(60.0*np.array(r))],[0,0],'--k')
    #plt.plot([0.05,2.0],[0,0],'--k')
    plt.gcf().set_size_inches(10,10)
    plt.xlim([0.05,2.0])
    plt.ylim([-0.05,0.05])
    plt.xscale('log')

    plt.xlabel(r'$\Delta\theta$ (arcmin)',fontsize=20)

    if corr_type == 'plus':
        plt.ylabel(r'$\delta\xi_+(\theta)$',fontsize=20)
    elif corr_type == 'minus':
        plt.ylabel(r'$\delta\xi_-(\theta)$',fontsize=20)
    elif corr_type == 'cross':
        plt.ylabel(r'$\delta\xi_\times(\theta)$',fontsize=20)
    elif corr_type == 'cross_prime':
        plt.ylabel(r'$\delta\xi_\times^\prime(\theta)$',fontsize=20)


    # Make axis ticks larger
    plt.gca().tick_params('both', length=5, width=1, which='major')
    plt.gca().tick_params('both', length=5, width=1, which='minor')

    plt.show()

    return


def plot_corr_color_demo(N=200):
    '''
    Runs the correlation color demo using the plot_lens_colors(), calc_corr_components(),
    and plot_corr_component() methods.
    '''

    # Load in the (0,0) Kappamap and Shearmap
    K = pangloss.Kappamap(PANGLOSS_DIR+'/data/GGL_los_8_0_0_N_4096_ang_4_rays_to_plane_37_f.kappa',FITS=False)
    S = pangloss.Shearmap([PANGLOSS_DIR+'/data/GGL_los_8_0_0_N_4096_ang_4_rays_to_plane_37_f.gamma_1',PANGLOSS_DIR+'/data/GGL_los_8_0_0_N_4096_ang_4_rays_to_plane_37_f.gamma_2'],FITS=False)

    # Create a background catalog and lens the galaxies
    B = pangloss.BackgroundCatalog(N=N,sigma_e=0.01,domain=[1.1175,1.0925,-1.54,-1.515],field=[0,0,0,0]) # High ellipticity used to highlight the intrinsic shape of background sources
    B.lens_by_map(K,S)

    # Plot background catalog with color scheme based upon galaxy distance from lens
    plot_lensed_colors(B)

    # Extract galaxy data
    ra = np.rad2deg(B.galaxies['RA'])
    dec = np.rad2deg(B.galaxies['Dec'])
    e1 = B.galaxies['e1']
    e2 = B.galaxies['e2']
    el = e1+1.0j*e2

    #
    points = []
    # Create set of galaxy points and ellipticities
    for i in range(len(ra)):
        # Don't calculate correlation for strongly-lensed sources
        if B.galaxies['strong_flag'][i] == 1: continue

        r = ra[i]
        d = dec[i]
        e = el[i]
        points.append((r,d,e))

    # Calculate the correlation components between galaxies
    r,del_xi_p,del_xi_m,del_xi_x,del_xi_xp,c = calc_corr_components(points)

    # Calculate the shear-shear correlation function
    #gg=B.calculate_corr()

    #
    plot_corr_component(r,del_xi_p,'plus',c)
    plot_corr_component(r,del_xi_m,'minus',c)
    plot_corr_component(r,del_xi_x,'cross',c)
    plot_corr_component(r,del_xi_xp,'cross_prime',c)

    # Set figure size for figure containing all subplots
    #plt.gcf().set_size_inches(50,10)

    return r,del_xi_p,del_xi_m,del_xi_x,del_xi_xp,c

#-------------------------------------------------------------------------------------------------------------------
