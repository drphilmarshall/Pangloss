import numpy as np
import matplotlib.pyplot as plt
import os, sys
from astropy.table import Table, Column
from matplotlib.patches import Ellipse
from matplotlib.collections import LineCollection
import pangloss

viewport = [0.1,0.1,0.8,0.8]

# ============================================================================

"""
NAME
    plot

PURPOSE
    A collection of functions related to coordinate conversions and plotting,
    useful for maps, catalogs, and lightcones.

COMMENTS

FUNCTIONS
    image2physical(self,i,j,mapfile=0):
    physical2image(self,x,y,mapfile=0):
    image2world(self,i,j,mapfile=0):
    world2image(self,a,d,mapfile=0):
    physical2world(self,x,y,mapfile=0):
    world2physical(self,a,d,mapfile=0):

BUGS
    ???

AUTHORS
  This file is part of the Pangloss project, distributed under the
  GPL v2, by Tom Collett (IoA) and  Phil Marshall (Oxford).
  Please cite: Collett et al 2013, http://arxiv.org/abs/1303.6564

HISTORY
  2015-07-5  Collected and extended by Everett (SLAC)
"""

# ----------------------------------------------------------------------------
# Conversions between image, physical, and world coordinate systems (world conversions are approximate)

'''
NOTE: these have been taken from WLMap and have not been modified to work more generally. This should be done
once its scope has been understood.
'''
'''
def image2physical(self,i,j,mapfile=0):
    x = (i - self.wcs[mapfile]['LTV1'])/self.wcs[mapfile]['LTM1_1'] # x in rad
    y = (j - self.wcs[mapfile]['LTV2'])/self.wcs[mapfile]['LTM2_2'] # y in rad
    return x,y

def physical2image(self,x,y,mapfile=0):
    i = self.wcs[mapfile]['LTV1'] + self.wcs[mapfile]['LTM1_1']*x # x in rad
    j = self.wcs[mapfile]['LTV2'] + self.wcs[mapfile]['LTM2_2']*y # y in rad
    return i,j

 # Only approximate WCS transformations - assumes dec=0.0 and small field
def image2world(self,i,j,mapfile=0):
    a = self.wcs[mapfile]['CRVAL1'] + self.wcs[mapfile]['CD1_1']*(i - self.wcs[mapfile]['CRPIX1'])
    #if a < 0.0: a += 360.0 :We are using nRA instead now
    d = self.wcs[mapfile]['CRVAL2'] + self.wcs[mapfile]['CD2_2']*(j - self.wcs[mapfile]['CRPIX2'])
    return a,d

def world2image(self,a,d,mapfile=0):
    i = (a - self.wcs[mapfile]['CRVAL1'])/self.wcs[mapfile]['CD1_1'] + self.wcs[mapfile]['CRPIX1']
    # if a negative pixel is returned for i, reinput a as a negative degree
    if i<0:
        a-=360
        i = (a - self.wcs[mapfile]['CRVAL1'])/self.wcs[mapfile]['CD1_1'] + self.wcs[mapfile]['CRPIX1']
    j = (d - self.wcs[mapfile]['CRVAL2'])/self.wcs[mapfile]['CD2_2'] + self.wcs[mapfile]['CRPIX2']
    return i,j

def physical2world(self,x,y,mapfile=0):
    a = -np.rad2deg(x) - self.field_x[mapfile]*self.field[mapfile]
    #if a < 0.0: a += 360.0 :we are using nRA instead now
    d = np.rad2deg(y) + self.field_y[mapfile]*self.field[mapfile]
    return a,d

def world2physical(self,a,d,mapfile=0):
    x = -np.deg2rad(a + self.field_x[mapfile]*self.field[mapfile])
    y = np.deg2rad(d - self.field_y[mapfile]*self.field[mapfile])
    return x,y
'''
# ----------------------------------------------------------------------------
#

def make_axes(fig,subplot,imsubplot=[0,1,0,1]):
    '''
    Creates axes for plotting maps and catalogs. Both sets of axes are
    contained in the current figure instance, and also returned for
    ease of use.
    '''
    
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

    return imshow,world


def set_figure_size(fig,fig_size,Lx,Ly):
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
    Write docstring
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

def plot_corr(corr,corr_type='gg',sep_units='arcmin',lensed=True,fig_size=10):
    '''
    '''
    
    # Get current figure (or make one if it doesn't exist)
    fig = plt.gcf()
    
    if fig._label != 'Correlation':
    
        # Set the figure name
        fig._label = 'Correlation'
    
        # Set max and min dtheta
        print 'sep_units: ',corr.sep_units
        if sep_units == 'arcmin':
            min_sep = np.rad2deg(corr.min_sep)*60
            max_sep = np.rad2deg(corr.max_sep)*60
            
        elif sep_units == 'deg':
            min_sep = np.rad2deg(corr.min_sep)
            max_sep = np.rad2deg(corr.max_sep)
        
        elif sep_units == 'rad':
            min_sep = corr.min_sep
            max_sep = corr.max_sep
            
        # Set figure size
        plt.gcf().set_size_inches(fig_size,fig_size)
            
        # Set plot settings
        plt.xscale('log')
        plt.gca().set_xlim(min_sep,max_sep)
        
        # Set axes labels
        if corr_type == 'gg':
            plt.xlabel(r'$\Delta\theta$ (arcmin)',fontsize=20)
            plt.ylabel(r'$\xi_+(\theta)\,\,$'+ '(blue), '+r'$\quad\xi_\times(\theta)\,\,$'+' (green)',fontsize=20)
            
        elif corr_type == 'ng':
            plt.ylabel(r'$\xi_{gm}(\theta)$',fontsize=20)
            plt.xlabel(r'$\Delta\theta$ (arcmin)',fontsize=20)
            
        else:
            # Can incorporate other correlation functions later if needed
            pass
        
        # Plot xi=0 for reference
        plt.plot([min(np.exp(corr.logr)),max(np.exp(corr.logr))],[0,0],c='k',linestyle='dashed')

        # Make axis ticks larger
        plt.gca().tick_params('both', length=5, width=1, which='major')
        plt.gca().tick_params('both', length=5, width=1, which='minor')
    
    if corr_type == 'gg':
        # For shear-shear (or ellipticity-ellipticity) correlation
        if lensed == True:
            # Plot lensed correlations as solid lines
            ls = 'solid'
            lab = 'Lensed '
        
        else:
            # Plot pre-lensed correlations as dashed lines
            ls = 'dashed'
            lab = 'Intrinsic '
        
        # Plot the shear-shear (or ellipticity-ellipticity) correlation function (xi_+ and xi_x)
        plt.errorbar(np.exp(corr.logr),corr.xip, np.sqrt(corr.varxi),c='b',linestyle=ls,label=lab+r'$\xi_+$')
        plt.errorbar(np.exp(corr.logr),0.5*(corr.xip_im+corr.xim_im),0.5*np.sqrt(corr.varxi+corr.varxi),c='g',linestyle=ls,label=lab+r'$\xi_\times$')
    
    elif corr_type == 'ng':
        # Plot the galaxi-mass correlation function (xi_gm)
        plt.errorbar(np.exp(corr.logr), corr.xi, np.sqrt(corr.varxi), c='b',label=r'$\operatorname{Re}\left(\xi_{gm}\right)$')
        plt.errorbar(np.exp(corr.logr), corr.xi_im, np.sqrt(corr.varxi), c='g',label=r'$\operatorname{Im}\left(\xi_{gm}\right)$')
    
    else:
        # Can incorporate other correlation functions here if needed
        pass
    
    # get legend handles
    handles, labels = plt.gca().get_legend_handles_labels()
    # remove the errorbars
    handles = [h[0] for h in handles]
    plt.legend(handles,labels,fontsize=18)
    
    return

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

# ----------------------------------------------------------------------------

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
