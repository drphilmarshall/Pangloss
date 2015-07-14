import numpy as np
import matplotlib.pyplot as plt
import os, sys
from astropy.table import Table, Column
from matplotlib.patches import Ellipse
from matplotlib.collections import LineCollection

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
    assert len(ra) == len(dec) and len(dec) == len(phi)
    
    # Convert the angles phi to rad
    phi = np.deg2rad(phi)

    # Preallocation
    pt1 = np.zeros(len(ra),dtype=tuple)
    pt2 = np.zeros(len(ra),dtype=tuple)
    lines = np.zeros(len(ra),dtype=tuple)
    
    # Set scale size of sticks
    
    xi, xf = axis.get_xlim()
    yi, yf = axis.get_ylim()
    Lx, Ly = abs(xf-xi), abs(yf-yi)
    L = np.mean([Lx,Ly]) 
    
    #scale = ((np.log10(mod)-0.5*np.mean(mod))/(max(mod)-0.5*np.mean(mod)))
    #scale = (mod-0.5*np.mean(mod))/(max(mod)-0.5*np.mean(mod))
    #floor = 0.01
    #size = 10.0*L*(scale*(scale > 0) + floor)
    
    # Need this to see sticks (weak lensing)
    scale = 1
    size = scale*mod*L
    #size=scale
    
    # Each stick has the same size
    #s = .001

    # For every object, create a line centered at (ra[i],dec[i]) with appropriate size
    # and orientation angle phi
    for i in range(len(ra)):
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
