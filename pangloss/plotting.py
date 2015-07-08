import numpy as np
import matplotlib.pyplot as plt
import os, sys
from astropy.table import Table, Column

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

def make_axes(fig,subplot):
    '''
    Creates axes for plotting maps and catalogs. Both sets of axes are
    contained in the current figure instance, and also returned for
    ease of use.
    '''
    # Make wcs axes for plotting catalogs in:
    world = fig.add_axes(viewport,label='world')
    world.set_xlabel('Right Ascension / deg')
    world.set_ylabel('Declination / deg')
    world.set_title('')
    world.set_xlim(subplot[0],subplot[1])
    world.set_ylim(subplot[2],subplot[3])

    # Make image axes for plotting maps in:
    image = fig.add_axes(viewport,label='image')
    image.get_xaxis().set_visible(False)
    image.get_yaxis().set_visible(False)

    return image,world

# ----------------------------------------------------------------------------

def set_figure_size(fig,fig_size,Lx,Ly):
    if Lx == Ly:
        fig.set_size_inches(fig_size,fig_size)
    elif Lx > Ly:
        fig.set_size_inches(fig_size,fig_size*(1.0*Ly/Lx))
    else:
        fig.set_size_inches(fig_size*(1.0*Lx/Ly),fig_size)
    return


# ----------------------------------------------------------------------------
'''
def plotCatalogOnMap(F,fig_size=10,subplot=None,mag_cutoff=[0,24],mass_cutoff=[0,10**20],z_cutoff=[0,1.3857]):
    ''
    Plot a foreground galaxy catalog field with filters over its corresponding convergence and shear map.
    NOTE: Currently only supports plotting one catalog over one map. Should expand in future.
    ''

    if subplot is None:
        # Default subplot is entire image
        ai, di = F.ra_max, F.dec_min
        af, df = F.ra_min, F.dec_max
        subplot = [ai,af,di,df]

    ai, af = subplot[0], subplot[1]    # RA limits for subplot
    di, df = subplot[2], subplot[3]    # DEC limits for subplot

    # Include Pangloss directory in path
    PANGLOSS_DIR = os.path.expandvars("$PANGLOSS_DIR")
    sys.path.append(PANGLOSS_DIR)

    # Create the corresponding convergence and shear maps
    K = Kappamap(PANGLOSS_DIR+'/data/GGL_los_8_'+str(F.map_x)+'_'+str(F.map_y)+'_N_4096_ang_4_rays_to_plane_37_f.kappa',FITS=False)
    S = Shearmap([PANGLOSS_DIR+'/data/GGL_los_8_'+str(F.map_x)+'_'+str(F.map_y)+'_N_4096_ang_4_rays_to_plane_37_f.gamma_1', \
                  PANGLOSS_DIR+'/data/GGL_los_8_'+str(F.map_x)+'_'+str(F.map_y)+'_N_4096_ang_4_rays_to_plane_37_f.gamma_2'],FITS=False)

    # Find world coordinates and masses of galaxies that meet cutoff criteria
    ra_cutoff, dec_cutoff = [ai, af], [di, df]     # RA flipped because RA is left-handed
    ra, dec, mass = F.findGalaxies(mag_cutoff,mass_cutoff,z_cutoff,ra_cutoff,dec_cutoff)
    print('max ra: ',max(ra),'min ra: ',min(ra),'max dec: ',max(dec),'min dec: ',min(dec))

    # Convert the galaxies' wc positions to pixel bins in the (x,y) map
    pix_ra = [K.world2image(a,0)[0] for a in ra]
    pix_dec = [K.world2image(0,d)[1] for d in dec]
    print('min pix_ra: ',min(pix_ra),'max pix_ra: ',max(pix_ra),'min pix_dec: ',min(pix_dec),'max pix_dec: ',max(pix_dec))

    # Scale galaxy plot size by its mass
    scale = ((np.log10(mass)-9.0)/(12.0-9.0))
    floor = 0.01
    size = 1000.0*(scale*(scale > 0) + floor)

    # Plot chosen galaxies over convergence and shear map
    K.plot(fig_size,subplot)
    S.plot(fig_size,subplot)
    #plt.scatter(pix_ra,pix_dec,s=size,color='orange',edgecolor=None,alpha=0.2)
    plt.xlabel('Right Ascension / deg')
    plt.ylabel('Declination / deg')
    #plt.gca().set_xlim((F.ra_min,F.ra_max))
    #plt.gca().set_ylim((F.dec_min,F.dec_max))

    # Set figure size to fig_size
    fig = plt.gcf()
    fig.set_size_inches(fig_size,fig_size)

    return
'''
