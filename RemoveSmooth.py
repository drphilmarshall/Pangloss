#!/usr/bin/env python
# ======================================================================

import pangloss

import sys,getopt,cPickle,numpy
import matplotlib.pyplot as plt

from math import pi


# ======================================================================

def RemoveSmooth(argv):
    """
    NAME
        RemoveSmooth.py

    PURPOSE


    COMMENTS

    INPUTS
        configfile    Plain text file containing Pangloss configuration

    OUTPUTS


    EXAMPLE

    BUGS

    AUTHORS

    HISTORY
      2014-04-23 started Mason (UCSB)
    """
    # --------------------------------------------------------------------

    try:
       opts, args = getopt.getopt(argv,"h",["help"])
    except getopt.GetoptError, err:
       print str(err) # will print something like "option -a not recognized"
       print RemoveSmooth.__doc__  # will print the big comment above.
       return

    for o,a in opts:
       if o in ("-h", "--help"):
          print RemoveSmooth.__doc__
          return
       else:
          assert False, "unhandled option"

    # Check for setup file in array args:
    if len(args) == 1:
        configfile = args[0]
        print pangloss.doubledashedline
        print pangloss.hello
        print pangloss.doubledashedline
        print "RemoveSmooth: finding the smooth component of convergence"
        print "RemoveSmooth: taking instructions from",configfile
    else:
        print RemoveSmooth.__doc__
        return



    # --------------------------------------------------------------------
    # Read in configuration
    # --------------------------------------------------------------------
    # --------------------------------------------------------------------
    # Read in the whole catalog
    
    experiment = pangloss.Configuration(configfile)

    catalog = experiment.parameters['ObservedCatalog'][0]
    
    table = pangloss.readCatalog(catalog,experiment)

    zd = experiment.parameters['StrongLensRedshift']
    zs = experiment.parameters['SourceRedshift']
    print 'Reading in the catalog...' 
    
    # Calculate the area of the patch in square radians
    xmax = table['nRA'].max()
    xmin = table['nRA'].min()
    ymax = table['Dec'].max()
    ymin = table['Dec'].min()     
    
    xrange_rad = numpy.abs(xmax - xmin)
    yrange_rad = numpy.abs(ymax - ymin)
    area_rad = xrange_rad * yrange_rad
    area_arcmin = area_rad * (pangloss.rad2arcmin**2)
    D = pangloss.Distance()

    print 'Calculating the area of the catalog...'    
    
    # --------------------------------------------------------------------    
    # Remove the subhalos if they remain (they have no stars...)
    try: 
        halos = table.where(table.Type != 2) 
    except AttributeError: pass    
    
    # Apply a magnitude cut
    nhalos_all = len(halos)
    halos = halos.where(halos.mag < 25)    
    nhalos = len(halos)   
    print 'There are',nhalos_all,'halos in this catalog'
    num_density = nhalos_all/area_arcmin
    print 'The number density of halos is',num_density,'in 1 square arcmin of sky...'
    # --------------------------------------------------------------------
    # Add galaxy property column, overwriting any values that already exist:

    def writeColumn(string,values):
        try:
           halos.add_column('%s'%string,values)
        except ValueError:
           halos["%s"%string]=values    
    
    # --------------------------------------------------------------------
    # Bin catalog into redshift slices      
    nplanes=100
    
    grid = pangloss.Grid(zd,zs,nplanes)

    halos = halos.where(halos.z_obs<zs+0.2)
    
    print 'Binning the catalog into',nplanes,'redshift slices...'
    
    # ----------------------------------------------------------------------------
    # Load the grid, find cosmology values for all the halos
    def loadGrid(Grid):
        if numpy.abs(zd-Grid.zltrue) > 0.05: print "Grid zl != lens zl" 
        if numpy.abs(zs-Grid.zs)     > 0.05: print "Grid zs != lens zs" 
        redshifts,dz = Grid.redshifts,Grid.dz
        Da_l,Da_s,Da_ls = Grid.Da_l,Grid.Da_s,Grid.Da_ls
        return
                 
    loadGrid(grid)
                           
    z=halos.z_obs
    sz,p = grid.snap(z)
    writeColumn('plane',p)
    
    # Rho_crit and sigma_crit at each redshift      
    writeColumn('rho_crit',grid.rho_crit[p])
    writeColumn('sigma_crit',grid.sigma_crit[p])
    
    # Angular diameter distance
    writeColumn('Da_p',grid.Da_p[p])
    
    
    # ----------------------------------------------------------------------------
                       
    # Get the parameters for the halos
    M200 = halos.Mhalo_obs        
    r200 = (3*M200/(800*pi*halos.rho_crit))**(1./3)
        
    c200 = pangloss.MCrelation(M200,scatter=True)
    r_s = r200/c200        
    writeColumn('rs',r_s)
    
    
    # ----------------------------------------------------------------------------
    # Now calculate total mass in each redshift slice, must be in solar mass
    # Divide by area of slice in Mpc^2
    
    def NFWtruncMass(truncationscale=10):
        """ Calculate all mass using formula from BMO09 within truncation scale """
        tau = truncationscale
        delta_c =  pangloss.delta_c(c200)
        
        M0 = 16*pi*(r_s**3)*delta_c*halos.rho_crit
        
        M_trunc = M0 * (tau**2/(tau**2 + 1.)**2) * ((tau**2 - 1.)*numpy.log(tau) + tau*pi - tau**2 - 1)
        writeColumn('M_trunc', M_trunc)
        
        return
        
    # --------------------------------------------------------------------

    NFWtruncMass(truncationscale=5)
    mass_total = numpy.sum(halos.M_trunc)
    log_mtot = numpy.log10(mass_total)
    print 'Calculating the truncated NFW profile mass of each halo...'
    print 'The total halo mass is 10^'+ '%.2f' % log_mtot +' M_sun if truncation radius is 10 R_vir...'

    kappa_slice = numpy.zeros(nplanes)
    
    for p in range(nplanes):
        
        # Get subset of halos at each redshift slice
        halos_slice = halos.where(halos.plane == p)
        
        if len(halos_slice) > 0:
        
            # Add up all the mass
            mass_tot = numpy.sum(halos_slice.M_trunc)
    
            # Proper area of each slice in Mpc^2
            area_prop = area_rad*(halos_slice.Da_p[0]**2)
            
            # Sigma of each slice in Msun/Mpc^2
            sigma_slice = mass_tot/area_prop
            
            # Kappa slice
            kappa_slice[p] = sigma_slice / halos_slice.sigma_crit[0]

        else:
            kappa_slice[p] = 0.
            
    kappa_smooth = numpy.sum(kappa_slice)
    
    print pangloss.doubledashedline       
    print 'The smooth component of the universe has kappa =',kappa_smooth
    print pangloss.doubledashedline 
    
    return
    
    

# ======================================================================

if __name__ == '__main__': 
    RemoveSmooth(sys.argv[1:])

# ======================================================================    