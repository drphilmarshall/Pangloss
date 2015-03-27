
import pangloss

import numpy
import distances
from scipy import interpolate,optimize

vb = False

# ============================================================================

class Grid():
    """
    NAME
        Grid

    PURPOSE
        Define a discrete redshift grid, and enable objects to be 
        snapped on to it.

    COMMENTS

    INITIALISATION
        zl            Strong lens redshift (needed for critical densities etc)
        zs            Source plane redshift
        nplanes       Number of redshift planes in grid (def=100)   
        cosmo         Cosmological parameters (def: [Om,Ol,h]=[0.25,0.75,0.73] 

    METHODS
        snap(self,z): Return redshift of nearest plane to z
    
    BUGS

    AUTHORS
      This file is part of the Pangloss project, distributed under the
      GPL v2, by Tom Collett (IoA) and  Phil Marshall (Oxford). 
      Please cite: Collett et al 2013, http://arxiv.org/abs/1303.6564

    HISTORY
      2013-03-23  Collett & Marshall (Cambridge)
    """

# ----------------------------------------------------------------------------

    def __init__(self,zl,zs,nplanes=100,cosmo=[0.25,0.75,0.73]): 

        assert zs > zl
        
        D = distances.Distance()
        self.name = '1D Redshift grid of lens planes, each containing precalculated quantities '
        self.zmax = zs*1.0
        self.zs = zs*1.0
        self.zltrue = zl
        self.nplanes = nplanes
        self.cosmo = cosmo
        
        # These are the plane redshifts:
        self.redshifts,self.dz = self.redshiftbins = numpy.linspace(0.0,self.zmax,self.nplanes,endpoint=True,retstep=True)
        self.redshifts += (self.dz/2.)
        self.nz = len(self.redshifts)
        
        # Snap lens to grid, and compute special distances:
        self.zl = self.snap([zl])[0][0]
        self.Da_l = D.Da(zl)
        self.Da_s = D.Da(zs)
        self.Da_ls = D.Da(zl,zs)
        self.plane = {}

        # Grid planes:
        self.Da_p = numpy.zeros(self.nz)
        self.rho_crit = numpy.zeros(self.nz)
        self.Da_ps = numpy.zeros(self.nz)
        self.Da_pl = numpy.zeros(self.nz)
        self.sigma_crit = numpy.zeros(self.nz)
        self.beta = numpy.zeros(self.nz)
        for i in range(self.nz):
            z = self.redshifts[i]
            self.Da_p[i] = D.Da(0,z)
            self.rho_crit[i] = D.rho_crit_univ(z)
            self.Da_ps[i] = D.Da(z,zs)
            self.Da_pl[i] = D.Da(z,zl)
            self.sigma_crit[i] = (1.663*10**18)*(self.Da_s/(self.Da_p[i]*self.Da_ps[i]))  # units M_sun/Mpc^2
            if z > zl: # 1 is lens, 2 is perturber
                D1s = self.Da_ls
                D2  = self.Da_p[i]
            else:      # 1 is perturber, 2 is lens
                D1s = self.Da_ps[i]
                D2  = self.Da_l
            D12 = self.Da_pl[i]
            self.beta[i] = (D12*self.Da_s)/(D2*D1s)

        return

# ---------------------------------------------------------------------------

    def snap(self,z):
        snapped_p = numpy.digitize(z,self.redshifts-(self.dz)/2.0)-1
        snapped_p[snapped_p < 0] = 0 # catalogs have some blue-shifted objects!
        snapped_z = self.redshifts[snapped_p]
        return snapped_z,snapped_p

# ---------------------------------------------------------------------------

    def __str__(self):
        return '1-D Grid of %i planes seperated in redshift by dz= %f' % (self.nplanes,self.dz)

# ============================================================================

if __name__ == '__main__':

    nplanes=100
    g=Grid(0.6,1.4,nplanes=nplanes)
    testfile = "/data/tcollett/Pangloss/grid%i.grid"%nplanes
    pangloss.writePickle(g,testfile)

# ============================================================================
