
import pangloss

import numpy
import distances
from scipy import interpolate,optimize,integrate

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

        # These are set in calculate_bin_volumes() if smooth-component speedup
        # is called in background.py
        self.volumes = None
        self.lc_radius = None

        return

    # ----------------------------------------------------------------------------

    def calculate_bin_volumes(self,lc_radius=None):
        '''
        Calculate all redshift bin proper volumes given a lightcone radius.
        '''

        assert lc_radius is not None, 'Must pass a lightcone radius!'

        self.volumes = numpy.zeros(self.nplanes)
        self.lc_radius = lc_radius

        for i in range(self.nplanes):
            z = self.redshifts[i]
            self.volumes[i] = self.bin_proper_volume(z,lc_radius)


    def bin_proper_volume(self,z,lc_radius):
        '''
        Calculate proper volume of a lightcone redshift bin at z in Mpc^3.
        Details in `pangloss/doc/void_correction.pdf`. Placed here instead of
        in `lightcone.py` so volumes can be pre-calculated in `background` to
        increase speed of smooth-component correction calculation.
        '''

        delta_z = self.dz/2.0
        omega_m0 = self.cosmo[0]
        omega_l0 = self.cosmo[1]
        omega_k = 1. - (omega_m0 + omega_l0)
        h = self.cosmo[2]
        H_0 = 100.*h # km/s/Mpc
        c = 299792.458 # km/s
        Dh = c/H_0 # Hubble distance in Mpc

        # Useful distance functions
        D = pangloss.Distance(self.cosmo)

        # z-range is z +/- delta_z (except at boundaries)
        assert z in self.redshifts
        if z == self.redshifts[0]: z1, z2 = z, z+delta_z
        elif z == self.redshifts[-1]: z1, z2 = z-delta_z, z
        else: z1, z2 = z-delta_z, z+delta_z

        # Solid angle of cone with apex angle of lightcone diameter:
        solid_angle = 2.*numpy.pi*( 1.-numpy.cos(lc_radius*pangloss.arcmin2rad) )

        # Integrand for proper volume
        f = lambda z,o_m,o_l,o_k: Dh * (D.angular_diameter_distance(0,z)**2) / ( (o_m*(1.+z)**3+o_k*(1.+z)**2+o_l)**0.5 * (1.+z) )

        # Return proper volume of redshift slice
        return solid_angle * integrate.romberg(f,z1,z2,(omega_m0,omega_l0,omega_k)) # M_sol/Mpc^3

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
