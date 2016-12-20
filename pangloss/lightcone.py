# ===========================================================================

import numpy as np
import pylab as plt
import scipy as sp
import pandas as pd
import pangloss

# ======================================================================

class Lightcone(object):
    """
    NAME
        Lightcone

    PURPOSE
        Define a conical region of space containing many galaxies, and
        enable mass to be assigned to those galaxies and consequent
        useful quantities to be calculated.

    COMMENTS
        The masses of galaxies in a lightcone are always uncertain.
        Methods are provided to characterise this uncertainty by drawing
        sample masses for each galaxy, from various specified relations.

    INITIALISATION
        catalog       astropy table of the parent galaxy catalog
        flavor        Is the catalog 'real' or 'simulated'?
        position      The sky position (J2000 deg) of the cone centre
        radius        The radius of the lightcone field of view (arcmin)
        ID            Lightcone ID number in a given background catalog
        maglimit      The depth of the galaxy selection (magnitudes)
        band          The band in which the selection is made

    METHODS
        galaxiesWithin(self,radius,cut=[18.5,24.5],band="F814W",radius_unit="arcsec"):

        numberWithin(self,radius,cut=[18.5,24.5],band="F814W",radius_unit="arcsec"):

        define_system(self,zl,zs,cosmo=[0.25,0.75,0.73]):

        loadGrid(self, Grid):

        mimicPhotozError(self,sigma=0.1):

        snapToGrid(self, Grid):

        drawMstars(self,model): Needs updating to take SHMR object

        drawMhalos(self,modelT):

        drawConcentrations(self,errors=False):

        makeKappas(self,errors=False,truncationscale=5,profile="BMO1"):

        combineKappas(self):

    BUGS

    HISTORY
      2013-03-23  Collett & Marshall (Cambridge)
      2015-07-29  Everett (SLAC)
    """

# ----------------------------------------------------------------------------

    def __init__(self, catalog, flavor, position, radius, ID=None, maglimit=99, band="r"):

        self.name = 'Lightcone through the Universe'
        self.flavor = flavor   # 'real' or 'simulated'
        self.catalog = catalog # table of galaxies
        self.ID = ID           # ID number from the cone's center background galaxy

        # Simulated lightcones have "true" (ray-traced) convergence:
        self.kappa_hilbert = None # until set!

        # Catalog limits:
        self.xmax = np.max(self.catalog['RA'])
        self.xmin = np.min(self.catalog['RA'])
        self.ymax = np.max(self.catalog['Dec'])
        self.ymin = np.min(self.catalog['Dec'])

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        # Cut out a square cone:
        self.rmax = radius
        self.xc = [position[0],position[1]]

        dx = self.rmax*pangloss.arcmin2rad
        # Ra is flipped as it is left-handed
        self.galaxies = pd.DataFrame.copy(self.catalog[(self.catalog['RA'] > (self.xc[0]-dx)) & \
                                              (self.catalog['RA'] < (self.xc[0]+dx)) & \
                                              (self.catalog['Dec'] > (self.xc[1]-dx)) & \
                                              (self.catalog['Dec'] < (self.xc[1]+dx)) & \
                                               self.catalog['Mh'] > 0
                                              ])

        # Trim it to a circle:
        x = -np.cos(self.galaxies['Dec'])*(self.galaxies['RA'] - self.xc[0])*pangloss.rad2arcmin
        y = (self.galaxies['Dec'] - self.xc[1])*pangloss.rad2arcmin
        r = np.sqrt(x*x + y*y)
        phi = np.arctan2(y,x)
        self.galaxies['x'] = x
        self.galaxies['y'] = y
        self.galaxies['r'] = r
        self.galaxies['phi'] = phi
        self.galaxies = pd.DataFrame.copy(self.galaxies[self.galaxies['r'] < self.rmax])

        try:
            self.galaxies = self.galaxies[self.galaxies['Type'] != 2]
        except AttributeError: pass

        # Why is this here? If we want to save all the galaxies, shouldn't it be placed above the galaxies[np.where()] line? - Spencer
        # self.allgalaxies = self.galaxies

        # Keep track of the number of galaxies in each lightcone
        self.galaxy_count = len(self.galaxies)

        # Now we have a small catalog, just for the requested lightcone.
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        # Halo parameters, to be varied during sampling analysis:
        self.galaxies['z'] = self.galaxies['z_obs'] * 1.0

        if self.flavor == 'simulated':
            # Take the log of the halo mass, and set up the parameter array:
            zero = np.sum(self.galaxies['Mhalo_obs'] == 0)
            if zero != 0: print 'There are {} galaxies with subhalo mass = 0'.format(zero)
            self.galaxies['Mh_obs'] = np.log10(self.galaxies['Mhalo_obs'])
            self.galaxies['Mh'] = self.galaxies['Mh_obs'] * 1.0
            # Stellar masses will be added by drawMstars
            # Halo masses will be replaced by drawMhalos

        elif self.flavor == 'real':
            # Mstar is already given as log M...
            self.galaxies['Mstar'] = self.galaxies['Mstar_obs'] * 1.0
            # Halo masses will be added by drawMhalos

        elif self.flavor == 'sampled':
            # Take the log of the halo mass, and set up the parameter array:
            zero = np.sum(self.galaxies['Mhalo_obs']==0)
            if zero != 0: print 'There are {} galaxies with subhalo mass = 0'.format(zero)
            self.galaxies['Mh_obs'] = np.log10(self.galaxies['Mhalo_obs'])

        if len(self.galaxies) == 0:
            print "Lightcone: WARNING: no galaxies here!"

        # Set lightcone

        # No smooth-component correction kappas until set!
        self.smooth_kappas = None

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        # Save memory!
        del self.catalog
        del catalog

        return

# ----------------------------------------------------------------------------

    def __str__(self):
        return 'Lightcone of radius %.2f arcmin, centred on (%.3f,%.3f) rad' % (self.rmax,self.xc[0],self.xc[1])

# ----------------------------------------------------------------------------
# Tell me the number of galaxies within a certain radius, that pass a
# certain magnitude cut.

    def galaxiesWithin(self,radius,cut=[18.5,24.5], band="F814W", radius_unit="arcmin"):

        if band == "u" or band ==  "g" or band == "r" or band ==  "i" or band == "z":
            col = "mag_SDSS_%s" % band
        elif band == "F814" or band == "F814W" or band == "814" or band == 814:
            col = "mag_F814W"
        elif band == "WFC125" or band == "F125W" or band == "F125" or band == "125" or band == 125:
            col = "WFC125"
        else:
            col = "mag_%s" % band
        if radius < 0.1:
            print "Warning: Default units for radius are arcmin!"
        if radius_unit == "arcmin":
            radius = radius
        if radius_unit == "arcsec":
            radius = radius/60.
        if col != "warning":
          #  self.N_cut=self.galaxies[np.where((self.galaxies.['r'] < radius)  & \
          #                                    (self.galaxies["%s"%col] < cut[1])& \
          #                                (self.galaxies["%s"%col] > cut[0]))]
            self.N_cut=self.galaxies[(self.galaxies['r'] < radius)  & \
                                           (self.galaxies["%s"%col] < cut[1])& \
                                           (self.galaxies["%s"%col] > cut[0])]

            return self.N_cut

    def numberWithin(self,radius,cut=[18.5,24.5],band="F125W",units="arcmin"):
        Ntable = self.galaxiesWithin(radius,cut,band,units)
        return len(Ntable.r)

# ----------------------------------------------------------------------------

    def defineSystem(self,zl,zs,cosmo=[0.25,0.75,0.73]):
        self.zl = zl
        self.zs = zs
        self.cosmo = cosmo
        self.galaxies = pd.DataFrame.copy(self.galaxies[self.galaxies['z_obs'] < zs])
        self.galaxy_count = len(self.galaxies)
        return

# ----------------------------------------------------------------------------

    def loadGrid(self, Grid):
        if np.abs(self.zl-Grid.zltrue) > 0.05: print "Grid zl != lens zl"
        if np.abs(self.zs-Grid.zs)     > 0.05: print "Grid zs != lens zs"
        self.redshifts,self.dz = Grid.redshifts,Grid.dz
        self.Da_l,self.Da_s,self.Da_ls = Grid.Da_l,Grid.Da_s,Grid.Da_ls
        self.volumes = Grid.volumes
        return

# ----------------------------------------------------------------------------
# Simple model for spectroscopic coverage - to be applied to the
# calibration lightcones to match the observed catalog.

    def configureForSurvey(self, experiment):

        PR = experiment.parameters['PhotometricRadius']
        PD = experiment.parameters['PhotometricDepth']
        assert len(PR) == len(PD)

        SR = experiment.parameters['SpectroscopicRadius']
        SD = experiment.parameters['SpectroscopicDepth']
        assert len(SR)==len(SD)

        band = experiment.parameters['LightconeDepthBand']

        if band == "u" or band ==  "g" or band == "r" or band ==  "i" or band == "z":
            col = "mag_SDSS_%s" % band
        elif band == "F814" or band == "F814W" or band == "814" or band == 814:
            col = "mag_F814W" #note that this isn't included atm
        elif band == "WFC125" or band == "F125" or band == "F125W" or band == "125" or band == 125:
            col = "WFC125" #this was the matching band for BoRG
        else:
            col = "mag_%s" % band

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        # Only include galaxies observed by photometry:

        self.galaxies['photo_flag'] = False
        self.galaxies['identifier'] = range(len(self.galaxies['x']))
        if PR != ['']:
            for i in range(len(PR)):
                R=PR[i]*60 # positions are stored in arcseconds
                D=PD[i]

                goodset = set(self.galaxies[(self.galaxies['r'] < R) & \
                   (self.galaxies["%s"%col] < D)['identifier']])

                self.galaxies['photo_flag'][np.array(\
                   [_ in goodset for _ in self.galaxies['identifier']])] = True

        self.galaxies = self.galaxies[self.galaxies['photo_flag']==True]

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        # Set spectroscopic flag of any galaxy that should have
        # spectroscopy, according to the given radius and depth:

        self.galaxies['spec_flag'] = False

        if SR!=['']:
            for i in range(len(SR)):
                R=SR[i]*60 # positions are stored in arcseconds
                D=SD[i]
                goodset=set(self.galaxies[(self.galaxies['r'] < R) & \
                   (self.galaxies["%s"%col] < D)]['identifier'])

                self.galaxies['spec_flag'][np.array(\
                   [_ in goodset for _ in self.galaxies['identifier']])]=True

        return

# ----------------------------------------------------------------------------
# The following methods are designed to be run multiple times
# (previous ones are single use per lightcone)
# ----------------------------------------------------------------------------

    def mimicPhotozError(self,sigma=0.1):

        # Start with the original, constant z_obs in the catalog:
        z_obs = self.galaxies['z_obs'].copy()

        # RMS error on this z is either sigma or 0.0:
        e = ( self.galaxies['spec_flag'] == False )*sigma

        # Add Gaussian noise to get this realisation's z:
        z = z_obs + e*(1+z_obs)*np.random.randn(len(z_obs))

        # Over-write the z parameter array:
        self.galaxies['z'] = z

        return

# ----------------------------------------------------------------------------
# Snap the parameters z onto the grid, to speed up calculations:

    def snapToGrid(self, Grid,foreground_kappas=None):
        z = self.galaxies['z']
        sz,p = Grid.snap(z)
        self.sigma_crit = Grid.sigma_crit
        self.galaxies['z_sz'] = sz
        self.galaxies['Da_p'] = Grid.Da_p[p]
        self.galaxies['rho_crit'] = Grid.rho_crit[p]
        self.galaxies['sigma_crit'] = Grid.sigma_crit[p]
        self.galaxies['beta'] = Grid.beta[p]
        rphys = self.galaxies['r']*pangloss.arcmin2rad*self.galaxies['Da_p']
        self.galaxies['rphys'] = rphys

        # Write foreground mean kappas for foreground correction (0 if foreground kappas not inputed)
        if foreground_kappas is None:
            f_kappa = np.zeros(np.size(z))
        else:
            self.foreground_kappas = foreground_kappas
            f_kappa = foreground_kappas[p]
        self.galaxies['f_kappa'] = f_kappa

        return

# ----------------------------------------------------------------------------
# Given Mhalo and z, draw an Mstar, and an identical Mstar_obs:

    def drawMstars(self,model):
        Mstar = model.drawMstars(self.galaxies['Mh'],self.galaxies['z'])
        self.galaxies['Mstar'] = Mstar
        # Copy this to Mstar_obs - mimicMstarError will deal with noise
        self.galaxies['Mstar_obs'] = Mstar
        return

# ----------------------------------------------------------------------------
# Given an observed Mstar_obs, what could the parameter Mstar be?

    def mimicMstarError(self,sigmaP,sigmaS):

        Mstar = self.galaxies['Mstar_obs'].copy()

        # Add offset due to photometric or spectroscopic error:
        Mstar[self.galaxies['spec_flag']==False] += np.random.randn(Mstar[self.galaxies['spec_flag']==False].size)*sigmaP
        Mstar[self.galaxies['spec_flag']==True] += np.random.randn(Mstar[self.galaxies['spec_flag']==True].size)*sigmaS

        # Over-write the Mstar parameter array:
        self.galaxies['Mstar'] = Mstar
        return


# ----------------------------------------------------------------------------
# Given an Mstar and z, what could the parameter Mh be?

    def drawMhalos(self,model):
        Mh = model.drawMhalos(self.galaxies['Mstar'],self.galaxies['z'])
        self.galaxies['Mh'] = Mh
        return

# ----------------------------------------------------------------------------
# Given an Mh, what could the halo concentration be?

    def drawConcentrations(self,errors=False):
        M200 = self.galaxies['Mh']
        r200 = (3*M200/(800*3.14159*self.galaxies['rho_crit']))**(1./3)
        self.galaxies["r200"] = r200
        c200 = pangloss.MCrelation(M200,scatter=errors)
        self.galaxies["c200"] = c200
        r_s = r200/c200
        self.galaxies['rs'] = r_s
        self.galaxies['X'] = self.galaxies['rphys']/r_s

        return

# ----------------------------------------------------------------------------
# Compute halos' contributions to the convergence:

    def makeKappas(self,errors=False,truncationscale=5,profile="BMO1",lensing_table=None):
        rho_s = pangloss.delta_c(self.galaxies['c200'])*self.galaxies['rho_crit']
        self.kappa_s = rho_s * self.galaxies['rs'] /self.galaxies['sigma_crit']  #kappa slice for each lightcone
        r_trunc = truncationscale*self.galaxies['r200']
        xtrunc = r_trunc/self.galaxies['rs']

        if profile=="BMO1":
            if lensing_table is None:
                # Compute the BMO truncated profile for each object
                F = pangloss.BMO1Ffunc(self.galaxies['X'].as_matrix(),xtrunc.as_matrix())
                G = pangloss.BMO1Gfunc(self.galaxies['X'].as_matrix(),xtrunc.as_matrix())
                #F = pangloss.BMO1FSpencerFunc(self.galaxies['X'],xtrunc)
                #G = pangloss.BMO1GSpencerFunc(self.galaxies['X'],xtrunc)
            else:
                # Use an interpolated lookup table for the profile for each object
                F = lensing_table.lookup_BMO1F(self.galaxies['X'].as_matrix(),xtrunc.as_matrix())
                G = lensing_table.lookup_BMO1G(self.galaxies['X'].as_matrix(),xtrunc.as_matrix())

        if profile=="BMO2":
            if lensing_table is None:
                # Compute the BMO truncated profile for each object
                F=pangloss.BMO2Ffunc(self.galaxies['X'].as_matrix(),xtrunc.as_matrix())
                G=pangloss.BMO2Gfunc(self.galaxies['X'].as_matrix(),xtrunc.as_matrix())
                #F = pangloss.BMO2FSpencerFunc(self.galaxies['X'],xtrunc)
                #G = pangloss.BMO2GSpencerFunc(self.galaxies['X'],xtrunc)
            else:
                # Use an interpolated lookup table for the profile for each object
                F = lensing_table.lookup_BMO2F(self.galaxies['X'].as_matrix(),xtrunc.as_matrix())
                G = lensing_table.lookup_BMO2G(self.galaxies['X'].as_matrix(),xtrunc.as_matrix())

        kappa = 1.0*self.kappa_s*F
        gamma = 1.0*self.kappa_s*(G-F)
        gamma1 = gamma*np.cos(2*self.galaxies['phi'])
        gamma2 = gamma*np.sin(2*self.galaxies['phi'])

        mu = 1.0/(((1.0 - kappa)**2.0) - (gamma**2.0))

        self.galaxies['kappa'] = kappa
        self.galaxies['gamma'] = gamma
        self.galaxies['gamma1'] = -gamma1
        self.galaxies['gamma2'] = -gamma2
        self.galaxies['mu'] = mu

        return

# ----------------------------------------------------------------------------

    def combineKappas(self,methods=['add','keeton','tom'],smooth_corr=False,foreground_kappas=None):

        # If a single string is passed for methods, place it in a list
        if type(methods) != list:
            methods = [methods]

        # Load in the convergence and shear from each galaxy in the lightcone
        K = self.galaxies['kappa']
        G1 = self.galaxies['gamma1']
        G2 = self.galaxies['gamma2']

        if 'keeton' in methods or 'tom' in methods:
                B = self.galaxies['beta']
                G = self.galaxies['gamma']
                D = K**2-G**2

        for method in methods:

            if method == 'keeton':
                # Calculate the 'keeton' convergence and shear for each galaxy in the lightcone
                kappa_keeton  = (1.-B) * (K- B*(D)) /  ( (1-B*K)**2   - (B*G)**2   )
                gamma1_keeton = (1.-B) * (G1) /  ( (1-B*K)**2   - (B*G)**2   )
                gamma2_keeton = (1.-B) * (G2) /  ( (1-B*K)**2   - (B*G)**2   )

                # Write the convergene and shear for each galaxy in the lightcone
                self.galaxies['kappa_keeton'] = kappa_keeton
                self.galaxies['gamma1_keeton'] = gamma1_keeton
                self.galaxies['gamma2_keeton'] = gamma2_keeton

                # Add the total convergence and shear from all galaxies
                self.kappa_keeton_total=np.sum(self.galaxies['kappa_keeton'])
                self.gamma1_keeton_total=np.sum(self.galaxies['gamma1_keeton'])
                self.gamma2_keeton_total=np.sum(self.galaxies['gamma2_keeton'])

                if foreground_kappas is not None:
                    # If foreground kappas are passed, implement foreground correction
                    self.kappa_keeton_total -= np.sum(foreground_kappas)

                if smooth_corr is True:
                    self.kappa_keeton_total += self.smooth_corr_kappas()

            if method == 'tom':
                # Calculate the 'tom' convergence and shear for each galaxy in the lightcone
                kappa_tom  = (1.-B) * K
                gamma1_tom = (1.-B) * G1
                gamma2_tom = (1.-B) * G2

                # Write the convergene and shear for each galaxy in the lightcone
                self.galaxies['kappa_tom'] = kappa_tom
                self.galaxies['gamma1_tom'] = gamma1_tom
                self.galaxies['gamma2_tom'] = gamma2_tom

                # Add the total convergence and shear from all galaxies
                self.kappa_tom_total=np.sum(self.galaxies['kappa_tom'])
                self.gamma1_tom_total=np.sum(self.galaxies['gamma1_tom'])
                self.gamma2_tom_total=np.sum(self.galaxies['gamma2_tom'])

                if foreground_kappas is not None:
                    # If foreground kappas are passed, implement foreground correction
                    self.kappa_tom_total -= np.sum(foreground_kappas)

                if smooth_corr is True:
                    self.kappa_tom_total += self.smooth_corr_kappas()


            if method == 'add':
                # Calculate the convergence and shear for each galaxy in the lightcone
                ## Why are these necessary ?? Should check if we can delete
                self.galaxies['kappa_add'] = K
                self.galaxies['gamma1_add'] = G1
                self.galaxies['gamma2_add'] = G2

                # Add the total convergence and shear from all galaxies
                self.kappa_add_total=np.sum(self.galaxies['kappa'])
                self.gamma1_add_total=np.sum(self.galaxies['gamma1'])
                self.gamma2_add_total=np.sum(self.galaxies['gamma2'])

                if foreground_kappas is not None:
                    # If foreground kappas are passed, implement foreground correction
                    self.kappa_add_total -= np.sum(foreground_kappas)

                if smooth_corr is True:
                    self.kappa_add_total += np.sum(self.smooth_corr_kappas())

        return

    #----------------------------------------------------------------------------
    # Functions used to calculate matter densities and the smooth-component correction

    def smooth_corr_kappas(self):
        '''
        Calculates the smooth-component correction kappas, which account for the additional mass
        in the Millennium Simulation that is not in the form of stellar mass, halos, or
        subhalos. For details, see void_correction.pdf in `pangloss/doc`.
        '''

        redshifts = self.redshifts
        delta_z = self.dz/2.0

        # Useful functions
        D = pangloss.Distance(self.cosmo)

        self.smooth_kappas = np.zeros(np.size(redshifts))

        # Calculate the smooth-component correction at each redshift slice
        for i in range(np.size(redshifts)):
            z = redshifts[i]
            # z-range is z +/- delta_z (except at boundaries)
            if z == redshifts[0]: z1, z2 = z, z+delta_z
            elif z == redshifts[-1]: z1, z2 = z-delta_z, z
            else: z1, z2 = z-delta_z, z+delta_z

            # Mean mass density at given z for assumed cosmology
            rho_m = D.mean_mass_density(z) # M_sol/Mpc^3

            # Halo mass density in this lightcone at this z
            rho_h = self.halo_density(z) # M_sol/Mpc^3

            # Smooth-component mass density is difference between mean mass and halo mass
            rho_s = rho_m - rho_h # M_sol/Mpc^3

            # Projected mass density found by integrating over z in redshift slice
            sigma_s = self.slice_sigma(rho_s,z)

            # Smooth-component kappa correction at redshift z:
            self.smooth_kappas[i] = 1.0*sigma_s / self.sigma_crit[i]

        return self.smooth_kappas

    def halo_density(self,z):
        '''
        Returns the halo mass density at the given redshift slice.
        '''
        assert z in self.redshifts
        return np.sum(self.galaxies['Mh'][self.galaxies['z_sz']==z] ) / self.slice_proper_volume(z) # M_sol/Mpc^3

    def slice_proper_volume(self,z):
        '''
        Calculate proper volume of redshift slice at z in Mpc^3. Details in `pangloss/doc/void_correction.pdf`.
        '''

        if self.volumes is not None:
            # If volumes have been pre-calculated, use them!
            return self.volumes[z==self.redshifts][0]

        # Important parameters
        redshifts = self.redshifts
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
        assert z in redshifts
        if z == redshifts[0]: z1, z2 = z, z+delta_z
        elif z == redshifts[-1]: z1, z2 = z-delta_z, z
        else: z1, z2 = z-delta_z, z+delta_z

        # Solid angle of cone with apex angle of lightcone diameter:
        solid_angle = 2.*np.pi*( 1.-np.cos(self.rmax*pangloss.arcmin2rad) )

        # Integrand for proper volume
        f = lambda z,o_m,o_l,o_k: Dh * (D.angular_diameter_distance(0,z)**2) / ( (o_m*(1.+z)**3+o_k*(1.+z)**2+o_l)**0.5 * (1.+z) )

        # Return proper volume of redshift slice
        return solid_angle * sp.integrate.romberg(f,z1,z2,(omega_m0,omega_l0,omega_k)) # M_sol/Mpc^3

    def slice_sigma(self,rho,z):
        '''
        Calculates the projected surface mass density of rho(z) for the z slice.
        '''

        # Important parameters
        redshifts = self.redshifts
        delta_z = self.dz/2.0
        omega_m0 = self.cosmo[0]
        omega_l0 = self.cosmo[1]
        omega_k = 1. - (omega_m0 + omega_l0)
        h = self.cosmo[2]
        H_0 = 100.*h # km/s/Mpc
        c = 299792.458 # km/s
        Dh = c/H_0 # Hubble distance in Mpc

        # z-range is z +/- delta_z (except at boundaries)
        assert z in redshifts
        if z == redshifts[0]: z1, z2 = z, z+delta_z
        elif z == redshifts[-1]: z1, z2 = z-delta_z, z
        else: z1, z2 = z-delta_z, z+delta_z

        f = lambda z,o_m,o_l,o_k: Dh / ( (o_m*(1.+z)**3+o_k*(1.+z)**2+o_l)**0.5 * (1.+z) )

        return rho * sp.integrate.romberg(f,z1,z2,(omega_m0,omega_l0,omega_k)) # M_sol/Mpc^2

# ----------------------------------------------------------------------------
# Calculate magnification along line of sight

    def combineMus(self,weakapprox=True):

        M=self.galaxies['mu']
        K=self.galaxies['kappa']
        G=self.galaxies['gamma']
        G1=self.galaxies['gamma1']
        G2=self.galaxies['gamma2']

        Ksum = self.kappa_add_total #np.sum(K)
        self.G1sum = np.sum(G1)
        self.G2sum = np.sum(G2)
        self.Gsum = np.sqrt(self.G1sum**2 + self.G2sum**2)

        if weakapprox is True:
            Msum = 1.0 + 2.0*Ksum
        else:
            inverseMsum = (1.0 - Ksum)**2.0 - self.Gsum**2.0
            Msum = 1.0/inverseMsum

        self.mu_add_total=Msum

        return self.mu_add_total

# ----------------------------------------------------------------------------
# Find contribution of various quantities along LoS at given z
# for plotting cumulative sums of parameters with z

    def findContributions(self,quantity):

       # Point positions:
       zmax = self.zs+0.1
       zbins = 15
       z = np.linspace(0.0,zmax,zbins)
       # Plot the points:
       if quantity == 'mass':
           contr = np.zeros(zbins)
           for i in range(zbins):
               galaxies = self.galaxies[self.galaxies['z'] <= z[i]]
               size = galaxies['Mhalo_obs']
               contr[i] = np.sum(size)

       elif quantity == 'kappa':
           contr = np.zeros(zbins)
           for i in range(zbins):
               galaxies = self.galaxies[self.galaxies['z'] <= z[i]]
               size = galaxies['kappa']
               contr[i] = np.sum(size)

       elif quantity == 'mu':
           contr = np.zeros(zbins)
           for i in range(zbins):
               galaxies = self.galaxies[self.galaxies['z'] <= z[i]]
               size = galaxies['mu']
               contr[i] = np.sum(size)

       elif quantity == 'stellarmass':
           contr = np.zeros(zbins)
           for i in range(zbins):
               galaxies = self.galaxies[self.galaxies['z'] <= z[i]]
               size = galaxies['Mstar_obs']
               contr[i] = np.sum(size)
       else:
           raise "Lightcone plotting error: unknown quantity "+quantity

       return contr

# ----------------------------------------------------------------------------
# Plotting
# ----------------------------------------------------------------------------

    def plot_kappas(self,fig_size=10,show=True):
        '''
        Compares the halo, foreground, and smooth-component kappas at each
        redshift slice.
        '''
        assert self.galaxies['kappa'] is not None
        if self.smooth_kappas is None:
            self.smooth_kappas = self.smooth_corr_kappas()

        kappa = self.galaxies['kappa']
        kappa_s = self.smooth_kappas
        redshifts = self.redshifts

        plt.plot(self.galaxies['z_sz'],kappa,'og',label='Individual Kappas')
        plt.plot(redshifts,kappa_s,'ob',label='Smooth-component Kappas')
        plt.gca().set_yscale("log")
        plt.xlabel(r'Redshift ($z$)',fontsize=18)
        plt.ylabel(r'Convergence $\kappa$',fontsize=18)
        plt.tick_params(axis='both', which='major', labelsize=16)
        plt.legend(fontsize=16)
        plt.gcf().set_size_inches(fig_size,0.75*fig_size)

        if show is True: plt.show()

        return

    def plotFieldOfView(self,quantity,AX):

       slicehalfwidth = self.rmax / 6.0

       if quantity == 'mass':
           self.galaxies['rtrunc_arcmin'] = (self.galaxies['r200']/self.galaxies['Da_p']) * \
                                             pangloss.rad2arcmin
           self.galaxies['rscale_arcmin'] = (self.galaxies['rs']/self.galaxies['Da_p']) * \
           pangloss.rad2arcmin
           for i in range(len(self.galaxies['x'])):
               trunc = plt.Circle([self.galaxies['x'][i], self.galaxies['y'][i]],radius=self.galaxies['rtrunc_arcmin'][i],fill=True,fc="b",alpha=0.05)
               AX.add_patch(trunc)
           for i in range(len(self.galaxies['x'])):
               core = plt.Circle([self.galaxies['x'][i], self.galaxies['y'][i]],radius=self.galaxies['rscale_arcmin'][i],fill=True,fc="r",alpha=0.5)
               AX.add_patch(core)
           plt.title('Halo Mass')

       elif quantity == 'kappa':
           plt.scatter(self.galaxies['x'], self.galaxies['y'], c='r', marker='o', s=(self.galaxies['kappa'])*30000)
           plt.title('Convergence')

       elif quantity == 'mu':
           plt.scatter(self.galaxies['x'], self.galaxies['y'], c='g', marker='o', s=((self.galaxies['mu']-1.0)*3E4))
           plt.title('Magnification')

       elif quantity == 'stellarmass':
           plt.scatter(self.galaxies['x'], self.galaxies['y'], c='y', marker='o', s=(np.log(self.galaxies['Mstar'])/2),edgecolor = 'none' )
           plt.title('Stellar Mass')

       elif quantity == 'light':
           plt.scatter(self.galaxies['x'], self.galaxies['y'], c='y', marker='o', s=(2**(25-self.galaxies['mag'])),edgecolor = 'none' )
           plt.title('Galaxy Light')

       else:
           raise "Lightcone plotting error: unknown quantity "+quantity

       # Lightcone boundary and centroid:
       circ = plt.Circle(self.xc,radius=self.rmax,fill=False,linestyle='dotted')
       AX.add_patch(circ)
       AX.plot([self.xc[0]],[self.xc[1]], c='k', marker='+',markersize=20,markeredgewidth=1)
       axlimits = [self.xc[0]-self.rmax-0.1,self.xc[0]+self.rmax+0.1,self.xc[1]-self.rmax-0.1,self.xc[1]+self.rmax+0.1]
       AX.axis(axlimits)

       # Show slice:
       plt.axvline(x=slicehalfwidth, ymin=axlimits[2], ymax=axlimits[3],color='black', ls='dotted')
       plt.axvline(x=-slicehalfwidth, ymin=axlimits[2], ymax=axlimits[3],color='black', ls='dotted')

       # Labels:
       plt.xlabel('x / arcmin')
       # plt.ylabel('y / arcmin')

       return

# ----------------------------------------------------------------------------

    def plotLineOfSight(self,quantity,AX):

       # Only plot a subset of points, in a slice down the middle of
       # the light cone:
       subset = self.galaxies[np.abs(self.galaxies['x']) < 0.3]

       # Point positions:
       z = subset['z']
       y = subset['y']

       # Plot the points:
       if quantity == 'mass':
           size = (10.0**(subset['Mh']-11.0))
           plt.scatter(z, y, c='k', marker='o', s=size, edgecolor='none' )
           plt.title('Line-of-sight Halo Mass')

       elif quantity == 'kappa':
           size = (subset['kappa']*30000.0)
           plt.scatter(z, y, c='r', marker='o', s=size, edgecolor='k' )
           plt.title('Line-of-sight Convergence')

       elif quantity == 'mu':
           size = ((subset['mu']-1.0)*3E4)
           plt.scatter(z, y, c='g', marker='o', s=size, edgecolor='k' )
           plt.title(r'Line-of-sight Magnification $(\mu - 1)$')

       elif quantity == 'stellarmass':
           size = (np.log(subset['Mstar'])/2.0)
           plt.scatter(z, y, c='y', marker='o', s=size, edgecolor='none' )
           plt.title('Line-of-sight Stellar Mass')

       elif quantity == 'light':
           size = (2**(25-subset['mag']))
           plt.scatter(z, y, c='y', marker='o', s=size, edgecolor='none' )
           plt.title('Line-of-sight Galaxy Light')

       else:
           raise "Lightcone plotting error: unknown quantity "+quantity

       # Axis limits:
       zmax = max(self.galaxies['z'].max(),self.zs)
       AX.axis([0,zmax+0.1,-self.rmax-0.1,self.rmax+0.1])

       # Labels:
       plt.xlabel('redshift z')
       plt.ylabel('y / arcmin')

       # Add lines marking source and lens plane, and optical axis:
       plt.axvline(x=self.zl, ymin=0, ymax=1,color='black', ls='dotted',label='bla')
       plt.axvline(x=self.zs, ymin=0, ymax=1,color='black', ls='dotted')
       plt.axhline(y=0.0, xmin=0.0, xmax=zmax, color='black', ls='dashed')

       return

# ----------------------------------------------------------------------------
# Plotting cumulative sum of parameters with redshift

    def plotContributions(self,quantity,output):

       plt.clf()

       # Point positions:
       zmax = self.zs+0.1
       z = np.linspace(0.0,zmax,100)
       # Plot the points:
       if quantity == 'mass':
           contr = np.zeros(len(z))
           for i in range(len(z)):
               galaxies = self.galaxies[self.galaxies['z'] <= z[i]]
               size = galaxies['Mhalo_obs']
               contr[i] = np.sum(size)
           plt.plot(z, contr)
           plt.title('Cumulative Sum of Halo Mass')

       elif quantity == 'kappa':
           contr = np.zeros(len(z))
           for i in range(len(z)):
               galaxies = self.galaxies[self.galaxies['z'] <= z[i]]
               size = galaxies['kappa']
               contr[i] = np.sum(size)
           plt.plot(z, contr)
           plt.title(r'Cumulative Sum of $\kappa_h$')

       elif quantity == 'mu':
           contr = np.zeros(len(z))
           for i in range(len(z)):
               galaxies = self.galaxies[self.galaxies['z'] <= z[i]]
               size = galaxies['mu']
               contr[i] = np.sum(size)
           plt.plot(z, contr)
           plt.title(r'Cumulative Sum of $\mu_h$')

       elif quantity == 'stellarmass':
           contr = np.zeros(len(z))
           for i in range(len(z)):
               galaxies = self.galaxies[self.galaxies['z'] <= z[i]]
               size = galaxies['Mstar_obs']
               contr[i] = np.sum(size)
           plt.plot(z, contr)
           plt.title('Cumulative Sum of Stellar Mass')

       else:
           raise "Lightcone plotting error: unknown quantity "+quantity

       # Axis limits:
       zmax = max(self.galaxies['z'].max(),self.zs+0.1)

       # Labels:
       plt.xlabel('Redshift z')

       plt.gcf().tight_layout()

       if output != None:
           pangloss.rm(output)
           plt.savefig(output,dpi=300)

       return
# ----------------------------------------------------------------------------

    def plot(self,var='kappa',output=None):

       plt.clf()

       # Panel 1: Galaxy positions:
       ax1 = plt.subplot(3,3,(1,4), aspect ='equal')
       self.plotFieldOfView('light',ax1)

       # Panel 2: Halo mass distributions:
       ax2 = plt.subplot(3,3,(2,5), aspect ='equal')
       self.plotFieldOfView('mass',ax2)

       # Panel 3: Kappa contributions:
       ax3 = plt.subplot(3,3,(3,6), aspect ='equal')
       self.plotFieldOfView(var,ax3)

       # Lower panel: View along redshift axis
       ax4 = plt.subplot(3,3,(7,9))
       self.plotLineOfSight(var,ax4)

       if output != None:
           pangloss.rm(output)
           plt.savefig(output,dpi=300)

       return None

# ----------------------------------------------------------------------------


#=============================================================================
