author: Stefan Hilbert
last modified:   2012 June 8


Galaxy catalogs
----------------

The files

GGL_los_8_x_y_i_j_N_4096_ang_4_SA_galaxies_on_plane_27_to_63.images.txt
GGL_los_8_x_y_i_j_N_4096_ang_4_Guo_galaxies_on_plane_27_to_63.images.txt
GGL_los_8_x_y_i_j_N_4096_ang_4_Bower_galaxies_on_plane_27_to_63.images.txt

where x and y can range from 0 to 7, and i and j currently can range from 0 to 3 contain tables of galaxy image properties in a simulated 1 deg^2 portion of the sky. The first row contains (hopefully) self-explanatory descriptions of the galaxy properties stored in the respective columns. Each other row contains the properties for one galaxy image.

List of currently stored properties:

GalID               unique ID of source galaxy of galaxy image                          
HaloID              unique ID of main halo (can be galaxy, group, or cluster halo) the source galaxy is located in
SubhaloID           unique ID of subhalo the galaxy is centre of (if there is a subhalo)
Type                0 = central galaxy of the main halo, 1 = non-central galaxy of the main halo, but with subhalo, 2 = non-central galaxy without subhalo (at least for the MPA Garching galaxy model)
PlaneNumber         ID of lens plane in ray-tracing the galaxy is associated with       
z_spec              "spectroscopic" redshift of the galaxy (includes cosmological redshift and peculiar motion Doppler shift)
pos_0[rad]          first sky position coordinate in radian (can be positive and negative)
pos_1[rad]          second sky position coordinate in radian (can be positive and negative)
Dc_los[Mpc/h]       comoving line-of-sight distance to source galaxy in Mpc/h
M_Halo[M_sol/h]     virial mass of the main halo of galaxy in solar mass /h (sorry, still have to check which overdensity threshold is used)
M_Subhalo[M_sol/h]  gravitationally bound mass of the subhalo of galaxy in solar mass/h (note that the Bower model does not really know subhalos)
M_Stellar[M_sol/h]  stellar mass of galaxy
mag_SDSS_u          SDSS u-band apparent observer-frame AB magnitude (also contains magnification effects)
mag_SDSS_g          SDSS g-band apparent observer-frame AB magnitude (also contains magnification effects)
mag_SDSS_r          SDSS r-band apparent observer-frame AB magnitude (also contains magnification effects)
mag_SDSS_i          SDSS i-band apparent observer-frame AB magnitude (also contains magnification effects)
mag_SDSS_z          SDSS z-band apparent observer-frame AB magnitude (also contains magnification effects)
mag_J               2MASS J-band apparent observer-frame magnitude (also contains magnification effects)
mag_H               2MASS H-band apparent observer-frame magnitude (also contains magnification effects)
mag_K               2MASS K-band apparent observer-frame magnitude (also contains magnification effects)

For the lightcone reconstruction, people should ideally only use the sky positions and the magnitudes (all else could be considered cheating).

The data was created by ray-tracing through the Millennium Simulation (see, e.g., Hilbert et al. 2009, A&A 499, 31). The galaxy properties are computed with various versions of semianalytic models:

*_SA_galaxies_*   : De Lucia & Blaizot 2007, MNRAS 375, 2,
*_Guo_galaxies_*  : Guo et al. 2011, MNRAS 413, 101,
*_Bower_galaxies_*: Bower et al. 2006, MNRAS 370, 645.

The simulated galaxy image catalogs are restricted to galaxies with mag_SDSS_i <= 26 (note that is higher than it used to be). Furthermore, the galaxy models become incomplete for galaxies with stellar masses M_Stellar ~ 1e9 Msolar and below due to limited resolution of the simulation. So, the calalog is stellar-mass limited at low redshifts, and magnitude-limited at higher redshift. The magnitude limit is supposed to represent some sensible limit for real observations (and can be changed if needed).

All 16 files with different i and j, but the same x and y put together give a catalog for a contiguous 4 x 4 deg^2 simulated field.


Convergence and shear maps
----------------

The files

GGL_los_8_x_y_N_4096_ang_4_rays_to_plane_37_f.kappa
GGL_los_8_x_y_N_4096_ang_4_rays_to_plane_37_f.gamma_1
GGL_los_8_x_y_N_4096_ang_4_rays_to_plane_37_f.gamma_2

where x and y can range from 0 to 7 contain the convergence and shear components for sources at redshift z_S = 1.3857 on a mesh of 4096 x 4096 pixels in a simulated 4 x 4 deg^2 portion of the sky. The files are raw binary data of a row-major 2-dimensional array of 4-byte (little-endian) floats. In C code, this would look like:

float kappa  [4096][4096];
float gamma_1[4096][4096];
float gamma_2[4096][4096];

The angular position of a pixel kappa[ix][iy] can be computed by something like:

void
get_pixel_position(const int ix, const int iy, double *pos_x,  double *pos_x)
{
  const double   degree        = M_PI / 180.;
  const double   L_field       = 4.0 * degree;
  const int      N_pix_per_dim = 4096;
  const double   L_pix         = L_field / N_pix_per_dim;

  *pos_x = -0.5 * L_field  + (ix + 0.5) * L_pix;
  *pos_y = -0.5 * L_field  + (iy + 0.5) * L_pix;
}

If you need to interpolate between pixels:

double
interpolated_kappa(const double pos_x, const double pos_y)
{
  const double   degree        = M_PI / 180.;
  const double   L_field       = 4.0 * degree;
  const int      N_pix_per_dim = 4096;
  const double   inv_L_pix     = N_pix_per_dim / L_field;

  const double x = inv_L_pix * (pos_x + 0.5 * L_field) - 0.5; /* positions in "units" of pixel indices */
  const double y = inv_L_pix * (pos_y + 0.5 * L_field) - 0.5;

  const int ix = floor(x); /* pixel indices for lower left pixel for bilinear interpolation */
  const int iy = floor(y);

  const double px = x - ix; /* weights for bilinear interpolation */
  const double py = y - iy;

  if((0 <= ix) && (ix < N_pix_per_dim - 1) && (0 <= iy) && (iy < N_pix_per_dim - 1))
      return kappa[ix][iy] * (1. - px) * (1. - py) + kappa[ix + 1][iy] * px * (1. - py) + kappa[ix][iy + 1] * (1. - px) * py + kappa[ix + 1][iy + 1] * px * py;
  else
    return 0.;
}




