Pangloss
--------

Pangloss is a code for reconstructing all the mass within a light cone through the Universe. Understanding such a mass distribution is important for accurate time delay lens cosmography, and also accurate lens magnification estimation. This code was used for Collett et al. 2013; MNRAS ????????? "Reconstructing the lensing mas in the universe from photometric catalogue data"


This README is aimed at the collaboration of people invetsigating this idea. If you are interested, but not one of these people, you should contact tcollett@ast.cam.ac.uk & dr.phil.marshall@gmail.com and read the wiki to find out more about our plans. We are keen to hear from other researchers in the field who would like to coordinate with us. 

License
-------

Pangloss is distributed under the Gnu Public License (GPL), a copy of which is delivered with this code. Copyright Tom Collett & Phil Marshall.


Installation
------------

Add the following lines (or similar) to your .login file:

  setenv PANGLOSS_DIR ${WORK_DIR}/Pangloss
  setenv PYTHONPATH ${PANGLOSS_DIR}/Pangloss:${PYTHONPATH}
  setenv PATH ${PATH}:${PANGLOSS_DIR}

Then "import Pangloss" from python should just work, and the command line scripts should too. You will need to have the following packages installed as well:

  atpy,asciitable,pyfits
  numpy,scipy,getopt,matplotlib

Most (if not all) are available via pip and/or homebrew.


Test Data:
----------
We include a 1x1 square degree of millenium simulation sky, and its associated reay traced convergence and shear maps.

Additional galaxy catalogs and ray-traced convergence maps from the Millenium Simulation are available from Stefan Hilbert on request.

Survey.params:
--------
Define your survey in survey.params so the code knows what you want it to do.


Classes:
--------

The primary "lightcone" class can be used for drilling out lines of sight from the Millenium catalogues of Hilbert et al, and returns estimates of convergence given various assumptions.

The "kappamap" class allows Hilbert's convergence maps, obtained by full-on ray-tracing, to be loaded and queried by position - this provides truth at each sky position.

distances.py is a class calculating cosmologically useful quantities i.e comoving distances. Written by Matt Auger.

The "grid" class creates a framework of redshift slices upon which useful quantities are pre-calculated/

Functionals: 
--------

The other files in pangloss contain useful functions.

LensingProfiles.py contains functions needed to calculate the convergence and shear caused by various halo profiles.

Relations.py contains the neto and maccio halo mass to concentration relation and the behroozi stellar mass to halo mass relation. 

LensingFunc.py contains functions to calculate the lensing scale factor, beta, and critical density, sigma_crit. # I think this is obsolete and never used. TC

MakeCones.py will drill out and save lightcones, this is rather slow and allows the process to be done only once, they can then be re-loaded as needed.

AnalyseCones.py will take cones and reconstruct their kappa_halo and mu_halo.

MakeCalibrationGuide.py takes the result of AnalyseCones applied to calibration lines of sight with known kappa and mu and then generates a joint distribution P(kappa_ext, mu_ext, kappa_halos, mu_halos)


lightcone:
----------

The user can either enter a position around which to draw a light cone, or can let Pangloss choose a lens, given redshift and stellar mass selections (currently hard coded, but can easily be changed), minor alterations could allow other lens selection proceedures. Feel free to add these, preferably as a function rather than hard coded.

Once a lens and lightcone has been selected the code calculates kappa_ext and gamma_ext for each halo in the catalogue. These are then weighted to create a kappa_keeton for each halo(see Keeton 2003, http://iopscience.iop.org/0004-637X/584/2/664/pdf/56371.web.pdf), which can be summed to give a kappa_keeton_total. Kappa_keeton is an approximation, but may well be good enough for our purposes.

Our code currently only assumes NFW halos, and ignores starlight information.


kappamap:
---------

Stefan's convergence maps can be read in and interpolated to a particular sky position. The map data are available as raw binary, but the kappamap class writes a map out in FITS format every chance it gets.
The maps are for convergence at a lens plane at zd=0.5, with source at redshift 1.6 (as in the case of B1608).





TC & PJM
