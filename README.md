# Pangloss

Pangloss is software module for reconstructing all the mass within a light cone
through the Universe.  Understanding complex mass distributions like
this is important for accurate time delay lens cosmography, and also for
accurate lens magnification estimation. It aspires to use all available
data in an  attempt to make the best of all mass maps. However, this
problem is sufficiently difficult that we are mostly spending our time
cultivating a garden of toy models and baby steps.

If you are interested in the ongoing investigation of this interesting
problem, you should:
* Read our first paper, [_"Reconstructing the lensing mass in the universe from photometric catalogue data;"_](http://arxiv.org/abs/1303.6564)
* [Browse the wiki](https://github.com/drphilmarshall/Pangloss/wiki) 
to see what our current plans are;
* [Fork the code](https://github.com/drphilmarshall/Pangloss) and play around with it;
* Contact us (Phil Marshall, dr.phil.marshall@gmail.com, and Tom
Collett, tcollett@ast.cam.ac.uk). We are keen to hear from other
researchers in the field who would like to coordinate and/or collaborate
with us!

## License

Pangloss is distributed under the Gnu Public License (GPL) v2, a copy of
which is delivered with this code. This basically means you can do
anything you like with it, as long as you keep propagating the license
information.

However, when reporting on your use of Pangloss, we do ask that you include:
* A citation to: **Collett et al. (2013), MNRAS accepted, arxiv/1303.6564**
* An acknowledgment including a link to this repo, such as:
_"This work made use of the Pangloss code, written by Tom Collett and
Phil Marshall, which is freely available at https://github.com/drphilmarshall/Pangloss."_

## Installation

After cloning or forking the repository, 
add the following lines (or similar) to your .login file:

    setenv PANGLOSS_DIR ${WORK_DIR}/Pangloss
    setenv PYTHONPATH ${PANGLOSS_DIR}/Pangloss:${PYTHONPATH}
    setenv PATH ${PATH}:${PANGLOSS_DIR}

Then "import pangloss" from python should just work, and the command
line scripts should too. However, you will need to have the following packages
installed as well:

    import atpy,asciitable,pyfits
    import numpy,scipy,getopt,matplotlib

Most (if not all) are available via pip and/or homebrew.


## Example use

Right now, there are three main scripts that carry out Pangloss's
functions:

* **[Drill.py](https://github.com/drphilmarshall/Pangloss/blob/master/Drill.py)**: 
from an input catalog of galaxies (either real or
simulated), drill out a narrow lightcone (or set of lightcones) centred
on some sky position of interest, and save the resulting smaller
catalog (or catalogs).
* **[Reconstruct.py](https://github.com/drphilmarshall/Pangloss/blob/master/Reconstruct.py)**:
read in a lightcone catalog (or set of catalogs), and
assign dark matter mass to the galaxies within, probabilistically. Each
lightcone yields its own Pr(kappah|D), where kappah is the total
gravitational lensing convergence at the centre of the lightcone - this
quantity is of great relevance to the strong lens time delay distance in
particular.
* **[Calibrate.py](https://github.com/drphilmarshall/Pangloss/blob/master/Calibrate.py)**:
 kappah is (at the moment) a fairly crude approximation
to the "true" convergence kappa - but we have shown that one can recover
a better approximation by considering Pr(kappah,kappa|C), where C is
large set of lightcone catalogs drawn from a cosmological simulation.
The resulting PDF Pr(kappah|D,C) contains the assumption that this
simulation is an accurate representation of our real Universe - but
we're working towards relaxing this.

They all take, as their sole input, the same configuration file, an
[example of which is given here](https://github.com/drphilmarshall/Pangloss/blob/master/pangloss/example_catalog.txt).

In the calib directory we include the means to obtain the halo catalog 
for a 1x1 square degree patch of Millennium Simulation sky, and its
associated ray traced convergence map, for making the calibration
lightcones. A small mock observed galaxy catalog is included in the
example directory, for testing.

To get the test data, please do the following:

    cd calib
    Fetch.csh

The Fetch script uses wget. 
Additional galaxy catalogs and ray-traced convergence maps from the
Millenium Simulation are available from Stefan Hilbert on request.

You should then be able to execute the following example analysis:

    cd example
    Drill.py example.config
    Reconstruct.py example.config
    Calibrate.py example.config
    
Analysing from start to finish will take some time. Be patient! 
(Drill ~2 mins, Reconstruct ~10 mins)

For more details of what the Pangloss scripts are doing, [start reading the code here.](https://github.com/drphilmarshall/Pangloss/wiki/Code-description)

## PDFs based on overdensities

The [original Pangloss](https://github.com/drphilmarshall/Pangloss) has been modified so that it can be used to find convergence and magnification pdfs for any given group of fields based on comparing overdensities of the observed fields with simulation lightcones.

