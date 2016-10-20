"""
Line of sight mass reconstruction in the Universe.
"""

# These are the static paths for data.
# They are fetched via data_fetcher.py

from os import environ

PANGLOSS_DIR = environ['PANGLOSS_DIR']
DATA_DIR = PANGLOSS_DIR + "/data"
CALIB_DIR = PANGLOSS_DIR + "/calib"

GAMMA_1_FILE = DATA_DIR + "/GGL_los_8_0_0_N_4096_ang_4_rays_to_plane_37_f.gamma_1"
GAMMA_2_FILE = DATA_DIR + "/GGL_los_8_0_0_N_4096_ang_4_rays_to_plane_37_f.gamma_2"
KAPPA_FILE = DATA_DIR + "/GGL_los_8_0_0_N_4096_ang_4_rays_to_plane_37_f.kappa"
GUO_FILE = DATA_DIR + "/GGL_los_8_0_0_0_0_N_4096_ang_4_Guo_galaxies_on_plane_27_to_63.images.txt"

CATALOG_EXAMPLE = CALIB_DIR + "/Millennium/catalog_example.txt"
KAPPA_EXAMPLE = CALIB_DIR + "/Millennium/kappa_example.fits"
HALO_MASS_REDSHIFT_CATALOG = CALIB_DIR + "/SHMR/HaloMassRedshiftCatalog.pickle"

# note: ordering of these imports is important because there is a dependency tree between them.

from lightcone import *
from wlmap import *
from kappamap import *
from shearmap import *
from grid import *
from pdf import *
from shmr import *
from catalog import *
from foreground import *
from background import *
from plotting import *
from config import *
from io import *
from lensing import *
from scalingrelations import *
from miscellaneous import *

# Matt Auger's classes.
from distances import *
from ndinterp import *

