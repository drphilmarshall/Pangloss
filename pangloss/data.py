from os import mkdir, path, remove
from traceback import print_exc
from requests import get

PANGLOSS_MODULE_DIR = path.dirname(path.realpath(__file__))
DATA_DIR = path.join(path.dirname(PANGLOSS_MODULE_DIR), 'data')
CALIB_DIR = path.join(path.dirname(PANGLOSS_MODULE_DIR), 'calib')
MILLENNIUM_DIR = path.join(CALIB_DIR, 'Millennium')
SHMR_DIR = path.join(CALIB_DIR, 'SHMR')
EXAMPLE_DIR = path.join(path.dirname(PANGLOSS_MODULE_DIR), 'example')

#files for demo notebooks
GAMMA_1_FILE = path.join(DATA_DIR, 'GGL_los_8_0_0_N_4096_ang_4_rays_to_plane_37_f.gamma_1')
GAMMA_2_FILE = path.join(DATA_DIR, 'GGL_los_8_0_0_N_4096_ang_4_rays_to_plane_37_f.gamma_2')
KAPPA_FILE = path.join(DATA_DIR, 'GGL_los_8_0_0_N_4096_ang_4_rays_to_plane_37_f.kappa')
GUO_FILE = path.join(DATA_DIR, 'GGL_los_8_0_0_0_0_N_4096_ang_4_Guo_galaxies_on_plane_27_to_63.images.txt')
CONFIG_FILE = path.join(EXAMPLE_DIR, 'example.config')

DEMO_FILES = [GAMMA_1_FILE, GAMMA_2_FILE, KAPPA_FILE, GUO_FILE, CONFIG_FILE]

#files for classic examples
CATALOG_EXAMPLE = path.join(MILLENNIUM_DIR, 'catalog_example.txt')
KAPPA_EXAMPLE = path.join(MILLENNIUM_DIR, 'kappa_example.fits')
HALO_MASS_REDSHIFT_CATALOG = path.join(SHMR_DIR, 'HaloMassRedshiftCatalog.pickle')

EXAMPLE_FILES = [CATALOG_EXAMPLE, KAPPA_EXAMPLE, HALO_MASS_REDSHIFT_CATALOG]

def fetch():
    """
    Downloads the data for both the demo notebooks and example configuration.
    """
    demo_data_url = 'http://www.slac.stanford.edu/~pjm/hilbert'
    calib_data_url = 'http://www.ast.cam.ac.uk/~tcollett/Pangloss/calib'

    for d in [DATA_DIR, CALIB_DIR, MILLENNIUM_DIR, SHMR_DIR]:
        if not path.exists(d):
            mkdir(d)

    for demo_file in DEMO_FILES:
        if not path.exists(demo_file):
            url = '{}/{}'.format(demo_data_url, path.basename(demo_file))
            download(url, demo_file)

    for example_file in EXAMPLE_FILES:
        if not path.exists(example_file):
            url = '{}/{}/{}'.format(calib_data_url, *example_file.split('/')[-2:])
            download(url, example_file)

def download(url, output):
    """
    Downloads the data for both the demo notebooks and example configuration.

    Note
    ----
    Does not provide progress feedback so be careful using this for HUGE downloads or SLOW network.

    Parameters
    ----------
    url : str
        Url to download from.
    output : str
        Path to write to.
    """
    print "Starting download \n" \
               "\t{}\n" \
               "\t >>> {}".format(url, output)
    try:
        r = get(url, stream=True)
        with open(output, 'wb') as f:
            for chunk in r.iter_content(chunk_size=1024):
                if chunk:
                    f.write(chunk)
        print "Done"
    except Exception:
        try:
            remove(output)
        except OSError:
            pass
        print_exc()
        msg = "Error downloading from '{}'. Please issue this on " \
              "https://github.com/drphilmarshall/Pangloss/issues".format(url)
        raise Exception(msg)
