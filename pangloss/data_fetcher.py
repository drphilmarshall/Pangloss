from os import mkdir, path, remove
from traceback import print_exc

from requests import get

from pangloss import CATALOG_EXAMPLE, DATA_DIR, GAMMA_1_FILE, GAMMA_2_FILE, GUO_FILE, \
    HALO_MASS_REDSHIFT_CATALOG, KAPPA_FILE, KAPPA_EXAMPLE, PANGLOSS_DIR

DEMO_DATA_URL = "http://www.slac.stanford.edu/~pjm/hilbert"
CALIB_DATA_URL = "http://www.ast.cam.ac.uk/~tcollett/Pangloss/calib"


def setup_data():
    """
    Downloads the data for both the demo notebooks and example configuration.
    """
    if not path.exists(DATA_DIR):
        try:
            mkdir(DATA_DIR)
        except OSError:
            print_exc()
            msg = "Error trying to create directory '{}', most likely your $PANGLOSS_DIR " \
                  "which is currently set to '{}' has not been set or is not correct. Please fix" \
                  " this and re-run.".format(DATA_DIR, PANGLOSS_DIR)
            raise Exception(msg)

    for demo_data in [GAMMA_1_FILE, GAMMA_2_FILE, GUO_FILE, KAPPA_FILE]:
        if not path.exists(demo_data):
            url = "{}/{}".format(DEMO_DATA_URL, path.basename(demo_data))
            download(url, demo_data)

    for calib_data in [CATALOG_EXAMPLE, HALO_MASS_REDSHIFT_CATALOG, KAPPA_EXAMPLE]:
        if not path.exists(calib_data):
            url = "{}/{}/{}".format(CALIB_DATA_URL, *calib_data.split('/')[-2:])
            download(url, calib_data)


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
               "\t --> {}"
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
        msg = "Error downloading from '{}'. Please email davidthomas5412@gmail.com if these " \
              "issues continue.".format(url)
        raise Exception(msg)
