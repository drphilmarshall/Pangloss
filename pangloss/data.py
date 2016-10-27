from os import mkdir, path, remove
from traceback import print_exc

from requests import get

pangloss_module_dir = path.dirname(path.realpath(__file__))
data_dir = path.join(path.dirname(pangloss_module_dir), 'data')
calib_dir = path.join(path.dirname(pangloss_module_dir), 'calib')
millennium_dir = path.join(calib_dir, 'Millennium')
shmr_dir = path.join(calib_dir, 'SHMR')

gamma_1_file = path.join(data_dir, 'GGL_los_8_0_0_N_4096_ang_4_rays_to_plane_37_f.gamma_1')
gamma_2_file = path.join(data_dir, 'GGL_los_8_0_0_N_4096_ang_4_rays_to_plane_37_f.gamma_2')
kappa_file = path.join(data_dir, 'GGL_los_8_0_0_N_4096_ang_4_rays_to_plane_37_f.kappa')
guo_file = path.join(data_dir, 'GGL_los_8_0_0_0_0_N_4096_ang_4_Guo_galaxies_on_plane_27_to_63'
                                '.images.txt')
demo_data_files = [gamma_1_file, gamma_2_file, kappa_file, guo_file]

catalog_example = path.join(millennium_dir, 'catalog_example.txt')
kappa_example = path.join(millennium_dir, 'kappa_example.fits')
halo_mass_redshift_catalog = path.join(shmr_dir, 'HaloMassRedshiftCatalog.pickle')

calib_data_files = [catalog_example, kappa_example, halo_mass_redshift_catalog]

def fetch():
    """
    Downloads the data for both the demo notebooks and example configuration.
    """
    demo_data_url = 'http://www.slac.stanford.edu/~pjm/hilbert'
    calib_data_url = 'http://www.ast.cam.ac.uk/~tcollett/Pangloss/calib'

    for d in [data_dir, calib_dir, millennium_dir, shmr_dir]:
        if not path.exists(d):
            mkdir(d)

    for demo_data in demo_data_files:
        if not path.exists(demo_data):
            url = '{}/{}'.format(demo_data_url, path.basename(demo_data))
            download(url, demo_data)

    for calib_data in calib_data_files:
        if not path.exists(calib_data):
            url = '{}/{}/{}'.format(calib_data_url, *calib_data.split('/')[-2:])
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
