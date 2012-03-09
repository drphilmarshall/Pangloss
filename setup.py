from distutils.core import setup, Extension
from numpy.distutils.misc_util import get_numpy_include_dirs

numpy_inc = get_numpy_include_dirs()

setup(
    name="Pangloss",
    version="n/a",
    author="Tom Collett (IoA) & Phil Marshall (Oxford)",
    author_email="",
    packages=["pangloss"],
    ext_modules = [],
    url="",
    license="GPLv2",
    description="line of sight mass reconstruction",
    long_description="TBD",
    classifiers=[
        "Development Status :: 2 - Pre-Alpha",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: GNU General Public License (GPL)",
        "Operating System :: OS Independent",
        "Programming Language :: Python",
    ],
)
