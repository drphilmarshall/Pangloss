from setuptools import setup

setup(
    name="pangloss",
    version='0.0.1',
    author="Phil Marshal",
    author_email="pjm@slac.stanford.edu",
    url="https://github.com/drphilmarshall/pangloss",
    packages=["pangloss"],
    description="Line of sight weak lensing and inference.",
    long_description=open("README.md").read(),
    package_data={"": ["README.md", "LICENSE"]},
    include_package_data=True,
    classifiers=[
        "Development Status :: 4 - Beta",
        "License :: OSI Approved :: GNU General Public License (GPL)",
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "Operating System :: OS Independent",
        "Programming Language :: Python",
    ],
    install_requires=["astropy", "numpy", "scipy", "matplotlib", "requests"],
)