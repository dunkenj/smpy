from os.path import join as pjoin

# Format expected by setup.py and doc/source/conf.py: string of form "X.Y.Z"
_version_major = 0
_version_minor = 1
_version_micro = ''  # use '' for first of series, number for 1 and above
_version_extra = 'dev'
#_version_extra = ''  # Uncomment this for full releases

# Construct full version string from these.
_ver = [_version_major, _version_minor]
if _version_micro:
    _ver.append(_version_micro)
if _version_extra:
    _ver.append(_version_extra)

__version__ = '.'.join(map(str, _ver))

CLASSIFIERS = ["Development Status :: 3 - Alpha",
               "Environment :: Console",
               "Intended Audience :: Science/Research",
               "License :: OSI Approved :: MIT License",
               "Operating System :: OS Independent",
               "Programming Language :: Python",
               "Topic :: Scientific/Engineering"]

# Description should be a one-liner:
description = "SMpy: Flexible stellar pops., masses and mock observations with python"
# Long description will go up on the pypi page
long_description = """
Stellar populations and Masses with Python
==========================================

This package contains Python software designed for building and processing composite stellar populations in a simple but flexible manner. It allows for easy synthetic photometry to be produced for single models or large suites of models.

The code makes use of the `Astropy <https://astropy.readthedocs.org>`_ module throughout and therefore allows for easy conversion of physical units and a wide range of allowed cosmologies.

Currently supported simple stellar population models are:

1. `Bruzual & Charlot 2003 <http://www.bruzual.org/bc03/Updated_version_2012/>`_
2. `BPASS V1 & V2 <http://bpass.auckland.ac.nz/>`_ 


License
=======
``SMpy`` is licensed under the terms of the MIT license. See the file
"LICENSE" for information on the history of this software, terms & conditions
for usage, and a DISCLAIMER OF ALL WARRANTIES.

All trademarks referenced herein are property of their respective holders.

Copyright (c) 2015--, Kenneth Duncan
"""

NAME = "astro-smpy"
MAINTAINER = "Kenneth Duncan"
MAINTAINER_EMAIL = "duncan@strw.leidenuniv.nl"
DESCRIPTION = description
LONG_DESCRIPTION = long_description
URL = "http://github.com/dunkenj/smpy"
DOWNLOAD_URL = ""
LICENSE = "MIT"
AUTHOR = "Kenneth Duncan"
AUTHOR_EMAIL = "duncan@strw.leidenuniv.nl"
PLATFORMS = "OS Independent"
MAJOR = _version_major
MINOR = _version_minor
MICRO = _version_micro
VERSION = __version__
PACKAGES = ['smpy',
            'smpy.tests']
PACKAGE_DATA = {'smpy': [pjoin('data', '*')]}
REQUIRES = ["numpy", "scipy", "h5py", "astropy", "six"]

