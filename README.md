## SMpy: Stellar populations and Masses with Python

This module is designed for building and processing composite stellar populations in a simple but flexible manner. It allows for easy synthetic photometry to be produced for single models or large suites of models.

The code makes use of the `Astropy <https://astropy.readthedocs.org>`_ module throughout and therefore allows for easy conversion of physical units and a wide range of allowed cosmologies.

Currently supported simple stellar population models are:

1. [Bruzual & Charlot 2003](http://www.bruzual.org/bc03/Updated_version_2012/)
2. [BPASS V1 & V2](http://bpass.auckland.ac.nz/)

These files are not included in the PyPI source distribution, however a set of unpacked Bruzual & Charlot models are included within the `scripts/data` folder on the github page.

The core code now includes only the code used to build and process the composite stellar populations. In the very near future, example code equivalent to the old `makeSEDS.py` and `processSEDs.py` scripts will be included in the `scripts` folder. Similarly, an updated version of `matchSEDS.py` will also be included shortly. _All functionality of the previous code version will therefore still be included_.


### Requirements

* Numpy version 1.6
* scipy
* Astropy version 1.0.1 or later
* h5py - for writing large model arrays to hdf files
* six - for Python 2 to 3 compatibility (_in progress_)

### Installation 
The code can now be installed as a module following standard python practice:

    git clone git://github.com/dunkenj/smpy.git
    python setup.py install    

