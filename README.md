# aMasSED Code - A Python based SED fitting code

## Pre-requisites:
* Python v2.5+
* scipy v0.9+
* numpy 
* Atpy:
	* ascii_table - *if using ascii format input/outputs*
	* PyFITS (For ) - *if using FITS format input/outputs*


Order of use:

makeSEDs -> processSEDs -> matchSEDs

To use non default parameter file from command line, use e.g.:

% python makeSEDs -p "parameter file"

