## aMasSED Code ##
__A Python based SED fitting code__

__Pre-requisites:__
* Python v2.5+
* scipy v0.9+
* numpy 
* Atpy:
	* ascii_table - *if using ascii format input/outputs*
	* PyFITS - *if using FITS format input/outputs*


__Order of use:__

makeSEDs -> processSEDs -> matchSEDs

To use non default parameter file from command line, use e.g.:

`% python makeSEDs -p "parameter file"`

