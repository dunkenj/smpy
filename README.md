## aMasSED Code 
__A Python based SED fitting code__

---

__Pre-requisites:__

* Python v2.5+  
* scipy v0.9+  
* numpy  
* Atpy + either/both of the following:  
	* ascii_table - *if using ascii format input/outputs*  
	* PyFITS - *if using FITS format input/outputs*  

__Basic usage:__

	% python makeSEDs.py
	% python processSED.py
	% python matchSEDs.py


To use non default parameter file from command line, use e.g.:  
`% python makeSEDs.py -p "parameter file"`

Quiet mode - _Suppress some of the screen printouts_:  
`% python processSEDs.py -q`

### Instructions ###

__Using the code__
The code is split into three seperate parts which must be run individually and in sequence in order to do the SED fitting. These three parts are:
1. makeSEDs.py:
    This part reads in and processes the Bruzual & Charlot SSP models before generating composite SEDs for a range of ages, star formation histories, dust extinction values and metallicities. For this grid of SEDs, a range of other parameters are calculated and stored for later use.
2. processSEDs.py:
    The second part takes the grid of SEDs and redshifts them to each of the desired redshift steps before convolving each of the spectra through the photometric filters described in the filter set.
3. matchSEDs.py:
    The final part of the code is the part which does the SED fitting in order to derive the stellar masses, stellar population parameters or rest-frame absolute magnitudes. Currently, this part of the code has two setttings:
    * Simple - what should be used in most circumstances and when testing the code. In this mode, the normalisations for each template are calculated analytically before the best-fitting template is selected. _Currently this does not include full error analysis._
    * Full PDF calculation - this mode calculates the full mass likelihood distribution for a range of stellar masses while marginalising over the other population parameters. It outputs a likelihood distribution for the mass, UV continuum slope and rest-frame UV magnitude. It doesn't include any improvements in the output of other rest-frame magnitudes so the simple mode will suffice.
    The desired mode is chosen by setting the parameter `calc_pdf = True` for the full PDF mode or `calc_pdf = False` for the simple mode. _See the explanation of the parameter file later._