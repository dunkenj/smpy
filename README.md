## aMasSED Code 
__A Python based SED fitting code__

---

__Pre-requisites:__

* Python v2.5+  
* scipy v0.9+  
* numpy  
* Atpy + either/both of the following:  
	* asciitable - *if using ascii format input/outputs*  
	* PyFITS - *if using FITS format input/outputs*  

__Basic usage:__

	% python makeSEDs.py
	% python processSED.py
	% python matchSEDs.py


To use non default parameter file from command line, use e.g.:  
`% python makeSEDs.py -p "parameter file"`

Quiet mode - _Suppress some of the screen printouts_:  
`% python processSEDs.py -q`


---


### Instructions ###

#### Using the code ####

The code is split into three seperate parts which must be run individually and in sequence in order to do the SED fitting. These three parts are:

1.  makeSEDs.py:

    This part reads in and processes the Bruzual & Charlot SSP models before generating composite SEDs for a range of ages, star formation histories, dust extinction values and metallicities. For this grid of SEDs, a range of other parameters are calculated and stored for later use.
        
2.  processSEDs.py:

    The second part takes the grid of SEDs and redshifts them to each of the desired redshift steps before convolving each of the spectra through the photometric filters described in the filter set.
        
3.  matchSEDs.py:

    The final part of the code is the part which does the SED fitting in order to derive the stellar masses, stellar population parameters or rest-frame absolute magnitudes. Currently, this part of the code has two setttings:
       
    * Simple - what should be used in most circumstances and when testing the code. In this mode, the normalisations for each template are calculated analytically before the best-fitting template is selected. _Currently this does not include full error analysis._
    
    * Full PDF calculation - this mode calculates the full mass likelihood distribution for a range of stellar masses while marginalising over the other population parameters. It outputs a likelihood distribution for the mass, UV continuum slope and rest-frame UV magnitude. It doesn't include any improvements in the output of other rest-frame magnitudes so the simple mode will suffice.
    
    The desired mode is chosen by setting the parameter `calc_pdf = True` for the full PDF mode or `calc_pdf = False` for the simple mode. _See also the explanation of the full parameter file below._
        
#### The parameter file ####

__Make Parameters__

* `ssp_input` - Path to the set of Bruzual & Charlot models being used. Can be either 'bc03' or 'cb07', Chabrier ('chab') or Salpeter ('salpeter') IMF, And either high or low resolution models (default is low). Explore the 'ssp/' to see which models are where and what is available.

* `ssp_output` - Path for the binary output.

* `tau` - Python list containing the desired star-formation history timescales _in Gyr_ for an exponentially declining or increasing taus (= negative taus).

* `tg` - Python list of the desired stellar population ages (since the onset of star-fomation) _in Gyr_. Limited to the 221 ages of the Bruzual & Charlot models and recommended to use ages of no less than 5 to 10 Myr. Only models less than the age of the universe at a given redshift will be included in the fitting.

* `tauv` - Dust extinction values. In units of _Av_ if using Calzetti or in units of optical depth `tau` if using Charlot and Fall.

* `dust_model` - 'calzetti' if using Calzetti 2000 dust extinction or 'charlot' if using Charlot & Fall.

* `mu` - Only used when `dust_model = 'charlot'`, fraction of dust extinction which comes from the wider ISM. Default = 0.3.

* `epsilon` - Fraction gas returned to the ISM which is recycled into new star-formation. Currently works correctly __only for declining models__. As this parameter is uncertain and not constrained in this context, default = 0. 

* `metallicities` - Python list of indices for the desired Bruzual & Charlot metallicities, between 0 and 5 (6 in total). __Placing a minus sign in front will turn on nebular emission, and can be either/or, e.g.__ `metallicities = [-1,1]`.

* `add_nebular` - Boolean determining the inclusion of nebular emission. May be deprecated but not sure so leave `= True` and follow the practice outlined above for the inclusion of nebular emission.

* `neb_file` - Path to the nebular continuum and line emission values for inclusion in the SEDs. Leave as default.

* `fesc` - Lyman continuum escape fraction determining the strength of nebular emission. `fesc = 1.` is equal to zero nebular emission. `fesc = 0.2` as default is consistent with measurements of escape fraction and nebular emission at high redshifts. 

__Process Parameters__

* `synmag_output` - Path for the binary output of the data cube.

* `filt_dir` - Directory containing the filter response curves through which the SEDs will be convolved.

* `filt_names` - Format for the filter file names.

*  `zspacing` - 'linear' or 'log' spacing for the desired redshift steps.

* `zmin` - If `zspacing = 'linear'`, then zmin must be = 0. If `zspacing = 'log'` then the lowest non-zero redshift step will be 10**zmin.

* `zmax` - If `zspacing = 'linear'`, then zmax is the largest desired redshift. If `zspacing = 'log'` then the highest redshift step will be 10**zmax.

* `zstep` - Only used if , the desired step sized between redshift bins.

* `n_zsteps` - Only used if `zspacing = 'log'`, the number redshift steps desired (total number will always be 1 larger than this value).

* `madau` - Boolean for the inclusion of Lyman continuum and lyman absorption due to hydrogen in the intergalactic medium. Turning off will not affect low redshifts but is crucial at high redshift. Leave as default (True) unless there is a specific need to turn off.

__Match Parameters__

* `input_catalog` - Path to the input galaxy catalog to which SEDs will be fit. See the Catalog section below for more details on the required formatting and data of the input catalogs.

* `input_format` - File format for the input catalog. As discussed in the next section, can be any format which the installed version of atpy can interpret.

* `ID_col` - Name of catalog column containing the ID numbers of the galaxies.

* `z_col` - Name of the catalog column containing the redshifts which should be used in the fitting.

* `flux_col_end` - String containing the suffix used to denote flux columns in the input catalog. E.g. `WFC3_F160W_FLUX` -> 'FLUX'. Not case sensitive. _Fluxes must be in units of micro Janskys_.

* `fluxerr_col_end` - Same as above but for the flux error columns. E.g. `WFC3_F160W_FLUXERR` -> 'FLUXERR'

* `filts_used` - List containing the indices of the filter columns which should included in the SED fitting.

* `flux_corr` - Correction required to convert input fluxes into total fluxes. Will depend on the input catalog, if the input fluxes are already total fluxes or similar then `flux_corr = 1.0`.

* `fit_mode` - _Deprecated: Fitting mode desired, to either fluxes or magnitudes. Now the code will always fit to fluxes._

* `calc_pdf` - Boolean which determines the whether the simple of full PDF fitting is done. Default = False. See the section earlier for the differences.

* `calc_mode` - Calculate the 'mode' mass when doing the fitting.

* `include_rest` - Include rest-frame absolute magnitudes in the output catalog.

* `mode_mass_percentage` - Only used if `calc_mode = True`. The top percentage of best-fitting templates to include when calculating the mode mass.

* `muv_max, muv_min, muv_bins` - Bright and faint limits and number of steps to use when constructing the M_UV likelihood distribution.

* `mass_max, mass_min, mass_bins` - The high and low mass limits to iterate over when constructing the marginalised stellar mass likelihood and the number of steps.

* `beta_max, beta_min, beta_bins` - The maximum and minimum betas and number of steps to use when constructing the Beta likelihood distribution.

* `output_name` - Path for the desired output catalog.

* `table_format` - Table format desired for output catalog. As with the input catalog, can be any format the installed version of Atpy is set up for.

#### Input Catalog ####

Because the code uses external modules for the input and output of catalogs, it is relatively flexible with regards to formatting. It will work with any of the formats the installed version of Atpy will work with.

Therefore, to use .fits format catalogs, the PyFITS package must also be installed. Similarly, to use ascii tables, the asciitables package must be installed.

The input catalog must contain the following:

* A column containing IDs for each of the galaxies so that the output can be easily matched. Identified with the `ID_col` parameter.

* A column containing a redshift for each of the galaxies, be it spectroscopic or photometric. Identified with the `z_col` parameter.

* At least one flux column with corresponding flux error in units of µJy.

__Please note: The order of the columns in the catalog MUST be the same as the order of the filters once the filter name array has been sorted. Hence it is best to number the filter files to match the order of the catalog.__



