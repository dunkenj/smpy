"""
Input Options
-------------

model_path : Input path to generated models - hdf5 file
input_catalog : Input catalog path
input_format : Input catalog format, see astropy.Table
               documentation for available formats

z_col : Column name for source redshifts
ID_col : Column name for source IDs
flux_col_end : Suffix corresponding to flux columns
fluxerr_col_end : Suffix corresponding to flux error columns

filts_used : Indices of filters in model file to be used for fitting.
            If 'None', assumes all filters to be used.

"""

model_path = 'candels.goodss.models.savetest.hdf'
output_name = 'candels_test.cat'
input_catalog = 'data/CANDELS.GOODSS.example.cat'
input_format = 'ascii.commented_header'

z_col = 'Photo_z'
ID_col = 'ID'
flux_col_end = '_FLUX'
fluxerr_col_end = '_FLUXERR'

filts_used = None

"""
Fitting Options
---------------

fitting_mode :
    Desired option for fitting mode/outputs,
        'simple' - Simple single best-fit model
        'hist' - Additional mass and sfr pdf outputs
include_rest : Calculate rest-frame magnitudes for best-fit model (True/False)
ncpus : Number of parallel processes to use when fitting
flux_corr : Correction to convert input fluxes to total
flux_err : Additional fractional flux error added to all bands in quadrature

"""
fitting_mode = 'hist'
include_rest = True
ncpus = 4

zp_offsets = 'test_offsets.txt'
temp_err = None #'TEMPLATE_ERROR.v2.0.zfourge.txt'

flux_corr = 1
flux_err = 0.
nmin_bands = 14.

"""
Output Options
--------------

output_catalog : Path for output catalog
output_format : Output catalog format, see astropy.Table
               documentation for available formats
output_hdf : Path for output hdf5 file ('hist' mode only)
"""

output_hdf_path = 'candels_test.hdf'
output_catalog_path = 'candels_test.cat'
output_format = 'ascii.commented_header'
