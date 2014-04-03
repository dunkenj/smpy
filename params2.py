"""
SSP Paramaters
make.py

"""
ssp_input = 'ssp/bc03/chab/lr/'               # Input SSP folder
#ssp_input = 'ssp/cb07/'
ssp_output = '../MassData/temp_fitting.npz' # Output filename
#ssp_output = '/data/stellarmasses/models/200912/ISM.jamie.npz' # Output filename

#tau = [0.01,0.025,0.05,0.1,      # Star formation timescales (Gyr)
#       0.25,0.5,0.8,1.0,1.5,
#       2.5,3.5,5.5,10.0,13.7]
tau = [2]
#tau = [-10.,-5.,-2.5,-1.0,-0.5,-0.25,0.05,0.25,0.5,1.0,2.5,5.,10.,1000]

tg = [0.005,0.0100,0.0158,0.0251,0.0398,0.0631,0.1000,0.1585,
      0.2512,0.3,0.4,0.5,0.6310,0.7,0.8,1.0,1.13,1.28,1.43,1.6]                 #Model ages (Gyr)

#tg = [0.0100,0.0158,0.0251,0.0398,0.0631,0.08,0.1000,0.1585,0.2,0.2512,0.3981,
#      0.5,0.6,0.7,0.8,0.9,1.0,1.13,1.28,1.43,1.6,2.0,2.3,2.5,2.7,3.0]                 #Model ages (Gyr)

#tauv = [ 0.  ,  0.05,  0.1 ,  0.15,  0.2 ,  0.25,  0.3 ,  0.35,  0.4 ,
#        0.45,  0.5 ,  0.55,  0.6 ,  0.65,  0.7 ,  0.75,  0.8 ,  0.85,
#        0.9 ,  0.95,  1.  ,  1.05,  1.1 ,  1.15,  1.2 ,  1.25,  1.3 ,
#        1.35,  1.4 ,  1.45,  1.5 ,  1.55,  1.6 ,  1.65,  1.7 ,  1.75,
#        1.8 ,  1.85,  1.9 ,  1.95,  2.  ,  2.05,  2.1 ,  2.15,  2.2 ,
#        2.25,  2.3 ,  2.35,  2.4 ,  2.45,  2.5 ,  2.55,  2.6 ,  2.65,
#        2.7 ,  2.75,  2.8 ,  2.85,  2.9 ,  2.95,  3.  ]

#tauv = [0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,
#        1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0]#4.5,5.0,5.5,6.0] # Dust Total Attenuation Values
#tauv = [0.,0.4,0.6,0.8,1.0]
tauv = [1.0] 

dust_model = 'calzetti'                # 'charlot' - Charlot & Fall or 'calzetti' - Caltzetti 2000
mu = 0.3                              # Fraction of tauv from ISM
epsilon = 0.                     # Gas recycling parameter
metallicities = [2]       # Metallicities to use

csf = False                           # Set to constant SFR (overrides, tau values)

#Nebular emission
add_nebular = True                    # Include nebular emission
neb_file = 'nebular_emission.dat'     # Nebular emission path
fesc = 0.2                            # Lyman continuum escape fraction

"""
Synthetic Magnitudes Parameters
process.py

"""

synmag_output = '/data/stellarmasses/synmags/200912/temp_fitting_neb.npz'

# Filter Directory and Naming Convention - Assumes Ordered Numerically
filt_dir = '/data/candels/catalogs/Filters/GOODS-S_18_FilterCurves/'
#filt_dir = '/data/candels/catalogs/Filters/ForNina/'
#filt_dir = '/home/ppxkd/astroraid/work/data/CANDELS_filters/test2/'
filt_names = 'Filter*.txt'
#filt_names = 'CANDELS_GOODS-S.filter*_test2.v1.txt'
#filt_dir = 'filters/GOODS-S_Alice/'
#filt_names = 'X*.txt'

tot = 11 # ID/No. of filter used for total magnitude measurement
mlr = 11 # ID/No. of filter used for mass-to-light ratio computation

# Redshift Grid 
zmin = 0.
zmax = 9.
zstep = 0.02

madau = True                           # Include HI absorption (Madau 1995)


"""
Template matching parameters
match.py

"""
#Observed magnitudes and redshifts to fit
#'input_catalog = '/data/candels/catalogs/mocks/gs_newmocks_full.fits'
#input_catalog = '/data/candels/catalogs/GOODS-S/gs_all_tf_h_120919a_multi.mag.highz.noAGN.fits'
#input_catalog = '/data/candels/photz/GOODS-S/130213/gs_all_tf_h_130213a_multi.sample.phot.fits'
#input_catalog = '/data/MRC0406_irac.mag.fits'
#input_catalog = '/data/candels/catalogs/mocks/regions/mockset1_wide.sample.fits'
input_catalog = '/data/candels/catalogs/GOODS-S/MCcatalogs/gs_all_tf_h_130213a_multi.sample.MCmaster.fits'

ID_col = 'ID'
z_col = 'z_a'

flux_col_end = 'flux'                 # Naming convention for input catalog
fluxerr_col_end ='fluxerr'            # E.g. ACS_F606W_*flux_col_end* = F606W Flux
mag_col_end = 'mag'                   # ACS_F606W_*fluxerr_col_end* = F606W Flux Error
magerr_col_end ='magerr'

filts_used = [0,1,2,3,4,5,6,7,        # Filter columns from input used
              8,9,10,11,12,13,14,15]              # (If not all)
#filts_used = [0,1,2,3]


flux_corr = 1.0                       # Flux correction to convert to Total Mags

fit_mode = 'flux'                     # Fit to Mag Colours - 'colours' or 'flux'
#fit_mode = 'colours'

calc_mode = False                     # Calculate Mode Mass
mode_mass_percentage = 10.            # Top percentage of fits to use in mode
                                      # mass calculations

muv_max, muv_min, muv_bins = -23., -16., 100 # PDF bins
mass_max, mass_min, mass_bins = 12., 6., 100
beta_max, beta_min, beta_bins = 3., -3., 50 

#Output table parameters:
output_name = '/data/candels/catalogs/mocks/MC/gs_all_tf_h_130213a_multi.chab.CB07.MC.neb.results'
#output_name = 'results/200912/gs_all_tf_h_130213a_multi.chab.MC.results'
#output_name = '/data/MRC0406_inc.SED.fits'
table_format = 'fits'                 # Available formats:'fits'/'ascii'/'IPAC'
                                      # (those available to local installation
                                      #  of AtPY)

