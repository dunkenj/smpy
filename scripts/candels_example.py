import numpy as np
from six.moves import cPickle
#import hickle
#import sys
#sys.path.append('/Users/ken/Documents/Astro/code/smpy/')

import smpy.smpy as S
import smpy.ssp as B

import matplotlib.pyplot as plt
import astropy.units as u
from astropy.cosmology import FlatLambdaCDM

cosmo = FlatLambdaCDM(H0=70, Om0=0.3)

bc03 = B.BC('data/ssp/bc03/chab/lr/', verbose=True)
models = S.CSP(bc03)


ages = np.logspace(7, 10, 10) * u.yr # Ages since onset of star-formation
sfhs = np.array([0.25, 0.5, 1., 3., 5.])*u.Gyr # SFH timescale (decreasing exponential)
metallicities = bc03.metallicities
dusts = np.linspace(0., 2., 11) # Av range

models.build(ages, sfhs, dusts, metallicities, verbose=True)

"""
Example: Save intermediate models
"""
with open('candels.goodss.csp.pkl', 'wb') as output:
    cPickle.dump(models, output, protocol=2)

"""
Example: Load saved models
"""
with open('candels.goodss.csp.pkl', 'rb') as input:
    models = cPickle.load(input)
    print('models_loaded')

# Load CANDELS Filter Set
filters = S.FilterSet('data/Filters/GS/*.txt')

zrange = np.linspace(0, 9, 101) # Decrease step size when doing full fits

# Convolve CSP models with filter set over redshift range
Obs = S.ObserveToFile()
Obs.build(models, filters, zrange, 'candels.goodss.models.test.hdf', verbose=True)
