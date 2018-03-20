import numpy as np
from six.moves import cPickle
import hickle
#import sys
#sys.path.append('/Users/ken/Documents/Astro/code/smpy/')

import smpy.smpy as S
import smpy.ssp as B

import matplotlib.pyplot as plt
import astropy.units as u
from astropy.cosmology import FlatLambdaCDM

cosmo = FlatLambdaCDM(H0=70, Om0=0.3)


bc03 = B.BC('/Users/ken/Astro/code/smpy/scripts/data/ssp/Miles_Atlas/Chabrier_IMF/bc2003_hr*')
models = S.CSP(bc03)


ages = np.logspace(7, 10, 5) * u.yr
sfhs = np.array([0.5, 1., 10.])*u.Gyr
metallicities = bc03.metallicities
dusts = np.linspace(0.,2.,5)
fesc = np.linspace(0, 1, 2)

models.build(ages, sfhs, dusts, metallicities, fesc=fesc, verbose=True)
hickle.dump(models, 'candels.goodss.csp.hkl', mode='w')

models = hickle.load('candels.goodss.csp.hkl')
print('models_loaded')


filters = S.FilterSet('data/Filters/*.res')

zrange = np.linspace(0, 9, 901)

#Obs = S.Observe(models, filters, zrange)
#Obs2 = S.ObserveToFile()
Obs2.build(models, filters, zrange, 'candels.goodss.models.test.hdf', verbose=True)

"""
for key in models.__dict__.keys():
    try:
        print key, hf(models.__dict__[key].nbytes)
    except:
        continue
"""
