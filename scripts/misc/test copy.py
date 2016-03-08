import numpy as np

import smpy as S
import smpy.ssp as B

import matplotlib.pyplot as plt
import astropy.units as u
from astropy.cosmology import FlatLambdaCDM

cosmo = FlatLambdaCDM(H0=70, Om0=0.3)


bc03 = B.BC('smpy/data/ssp/bc03/chab/lr/bc2003_lr_m*')
models = S.CSP2(bc03)

EAZYfilters = S.LoadEAZYFilters('smpy/FILTER.RES.V8')

filters = S.FilterSet()
filters.addEAZYFilter(EAZYfilters, EAZYfilters.search('Bessel_B'))
filters.addEAZYFilter(EAZYfilters, EAZYfilters.search('Bessel_R'))
filters.addEAZYFilter(EAZYfilters, EAZYfilters.search('Bessel_I'))
filters.addEAZYFilter(EAZYfilters, EAZYfilters.search('SDSS/z'))
filters.addEAZYFilter(EAZYfilters, EAZYfilters.search('KPNO/IRIMJ'))
filters.addEAZYFilter(EAZYfilters, EAZYfilters.search('KPNO/IRIMH'))
filters.addEAZYFilter(EAZYfilters, EAZYfilters.search('irac'))

obs_fluxes = []

ages = [1e7] * u.yr
sfhs = np.array([10.])*u.Gyr
metallicities = [0.5, 1.] 
dusts = np.linspace(0.,2.,20)
fesc = [1.]
#models.build(ages, sfhs, dusts, metallicities, verbose=True)
models.build(ages, sfhs, dusts, metallicities, fesc=fesc, verbose=True)

zrange = np.linspace(0, 6, 30)

Obs = S.Observe(models, filters, zrange)
Obs2 = S.FileObserve()
Obs2.build(models, filters, zrange, 'test_output.hdf',verbose=True)


#models = S.CSP(bc03)
#models.build(1*u.Gyr, 0.5*u.Gyr, 0, 1.)
#Obs = S.Observe(models, Filts, zrange, verbose=True)
"""
EAZY
Filts = S.FilterSet('smpy/data/filters/PanStarrs*.txt')

zrange = np.linspace(0, 7, 11)

Obs2 = S.FileObserve()
Obs2.build(models2, Filts, zrange, savepath='test.hdf', verbose=True)
#print((models.SED == models2.SED).all())

"""
