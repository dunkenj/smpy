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

stop
obs_fluxes = []

ages = [1e7] * u.yr
sfhs = np.array([10.])*u.Gyr
metallicities = [0.5, 1.] #[0.01, 0.1, 0.3, 0.5, 0.7, 1., 2.]
dusts = np.linspace(2.,2.2,30)
fesc = [1.]
#models.build(ages, sfhs, dusts, metallicities, verbose=True)
models.build(ages, sfhs, dusts, metallicities, fesc=fesc, verbose=True)

Obs = S.Observe(models, filters, 1.7)

norm_sfr = 1200

Ob_Mags = [26.04, 24.9, 24.144, 23.12, 22.07, 22.17, 21.29, 20.89, 20.67, 20.84]


Norm_flux = (10**(-0.4*(Ob_Mags[6]-23.9))/Obs.fluxes[:,6].value)[None, None, :]

Norm = (norm_sfr/models.SFR.value)[None, None, :]
#
Mags = 23.9 - 2.5*np.log10(Norm * Obs.fluxes.value).reshape((10, np.product(Norm.shape)))


plt.semilogx(Obs.wl, Mags, color='firebrick', alpha=0.2)
plt.semilogx(Obs.wl, Ob_Mags, 'o')
plt.xlim([5000, 200000])
plt.ylim([23, 19])
plt.show()


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
