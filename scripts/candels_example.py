import numpy as np
import sys
sys.path.append('/Users/ken/Documents/Astro/code/smpy/')
import smpy as S
import smpy.ssp as B

import matplotlib.pyplot as plt
import astropy.units as u
from astropy.cosmology import FlatLambdaCDM

cosmo = FlatLambdaCDM(H0=70, Om0=0.3)

bc03 = B.BC('../smpy/data/ssp/bc03/chab/lr/bc2003_lr_m*')
models = S.CSP2(bc03)

filters = S.FilterSet('data/Filters/Filter*.txt')

obs_fluxes = []

ages = np.logspace(7, 10, 1) * u.yr
sfhs = np.array([0.1, 0.25, 0.5, 1., 2.5,  5., 7.5, 10.])*u.Gyr
metallicities = [0.2, 0.4, 0.6, 0.8, 1.] 
dusts = np.linspace(0.,2.,10)
fesc = [0., 0.2, 0.4, 0.6, 0.8, 1.]

models.build(ages, sfhs, dusts, metallicities, fesc=fesc, verbose=True)

zrange = np.linspace(0, 6, 1)

#Obs = S.Observe(models, filters, zrange)
Obs2 = S.FileObserve()
Obs2.build(models, filters, zrange, 'candels.goodss.models.hdf',verbose=True)


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
