import numpy as np
from six.moves import cPickle
#import sys
#sys.path.append('/Users/ken/Documents/Astro/code/smpy/')

import smpy.smpy as S
import smpy.ssp as B

import matplotlib.pyplot as plt
import astropy.units as u
from astropy.cosmology import FlatLambdaCDM

cosmo = FlatLambdaCDM(H0=70, Om0=0.3)

bc03 = B.BC('data/ssp/bc03/chab/lr/bc2003_lr_m*')
models = S.CSP(bc03)


ages = np.logspace(7, 10, 5) * u.yr
sfhs = np.array([0.1, 0.25, 0.5, 1., 2.5,  5., 7.5, 10.])*u.Gyr
metallicities = [0.2, 1., 2.5] 
dusts = np.linspace(0.,2.,5)
fesc = [0., 1.]

models.build(ages, sfhs, dusts, metallicities, fesc=fesc, verbose=True)
cPickle.dump(models, file('candels.goodss.csp.pickle','wb'), protocol=-1)

models = cPickle.load(file('candels.goodss.csp.pickle', 'rb'))
print('models_loaded')


filters = S.FilterSet('data/Filters/GS/Filter*.txt')

zrange = np.linspace(0, 9, 45)

#Obs = S.Observe(models, filters, zrange)
Obs2 = S.ObserveToFile()
Obs2.build(models, filters, zrange, 'candels.goodss.models.test.hdf',verbose=True)

"""
for key in models.__dict__.keys():
    try:
        print key, hf(models.__dict__[key].nbytes)
    except:
        continue
"""

