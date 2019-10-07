import numpy as np
from six.moves import cPickle
#import hickle
#import sys
#sys.path.append('/Users/ken/Documents/Astro/code/smpy/')

import smpy.smpy as S
import smpy.ssp as B
from smpy.sfh import dblpower
from smpy.dust import Charlot

import matplotlib.pyplot as plt
import astropy.units as u
from astropy.cosmology import FlatLambdaCDM

cosmo = FlatLambdaCDM(H0=70, Om0=0.3)

bc03 = B.BC('data/ssp/bc03/chab/lr/', verbose=True)
models = S.CSP(bc03)
#BC = B.BC('')

def make_dbl_sfhs(nalpha, nbeta, ntau,
                  lalpha_min=-1., lalpha_max=3.,
                  lbeta_min=-1., lbeta_max=3.,
                  tau_min=0.1, tau_max=1.0):

    alphas, betas, taus = np.meshgrid(np.logspace(lalpha_min, lalpha_max, nalpha),
                                      -1*np.logspace(lbeta_min, lbeta_max, nbeta),
                                      np.linspace(tau_min, tau_max, ntau))

    sfh_params = np.vstack([alphas.flatten(), betas.flatten(), taus.flatten()]).T
    return sfh_params

sfhs = make_dbl_sfhs(7, 7, 7)

zrange = np.linspace(0., 2., 21.)
zform = 20.

ages = cosmo.age(zrange) - cosmo.age(zform)

metallicities = 1.
dusts = [0., 0.5, 1.] #np.linspace(0., 2., 11) # Av range

models.build(ages, sfhs, dusts, metallicities, sfh_law=dblpower, verbose=True)

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

# Convolve CSP models with filter set over redshift range
Obs = S.ObserveToFile()
Obs.build(models, filters, zrange, 'candels.goodss.models.test.hdf', verbose=True)
