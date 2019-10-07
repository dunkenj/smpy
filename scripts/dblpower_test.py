import numpy as np
import astropy.units as u
from smpy import smpy as S
from smpy import ssp as B
from smpy.sfh import dblpower
import matplotlib.pyplot as plt

BC = B.BC('data/ssp/bc03/chab/lr/')
models = S.CSP(BC)

tau = 1*u.Gyr

age = 3*u.Gyr
sfh = [(1., -1., 0.3), (1., -1., 0.5)]
dust = 0.
metal = 1.

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

sfh = make_dbl_sfhs(9, 9, 7)

#sfh = [(31.6227766,-31.6227766,0.7)]
#sfh = [(31.62,-10.,0.7)]

models.build(age, sfh, dust, metal, fesc=[0., 0.1, 0.2], sfh_law=dblpower, verbose=True)
