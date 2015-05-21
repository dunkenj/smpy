import smpy as S
from ssp import BC, SSP
from dust import Charlot
import numpy as np
import matplotlib.pyplot as plt
reload(S)
import time
from astropy.utils.console import human_time
import astropy.units as u

path = '/Users/ken/Documents/PhD/code/bc03/models/Padova1994/chabrier/bc2003_lr_BaSeL_m'

ssp = BC(path)

gal = S.CSP(ssp)

#start = time.time()
#ages = np.logspace(-2, 1.1, 25) * u.Gyr
#taus =  np.array([0.25, 0.5, 1.,]) * u.Gyr
fescs = [0, 0.5, 1.]
#gal.build(ages,taus,0.,1., fesc=fescs, verbose=True)
#finish = (time.time()-start)
#print finish, finish/np.product(gal.SFR.shape), np.product(gal.SFR.shape)

start = time.time()
#ages2 = np.array([0.1, 0.25, 0.5, 1., 1.5, 2., 5., 10.]) * u.Gyr
ages2 = np.logspace(-2, 1.1, 25) * u.Gyr
taus2 =  np.array([0.05, 0.5, 1., 2.5, 10.]) * u.Gyr
tauv = np.array([0, 0.25, 0.5, 0.75, 1., 1.5, 2.])
metals = np.array([0.2, 1.])
gal.build(ages2,taus2,tauv,metals, fesc=fescs, verbose=True)

finish = (time.time()-start)
print finish, finish/np.product(gal.SFR.shape), np.product(gal.SFR.shape)
