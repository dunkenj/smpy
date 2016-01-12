import numpy as np

import smpy as S
import smpy.ssp as B

import astropy.units as u

bc03 = B.BC('smpy/data/ssp/bc03/chab/lr/bc2003_lr_m*')
models = S.CSP(bc03)


ages = np.logspace(7, 10, 20) * u.yr
sfhs = np.array([0.05,0.1,0.25,0.5,1.,2.,3.,5.,10.,-0.25,-0.5,-1.,2.,-3.,-5.,-10.,1000])*u.Gyr
metallicities = [0.1, 0.3, 0.5, 0.7, 1.5]
dusts = [0.,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]

models.build(ages, sfhs, dusts, metallicities, verbose=True)

