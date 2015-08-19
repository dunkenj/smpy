import ssp as B
import smpy as S
import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np

#single = B.BPASS(path='ssp/bpass_v1.1/SEDS/sed.bpass.instant.cloudy.sin.*')
single = B.BPASS2(path='ssp/BPASSv2_imf135all_100/OUTPUT_POP/spectra.z*.dat')


#binary = B.BPASS(path='ssp/bpass_v1.1/SEDS/sed.bpass.instant.cloudy.bin.*')
binary = B.BPASS2(path='ssp/BPASSv2_imf135all_100/OUTPUT_POP/spectra-bin.z*.dat')

#exit()

path = '/Users/ken/Documents/Astro/code/bc03/models/Padova1994/salpeter/bc2003_lr_BaSeL_m'
BC03 = B.BC(path)

bc03_gal = S.CSP(BC03)
bc03_gal.build(0.5*u.Gyr, 0.5*u.Gyr, 0., 1.)

bc03_gal *= bc03_gal.STR

bpass_gal = S.CSP(single)
bpass_gal.build(0.5*u.Gyr, 0.5*u.Gyr, 0., 1.)

bpass_bin = S.CSP(binary)
bpass_bin.build(0.5*u.Gyr, 0.5*u.Gyr, 0., 1.)


Filts = S.FilterSet()
Filts.addTophatFilter(6000*u.AA, 600*u.AA)
Filts.addTophatFilter(20000*u.AA, 1000*u.AA)

Obs1 = S.Observe(bc03_gal, Filts, 0.)
Obs2 = S.Observe(bpass_gal, Filts, 0.)
Obs3 = S.Observe(bpass_bin, Filts, 0.)


Fig, Ax = plt.subplots(1,figsize=(4.5,3))

Ax.loglog(bc03_gal.wave, np.squeeze(bc03_gal.SED))
Ax.loglog(bpass_gal.wave, np.squeeze(bpass_gal.SED))
Ax.loglog(bpass_bin.wave, np.squeeze(bpass_bin.SED))

#Ax.set_xlim([200, 30000])
#Ax.set_ylim([1e-5, 1e-1])
plt.show()