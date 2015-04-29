import mkCSPs_dev as d
import mkCSPs as C
import matplotlib.pyplot as plt
import numpy as np
import astropy.units as u
import astropy.constants as c
reload(d)
path = '/Users/ken/Documents/PhD/code/bc03/models/Padova1994/chabrier/bc2003_lr_BaSeL_m'

bcflux = np.loadtxt('/Users/ken/Documents/PhD/code/bc03/test3_lr_BaSeL_ssp.136')
bcdata = np.loadtxt('/Users/ken/Documents/PhD/code/bc03/test3_lr_BaSeL_ssp.4color')
bcdata2 = np.loadtxt('/Users/ken/Documents/PhD/code/bc03/test3_lr_BaSeL_ssp.3color')

bcage = 10**bcdata[134,0]
print bcage
bcmass = bcdata[134,5]
print bcmass


bc_vmag = bcdata[134,2]
print bc_vmag

ssp = d.BC(path)

gal = d.CSP(ssp)

gal.build(bcage*u.yr,0.5*u.Gyr,0.,1.)
gal *= bcmass

print bcdata2[134, 5], gal.Nly

gal2 = C.CSP(path,bcage/1e9,0.5,0,4)
gal2 *= bcmass

Fi = d.FilterSet('../smpy-fit/GOODS-S_18_FilterCurves/Filter*.txt')
#Fi.addFileFilter

Fi2 = C.FilterSet('../smpy-fit/GOODS-S_18_FilterCurves/Filter03*.txt')

Obs = d.Observe(gal,Fi,0.,v=1)
print Obs.AB
Obs = d.Observe(gal,Fi,0.,v=2)
print Obs.AB

#Obs2 = C.Observe(gal2,Fi,0.)

#print Obs.fluxes
#print Obs.AB

#print Obs2.AB

#Obs = d.Observe(gal,Fi,[0,1,2,3,4],v=2)
#print Obs.fluxes
#print Obs.AB


Fig, Ax = plt.subplots(1)

lbc = len(bcflux[:,0])
#Ax.loglog(bcflux[:,0],bcflux[:,1])
#Ax.loglog(gal.wave[:lbc],gal.SED[:lbc])

Ax.semilogx(bcflux[:,0], (bcflux[:,1] - gal.SED.value[:lbc]) / bcflux[:,1])
Ax.semilogx(bcflux[:,0], (bcflux[:,1] - gal2.SED[:lbc]) / bcflux[:,1])
Ax.semilogx(gal.wave[:lbc], (gal2.SED[:lbc] - gal.SED.value[:lbc]) / bcflux[:,1],'--')
#Ax.set_ylim([-0.3,0.3])
Ax.set_xlim([100,50000])

Fig2, Ax = plt.subplots(1)

lbc = len(bcflux[:,0])
#Ax.loglog(bcflux[:,0],bcflux[:,1])
#Ax.loglog(gal.wave[:lbc],gal.SED[:lbc])
Ax.semilogx(bcflux[:,0], 2.5*np.log10(bcflux[:,1]))
Ax.semilogx(bcflux[:,0], 2.5*np.log10(gal.SED.value[:lbc]))
Ax.semilogx(bcflux[:,0], 2.5*np.log10(gal2.SED[:lbc]))

#Ax.set_ylim([-0.3,0.3])
Ax.set_xlim([100,10000])

plt.show()
