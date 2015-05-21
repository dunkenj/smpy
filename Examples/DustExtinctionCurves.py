
# coding: utf-8

from pylab import *
from mkCSPs import *


# In[18]:

# In[2]:

SHARDS = FilterSet('../../HFF/transmCurves_SHARDS/shards_f*.res')
gs = glob('../smpy-fit/GOODS-S_18_FilterCurves/Filter*.txt')
print gs
SHARDS.addFileFilter(gs[2])
SHARDS.addFileFilter(gs[3])
SHARDS.addFileFilter(gs[4])
SHARDS.addFileFilter(gs[5])
SHARDS.addFileFilter(gs[6])
SHARDS.addFileFilter(gs[8])
SHARDS.addFileFilter(gs[9])
SHARDS.addFileFilter(gs[10])
SHARDS.addFileFilter(gs[11])

for filt in gs[-5:]:
    SHARDS.addFileFilter(filt)

FStar = FilterSet('*_cam_optics.txt')

# In[4]:

NoDust = CSP('../ssp/bc03/chab/lr/')
NoDust.build(1.,1.4,0.,2,sfh_law='pow')
#NoDust.addEmissionLine(1216,100)
NoDust /= NoDust.SFR
Av = 0.5

Calz = CSP('../ssp/bc03/chab/lr/',1.,1.4,Av,2,sfh_law='pow',dustmodel='calzetti')
Calz /= Calz.SFR
#Calz.addEmissionLine(1216,100)

SMC = CSP('../ssp/bc03/chab/lr/')
SMC.build(1.,1.4,Av,2,sfh_law='pow',dustmodel='smc')
SMC /= SMC.SFR
#SMC.addEmissionLine(1216,100)


LMC = CSP('../ssp/bc03/chab/lr/')
LMC.build(1.,1.4,Av,2,sfh_law='pow',dustmodel='lmc')
LMC /= LMC.SFR
#LMC.addEmissionLine(1216,100)


MW = CSP('../ssp/bc03/chab/lr/')
MW.build(1.,1.4,Av,2,sfh_law='pow',dustmodel='mw')
#MW.addEmissionLine(1216,100)
MW /= MW.SFR

# In[12]:

M = 30.
z = 2.2
O_N = Observe(NoDust*M,SHARDS,z)
O_Calz = Observe(Calz*M,SHARDS,z)
O_SMC = Observe(SMC*M,SHARDS,z)
O_LMC = Observe(LMC*M,SHARDS,z)
O_MW = Observe(MW*M,SHARDS,z)

FS = Observe(MW*M,FStar,z)

# In[15]:

wavelengths = numpy.array([filt.lambda_c for filt in SHARDS.filters])

ww = numpy.argmin(numpy.abs(NoDust.wave*(1+z) - 35573))
N_sed_AB = -2.5*numpy.log10((NoDust.SED / (2.997925e18/((NoDust.wave/(1+z))**2))))
N_sed_AB += (O_N.AB[-4].value - N_sed_AB[ww])

ww = numpy.argmin(numpy.abs(SMC.wave*(1+z) - 35573))
SMC_sed_AB = -2.5*numpy.log10((SMC.SED / (2.997925e18/((SMC.wave/(1+z))**2))))
SMC_sed_AB += (O_SMC.AB[-4].value - SMC_sed_AB[ww])
#semilogx(wavelengths,O_Calz.AB,'o',label='Calzetti')

ww = numpy.argmin(numpy.abs(MW.wave*(1+z) - 35573))
MW_sed_AB = -2.5*numpy.log10((MW.SED / (2.997925e18/((MW.wave/(1+z))**2))))
MW_sed_AB += (O_MW.AB[-4].value - MW_sed_AB[ww])

fig = figure(figsize=(5.5,3.))
Ax = fig.add_subplot(111)

noise = normal(scale=0.1,size=len(SHARDS.filters))
Ax.plot(SMC.wave*(1+z)/1e4,SMC_sed_AB,'-',lw=2,color='dodgerblue',alpha=0.5)
Ax.plot(MW.wave*(1+z)/1e4,MW_sed_AB,'-',lw=2,color='olivedrab',alpha=0.5)

Ax.plot(NoDust.wave*(1+z)/1e4,N_sed_AB,'-',lw=2,color='0.5',alpha=0.5,label='Intrinsic')

#Ax.errorbar(wavelengths/1e4,O_N.AB.value,yerr=0,fmt='s',label='Intrinsic',mew=0,ms=7,color='0.5',capsize=0)
Ax.errorbar(wavelengths/1e4,O_SMC.AB.value+noise,xerr = O_SMC.fwhm.value/2e4, yerr=0.1,fmt='o',label='SMC',mew=0,ms=7,capsize=0)
#Ax.plot(wavelengths,O_LMC.AB.value+normal(scale=0.001,size=len(SHARDS.filters)),'p',label='LMC')
Ax.errorbar(wavelengths/1e4,O_MW.AB.value+noise,xerr = O_MW.fwhm.value/2e4, yerr=0.1,fmt='d',label='MW',mew=0,ms=7,capsize=0)

#Ax.errorbar(FS.wl.value/1e4,FS.AB.value,fmt='x',xerr = FS.fwhm.value/2e4,ms=8,mew=2,color='maroon',capsize=0,zorder=10000,label = 'FourStar')

#Ax.fill_between([0.996,1.1123],21,28,color = 'orange',alpha=0.5)
#Ax.fill_between([1.066,1.227],21,28,color  = 'orange',alpha=0.5)
#Ax.fill_between([1.209,1.368],21,28,color  = 'orange',alpha=0.5)

#Ax.fill_between([1.46,1.643],21,28,color   = 'orange',alpha=0.5)
#Ax.fill_between([1.612,1.7960],21,28,color = 'orange',alpha=0.5)




#Ax.text(0.4000, 24, 'HETDEX Coverage', color='maroon',size=10,horizontalalignment='center',
#        verticalalignment='center',rotation=90)
Ax.set_ylim([26,21.5])
Ax.set_yticks([26,25,24,23,22])
Ax.set_xscale('log')
Ax.set_xticks([0.5,1,2,3,4,5])
Ax.set_xticklabels(['0.5','1','2','3','4','5'])
Ax.set_xlim([0.35,2.3])
Ax.legend(loc='lower right',prop={'size':9})
Ax.set_xlabel(r'Wavelength, [$\mu m$]')
Ax.set_ylabel('AB Mag')
Ax.text(0.8, 25.8,r'$\rm{SFR}\/=\/%.0f\/\rm{M}_{\odot}\/yr^{-1}$' % M,size=10,horizontalalignment='left')
Ax.text(0.80, 25.45,r'$A_{v}\/=\/%.1f$' % Av,size=10,horizontalalignment='left')
Ax.text(0.8,25.1,r'$z\/=\/%.1f$' % z,size=10,horizontalalignment='left')
fig.subplots_adjust(bottom=0.18)
fig.savefig('/Users/ken/Documents/jorbz/Michigan/plots/DustShape2.pdf',fmt='pdf',bbox_inches='tight')

fig.show()

