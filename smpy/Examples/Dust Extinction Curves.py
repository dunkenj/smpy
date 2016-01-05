
# coding: utf-8

# In[1]:

get_ipython().magic(u'pylab inline')
from pylab import *
from mkCSPs import *


# In[18]:

# In[2]:

SHARDS = FilterSet('../../HFF/transmCurves_SHARDS/shards_f*.res')
gs = glob('../smpy-fit/GOODS-S_18_FilterCurves/Filter*.txt')
print gs
#for filt in gs:
#    SHARDS.addFileFilter(filt)


# In[4]:

NoDust = CSP('../ssp/bc03/chab/lr/')
NoDust.build(1.,1.4,0.,2,sfh_law='pow')
NoDust /= NoDust.SFR
Av = 0.5

Calz = CSP('../ssp/bc03/chab/lr/',0.3,1.4,Av,2,sfh_law='pow',dustmodel='calzetti')
Calz /= Calz.SFR
Calz.addEmissionLine(1216,100)

SMC = CSP('../ssp/bc03/chab/lr/')
SMC.build(1.,1.4,Av,2,sfh_law='pow',dustmodel='smc')
SMC /= SMC.SFR
SMC.addEmissionLine(1216,100)


LMC = CSP('../ssp/bc03/chab/lr/')
LMC.build(1.,1.4,Av,2,sfh_law='pow',dustmodel='lmc')
LMC /= LMC.SFR
LMC.addEmissionLine(1216,100)


MW = CSP('../ssp/bc03/chab/lr/')
MW.build(1.,1.4,Av,2,sfh_law='pow',dustmodel='mw')
MW.addEmissionLine(1216,100)
MW /= MW.SFR

# In[12]:

M = 30.
z = 2.7
O_Calz = Observe(Calz*M,SHARDS,z)
O_SMC = Observe(SMC*M,SHARDS,z)
O_LMC = Observe(LMC*M,SHARDS,z)
O_MW = Observe(MW*M,SHARDS,z)


# In[15]:

wavelengths = numpy.array([filt.lambda_c for filt in SHARDS.filters])

#semilogx(wavelengths,O_Calz.AB,'o',label='Calzetti')

fig = figure(figsize=(5,3))
Ax = fig.add_subplot(111)

Ax.plot(wavelengths,O_SMC.AB.value+normal(scale=0.001,size=25),'s')
Ax.plot(wavelengths,O_LMC.AB.value+normal(scale=0.001,size=25),'p',label='LMC')
Ax.plot(wavelengths,O_MW.AB.value+normal(scale=0.001,size=25),'d',label='MW')
Ax.fill_between([3500,5800],22,28,color='orange',alpha=0.5)
Ax.set_ylim([27,23])
Ax.set_xlim([3000,13000])
Ax.legend(loc='lower right')
Ax.set_xlabel(r'Wavelength, [$\AA$]')
Ax.set_ylabel('AB Mag')
Ax.text(10000, 24.75,r'SFR = %.0f M$_{\odot}$/year' % M,size=10)
Ax.text(10000, 25.0,'Av = %.1f' % Av,size=10)
Ax.text(10000,24.5,'z = %.1f' % z,size=10)
Ax.savefig('DustShape2.pdf',fmt='pdf')

show()

