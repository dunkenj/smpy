"""
Say stuff
"""
import numpy as np
import array, time, os, sys, re
from scipy.interpolate import griddata
from scipy.integrate import quad
from sm_functions import tl,t0,dm,dist,dec_a_func,dec_b_func
from glob import glob
import sm_params as p

from astropy import units as u
from astropy import constants as c
from astropy import cosmology as cos
from astropy.utils.console import ProgressBar

cosmo = cos.FlatLambdaCDM(H0=70, Om0=0.3)

"""
Command line arguments and parameter import section
"""
version = sys.version_info[1]
if version == 7:
	import argparse
	parser = argparse.ArgumentParser()
	parser.add_argument("-p","--params", type=str, default="sm_params",
						help = "Parameter file, default = sm_params")
	parser.add_argument("-q", "--quiet", help = "Suppress extra outputs",
						action = "store_true")
	args = parser.parse_args()
	quiet = args.quiet

	params_root = re.split(".py",args.params)[0]
	if os.path.isfile(params_root+".pyc"):
		os.remove(params_root+".pyc")

	import importlib
	try:
		params = importlib.import_module(params_root)
		print 'Loaded '+args.params+' as params'
	except:
		print 'Failed to load "'+args.params+'" as params, loading default instead'
		import sm_params as params
		
if version == 6:
	import optparse
	parser = optparse.OptionParser()
	parser.add_option("-p","--params", type=str, default="sm_params",
						dest="params",help = "Parameter file, default = sm_params")
	parser.add_option("-q", "--quiet", help = "Suppress extra outputs",
						dest="quiet", default=False,
						action = "store_true")
	args, dump = parser.parse_args()
	quiet = args.quiet

	params_root = re.split(".py",args.params)[0]
	if os.path.isfile(params_root+".pyc"):
		os.remove(params_root+".pyc")
	import imp
	try:
		fp, pathname, description = imp.find_module(params_root)
		params = imp.load_module(params_root, fp, pathname, description)
	except:
		print 'Failed to load "'+args.params+'" as params, loading default instead'
		import sm_params as params
	
elif version not in [6, 7]:		
	print 'Import option only coded for python versions 2.6 and 2.7... \n',
	print 'Loading default instead'
	import sm_params as params
	quiet = False

if quiet:
	print "Shhhh"

if quiet:
	quietprint = lambda *a: None
else:
	def quietprint(*args):
		for arg in args:
			print arg,
		print

f = open("error.log", "w")
original_stderr = sys.stderr
sys.stderr = f
 
start_time = time.clock()


def tau_madau(wave, z):
    """ Lyman limit system absorption following Madau 1995
    
    This version is in part based on my previous python implementation
    and the simpler version shown here
    https://github.com/spacetelescope/pysynphot/issues/77
    
    This version has been extended to include the additional
    higher order lines calculated in EAZY. 
    
    Parameters:
        wave (array): Array with astropy.units of length
        z (float): redshift at which the source is being observed
        
    Returns:
        LyL_absorption (array): Throughput along the line of sight
            for wavelengths in wave from Lyman-alpha forest and 
            Lyman-limit systems
    
    """

    lymanlim = 912. * u.AA

    lyn_w = np.array([1216, 1026., 972.8, 950., 938.1,
                      931.0, 926.5, 923.4, 921.2,
                      919.6, 918.4, 917.5, 916.7,
                      916.1, 915.6, 915.2]) * u.AA
    lyn_c = np.array([3.6e-3, 1.7e-3, 1.2e-3, 9.3e-4, 8.2e-4,
                      7.5e-4, 7.1e-4, 6.8e-4, 6.6e-4,
                      6.4e-4, 6.3e-4, 6.2e-4, 6.1e-4,
                      6.0e-4, 6.0e-4, 6.0e-4])
    tau = np.zeros(len(wave))

    for i in range(len(lyn_w)):
        tau = np.where(wave <= lyn_w[i] * (1 + z), tau + lyn_c[i] * (wave / lyn_w[i]) ** 3.46, tau)

    wly = wave / lymanlim
    wly3 = wly ** 3
    pe_abs = (0.25 * wly3 * ((1 + z) ** 0.46 - wly ** 0.46) +
              9.4 * wly ** 1.5 * ((1 + z) ** 0.18 - wly ** 0.18) -
              0.7 * wly3 * (wly ** (-1.32) - (1 + z) ** (-1.32)) -
              0.023 * ((1 + z) ** 1.68 - wly ** 1.68))
    tau = np.where(wave <= lymanlim * (1 + z), tau + pe_abs, tau)
    lyl_absorption = np.where(tau > 700., 0., np.exp(-tau))
    return np.clip(lyl_absorption, 0, 1)

input_binary = params.ssp_output
output_binary = params.synmag_output
input_head,input_tail = os.path.split(input_binary)

files = glob(params.filt_dir+params.filt_names) # Array of filter paths
files.sort()

if params.zspacing == 'linear':
    z = np.arange(params.zmin,params.zmax,params.zstep)
elif params.zspacing == 'log':
    z = np.logspace(params.zmin,params.zmax,params.n_zsteps)
    z = np.insert(z,[0],[0.])
z1= z+1

quietprint('{0:s} {1:.1f} {2:s} {3:.1f} {4:s} {5} {6:s} {7}'.format('Redshifts range from',z.min(),'to',z.max(),'in',len(z),'steps','\n')) 

quietprint('{0:s} {1:s}{2:s}'.format('Loading input SEDs from -',input_tail,':'))

parameters = np.load(input_binary+'.par.npy').item()

with np.load(input_binary) as data:
    #quietprint('{0:30s}'.format('Parameters')),
    #parameters = data['parameters'][()]

    quietprint('{0:30s}'.format('Wavelength Array')),
    wave = parameters['wavelengths'] * u.AA

    quietprint('{0:30s}'.format('Age Array')),
    tg = parameters['tg'] * u.yr

    quietprint('{0:30s}'.format('SEDs')),
    SED = data['SED'] * u.solLum / u.AA

    quietprint('{0:30s}'.format('Stellar Masses')),
    STR = data['STR']

    quietprint('{0:30s}'.format('Star Formation Rates')),
    SFR = data['SFR']
    quietprint('Done \n')


S = SED.shape
Flux = np.zeros((len(files),len(z),S[1],S[2],S[3],S[4])) * (u.erg / u.cm ** 2 / u.s / u.Hz)
Flux_UV = np.zeros((len(z),S[1],S[2],S[3],S[4])) * (u.erg / u.cm ** 2 / u.s / u.Hz)

"""
SECTION 1 
Calculate cosmological distance modulus 'dm' and the age of the universe
at each redshift 'z' - i.e. those variables which are functions of z
only
"""
quietprint("Setting up redshift dependent arrays:")

#dl = [dm(z[i]) for i in range(len(z))]
#dl = np.array(dl)
age = cosmo.age(z).to(u.yr)
quietprint('{0:30s}'.format('Ages')),

dl = cosmo.luminosity_distance(z).cgs
dl[0] = 10*c.pc.cgs

dmo = cosmo.distmod(z)
dmo[0] = 0.

quietprint('{0:30s}'.format('Distance Moduli')),

ai = np.zeros(len(z))
lyman_abs = np.ones((len(wave),len(z)))
if params.madau:
    for zi in range(len(z)):
        ai[zi] = np.array(np.where(tg<age[zi]))[0][-1]
        lyman_abs[:,zi] = tau_madau(wave, z[zi])
else:
    for zi in range(len(z)):
        ai[zi] = np.array(np.where(tg<age[zi]))[0][-1]
quietprint('{0:30s}'.format('Madau Absorption Decrements')),


"""
SECTION 2

Compute fluxes for each filter in turn
"""

print('{0:30s} {1:15s} {2:15s}').format('Filter', 'Filter length','Prep Time') 
print('-'*60)
print('\n')

for filt in range(len(files)):
    
    head, tail = os.path.split(files[filt])
    print '{0:30s}'.format(tail),

    start_filtinterp = time.clock()
    wf = np.loadtxt(files[filt], usecols=[0]) * u.AA
    tp = np.loadtxt(files[filt], usecols=[1])
    print('{0:<15d}').format(len(wf)),
    
    if len(wf) > 1000: # Re-sample large filters for performance
        wfx = np.linspace(wf[0],wf[-1],1000)
        tpx = griddata(wf,tp,wfx)

        wf = wfx
        tp = tpx


    # Find SED wavelength entries within filter range
    wff = np.array([wf[0] < wave[i] < wf[-1]
                    for i in range(len(wave))])
    wft = wave[wff]

    # Interpolate to find throughput values at new wavelength points
    tpt = griddata(wf, tp, wft)

    # Join arrays and sort w.r.t to wf
    # Also replace units stripped by concatenate
    wf = np.array(np.concatenate((wf, wft))) * u.AA
    tp = np.concatenate((tp, tpt))

    order = np.argsort(wf)
    wf = wf[order]
    tp = tp[order]
       
    print '{0:<15.2f}'.format(time.clock()-start_filtinterp)

    start_conv = time.clock()
    with ProgressBar(len(z)) as bar:
        for zi in range(len(z)):
            # Interpolate redshifted SED and LyAbs at new wavelength points
            sed = griddata(wave * (1 + z[zi]), SED[:,:ai[zi]+1,:,:], wf) * SED.unit
            lyabs = griddata(wave, lyman_abs[:,zi], wf)

            # Calculate f_nu mean
            # Integrate SED through filter, as per BC03 Fortran
            # As: f_nu=int(dnu Fnu Rnu/h*nu)/int(dnu Rnu/h*nu)
            # ie: f_nu=int(dlm Flm Rlm lm / c)/int(dlm Rlm/lm)

            top = np.trapz(sed * (lyabs * tp * wf)[:,None,None,None,None] / c.c.to(u.AA / u.s), 
                           wf[:,None, None,None,None], axis=0)
            bottom = np.trapz(tp / wf, wf)

            area = (4 * np.pi * (dl[zi] ** 2))

            F_mean = top / bottom / (1 + z[zi]) / area

            # Set flux to appropriate units and calculate AB magnitude
            Flux[filt,zi,:ai[zi]+1,:,:] = F_mean.to(u.erg / u.cm ** 2 / u.s / u.Hz)
            bar.update()
    print('\n')

"""
Compute rest-frame UV (1500A) flux
"""

compute_MUV = True
if compute_MUV:
    print('{0}{1}').format('\n','Calculating Rest-Frame UV (1500AA) Fluxes: ')
    print('-'*60)
    start_filtinterp = time.clock()
    wf = np.arange(1445, 1555) * u.AA
    tp = np.zeros(len(wf))
    tp[(wf >= 1450*u.AA) & (wf < 1551*u.AA)] = 1.0

    print '{0:30s}'.format('1500A Tophat filter'),
    print('{0:<15d}').format(len(wf)),
    
    # Find SED wavelength entries within filter range
    wff = np.array([wf[0] < wave[i] < wf[-1]
                    for i in range(len(wave))])
    wft = wave[wff]

    # Interpolate to find throughput values at new wavelength points
    tpt = griddata(wf, tp, wft)

    # Join arrays and sort w.r.t to wf
    # Also replace units stripped by concatenate
    wf = np.array(np.concatenate((wf, wft))) * u.AA
    tp = np.concatenate((tp, tpt))

    order = np.argsort(wf)
    wf = wf[order]
    tp = tp[order]
       
    print '{0:<15.2f}'.format(time.clock()-start_filtinterp)

    start_conv = time.clock()
    with ProgressBar(len(z)) as bar:
        for zi in range(len(z)):
            # Interpolate redshifted SED and LyAbs at new wavelength points
            sed = griddata(wave * (1 + z[zi]), SED[:,:ai[zi]+1,:,:], wf) * SED.unit
            lyabs = griddata(wave, lyman_abs[:,zi], wf)

            """
            Old Method:
    
            WR = (np.trapz(sed * lyabs * tp * wf, wf) /c.c.to(u.AA/u.s)).value
            F_mean = WR/f_mean2.value/z1#/c.c.to(u.AA / u.s)
    
            #Convert fluxes to AB magnitudes
            Mag = -2.5*np.log10( F_mean * c.L_sun.cgs.value/ (4* np.pi* (c.pc.cgs.value*10)**2) ) - 48.6
            Mag += self.dm
    
            # Coded as this in BC03
            #
            # AB0 = -2.5*np.log10( c.L_sun.cgs.value/ (4* np.pi* (10 * c.pc.cgs.value)**2) )
            # Mag = AB0 - 2.5*np.log10(F_mean) - 48.6
            
            Flux = 10**((23.9 - Mag)/2.5) * u.uJy #uJy
            #print Flux/Flux2.to(u.mJy)
            print Flux, self.tmp.to(u.uJy)
    
            """

            # Calculate f_nu mean
            # Integrate SED through filter, as per BC03 Fortran
            # As: f_nu=int(dnu Fnu Rnu/h*nu)/int(dnu Rnu/h*nu)
            # ie: f_nu=int(dlm Flm Rlm lm / c)/int(dlm Rlm/lm)

            top = np.trapz(sed * (lyabs * tp * wf)[:,None,None,None,None] / c.c.to(u.AA / u.s), 
                           wf[:,None, None,None,None], axis=0)
            bottom = np.trapz(tp / wf, wf)

            area = (4 * np.pi * (dl[zi] ** 2))

            F_mean = top / bottom / (1 + z[zi]) / area

            # Set flux to appropriate units and calculate AB magnitude
            Flux_UV[zi,:ai[zi]+1,:,:] = F_mean.to(u.erg / u.cm ** 2 / u.s / u.Hz)
            bar.update()

print('\n')
print('{0}').format('Converting Flux arrays to AB Magnitudes: '),

Mags = -2.5 * np.log10(Flux.to(u.Jy) / (3631 * u.Jy))
MUV = -2.5 * np.log10(Flux_UV.to(u.Jy) / (3631 * u.Jy))

#Store mass-to-light ratio for selected filter


print('Done')


print('{0} {1}{2}').format('Saving output binaries to',output_binary,':'),
if os.path.isfile(output_binary+'.main'+'.npz'):
    os.remove(output_binary+'.main'+'.npz')
if os.path.isfile(output_binary+'.mags'+'.npy'):
    os.remove(output_binary+'.mags'+'.npy')
if os.path.isfile(output_binary+'.fluxes'+'.npy'):
    os.remove(output_binary+'.fluxes'+'.npy')

np.savez(output_binary+'.main',parameters=parameters,z=z,filters=files,SFR=SFR,Mshape=Flux.shape,MUV=MUV.value,UV_flux=Flux_UV.to(u.uJy).value)
np.save(output_binary+'.mags',Mags.value)
np.save(output_binary+'.fluxes',Flux.to(u.uJy).value)
print('Done')

print('{0}{1}{2}{3:s} {4:.1f}').format('\n','All finished!','\n','Time elapsed:',(time.clock()-start_time))


sys.stderr = original_stderr
f.close() 
