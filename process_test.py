"""
Say stuff
"""
import argparse
import importlib
import numpy, array, time, os, sys, re
from scipy.interpolate import griddata
from scipy.integrate import quad
from sm_functions import tl,t0,dm,dist,dec_a_func,dec_b_func
from glob import glob
import sm_params as p

"""
Command line arguments and parameter import section
"""

parser = argparse.ArgumentParser()
parser.add_argument("-p","--params", type=str, default="sm_params",
					help = "Parameter file, default = sm_params")
parser.add_argument("-q", "--quiet", help = "Suppress extra outputs",
					action = "store_true")
args = parser.parse_args()

params_root = re.split(".py",args.params)[0]
if os.path.isfile(params_root+".pyc"):
	os.remove(params_root+".pyc")

params = importlib.import_module(params_root) 
print 'Loaded '+args.params+' as params'

if args.quiet:
	print "Shhhh"

if args.quiet:
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

input_binary = params.ssp_output
output_binary = params.synmag_output
input_head,input_tail = os.path.split(input_binary)

files = glob(params.filt_dir+params.filt_names) #Array of filter paths
files.sort()

tot = params.tot #Filter to use as total magnitude
mlr = params.mlr #Filter to use for mass-to-light ratio

z = numpy.arange(params.zmin,params.zmax,params.zstep)
z1= z+1

quietprint('{0:s} {1:.1f} {2:s} {3:.1f} {4:s} {5} {6:s} {7}'.format('Redshifts range from',params.zmin,'to',params.zmax,'in',len(z),'steps','\n')) 

quietprint('{0:s} {1:s}{2:s}'.format('Loading input SEDs from -',input_tail,':'))

with numpy.load(input_binary) as data:
    quietprint('{0:30s}'.format('Parameters')),
    parameters = data['parameters'][0]
    quietprint('Done')

    quietprint('{0:30s}'.format('Wavelength Array')),
    wave = parameters[0]
    quietprint('Done')

    quietprint('{0:30s}'.format('Age Array')),
    tg = parameters[1]
    quietprint('Done')

    quietprint('{0:30s}'.format('SEDs')),
    SED = data['SED']
    quietprint('Done')

    quietprint('{0:30s}'.format('Stellar Masses')),
    STR = data['STR']
    quietprint('Done \n')

    quietprint('{0:30s}'.format('Star Formation Rates')),
    SFR = data['SFR']
    quietprint('Done \n')


S = SED.shape
F_mean = numpy.zeros((len(files),len(z),S[1],S[2],S[3],S[4]))
F_mean_UV = numpy.zeros((len(z),S[1],S[2],S[3],S[4]))

AB0 = 5*numpy.log10(1.7684e8*1e-5)


"""
SECTION 1 
Calculate cosmological distance modulus 'dm' and the age of the universe
at each redshift 'z' - i.e. those variables which are functions of z
only
"""
quietprint("Setting up redshift dependent arrays:")

#dl = [dm(z[i]) for i in range(len(z))]
#dl = numpy.array(dl)
age = [tl(z[i]) for i in range(len(z))]
age = (t0()-numpy.array(age))*1e9
quietprint('{0:30s}'.format('Ages')),

dmo = [dist(z[i]) for i in range(len(z))]
dmo = numpy.array(dmo)
dmo[0] = 0.
quietprint('{0:30s}'.format('Distance Moduli')),

ai = numpy.zeros(len(z))
dec_a = numpy.zeros(len(z))
dec_b = numpy.zeros(len(z))

lyman_abs = numpy.ones((len(wave),len(z)))
ly_cont_w = numpy.array([(wave<=912.)][0])
ly_b_w = numpy.array([(wave > 912.) & (wave <= 1026.)][0])
ly_a_w = numpy.array([(wave > 1026.) & (wave <= 1216.)][0])

if params.madau:
    for zi in range(len(z)):
        ai[zi] = numpy.array(numpy.where(tg<age[zi]))[0][-1]
        dec_a[zi] = (1/(120*(1+z[zi])))*quad(dec_a_func,1050*(1+z[zi]),1170*(1+z[zi]))[0]
        dec_b[zi]= (1/(95*(1+z[zi])))*quad(dec_b_func,920*(1+z[zi]),1015*(1+z[zi]))[0]
        
        lyman_abs[ly_cont_w,zi] = 0.
        lyman_abs[ly_b_w,zi] = dec_b[zi]
        lyman_abs[ly_a_w,zi] = dec_a[zi]
else:
    for zi in range(len(z)):
        ai[zi] = numpy.array(numpy.where(tg<age[zi]))[0][-1]
quietprint('{0:30s}'.format('Madau Absorption Decrements')),


"""
SECTION 2

Compute fluxes for each filter in turn
"""

print('{0:30s} {1:15s} {2:15s} {3:20s}').format('Filter', 'Filter length','Prep Time', 'Calc Time (per z)') 

for filt in range(len(files)):
    head, tail = os.path.split(files[filt])
    print '{0:30s}'.format(tail),
    sys.stdout.flush()

    start_filtinterp = time.clock()
    wf = numpy.loadtxt(files[filt], usecols=[0])
    tp = numpy.loadtxt(files[filt], usecols=[1])
    print('{0:<15d}').format(len(wf)),
    if len(wf) > 1000: #Re-sample large filters for performance
        wfx = numpy.linspace(wf[0],wf[-1],1000)
        tpx = griddata(wf,tp,wfx)

        wf = wfx
        tp = tpx
        

    #Find SED wavelength entries within filter range
    wff = numpy.array([wf[0] < wave[i] < wf[-1] for i in range(len(wave))])
    wft = wave[wff]
    
    #Interpolate to find throughput values at new wavelength points
    tpt = griddata(wf,tp,wft)
    
    #Join arrays and sort w.r.t to wf
    wf = numpy.concatenate((wf,wft))
    tp = numpy.concatenate((tp,tpt))
    
    order = numpy.argsort(wf)
    wf = wf[order]
    tp = tp[order]

    dwf = numpy.diff(wf)
    nwf = len(wf)

    tpwf = tp/wf
    
    f_mean2 = numpy.dot(dwf,(tpwf[:nwf-1]+tpwf[1:])/2)

    tpwf = tp*wf #Reassign tpwf as product

    print '{0:<15.2f}'.format(time.clock()-start_filtinterp),
    sys.stdout.flush()

    start_conv = time.clock()
    for zi in range(len(z)):
        wf1 = wf/z1[zi]

        WR = 0.

        for i in range(nwf):

            #Interpolation indices
            j = numpy.where(wave<wf1[i])[0][-1]
            #print j

            a = (wf1[i] - wave[j])/(wave[j+1]-wave[j])
            tpa = tpwf[i]*((1-a)*(SED[j,:ai[zi]+1,:]*lyman_abs[j,zi]) + a*SED[j+1,:ai[zi]+1,:]*lyman_abs[j+1,zi])
            
            if i != 0:
                WR += dwf[i-1]*(tpb+tpa)

            tpb = tpa
    
        F_mean[filt,zi,:ai[zi]+1,:] = WR/2/z1[zi]/f_mean2/2.997925e18

    print '{0:<20.2f}'.format((time.clock()-start_conv)/len(z))


"""
Compute rest-frame UV (1500A) flux
"""

compute_MUV = True
if compute_MUV:
    print('{0}{1}').format('\n','Calculating Rest-Frame UV (1500AA) Fluxes: '),
    start_filtinterp = time.clock()
    wf = numpy.arange(1445,1551)
    tp = numpy.zeros(len(wf))
    tp[(wf>=1450) & (wf<1550)] = 1.0

    #Find SED wavelength entries within filter range
    wff = numpy.array([wf[0] < wave[i] < wf[-1] for i in range(len(wave))])
    wft = wave[wff]

    #Interpolate to find throughput values at new wavelength points
    tpt = griddata(wf,tp,wft)

    #Join arrays and sort w.r.t to wf
    wf = numpy.concatenate((wf,wft))
    tp = numpy.concatenate((tp,tpt))

    order = numpy.argsort(wf)
    wf = wf[order]
    tp = tp[order]

    dwf = numpy.diff(wf)
    nwf = len(wf)

    tpwf = tp/wf

    f_mean2 = numpy.dot(dwf,(tpwf[:nwf-1]+tpwf[1:])/2)

    tpwf = tp*wf #Reassign tpwf as product

    

    start_conv = time.clock()
    for zi in range(len(z)):
        wf1 = wf/z1[zi]

        WR = 0.

        for i in range(nwf):

            #Interpolation indices
            j = numpy.where(wave<wf1[i])[0][-1]

            a = (wf1[i] - wave[j])/(wave[j+1]-wave[j])
            tpa = tpwf[i]*((1-a)*(SED[j,:ai[zi]+1,:]*lyman_abs[j,zi]) + a*SED[j+1,:ai[zi]+1,:]*lyman_abs[j+1,zi])

            if i != 0:
                WR += dwf[i-1]*(tpb+tpa)

            tpb = tpa


        #print WR

        F_mean_UV[zi,:ai[zi]+1,:] = WR/2/z1[zi]/f_mean2/2.997925e18

    print('Done')


print('{0}').format('Converting Flux arrays to AB Magnitudes: '),
#Convert fluxes to AB magnitudes
Mags = numpy.empty(F_mean.shape)
Mags = AB0 - 2.5*numpy.log10(F_mean) - 48.6

MUV = numpy.empty(F_mean_UV.shape)
MUV = AB0 - 2.5*numpy.log10(F_mean_UV) - 48.6

#Store mass-to-light ratio for selected filter

MLRb = F_mean[mlr,0,:,:,:,:]**-1
STRx = numpy.empty(F_mean[mlr,0,:,:,:,:].shape)
for i in range(S[2]):
    STRx[:,i,:,:] = STR
MLR = STRx/F_mean[mlr,0,:,:,:,:]

print('Done')

"""
SECTION 3

"""

S = Mags.shape
MS = numpy.copy(Mags[tot,:,:,:,:,:])
MB = numpy.copy(Mags[tot,:,:,:,:,:])
SFRb = numpy.copy(Mags[tot,:,:,:,:,:])

for r in range(S[1]):
    MS[r,:] = numpy.power(10,((MS[r,:] + dmo[r])/2.5))
    MB[r,:] = numpy.power(10,((MB[r,:] + dmo[r])/2.5))
    SFRb[r,:] = numpy.power(10,((SFRb[r,:] + dmo[r])/2.5))
   
    for t in range(S[3]):
        MS[r,:,t,:,:]=(MS[r,:,t,:,:])*STR
        SFRb[r,:,t,:,:]=(SFRb[r,:,t,:,:])*SFR


print('{0} {1}{2}').format('Saving output binaries to',output_binary,':'),
if os.path.isfile(output_binary+'.main'+'.npz'):
    os.remove(output_binary+'.main'+'.npz')
if os.path.isfile(output_binary+'.mags'+'.npy'):
    os.remove(output_binary+'.mags'+'.npy')
if os.path.isfile(output_binary+'.fluxes'+'.npy'):
    os.remove(output_binary+'.fluxes'+'.npy')

numpy.savez(output_binary+'.main',parameters=parameters,z=z,filters=files,MS=MS,MB=MB,MLR=MLR,SFR=SFRb,Mshape=S,MUV=MUV)
numpy.save(output_binary+'.mags',Mags)
numpy.save(output_binary+'.fluxes',F_mean)
print('Done')

print('{0}{1}{2}{3:s} {4:.1f}').format('\n','All finished!','\n','Time elapsed:',(time.clock()-start_time))


sys.stderr = original_stderr
f.close() 
