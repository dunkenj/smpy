"""
Kenneth Duncan, 2012
Based on Alice Mortlock, Asa Bluck and Steven Bevan, (2011,2009,2007)
This code produces synthetic spectral energy distributions for model galaxies. The process precisely mimmicks that of Bruzual & Charlot's csp_galaxev code - except that only a selection of model ages are output.

Data is read from the Bruzual & Charlot S.S.params. files for each metallicity in turn. An array of star formation timescales 'tau' can be specified in units of Gyr, as well as an array of attenuation 'tauv' values.Only single values of the gas recycling (epsilon) and ISM attenuation (mu) parameters can be specified.
"""

import numpy
import array
import re
import os
import sys
from glob import glob
from scipy.interpolate import griddata
from scipy.integrate import simps
from sm_functions import read_ised,read_ised2,calc_lyman,calc_beta


from astropy import units as u
from astropy import constants as c
from astropy import cosmology as cos
from astropy.utils.console import ProgressBar
cosmo = cos.FlatLambdaCDM(H0=70, Om0=0.3)


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

output = params.ssp_output
ised_input = params.ssp_input

files = glob(ised_input+'*.ised')
files.sort()
quietprint('SSP binary files found:')
for file in files:
    quietprint(file)
quietprint('')

tau = numpy.array(params.tau)*1e9
tg = numpy.array(params.tg)*1e9
tauv = numpy.array(params.tauv)
mu = params.mu
epsilon = params.epsilon

#files = [files[a] for a in params.metallicities]

add_nebular = False
if params.add_nebular:
    add_nebular = True

# Load first SSP file to obtain age array 'ta' which is required
# for section 1. It is assumed that 'ta' and 'wave' are the same
# for all files.

data = read_ised2(files[abs(params.metallicities[0])-1])[0]

ta, metal, iw, wave, sed, strm, rmtm = data

#Find closest match for each tg value in ta - set tg to these values
[T1,T2] = numpy.meshgrid(tg,ta)
tdiff, tgi = numpy.min(abs(T1-T2),0),numpy.argmin(abs(T1-T2),0)
tg = ta[tgi]


#Calculate coefficients for the attenuation of the SEDs using
#2-component model of Charlot & Fall (2000). Used in Section 3

if params.dust_model == "charlot":
    ATT = numpy.empty([len(tauv),len(wave),len(ta)])
    for tvi in range(0,len(tauv)):
        tv = (tauv[tvi]*numpy.ones(len(ta)))
        tv[ta>1e7] = mu*tauv[tvi]
        lam = numpy.array((5500/wave)**0.7)
        ATT[tvi,:,:] = (numpy.exp(-1*numpy.outer(lam,tv)))

elif params.dust_model == "calzetti":
    ATT = numpy.ones([len(tauv),len(wave),len(ta)])

    k = numpy.zeros_like(wave)
    w1 = [wave < 6300]
    w2 = [wave >= 6300]
    w_u = wave/1e4

    k[w2] = 2.659*(-1.857 + 1.040/w_u[w2])
    k[w1] = 2.659*(-2.156 + (1.509/w_u[w1]) - (0.198/w_u[w1]**2) + (0.011/w_u[w1]**3))
    k += 4.05

    k[k < 0.] = 0.

    for tvi in range(0,len(tauv)):
        tv = tauv[tvi]*k/4.05
        for ti in range(0,len(ta)):
            ATT[tvi,:,ti] *= numpy.power(10,-0.4*tv)
            
#def Pei(wave):
#    ai = [185., 27., 0.005, 0.01, 0.12, 0.03]
#    lambda_i = [0.042, 0.08, 0.22, 9.7, 18., 25.]


    
"""
def Calzetti(wave):
    k = numpy.zeros_like(wave)
    w1 = [wave < 6300]
    w2 = [wave >= 6300]
    w_u = wave/1e4

    k[w2] = 2.659*(-1.857 + 1.040/w_u[w2])
    k[w1] = 2.659*(-2.156 + (1.509/w_u[w1]) - (0.198/w_u[w1]**2) + (0.011/w_u[w1]**3))
    k += 4.05

    k[k < 0.] = 0.
    tv = k/4.05   
    return tv

def red(wave,C1,C2,C3,C4,x0,gamma,Rv):
    x = 1e4/wave
    xsq = x**2
    #SMC Coefficients   
    #C1 = C1/Rv
    #C2 /= Rv
    #C3 /= Rv
    #C4 /= Rv
    
    #LMC
    C1 = -0.890
    C2 = 0.998
    C3 = 2.719
    C4 = 0.400
    x0 = 4.579
    gamma = 0.934
    Rv = 3.41
    
    D = xsq / ((xsq - x0**2)**2 + xsq*(gamma**2))
    #D=xsq*((xsq-x0*x0)**2+xsq*gamma*gamma)**-2
    Cx4 = (x >= 5.9)
    F = 0.5392 *((x[Cx4]-5.9)**2) + 0.05644*((x[Cx4]-5.9)**3)
    Elam = C1 + C2*x + C3*D
    Elam[Cx4] += C4*F
    
    # Elam = E(lam - V)/E(B-V)
    Alam = (Elam/Rv) +1    
    return Alam

def SMC(wave):            
#elif params.dust_model == 'smc':
    #ATT = numpy.ones([len(tauv),len(wave),len(ta)])
    x = 1e4/wave
    xsq = x**2
    #SMC Coefficients   
    C1 = -4.959
    C2 = 2.264
    C3 = 0.389
    C4 = 0.461
    x0 = 4.6
    gamma = 1
    Rv = 2.74
    #C1 = C1/Rv
    #C2 /= Rv
    #C3 /= Rv
    #C4 /= Rv
    
    #LMC
    C1 = -0.890
    C2 = 0.998
    C3 = 2.719
    C4 = 0.400
    x0 = 4.579
    gamma = 0.934
    Rv = 3.41
    
    D = xsq / ((xsq - x0**2)**2 + xsq*(gamma**2))
    #D=xsq*((xsq-x0*x0)**2+xsq*gamma*gamma)**-2
    Cx4 = (x >= 5.9)
    F = 0.5392 *((x[Cx4]-5.9)**2) + 0.05644*((x[Cx4]-5.9)**3)
    Elam = C1 + C2*x + C3*D
    Elam[Cx4] += C4*F
    
    # Elam = E(lam - V)/E(B-V)
    Alam = (Elam/Rv) +1    
    return Alam
    
def LMC(wave):            
#elif params.dust_model == 'smc':
    #ATT = numpy.ones([len(tauv),len(wave),len(ta)])
    x = 1e4/wave
    xsq = x**2
    
    #SMC Coefficients   
    C1 = -4.959
    C2 = 2.264
    C3 = 0.389
    C4 = 0.461
    x0 = 4.6
    gamma = 1
    Rv = 2.74
    
    #LMC
    C1 = -0.890
    C2 = 0.998
    C3 = 2.719
    C4 = 0.400
    x0 = 4.579
    gamma = 0.934
    Rv = 3.41
    
    #C1 = C1/Rv +1
    #C2 /= Rv
    #C3 /= Rv
    #C4 /= Rv
    
    D = xsq / ((xsq - x0**2)**2 + xsq*(gamma**2))
    #D=xsq*((xsq-x0*x0)**2+xsq*gamma*gamma)**-2
    Cx4 = (x >= 5.9)
    F = (0.5392 * (x[Cx4]-5.9)**2) + (0.05644 * (x[Cx4]-5.9)**3)
    
    Elam = C1 + C2*x + C3*D
    Elam[Cx4] += C4*F
    
    # Elam = E(lam - V)/E(B-V)
    Alam = (Elam/Rv) +1  
    return Alam
        
    ext = numpy.maximum(ext,0)
    
    for tvi in range(0,len(tauv)):
        tv = tauv[tvi]*ext
        for ti in range(0,len(ta)):
            ATT[tvi,:,ti] *= numpy.power(10,-0.4*tv)
"""

if params.add_nebular:
    nebular = numpy.loadtxt(params.neb_file,skiprows=1)
    neb_cont = nebular[:,1]
    neb_hlines = nebular[:,2]
    neb_metal = nebular[:,3:]
    neb_wave = nebular[:,0]

    if len(neb_wave) != len(wave):
        neb_cont = griddata(neb_wave,neb_cont,wave)
        neb_hlines = griddata(neb_wave,neb_hlines,wave)
        neb_metaln = numpy.zeros((len(wave),3))
        neb_metaln[:,0] = griddata(neb_wave,neb_metal[:,0],wave)
        neb_metal = neb_metaln

"""
SECTION 1
First calculate and store those parameters that are functions of the age array 'ta' only - these are the same for every model to be made. The parameters are the age array TP, the time interval array DT, the interpolation coefficient 'a' and the interpolation indices J. Each are stored in cell arrays of size ks, with the data corresponding to the original age array first, and the interpolated data second.
"""
TP = {}
A = {}
J = {}
DT = {}

for ai in range(max(tgi)+1):
    #Calculate taux2: the reverse age array; remove those values which
    #are less than the first non-zero entry of taux1 - these values
    #are treated differently in the original BC code
    taux1 = ta[:ai+1]
    taux2 = ta[ai]-ta[ai::-1]
    if max(taux1) > 0.:
        taux2 = numpy.delete(taux2,numpy.where(taux2<taux1[numpy.flatnonzero(taux1)[0]]))
    #Remove values common to taux1 and taux2; calulate array TP


    [T1,T2] = numpy.meshgrid(taux1,taux2)
    [i,j] = numpy.where(T1-T2==0)
    taux2 = numpy.delete(taux2, i)
    TP[ai] = ta[ai]-numpy.concatenate((taux1,taux2),axis=0)
    l = len(taux2)

    #If taux2 has entries, calculate the interpolation parameters a and J.
    #The indicies correspond to those values of 'ta' which are just below
    #the entries in taux2. They are calculated by taking the difference
    #between the two arrays, then finding the last negative entry in the
    #resulting array.

    if l == 0:
        J[ai] = numpy.array([])
        A[ai] = numpy.array([])
    if l>0:
        [T1,T2] = numpy.meshgrid(ta,taux2)
        T = T1-T2
        T[numpy.where(T<=0)] = 0
        T[numpy.where(T!=0)] = 1
        T = numpy.diff(T,1,1)
        (i,J[ai]) = T.nonzero()

        A[ai] = numpy.log10(taux2/ta[J[ai]])/numpy.log10(ta[J[ai]+1]/ta[J[ai]])

    #Calculate age difference array: the taux arrays are joined and
    #sorted, the differences calculated, then rearranged back to the order
    #of the original taux values.
    taux = numpy.concatenate((taux1,taux2),axis=0)
    taux.sort()

    b = numpy.searchsorted(taux,taux1)
    c = numpy.searchsorted(taux,taux2)
    order = numpy.concatenate((b,c))

    d = numpy.diff(taux)
    dt = numpy.append(d,0) + numpy.append(0,d)
    DT[ai] = numpy.copy(dt[order])


SED = numpy.empty([len(wave),len(tgi),len(tauv),len(tau),len(params.metallicities)])
MF = numpy.empty([len(wave),len(tgi),len(tauv),len(tau),len(params.metallicities)])
Nlyman = numpy.empty([len(tgi),len(tauv),len(tau),len(params.metallicities)])
Nlyman_final = numpy.empty([len(tgi),len(tauv),len(tau),len(params.metallicities)])
beta = numpy.empty([len(tgi),len(tauv),len(tau),len(params.metallicities)])
norm = numpy.empty([len(tgi),len(tauv),len(tau),len(params.metallicities)])
STR = numpy.empty([max(tgi)+1,len(tau),len(params.metallicities)])
SFR = numpy.empty([max(tgi)+1,len(tau),len(params.metallicities)])
W = {}
metal=[str((data[1]))[12:-3].strip()]*len(params.metallicities)

RMr = numpy.empty([max(tgi)+1])
PRr = numpy.empty([max(tgi)+1])
URr = numpy.empty([max(tgi)+1])
Tr = numpy.empty([max(tgi)+1])
for mi in range(len(params.metallicities)):
    SSP = abs(params.metallicities[mi])-1
    if params.metallicities[mi] < 0:
        add_nebular = True
    else:
        add_nebular = False

        quietprint("Metallicity "+str(mi+1)+":")
    #print ".ised file: "+files[abs(SSP)]
 
    data = read_ised2(files[SSP])[0]
    sed = data[4]
    strm = data[5]
    rmtm = data[6]
    metal[mi]=str((data[1]))[12:-3].strip()
    quietprint(metal[mi] + "\nInclude nebular emission: " + str(add_nebular))

    SSP_Z = float(re.split("Z=?",metal[mi])[1])
    if SSP_Z <= 0.0004: neb_z = 0
    elif SSP_Z > 0.0004 and SSP_Z <= 0.004: neb_z = 1
    elif SSP_Z > 0.004: neb_z = 2

    """
    SECTION 2
    Now calculate the integration coefficients w, and store them in the
    cell array W. Also calculate the stellar mass fraction str. The so
    array is expanded and used by each successive iteration of the inner
    loop (ai). The outer loop repeats the operation for each tau value.

    """

    for ti in range(len(tau)):
        prgas = numpy.zeros(max(tgi)+1)

        for ai in xrange(max(tgi)+1):
            j = J[ai]   #Interpolation indices
            tp = TP[ai] #Integration timescale

            pgas = numpy.zeros_like(tp)
            if ai ==0:
                prgas = numpy.zeros_like(ta)
            else:
                i = numpy.where(tp<=ta[ai-1])
                ii = numpy.where(tp>ta[ai-1])
                pgas[i] = griddata(ta,prgas,tp[i])
                pgas[ii] = prgas[ai-1]
            #print prgas[ai]
            
            tbins = numpy.logspace(0,numpy.log10(max(tp)),10000)
            npgas = numpy.zeros_like(tbins)
            if tau[ti] > 0.:
                sr = (1 + epsilon*pgas)*numpy.exp(-1*tp/tau[ti])/abs(tau[ti])
                norma = 1
                if len(sr) > 1:
                    i = numpy.where(tbins <= ta[ai-1])
                    ii = numpy.where(tbins > ta[ai-1])
                    npgas[i] = griddata(ta,prgas,tbins[i])
                    npgas[ii] = prgas[ai-1]
                    norma = simps((1+ epsilon*npgas)*numpy.exp(-1*tbins/tau[ti])/abs(tau[ti]),tbins)
                    sr /= norma

            elif tau[ti] < 0.:
                sr = numpy.exp(-1*tp/tau[ti])/abs(tau[ti])
                norma = 1
                if len(sr) > 1:
                    norma = simps(numpy.exp(-1*tbins/tau[ti])/abs(tau[ti]),tbins)
                    sr /= norma
                #print sr[0]

            w = sr*DT[ai]/2
            w1 = numpy.array(w[:ai+1])
            W[0,ai,ti] = w1

            strr = numpy.array(numpy.dot(w1,strm[:ai+1]))
            rm = numpy.array(numpy.dot(w1,rmtm[:ai+1]))

            l = len(A[ai])
            if l>0:

                w2 = w[ai+1:ai+l+1]
                wa = w2*A[ai]
                wb = w2-wa

                W[1,ai,ti] = wa
                W[2,ai,ti] = wb
                strr += (numpy.dot(wb,strm[j]) + numpy.dot(wa,strm[j+1]))
                rm += (numpy.dot(wb,rmtm[j]) + numpy.dot(wa,rmtm[j+1]))


            if strr > 1: strr= 1

            if tau[ti] > 0.:
                ugas = numpy.exp(-1*ta[ai]/tau[ti])
            elif tau[ti] < 0.:
                ugas = numpy.exp(-1*ta[ai]/tau[ti])/numpy.exp(-1*max(ta)/tau[ti])
                #ugas = 1.
            #Processed gas = gas formed into stars - mass in stars - remnants
            prgas[ai] = 1 - ugas - strr -rm
            if prgas[ai] < 0.: prgas[ai] = 0

            #print prgas[ai]
            URr[ai] = ugas
            PRr[ai] = prgas[ai]
            RMr[ai] = rm
            Tr[ai] = simps(numpy.exp(-1*numpy.sort(tp)/tau[ti])/tau[ti],numpy.sort(tp))

            STR[ai,ti,mi] = strr
            if tau[ti] > 0:
                SFR[ai,ti,mi] = (1 + epsilon*prgas[ai])*numpy.exp(-ta[ai]/tau[ti])/abs(tau[ti])/norma
            elif tau[ti] < 0:
                SFR[ai,ti,mi] = numpy.exp(-ta[ai]/tau[ti])/abs(tau[ti])/norma
            #print SFR[ai,ti,mi]
            
            SFR[ai,ti,mi] /= STR[ai,ti,mi]


    """
    SECTION 3
    Finally, for each tauv/tau/tg combination, perform a weighted
    sum of the S.S.params. spectral energy distribution 'sed1' to obtain the
    model S.E.D. 'y'. Add each record to the SED array.
    """
    with ProgressBar(len(tauv)*len(tgi)*len(tau)) as bar:
        for tvi in range(len(tauv)):
            sed1 = sed*ATT[tvi] #dust-attenuated SED
            for ti in range(len(tau)):
                for ai1 in range(len(tgi)):
                    ai = tgi[ai1]
                    cf = numpy.zeros([1,iw])
                    y = numpy.zeros([1,iw])
                    y_nodust = numpy.zeros([1,iw])
                    j = J[ai]


                    w1 = W[0,ai,ti]
                    wa = W[1,ai,ti]
                    wb = W[2,ai,ti]

                    for i in range(ai):
                        y += (w1[i]*sed1[:,i])
                        y_nodust += (w1[i]*sed[:,i])

                    for i in range(len(wb)):
                        y += (wb[i]*sed1[:,j[i]] + wa[i]*sed1[:,j[i]+1])
                        y_nodust += (wb[i]*sed[:,j[i]] + wa[i]*sed[:,j[i]+1])


                    Nly = calc_lyman(wave,y_nodust[0])
                    if Nly > 0.:
                        Nlyman[ai1,tvi,ti,mi] = numpy.log10(Nly)
                    else:
                        Nlyman[ai1,tvi,ti,mi] = 0.

                    if add_nebular:
                        total = neb_cont + neb_hlines + neb_metal[:,neb_z]
                        total *= 2.997925e18/(wave**2) #Convert to Flambda
                        total *= (Nly*(1-params.fesc))

                        y += total
                    
                    Nly = calc_lyman(wave,y[0])
                    if Nly > 0.:
                        Nlyman_final[ai1,tvi,ti,mi] = numpy.log10(Nly) + 33. + numpy.log10(3.826)
                    else:
                        Nlyman_final[ai1,tvi,ti,mi] = 0.


                    beta[ai1,tvi,ti,mi] = calc_beta(wave,y[0])
                    #print ai,ai1
                    #print STR[ai1,ti,mi]
                    SED[:,ai1,tvi,ti,mi] = y/STR[ai,ti,mi] #normalised to 1 solar mass
                    norm[ai1,tvi,ti,mi] = simps(numpy.exp(-1*numpy.logspace(0,numpy.log10(ta[tgi[ai1]]),10000)/tau[ti]),numpy.logspace(0,numpy.log10(ta[tgi[ai1]]),10000))
                    bar.update()

"""
Section 4

"""
b = SFR
print "Saving to numpy binaries:",
#print(SED[:,0,0,0,0])
#print(STR[0,0,0])

STR = STR[tgi,:,:]
SFR = SFR[tgi,:,:]

"""
Tr = Tr[tgi]
RMr = RMr[tgi]/Tr
PRr = PRr[tgi]/Tr
SRr = STR[:,0,0]/Tr

os.remove("temp_output.txt")
temp_file = open("temp_output.txt","w")

temp_file.write('{0} {1} {2} {3} {4}'.format('Time','Stellar Mass','Processed Gas','Remnant Mass','\n'))
for i in range(len(tgi)):
    temp_file.write('{0} {1} {2} {3} {4}'.format(ta[tgi][i],SRr[i],PRr[i],RMr[i],'\n'))

temp_file.close()
"""

parameters = dict(wavelengths = wave,
                  tg = tg,
                  Av = tauv,
                  SFH = tau,
                  SSPs = metal,
                  ISM_attenuation_fraction = mu,
                  Gas_recycling_parameter = epsilon)

numpy.savez(params.ssp_output,parameters=parameters,SED=SED,STR=STR,SFR=SFR,Nlyman=Nlyman,Nlyman_final = Nlyman_final, beta=beta)
numpy.save(params.ssp_output+'.par',parameters)
print('Done')


sys.stderr = original_stderr
f.close()
