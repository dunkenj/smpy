import numpy
import array
import re,os,sys, copy
from glob import glob
from scipy.interpolate import griddata
from scipy.integrate import simps,quad
from sm_functions import read_ised,read_ised2,calc_lyman,calc_beta

from astropy import cosmology as C

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


class SSP:
    def __init__(self):
        self.files = glob(ised_input+'*.ised')
        self.files.sort()
        self.ta_arr = []
        self.metal_arr = []
        self.iw_arr = []
        self.wave_arr = []
        self.sed_arr = []
        self.strm_arr = []
        self.rmtm_arr = []
        #Set up 
        for file in self.files:
            path = self.files[0]
            print path
            data = read_ised2(path)[0]
            ta, metal, iw, wave, sed, strm, rmtm = data
            self.ta_arr.append(ta)
            self.metal_arr.append(metal)
            self.iw_arr.append(iw)
            self.wave_arr.append(wave)
            self.sed_arr.append(sed)
            self.strm_arr.append(strm)
            self.rmtm_arr.append(rmtm)

        #Find closest match for each tg value in ta - set tg to these values
        
        nebular = numpy.loadtxt(params.neb_file,skiprows=1)
        self.neb_cont = nebular[:,1]
        self.neb_hlines = nebular[:,2]
        self.neb_metal = nebular[:,3:]
        self.neb_wave = nebular[:,0]

    def build(self,age,sfh,dust,metal,fesc=1.):
        self.tg = age*1e9
        self.tau = sfh*1e9
        self.tauv = dust
        self.mi = int(abs(metal))
        self.fesc = fesc
        
        mu = params.mu
        epsilon = params.epsilon
        
        self.ta = self.ta_arr[self.mi]
        self.wave = self.wave_arr[self.mi]
        
        [T1,T2] = numpy.meshgrid(self.tg,self.ta)
        tgi = numpy.argmin(numpy.abs(self.tg-self.ta))
        self.tg = self.ta[tgi]
        
        add_nebular = True
        
        if len(self.neb_wave) != len(self.wave):
            self.neb_cont = griddata(self.neb_wave,self.neb_cont,self.wave)
            self.neb_hlines = griddata(self.neb_wave,self.neb_hlines,self.wave)
            neb_metaln = numpy.zeros((len(self.wave),3))
            for i in range(3):
                neb_metaln[:,i] = griddata(self.neb_wave,self.neb_metal[:,i],self.wave)
            self.neb_metal = neb_metaln
            self.neb_wave = self.wave
        #quietprint("Metallicity "+str(self.mi+1)+":")
    #print ".ised file: "+files[abs(SSP)]
        sed = self.sed_arr[self.mi]
        strm = self.strm_arr[self.mi]
        rmtm = self.rmtm_arr[self.mi]
        self.iw = self.iw_arr[self.mi]
        metal=str((self.metal_arr[self.mi]))[12:-3].strip()
        #quietprint(metal[self.mi] + "\nInclude nebular emission: " + str(add_nebular))

        SSP_Z = float(re.split("Z=?",metal)[1])
        if SSP_Z <= 0.0004: neb_z = 0
        elif SSP_Z > 0.0004 and SSP_Z <= 0.004: neb_z = 1
        elif SSP_Z > 0.004: neb_z = 2
        

        if params.dust_model == "charlot":
            ATT = numpy.empty([len(self.wave),len(self.ta)])
            tv = (self.tauv*numpy.ones(len(self.ta)))
            tv[self.ta>1e7] = mu*self.tauv
            lam = numpy.array((5500/self.wave)**0.7)
            ATT[:,:] = (numpy.exp(-1*numpy.outer(lam,tv)))

        elif params.dust_model == "calzetti":
            ATT = numpy.ones([len(self.wave),len(self.ta)])

            k = numpy.zeros_like(self.wave)
            w1 = [self.wave < 6300]
            w2 = [self.wave >= 6300]
            w_u = self.wave/1e4

            k[w2] = 2.659*(-1.857 + 1.040/w_u[w2])
            k[w1] = 2.659*(-2.156 + (1.509/w_u[w1]) - (0.198/w_u[w1]**2) + (0.011/w_u[w1]**3))
            k += 4.05

            k[k < 0.] = 0.
            tv = self.tauv*k/4.05
            for ti in range(0,len(self.ta)):
                ATT[:,ti] *= numpy.power(10,-0.4*tv)
                
        
        """
        SECTION 1
        First calculate and store those parameters that are functions of the age array 
        'ta' only - these are the same for every model to be made. The parameters are 
        the age array TP, the time interval array DT, the interpolation coefficient 
        'a' and the interpolation indices J. Each are stored in cell arrays of size ks,
        with the data corresponding to the original age array first, and the 
        interpolated data second.
        """
        TP = {}
        A = {}
        J = {}
        DT = {}

        for ai in range(tgi+1):
            #Calculate taux2: the reverse age array; remove those values which
            #are less than the first non-zero entry of taux1 - these values
            #are treated differently in the original BC code
            taux1 = self.ta[:ai+1]
            taux2 = self.ta[ai]-self.ta[ai::-1]
            if max(taux1) > 0.:
                taux2 = numpy.delete(taux2,numpy.where(taux2<taux1[numpy.flatnonzero(taux1)[0]]))
            #Remove values common to taux1 and taux2; calulate array TP


            [T1,T2] = numpy.meshgrid(taux1,taux2)
            [i,j] = numpy.where(T1-T2==0)
            taux2 = numpy.delete(taux2, i)
            TP[ai] = self.ta[ai]-numpy.concatenate((taux1,taux2),axis=0)
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
                [T1,T2] = numpy.meshgrid(self.ta,taux2)
                T = T1-T2
                T[numpy.where(T<=0)] = 0
                T[numpy.where(T!=0)] = 1
                T = numpy.diff(T,1,1)
                (i,J[ai]) = T.nonzero()

                A[ai] = numpy.log10(taux2/self.ta[J[ai]])/numpy.log10(self.ta[J[ai]+1]/self.ta[J[ai]])

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


        SED = numpy.empty([len(self.wave)])
        Nlyman = numpy.empty([1])
        Nlyman_final = numpy.empty([1])
        beta = numpy.empty([1])
        norm = numpy.empty([1])
        STR = numpy.empty([tgi+1])
        SFR = numpy.empty([tgi+1])
        W = {}
       # metal=[str((self.data[1]))[12:-3].strip()]*len(params.metallicities)

        RMr = numpy.empty([tgi+1])
        PRr = numpy.empty([tgi+1])
        URr = numpy.empty([tgi+1])
        Tr = numpy.empty([tgi+1])
        

        """
        SECTION 2
        Now calculate the integration coefficients w, and store them in the
        cell array W. Also calculate the stellar mass fraction str. The so
        array is expanded and used by each successive iteration of the inner
        loop (ai). The outer loop repeats the operation for each tau value.

        """

        prgas = numpy.zeros(tgi+1)

        for ai in xrange(tgi+1):
            j = J[ai]   #Interpolation indices
            tp = TP[ai] #Integration timescale

            pgas = numpy.zeros_like(tp)
            if ai ==0:
                prgas = numpy.zeros_like(self.ta)
            else:
                i = numpy.where(tp<=self.ta[ai-1])
                ii = numpy.where(tp>self.ta[ai-1])
                pgas[i] = griddata(self.ta,prgas,tp[i])
                pgas[ii] = prgas[ai-1]
            #print prgas[ai]
    
            tbins = numpy.logspace(0,numpy.log10(max(tp)),10000)
            npgas = numpy.zeros_like(tbins)
            if self.tau > 0.:
                sr = (1 + epsilon*pgas)*numpy.exp(-1*tp/self.tau)/abs(self.tau)
                norma = 1
                if len(sr) > 1:
                    i = numpy.where(tbins <= self.ta[ai-1])
                    ii = numpy.where(tbins > self.ta[ai-1])
                    npgas[i] = griddata(self.ta,prgas,tbins[i])
                    npgas[ii] = prgas[ai-1]
                    norma = simps((1+ epsilon*npgas)*numpy.exp(-1*tbins/self.tau)/abs(self.tau),tbins)
                    sr /= norma

            elif self.tau < 0.:
                sr = numpy.exp(-1*tp/self.tau)/abs(self.tau)
                norma = 1
                if len(sr) > 1:
                    norma = simps(numpy.exp(-1*tbins/self.tau)/abs(self.tau),tbins)
                    sr /= norma
                #print sr[0]

            w = sr*DT[ai]/2
            w1 = numpy.array(w[:ai+1])
            W[0,ai] = w1

            strr = numpy.array(numpy.dot(w1,strm[:ai+1]))
            rm = numpy.array(numpy.dot(w1,rmtm[:ai+1]))

            l = len(A[ai])
            if l>0:

                w2 = w[ai+1:ai+l+1]
                wa = w2*A[ai]
                wb = w2-wa

                W[1,ai] = wa
                W[2,ai] = wb
                strr += (numpy.dot(wb,strm[j]) + numpy.dot(wa,strm[j+1]))
                rm += (numpy.dot(wb,rmtm[j]) + numpy.dot(wa,rmtm[j+1]))


            if strr > 1: strr= 1

            if self.tau > 0.:
                ugas = numpy.exp(-1*self.ta[ai]/self.tau)
            elif self.tau < 0.:
                ugas = numpy.exp(-1*self.ta[ai]/self.tau)/numpy.exp(-1*max(self.ta)/self.tau)
                #ugas = 1.
            #Processed gas = gas formed into stars - mass in stars - remnants
            prgas[ai] = 1 - ugas - strr -rm
            if prgas[ai] < 0.: prgas[ai] = 0

            #print prgas[ai]
            URr[ai] = ugas
            PRr[ai] = prgas[ai]
            RMr[ai] = rm
            Tr[ai] = simps(numpy.exp(-1*numpy.sort(tp)/self.tau)/self.tau,numpy.sort(tp))

            STR[ai] = strr
            if self.tau > 0:
                SFR[ai] = (1 + epsilon*prgas[ai])*numpy.exp(-self.ta[ai]/self.tau)/abs(self.tau)/norma
            elif self.tau < 0:
                SFR[ai] = numpy.exp(-ta[ai]/self.tau)/abs(self.tau)/norma
            #print SFR[ai,ti,mi]
    
            SFR[ai] /= STR[ai]


        """
        SECTION 3
        Finally, for each tauv/tau/tg combination, perform a weighted
        sum of the S.S.params. spectral energy distribution 'sed1' to obtain the
        model S.E.D. 'y'. Add each record to the SED array.
        """

        sed1 = sed*ATT #dust-attenuated SED

        ai = tgi
        y = numpy.zeros([1,self.iw])
        y_nodust = numpy.zeros([1,self.iw])
        j = J[ai]


        w1 = W[0,ai]
        wa = W[1,ai]
        wb = W[2,ai]

        for i in range(ai):
            y += (w1[i]*sed1[:,i])
            y_nodust += (w1[i]*sed[:,i])

        for i in range(len(wb)):
            y += (wb[i]*sed1[:,j[i]] + wa[i]*sed1[:,j[i]+1])
            y_nodust += (wb[i]*sed[:,j[i]] + wa[i]*sed[:,j[i]+1])


        Nly = calc_lyman(self.wave,y_nodust[0])
        if Nly > 0.:
            Nlyman = numpy.log10(Nly)
        else:
            Nlyman = 0.

        if add_nebular:
            total = self.neb_cont + self.neb_hlines + self.neb_metal[:,neb_z]
            total *= 2.997925e18/(self.wave**2) #Convert to Flambda
            total *= (Nly*(1-self.fesc))

            y += total
    
        Nly = calc_lyman(self.wave,y[0])
        if Nly > 0.:
            Nlyman_final = numpy.log10(Nly) + 33. + numpy.log10(3.826)
        else:
            Nlyman_final = 0.


        beta = calc_beta(self.wave,y[0])
        #print ai,ai1
        #print STR[ai1,ti,mi]
        SED[:] = y/STR[ai] #normalised to 1 solar mass
        norm = simps(numpy.exp(-1*numpy.logspace(0,numpy.log10(self.ta[tgi]),10000)/self.tau),numpy.logspace(0,numpy.log10(self.ta[tgi]),10000))

        STR = STR[tgi]
        SFR = SFR[tgi]
        
        self.SED = SED
        self.SFR = SFR
        self.STR = STR
        self.beta = beta
        self.Nly = Nlyman_final
        
    def __str__(self):
        AA = ['Age', 'SFH Tau', 'Dust Tau', 'SFR',  'Beta']
        BB = [self.tg/1e9, self.tau/1e9, self.tauv, self.SFR, self.beta]
        CC = ['Gyr', 'Gyr', 'Av', 'Ms/yr','']
        output = ['{:>10s}'.format(AA[i]) +': '+ '{:<.3g}'.format(BB[i]) +' '+CC[i] for i in range(len(AA))]
        
        return '\n'.join(output)
    
    def __add__(self,other):
        if isinstance(other,SSP):
            self.SED += other.SED
            self.SFR += other.SFR 
            self.Nly = numpy.log10( 10**self.Nly + 10**other.Nly )
            self.beta = calc_beta(self.wave,self.SED)
        return self
    
    def __mul__(self,other):
        self.SED *= other
        self.SFR *= other
        self.Nly += numpy.log10(other)
        
    def __rmul__(self,other):
        self.SED *= other
        self.SFR *= other
        self.Nly += numpy.log10(other)
        
class Filter(object):
    def __init__(self):
        self.wave = []
        self.response = []
        
class FileFilter(Filter):
    def __init__(self,filepath):
        self.path = filepath
        
        try:
            data = numpy.loadtxt(self.path)
            self.wave = data[:,0]
            self.response = data[:,1]
        
        except:
            print 'Ohhhhh dear.'
        
        
class TophatFilter(Filter):
    def __init__(self, centre, width, steps = 200):
        self.centre = centre
        self.width = width
        self.steps = steps

        upper, lower = centre+width, centre-width
        resp_upper, resp_lower = centre+(width*0.5), centre-(width*0.5)

        self.wave = numpy.linspace(lower,upper,steps)
        self.response = numpy.zeros_like(self.wave)

        tophat = (self.wave >= resp_lower)*(self.wave < resp_upper)

        self.response[tophat] = 1
        
class FilterSet:
    def __init__(self,path):
        self.directory = path
        self.files = glob(self.directory)
        
        self.filters = []
        
        for file in self.files:
            self.filters.append(FileFilter(file))

    def addFileFilter(self,path):
        self.filters.append(FileFilter(path))
    
    def addTophatFilter(self,centre, width, steps = 200):
        self.filters.append(TophatFilter(centre, width, steps))


class Observe:
    def __init__(self,SED,Filters,redshift,force_age = True):
        self.SED = SED
        self.F = Filters
        self.z = redshift
        self.wave = self.SED.wave
        
        self.lyman_abs = numpy.ones(len(self.wave))
        ly_cont_w = numpy.array([(self.wave<=912.)][0])
        ly_b_w = numpy.array([(self.wave > 912.) & (self.wave <= 1026.)][0])
        ly_a_w = numpy.array([(self.wave > 1026.) & (self.wave <= 1216.)][0])
        
        dec_a = (1/(120*(1+self.z)))*quad(self.dec_a_func,
                1050*(1+self.z),1170*(1+self.z))[0]
        dec_b= (1/(95*(1+self.z)))*quad(self.dec_b_func,
                920*(1+self.z),1015*(1+self.z))[0]
        
        self.lyman_abs[ly_cont_w] = 0.
        self.lyman_abs[ly_b_w] = dec_b
        self.lyman_abs[ly_a_w] = dec_a
        
        if self.z > 0:
            self.dm = cosmo.distmod(self.z)
        else:
            self.dm = 0.
        
        if (self.SED.tg/1e9 > cosmo.age(self.z)) and force_age:
            print 'SSP age older than universe...stopping.'
        else:            
            self.fluxes = []
            self.mags = []
            
            for filt in self.F.filters:
                print filt.wave[0]
                flux, mag = self.calcFlux(filt)
                self.fluxes.append(flux)
                self.mags.append(mag)
    
    def dec_a_func(self,wave_obs):
        return numpy.exp(-1*0.0036*(numpy.power(wave_obs/1216.,3.46)))

    def dec_b_func(self,wave_obs):
        teff_beta=1.7e-3*(numpy.power(wave_obs/1026.,3.46))
        teff_gamma=1.2e-3*(numpy.power(wave_obs/973.,3.46))
        teff_delta=9.3e-4*(numpy.power(wave_obs/950.,3.46))
        #The above absorption lines dominate the average absorption over the range 
        #but higher order lines should be added in future
    
        teff_total=teff_beta+teff_gamma+teff_delta

        return numpy.exp(-1*teff_total)
        
    def calcFlux(self,filt):
        wf = filt.wave
        tp = filt.response
        z1 = self.z+1

        if len(wf) > 1000: #Re-sample large filters for performance
            wfx = numpy.linspace(wf[0],wf[-1],1000)
            tpx = griddata(wf,tp,wfx)

            wf = wfx
            tp = tpx
        

        #Find SED wavelength entries within filter range
        wff = numpy.array([wf[0] < self.wave[i] < wf[-1] 
                            for i in range(len(self.wave))])
        wft = self.wave[wff]
    
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

        wf1 = wf/z1

        WR = 0.
        for i in range(nwf):
            #Interpolation indices
            j = numpy.where(self.wave<wf1[i])[0][-1]
            
            a = (wf1[i] - self.wave[j])/(self.wave[j+1]-self.wave[j])
            tpa = (tpwf[i]*((1-a)*(self.SED.SED[j]*self.lyman_abs[j]) + 
                    a*self.SED.SED[j+1]*self.lyman_abs[j+1]))
            if i != 0:
                WR += dwf[i-1]*(tpb+tpa)

            tpb = tpa
            
        F_mean = WR/2/z1/f_mean2/2.997925e18
        #print 'F_mean shape '+ len(F_mean)
        AB0 = 5*numpy.log10(1.7684e8*1e-5)
        # dl = 10pc in Mpc
        # this factor is sqrt(4*pi*(3.0856e24)^2 Lsun)

        #Convert fluxes to AB magnitudes
        Mag = AB0 - 2.5*numpy.log10(F_mean) - 48.6
        Mag += self.dm
        
        Flux = 10**((23.9 - Mag)/2.5) #uJy
        
        return Flux, Mag 
        
#must define cosmo before calling an Observe
cosmo = C.FlatLambdaCDM(H0=70,Om0=0.3)

a = SSP()
a.build(0.5,0.5,1,2) #(Age, Tau, Dust, Metallicity Index)
#print a

b = SSP()
b.build(1.5,0.01,3,4)
#print b

a*1e9 #Multiply by a factor, equivalent to *=
b*5e9 
a+b #Add b to a (a += b)
print a

filt_dir = 'GOODS-S_18_FilterCurves/Filter*.txt'

Filts = FilterSet(filt_dir)
Filts.addTophatFilter(1500,100)

AA = Observe(a,Filts,2) # Observe a (built from two initial SEDs) 
                        # through the filters in Filts at redshift of z = 2 
BB = Observe(b,Filts,2) 
