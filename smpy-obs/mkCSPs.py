import numpy
import array
import copy
import re,os,sys,copy
from glob import glob
from scipy.interpolate import griddata
from scipy.integrate import simps,quad
from scipy.optimize import leastsq, fsolve

#from sm_functions import read_ised,read_ised2,calc_lyman,calc_beta
from astropy import units as U
from astropy import constants as C
from astropy import cosmology as cos
cosmo = cos.FlatLambdaCDM(H0=70,Om0=0.3)

f = open("error.log", "w")
original_stderr = sys.stderr
sys.stderr = f


class ised(object):
    def __init__(self,path):
        self.file = path
        
        self.read_ised(self.file)
        
    def read_ised(self,filename):
        """
        This function reads data from Bruzual & Charlot binary format
        SSP files and returns the necessary data in an array The input files
        should be '.ised' files, either 2003 or 2007.
    
        'ks' in the binary files is slightly different between 03/07 files
        so the read length and index should be set appropriately, therefore 
        the function tries '03 format first and retries with the '07 format
        if the returned number of ages isn't as expected (e.g. 221 ages)
        """

        with open(filename,'rb') as f:
            check = array.array('i')
            check.fromfile(f,2)
        
        if check[1] == 221:
            ksl, ksi = 2, 1
            F_l, F_i = 3, 2
        else:
            ksl, ksi = 3, 2
            F_l, F_i = 5, 4
            
        with open(filename,'rb') as f:
            ks = array.array('i')
            ks.fromfile(f,ksl)

            ta = array.array('f')
            ta.fromfile(f,ks[ksi])
            self.ta = numpy.array(ta)

            tmp = array.array('i')
            tmp.fromfile(f,3)
            self.ml,self.mul,iseg = tmp

            if iseg > 0:
                tmp = array.array('f')
                tmp.fromfile(f,iseg*6)

            tmp = array.array('f')
            tmp.fromfile(f,5)
            self.totm, self.totn, self.avs, self.jo, self.tauo = tmp


            self.ids= array.array('c')
            self.ids.fromfile(f,80)

            tmp = array.array('f')
            tmp.fromfile(f,4)
            self.tcut = tmp[0]
            self.ttt = tmp[1:]

            ids = array.array('c')
            ids.fromfile(f,80)

            self.ids = array.array('c')
            self.ids.fromfile(f,80)

            self.igw = array.array('i')
            self.igw.fromfile(f,1)

            tmp = array.array('i')
            tmp.fromfile(f,F_l)

            self.iw = array.array('i')
            self.iw.fromfile(f,1)

            wave = array.array('f')
            wave.fromfile(f,self.iw[0])
            self.wave = numpy.array(wave)

            #SED Section
            self.F = array.array('i')
            self.F.fromfile(f,F_l)
            self.iw = self.F[F_i] #Number of wavelength elements

            self.sed = numpy.zeros((self.iw,ks[ksi]),dtype=numpy.float32)
            G = array.array('f')
            G.fromfile(f,self.iw)
            self.sed[:,0] = G
            ik = array.array('i')
            ik.fromfile(f,1)

            self.h = numpy.empty((ik[0],ks[ksi]),'f')
            H = array.array('f')
            H.fromfile(f,ik[0])
            self.h[:,0] = H

            for i in range(1,ks[ksi]): #Fill rest of array with SEDs
                F = array.array('i')
                F.fromfile(f,F_l)
                iw = F[F_i]

                G = array.array('f')
                G.fromfile(f,iw)
                self.sed[:,i] = G
                ik = array.array('i')
                ik.fromfile(f,1)

                H = array.array('f')
                H.fromfile(f,ik[0])
                self.h[:,i] = H

            tmp = array.array('i')
            tmp.fromfile(f,F_l)

            self.bflx = array.array('f')
            self.bflx.fromfile(f,tmp[F_i])

            tmp = array.array('i')
            tmp.fromfile(f,F_l)

            strm = array.array('f')
            strm.fromfile(f,tmp[F_i])
            self.strm = numpy.array(strm)

            tmp = array.array('i')
            tmp.fromfile(f,F_l)

            self.evf = array.array('f')
            self.evf.fromfile(f,tmp[F_i])

            tmp = array.array('i')
            tmp.fromfile(f,F_l)

            self.evf = array.array('f')
            self.evf.fromfile(f,tmp[F_i])

            tmp = array.array('i')
            tmp.fromfile(f,F_l)

            self.snr = array.array('f')
            self.snr.fromfile(f,tmp[F_i])

            tmp = array.array('i')
            tmp.fromfile(f,F_l)

            self.pnr = array.array('f')
            self.pnr.fromfile(f,tmp[F_i])

            tmp = array.array('i')
            tmp.fromfile(f,F_l)

            self.sn = array.array('f')
            self.sn.fromfile(f,tmp[F_i])

            tmp = array.array('i')
            tmp.fromfile(f,F_l)

            self.bh = array.array('f')
            self.bh.fromfile(f,tmp[F_i])

            tmp = array.array('i')
            tmp.fromfile(f,F_l)

            self.wd = array.array('f')
            self.wd.fromfile(f,tmp[F_i])

            tmp = array.array('i')
            tmp.fromfile(f,F_l)

            rmtm = array.array('f')
            rmtm.fromfile(f,tmp[F_i])
            self.rmtm = numpy.array(rmtm)

class CSP:
    def __init__(self,SSPpath = '../ssp/bc03/salpeter/lr/',
                 age=None,sfh=None,dust=None,metal_ind=None,fesc=None,
                 sfh_law='exp',dustmodel = 'calzetti',neb_cont=True,neb_met=True):
        self.SSPpath = SSPpath
        self.files = glob(self.SSPpath + '*.ised')
        self.files.sort()
        self.iseds = []
        self.ta_arr = []
        self.metal_arr = []
        self.iw_arr = []
        self.wave_arr = []
        self.sed_arr = []
        self.strm_arr = []
        self.rmtm_arr = []
        #Set up 
        for file in self.files:
            ised_binary = ised(file)
             
            self.ta_arr.append(ised_binary.ta)
            self.metal_arr.append(ised_binary.ids)
            self.iw_arr.append(ised_binary.iw)
            self.wave_arr.append(ised_binary.wave)
            self.sed_arr.append(ised_binary.sed)
            self.strm_arr.append(ised_binary.strm)
            self.rmtm_arr.append(ised_binary.rmtm)
            self.iseds.append(ised_binary)

        #Find closest match for each tg value in ta - set tg to these values
        
        nebular = numpy.loadtxt('nebular_emission.dat',skiprows=1)
        self.neb_cont = nebular[:,1]
        self.neb_hlines = nebular[:,2]
        self.neb_metal = nebular[:,3:]
        self.neb_wave = nebular[:,0]
        
        if None not in (age,sfh,dust,metal_ind):
            if fesc == None:
                self.build(age,sfh,dust,metal_ind,sfh_law=sfh_law,dustmodel=dustmodel,
                           neb_cont=neb_cont,neb_met=neb_met)
                
            else:
                self.build(age,sfh,dust,metal_ind,fesc,sfh_law,dustmodel,neb_cont,neb_met)
    
    def _sfh_exp(self,t,tau):
        sfh = numpy.exp(-1*t/tau)/abs(tau)
        return sfh

    def _sfh_pow(self,t,alpha):
        sfh = numpy.power(t/1.e9,alpha)
        return sfh
        
    def _sfh_del(self,t,tau):
        sfh = t/(tau**2)*numpy.exp(-t/tau)
        return sfh
        
    def _sfh_tru(self,t,tstop):
        sfh = numpy.ones_like(t)
        sfh[t > tstop*numpy.max(t)] = 0.
    
        sfh /= numpy.trapz(sfh,t)
        return sfh
    
    def dust_func(self,lam,ai,bi,ni,li):
        """
        Functional form for SMC, LMC and MW extinction curves of
        Pei et al. 1992
        """
        lam = numpy.array(lam) / 1e4 
        ki = numpy.power((lam / li),ni) + numpy.power((li / lam),ni) + bi
        eta_i = ai / ki
        return eta_i

    def build(self,age,sfh,dust,metal,fesc=1.,sfh_law='exp',dustmodel = 'calzetti',
              neb_cont=True,neb_met=True):
        """
        
        
        """
        self.tg = age*1.e9
        if sfh_law == 'exp':
            self.tau = sfh*1.e9
        elif sfh_law == 'del':
            self.tau = sfh*1.e9
        else:
            self.tau = sfh
        self.tauv = dust
        self.mi = int(abs(metal))
        self.fesc = fesc
        self.sfh_law = sfh_law
        self.inc_cont= neb_cont
        self.inc_met = neb_met
        self.dust_model = dustmodel
        
        mu = 0.3
        epsilon = 0.
        
        self.ta = self.ta_arr[self.mi]
        self.wave = self.wave_arr[self.mi]
        
        [T1,T2] = numpy.meshgrid(self.tg,self.ta)
        tgi = numpy.argmin(numpy.abs(self.tg-self.ta))
        self.tg = self.ta[tgi]
                
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
        #print SSP_Z,
        if SSP_Z <= 0.0004: neb_z = 0
        elif SSP_Z > 0.0004 and SSP_Z <= 0.004: neb_z = 1
        elif SSP_Z > 0.004: neb_z = 2
        #print neb_z

        if self.dust_model == "charlot":
            ATT = numpy.empty([len(self.wave),len(self.ta)])
            tv = ((self.tauv/1.0857)*numpy.ones(len(self.ta)))
            tv[self.ta>1e7] = mu*self.tauv
            lam = numpy.array((5500/self.wave)**0.7)
            ATT[:,:] = (numpy.exp(-1*numpy.outer(lam,tv)))

        elif self.dust_model == "calzetti":
            ATT = numpy.ones([len(self.wave),len(self.ta)])
            k = numpy.zeros_like(self.wave)

            w0 = [self.wave <= 1200]
            w1 = [self.wave < 6300]
            w2 = [self.wave >= 6300]
            w_u = self.wave/1e4

            x1 = numpy.argmin(numpy.abs(self.wave-1200))
            x2 = numpy.argmin(numpy.abs(self.wave-1250))
    
            k[w2] = 2.659*(-1.857 + 1.040/w_u[w2])
            k[w1] = 2.659*(-2.156 + (1.509/w_u[w1]) - (0.198/w_u[w1]**2) + (0.011/w_u[w1]**3))
            k[w0] = k[x1] + ((self.wave[w0]-1200.)  * (k[x1]-k[x2]) / (self.wave[x1]-self.wave[x2]))

            k += 4.05
            k[k < 0.] = 0.

            tv = self.tauv*k/4.05
            for ti in range(0,len(self.ta)):
                ATT[:,ti] *= numpy.power(10,-0.4*tv)

        elif self.dust_model == "calzetti2":
            ATT = numpy.ones([len(self.wave),len(self.ta)])
            k = numpy.zeros_like(self.wave)

            w0 = [self.wave <= 1000]
            w1 = [(self.wave > 1000)*(self.wave < 6300)]
            w2 = [self.wave >= 6300]
            w_u = self.wave/1e4
    
            k[w2] = 2.659*(-1.857 + 1.040/w_u[w2])
            k[w1] = 2.659*(-2.156 + (1.509/w_u[w1]) - (0.198/w_u[w1]**2) + (0.011/w_u[w1]**3))
    
            p1 = self.dust_func(self.wave,27,4,5.5,0.08) + self.dust_func(self.wave,185,90,2,0.042)
    
            k[w0] = p1[w0] / (p1[w1][0]/k[w1][0])
            k += 4.05
            k[k < 0.] = 0.
            tv = self.tauv*k/4.05
            for ti in range(0,len(self.ta)):
                ATT[:,ti] *= numpy.power(10,-0.4*tv)

        elif self.dust_model == "smc":
            ai = [185., 27., 0.005, 0.01, 0.012, 0.03]
            bi = [90., 5.5, -1.95, -1.95, -1.8, 0.]
            ni = [2., 4., 2., 2., 2., 2.]
            li = [0.042, 0.08, 0.22, 9.7, 18., 25.]
            
            eta = numpy.zeros_like(self.wave)
            for i in xrange(len(ai)):
                eta += self.dust_func(self.wave, ai[i], bi[i], ni[i], li[i])
            
            Rv = 2.93
            Ab = self.tauv * (1 + (1/Rv))
            
            print numpy.exp(self.tauv*eta)
            ATT = numpy.ones([len(self.wave),len(self.ta)])
            for ti in range(0,len(self.ta)):
                ATT[:,ti] *= numpy.power(10,-0.4*(Ab*eta))
                #Offset added to renormalise from B to V band
                #ATT[:,ti] *= numpy.exp(-1*self.tauv*eta)
        
        elif self.dust_model == "lmc":
            ai = [175., 19., 0.023, 0.005, 0.006, 0.02]
            bi = [90., 4.0, -1.95, -1.95, -1.8, 0.]
            ni = [2., 4.5, 2., 2., 2., 2.]
            li = [0.046, 0.08, 0.22, 9.7, 18., 25.]
            
            eta = numpy.zeros_like(self.wave)
            for i in xrange(len(ai)):
                eta += self.dust_func(self.wave, ai[i], bi[i], ni[i], li[i])

            Rv = 3.16
            Ab = self.tauv * (1 + (1/Rv))
            
            ATT = numpy.ones([len(self.wave),len(self.ta)])
            for ti in range(0,len(self.ta)):
                ATT[:,ti] *= numpy.power(10,-0.4*(Ab*eta))
                #Offset added to renormalise from B to V band
                #ATT[:,ti] *= numpy.exp(-1*self.tauv*eta)
                 
        elif self.dust_model == "mw":
            ai = [165., 14., 0.045, 0.002, 0.002, 0.012]
            bi = [90., 4., -1.95, -1.95, -1.8, 0.]
            ni = [2., 6.5, 2., 2., 2., 2.]
            li = [0.047, 0.08, 0.22, 9.7, 18., 25.]
            
            eta = numpy.zeros_like(self.wave)
            for i in xrange(len(ai)):
                eta += self.dust_func(self.wave, ai[i], bi[i], ni[i], li[i])
            
            Rv = 3.08
            Ab = self.tauv * (1 + (1/Rv))

            ATT = numpy.ones([len(self.wave),len(self.ta)])
            for ti in range(0,len(self.ta)):
                ATT[:,ti] *= numpy.power(10,-0.4*(Ab*eta))
                #Offset added to renormalise from B to V band
                #ATT[:,ti] *= numpy.exp(-1*self.tauv*eta)
                 
        

        
        """
        SECTION 1
        First calculate and store those parameters that are functions of the age array 
        'ta' only - these are the same for every model to be made. The parameters are 
        the age array TP, the time interval array DT, the interpolation coefficient 
        'a' and the interpolation indices J. Each are stored in cell arrays of size ks,
        with the data corresponding to the original age array first, and the 
        interpolated data second.
        """
        self.TP = {}
        self.A = {}
        self.J = {}
        self.DT = {}

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
            self.TP[ai] = self.ta[ai]-numpy.concatenate((taux1,taux2),axis=0)
            l = len(taux2)

            #If taux2 has entries, calculate the interpolation parameters a and J.
            #The indicies correspond to those values of 'ta' which are just below
            #the entries in taux2. They are calculated by taking the difference
            #between the two arrays, then finding the last negative entry in the
            #resulting array.

            if l == 0:
                self.J[ai] = numpy.array([])
                self.A[ai] = numpy.array([])
            if l>0:
                [T1,T2] = numpy.meshgrid(self.ta,taux2)
                T = T1-T2
                T[numpy.where(T<=0)] = 0
                T[numpy.where(T!=0)] = 1
                T = numpy.diff(T,1,1)
                (i,self.J[ai]) = T.nonzero()

                self.A[ai] = (numpy.log10(taux2/self.ta[self.J[ai]]) / 
                              numpy.log10(self.ta[self.J[ai]+1]/self.ta[self.J[ai]]))

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
            self.DT[ai] = numpy.copy(dt[order])


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
            j = self.J[ai]   #Interpolation indices
            tp = self.TP[ai] #Integration timescale

            pgas = numpy.zeros_like(tp)
            if ai ==0:
                prgas = numpy.zeros_like(self.ta)
            else:
                i = numpy.where(tp<=self.ta[ai-1])
                ii = numpy.where(tp>self.ta[ai-1])
                pgas[i] = griddata(self.ta,prgas,tp[i])
                pgas[ii] = prgas[ai-1]
            #print prgas[ai]
    
            tbins = numpy.logspace(0,numpy.log10(max(tp)),1000)
            npgas = numpy.zeros_like(tbins)
            
            if self.sfh_law == 'exp':
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
                    self.sr = sr
                    if len(sr) > 1:
                        norma = simps(numpy.exp(-1*tbins/self.tau)/abs(self.tau),tbins)
                        sr /= norma
                    #print sr[0]
                    self.norma = norma

                w = sr*self.DT[ai]/2
                w1 = numpy.array(w[:ai+1])
                W[0,ai] = w1

                strr = numpy.array(numpy.dot(w1,strm[:ai+1]))
                rm = numpy.array(numpy.dot(w1,rmtm[:ai+1]))

                l = len(self.A[ai])
                if l>0:

                    w2 = w[ai+1:ai+l+1]
                    wa = w2*self.A[ai]
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
                    SFR[ai] = numpy.exp(-self.ta[ai]/self.tau)/abs(self.tau)/norma
                #print SFR[ai,ti,mi]
                #self.SN = float(snn)    
                SFR[ai] /= STR[ai]
            else:
                if self.sfh_law == 'pow':
                    sfr = self._sfh_pow
                elif self.sfh_law == 'del':
                    sfr = self._sfh_del
                elif self.sfh_law == 'tru':
                    sfr = self._sfh_tru
                
                sr = sfr(tp,self.tau)
                self.tp=tp
                norma = 1
                self.sr = sr
                if len(sr) > 1:
                    norma = simps(sfr(tbins,self.tau),tbins)
                    sr /= norma
                self.norma = norma
                #print sr[0]

                w = sr*self.DT[ai]/2
                w1 = numpy.array(w[:ai+1])
                W[0,ai] = w1

                strr = numpy.array(numpy.dot(w1,strm[:ai+1]))
                rm = numpy.array(numpy.dot(w1,rmtm[:ai+1]))

                l = len(self.A[ai])
                if l>0:

                    w2 = w[ai+1:ai+l+1]
                    wa = w2*self.A[ai]
                    wb = w2-wa

                    W[1,ai] = wa
                    W[2,ai] = wb
                    strr += (numpy.dot(wb,strm[j]) + numpy.dot(wa,strm[j+1]))
                    rm += (numpy.dot(wb,rmtm[j]) + numpy.dot(wa,rmtm[j+1]))


                if strr > 1: strr= 1

                if self.tau > 0.:
                    ugas = sfr(self.ta,self.tau)[ai]
                elif self.tau < 0.:
                    ugas = sfr(self.ta,self.tau)[ai]/sfr(max(self.ta),self.tau)
                    #ugas = 1.
                #Processed gas = gas formed into stars - mass in stars - remnants
                prgas[ai] = 1 - ugas - strr -rm
                if prgas[ai] < 0.: prgas[ai] = 0

                #print prgas[ai]
                URr[ai] = ugas
                PRr[ai] = prgas[ai]
                RMr[ai] = rm
                Tr[ai] = simps(sfr(numpy.sort(tp)/1.e9,self.tau),numpy.sort(tp))

                STR[ai] = strr
                if self.tau > 0:
                    SFR[ai] = (1 + epsilon*prgas[ai])*sfr(self.ta,self.tau)[ai]/norma
                elif self.tau < 0:
                    SFR[ai] = sfr(self.ta[ai],self.tau)/norma
                #print SFR[ai,ti,mi]
                #self.SN = float(snn)    
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
        j = self.J[ai]


        w1 = W[0,ai]
        wa = W[1,ai]
        wb = W[2,ai]

        for i in range(ai):
            y += (w1[i]*sed1[:,i])
            y_nodust += (w1[i]*sed[:,i])

        for i in range(len(wb)):
            y += (wb[i]*sed1[:,j[i]] + wa[i]*sed1[:,j[i]+1])
            y_nodust += (wb[i]*sed[:,j[i]] + wa[i]*sed[:,j[i]+1])


        Nly = self.calc_lyman(self.wave,numpy.nan_to_num(y_nodust[0]))
        #print Nly
        if Nly > 0.:
            Nlyman = numpy.log10(Nly)
        else:
            Nlyman = 0.

        total = (self.neb_cont*self.inc_cont) + self.neb_hlines + (self.neb_metal[:,neb_z]*self.inc_met)
        total *= 2.997925e18/(self.wave**2) #Convert to Flambda
        total *= (Nly*(1-self.fesc))

        y += total
    
        Nly = self.calc_lyman(self.wave,numpy.nan_to_num(y[0] / STR[ai]))
        #print Nly
        self.fesc_tot = (self.fesc*Nly) / 10**Nlyman
        if Nly > 0.:
            Nlyman_final = numpy.log10(Nly) + 33. + numpy.log10(3.826)
            if self.fesc > 0.:
                Nlyman_final = numpy.log10(10**Nlyman_final * self.fesc)
            elif self.fesc == 0:
                Nlyman_final = 0.
        else:
            Nlyman_final = 0.

        beta = self.calc_beta(self.wave,y[0])
        #print ai,ai1
        #print STR[ai1,ti,mi]
        SED[:] = y/STR[ai] #normalised to 1 solar mass
        norm = simps(numpy.exp(-1*numpy.logspace(0,numpy.log10(self.ta[tgi]),10000)/self.tau),numpy.logspace(0,numpy.log10(self.ta[tgi]),10000))

        STR = STR[tgi]
        SFR = SFR[tgi]
        
        self.SED = SED
        self.SFR = SFR / STR
        self.STR = STR
        self.beta = beta
        self.Nly = Nlyman_final
        self.Ms = 1.
        
    def calc_beta(self,wave, SED):
        """
        wave = wavelength array
        SED = Rest-frame flux

        Returns UV slope index and the error on that fit
        """
        #new_wave = numpy.arange(1200,2650)
        #new_SED = griddata(wave,SED,new_wave)
        #wave = new_wave
        #SED = new_SED

        #window_lower = numpy.array([1268.,1309.,1342.,1407.,1562.,1677.,1760.,1866.,1930.,2400.])
        #window_upper = numpy.array([1284.,1316.,1371.,1515.,1583.,1740.,1833.,1890.,1950.,2580.])

        #window_lower = numpy.array([1407,1562,1677.,1760.,1866.,1930.,2400.])
        #window_upper = numpy.array([1515,1583,1740.,1833.,1890.,1950.,2580.])

        window_lower = numpy.array([1600,2400])
        window_upper = numpy.array([1950,2600])
        
        ww = numpy.zeros_like(wave,dtype=bool)
        for w in numpy.arange(len(window_lower)):
            ww[(wave >= window_lower[w])*(wave < window_upper[w])] = True

        #window_mean = (window_lower+window_upper)/2 #midpoint for power-law fitting

        fluxes = numpy.zeros_like(window_lower)

        #for w, window_lower in enumerate(window_lower):
        #    wf = numpy.where((wave > window_lower) & (wave <= window_upper[w]))
        #    fluxes[w] = numpy.mean(SED[wf])

        #fluxes *= 2.997925e18/(window_mean**2)
        fluxerr = numpy.sqrt(fluxes)

        logx = numpy.log10(wave[ww])
        logy = numpy.log10(SED[ww])
        logyerr = 1.#fluxerr/fluxes
    
        fitfunc = lambda p, x: p[0] + (x*p[1])
        errfunc = lambda p, x, y, err: (y - fitfunc(p,x))/err
    
        pinit = [numpy.max(SED[ww]), -2.0]
        out = leastsq(errfunc, pinit, args=(logx,logy,logyerr))
        #out = leastsq(errfunc, pinit, args=(log,fluxes,fluxerr))

        pfinal = out[0]
        covar = out[1]
        #print pfinal
        index = pfinal[1]
        #indexerr = numpy.sqrt(covar[0])
    
        return (index)#, indexerr)
    
    def calc_lyman(self,x,s):
        wly = 912.
        const_c = 2.997925e10
        const_h = 6.6262e-27
        const = 1e-8/const_h/const_c

        n = int(sum([x < wly][0]))
        f = numpy.zeros(n+1)
        w = numpy.zeros(n+1)

        for i in range(n+1):
            if x[i] < wly:
                w[i] = x[i]
                f[i] = w[i]*s[i]
            elif x[i] == wly:
                w[i] = x[i]
                f[i] = w[i]*s[i]
            elif x[i] > wly:
                w[i] = wly
                f[i] = s[i-1] + ((w[i]-x[i-1])*(s[i]-s[i-1])/(x[i]-x[i-1]))
                f[i] = w[i]*f[i]

        nlyman = const*numpy.trapz(f,w)
        #print numpy.log10(N_lyman)
        return nlyman
        
    def __str__(self):
        params = ['Age', 'SFH Tau', 'Dust Tau', 'SFR', 'Stellar Mass','Beta']
        values = [self.tg/1e9, self.tau/1e9, self.tauv, self.SFR, self.Ms, self.beta]
        units = ['Gyr', 'Gyr', 'Av', 'Ms/yr','Msol','']
        output = ['{:>14s}'.format(params[i]) +': '+ '{:<.3g}'.format(values[i]) + ' ' +units[i] for i in range(len(params))]
        
        return '\n'.join(output)
    
    def __add__(self,other):
        if isinstance(other,CSP):
            new = copy.deepcopy(self)
            new.SED += other.SED
            new.SFR += other.SFR
            new.Ms += other.Ms
            if new.Nly == 0:
                new.Nly = other.Nly
            if other.Nly == 0:
                new.Nly = new.Nly
            if (new.Nly > 0.) and (other.Nly > 0.):
                new.Nly = numpy.log10(10**self.Nly + 10**other.Nly )
            new.beta = new.calc_beta(new.wave,new.SED)
        return new

    def __iadd__(self,other):
        if isinstance(other,CSP):
            new = copy.deepcopy(self)
            new.SED += other.SED
            new.SFR += other.SFR
            new.Ms += other.Ms
            if new.Nly == 0:
                new.Nly = other.Nly
            if other.Nly == 0:
                new.Nly = new.Nly
            if (new.Nly > 0.) and (other.Nly > 0.):
                new.Nly = numpy.log10(10**self.Nly + 10**other.Nly )
            new.beta = new.calc_beta(new.wave,new.SED)
        return new


    def __mul__(self,other):
        new = copy.deepcopy(self)
        new.SED *= other
        new.SFR *= other
        new.Ms *= other
        if other == 0.:
            new.Nly = 0.
        else:
            new.Nly += numpy.log10(other)
            new.Nly = numpy.maximum(new.Nly,0)
        return new

    def __imul__(self,other):
        new = copy.deepcopy(self)
        new.SED *= other
        new.SFR *= other
        new.Ms *= other
        if other == 0.:
            new.Nly = 0.
        else:
            new.Nly += numpy.log10(other)
            new.Nly = numpy.maximum(new.Nly,0)
        return new

    def __div__(self,other):
        new = copy.deepcopy(self)
        new.SED /= other
        new.SFR /= other
        new.Ms /= other
        if other == 0.:
            new.Nly = 0.
        else:
            new.Nly -= numpy.log10(other)
            new.Nly = numpy.maximum(new.Nly,0)
        return new

    def __idiv__(self,other):
        new = copy.deepcopy(self)
        new.SED /= other
        new.SFR /= other
        new.Ms /= other
        if other == 0.:
            new.Nly = 0.
        else:
            new.Nly -= numpy.log10(other)
            new.Nly = numpy.maximum(new.Nly,0)
        return new
        
    def __rmul__(self,other):
        new = copy.deepcopy(self)
        new.SED *= other
        new.SFR *= other
        new.Ms *= other
        if other == 0.:
            new.Nly = 0.
        else:
            new.Nly += numpy.log10(other)
            new.Nly = numpy.maximum(new.Nly,0)
        return new
    
    def addEmissionLine(self,wavelength,EqW):
        wbin = numpy.argmin(numpy.abs(self.wave-wavelength))
        #print wbin
        binwidth = numpy.mean(numpy.diff(self.wave)[wbin-1:wbin+1])
        #print binwidth
        continuum = numpy.mean(self.SED[wbin:wbin+1])
        #print continuum
        lineluminosity = continuum * EqW
        #print lineluminosity, lineluminosity/binwidth
        self.Lalpha = lineluminosity
        self.SED[wbin] += lineluminosity/binwidth


        
class Filter(object):
    def __init__(self):
        self.wave = []
        self.freq = []
        self.response = []

        self.lambda_c = []
        self.nu_c = []
        
class FileFilter(Filter):
    def __init__(self,filepath,minbins=200):
        self.path = filepath
        
#        try:
        data = numpy.loadtxt(self.path)
        wf = data[:,0]
        tp = data[:,1]
        if len(data[:,0]) < minbins: #Re-sample large filters for performance
            wfx = numpy.linspace(wf[0],wf[-1],minbins)
            tpx = griddata(wf,tp,wfx)

            wf = wfx
            tp = tpx
            
        self.wave = wf * U.angstrom
        self.response = tp
        
        self.freq = (C.c/self.wave).to(U.Hz)
        
        nmax = numpy.argmax(self.response)
        halfmax_low = self.wave[:nmax][numpy.argmin(numpy.abs(self.response[nmax] - 2*self.response[:nmax]))]
        halfmax_hi = self.wave[nmax:][numpy.argmin(numpy.abs(self.response[nmax] - 2*self.response[nmax:]))]
        print self.wave[nmax],halfmax_low, halfmax_hi
        self.fwhm = halfmax_hi-halfmax_low
        
        self.lambda_c = (simps(self.wave*self.response,self.wave) / 
                         simps(self.response,self.wave))

        self.nu_c = (simps(self.freq*self.response,self.freq) / 
                     simps(self.response,self.freq))
        
        
#        except:
#            print 'Ohhhhh dear.'
        
        
class TophatFilter(Filter):
    def __init__(self, centre, width, steps = 200):
        self.centre = centre * U.angstrom
        self.width = width * U.angstrom
        self.steps = steps

        upper, lower = self.centre+self.width, self.centre-self.width
        resp_upper, resp_lower = self.centre+(self.width*0.5), self.centre-(self.width*0.5)

        self.wave = numpy.linspace(lower,upper,steps)
        self.response = numpy.zeros_like(self.wave.value)

        tophat = (self.wave >= resp_lower)*(self.wave < resp_upper)
        self.response[tophat] = 1

        self.freq = (C.c/self.wave).to(U.Hz)
        
        self.lambda_c = (simps(self.wave*self.response,self.wave) / 
                         simps(self.response,self.wave))

        self.nu_c = (simps(self.freq*self.response,self.freq) / 
                     simps(self.response,self.freq))

class LoadEAZYFilters(object):
    def __init__(self,path):
        self.path = path
        
        self.filters = []
        self.filternames = []
        self.central_wlengths = []
        
        with open(self.path) as file:
            for f in xrange(1000):
                x = file.readline().split()
                if len(x) < 1:
                    break
                nwave, name, lambda_c = x[0],x[1],x[4] 
                nwave = int(nwave)
            
                wavelength = []
                response = []
                #print nwave
                for i in xrange(nwave):
                
                    N, w, r = numpy.array(file.readline().split()).astype('float')
                    wavelength.append(w)
                    response.append(r)
                
                wavelength *= U.angstrom
                freq = (C.c/wavelength).to(U.Hz)
        
                lambda_c = (simps(wavelength*response,wavelength) / 
                            simps(response,wavelength))

                nu_c = (simps(freq*response,freq) / 
                        simps(response,freq))
                
                new_filter = Filter()
                new_filter.wave = numpy.array(wavelength) * U.angstrom
                new_filter.response = numpy.array(response)
                new_filter.freq = numpy.array(freq) * U.Hz
                new_filter.lambda_c = lambda_c
                new_filter.nu_c = nu_c
                nmax = numpy.argmax(new_filter.response)
                halfmax_low = new_filter.wave[:nmax][numpy.argmin(numpy.abs(new_filter.response[nmax] - 2*new_filter.response[:nmax]))]
                halfmax_hi = new_filter.wave[nmax:][numpy.argmin(numpy.abs(new_filter.response[nmax] - 2*new_filter.response[nmax:]))]
                #print new_filter.wave[nmax],halfmax_low, halfmax_hi
                new_filter.fwhm = halfmax_hi-halfmax_low                
                
                self.filters.append(new_filter)
                self.filternames.append(name)
                self.central_wlengths.append(lambda_c)
                
                
        self.central_wlengths *= U.angstrom
  
  
class FilterSet:
    def __init__(self,path=None):
        self.directory = path
        self.filters = []
        if type(self.directory) == str:
            try:
                self.files = glob(self.directory)
                self.files.sort()
                for file in self.files:
                    self.filters.append(FileFilter(file))
            except:
                a = ''

    def addFileFilter(self,path):
        self.filters.append(FileFilter(path))
    
    def addTophatFilter(self,centre, width, steps = 200):
        self.filters.append(TophatFilter(centre, width, steps))

    def addEAZYFilter(self,EAZYset,N):
        if type(N) == int:
            self.filters.append(EAZYset.filters[N])
        elif type(N) == list:
            for x in N:
                self.filters.append(EAZYset.filters[x])

class Observe:
    def __init__(self,SED,Filters,redshift,force_age = True,madau=True):
        self.SED = SED
        self.F = Filters
        self.redshifts = numpy.array(redshift,ndmin=1)
        self.wave = self.SED.wave

        self.fluxes = []
        self.AB = []
        self.wl = []
        self.fwhm = []
        for z in self.redshifts:
            self.lyman_abs = numpy.ones(len(self.wave))
            if madau:
                ly_cont_w = numpy.array([(self.wave<=912.)][0])
                ly_b_w = numpy.array([(self.wave > 912.) & (self.wave <= 1026.)][0])
                ly_a_w = numpy.array([(self.wave > 1026.) & (self.wave <= 1216.)][0])
        
                dec_a = (1/(120*(1+z)))*quad(self.dec_a_func,
                        1050*(1+z),1170*(1+z))[0]
                dec_b= (1/(95*(1+z)))*quad(self.dec_b_func,
                        920*(1+z),1015*(1+z))[0]
        
                self.lyman_abs[ly_cont_w] *= 0.1
                self.lyman_abs[ly_b_w] = dec_b
                self.lyman_abs[ly_a_w] = dec_a
        
            if z > 0:
                self.dm = cosmo.distmod(z).value
            else:
                self.dm = 0.
        
            if (self.SED.tg/1e9 > cosmo.age(z).value) and force_age:
                print 'SSP age older than universe...stopping.'
            else:            
                tfluxes = []
                tAB = []
                twl = []
                tfwhm = []
            
                for filt in self.F.filters:
                    #print filt.wave[0]
                    flux, mag = self.calcFlux(filt,z)

                    tfluxes.append(flux)
                    tAB.append(mag)
                    twl.append(filt.lambda_c)
                    tfwhm.append(filt.fwhm.value)
                self.fluxes.append(tfluxes)
                self.AB.append(tAB)
                self.wl = twl
                self.fwhm = tfwhm
        
        self.wl *= U.angstrom
        self.fwhm *= U.angstrom
        self.fluxes = (numpy.squeeze(self.AB) * 1e-6*U.Jy)
        self.AB = (numpy.squeeze(self.AB) * U.mag)
    
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
        
    def calcFlux(self,filt,z):
        wf = filt.wave.value
        tp = filt.response
        z1 = z+1

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
        
        return Flux , Mag


     
#must define cosmo before calling an Observe
"""
a = SSP()
a.build(0.8,1000,0,4,fesc=0.) #(Age, Tau, Dust, Metallicity Index)
#print a

b = SSP()
b.build(1.5,0.5,1,2,fesc=1.)
#print b

a = a*1e9 #Multiply by a factor, equivalent to *=
b = b*5e9 
c = a+b
#a+b #Add b to a (a += b)
print a

filt_dir = '../aMasSED-code/GOODS-S_18_FilterCurves/Filter*.txt'

Filts = FilterSet(filt_dir)
Filts.addTophatFilter(1500,100)

AA = Observe(a,Filts,2) # Observe a (built from two initial SEDs) 
                        # through the filters in Filts at redshift of z = 2 
#BB = Observe(b,Filts,2) 
print AA.AB
"""