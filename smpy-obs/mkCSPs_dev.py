import numpy as np
import array
import copy
import re
import sys
from glob import glob
from scipy.interpolate import griddata
from scipy.integrate import simps
from scipy.optimize import leastsq
from scipy.spatial import Delaunay

from astropy import units as u
from astropy import constants as c
from astropy import cosmology as cos

cosmo = cos.FlatLambdaCDM(H0=70, Om0=0.3)

f = open("error.log", "w")
original_stderr = sys.stderr
sys.stderr = f


class Ised(object):
    def __init__(self, path):
        self.file = path

        self.read_ised(self.file)

    def read_ised(self, filename):
        """
        This function reads data from Bruzual & Charlot binary format
        SSP files and returns the necessary data in an array The input files
        should be '.ised' files, either 2003 or 2007.
    
        'ks' in the binary files is slightly different between 03/07 files
        so the read length and index should be set appropriately, therefore 
        the function tries '03 format first and retries with the '07 format
        if the returned number of ages isn't as expected (e.g. 221 ages)

        :type filename: string
        """

        with open(filename, 'rb') as file:
            check = array.array('i')
            check.fromfile(file, 2)

        if check[1] == 221:
            ksl, ksi = 2, 1
            F_l, F_i = 3, 2
        else:
            ksl, ksi = 3, 2
            F_l, F_i = 5, 4

        with open(filename, 'rb') as file:
            ks = array.array('i')
            ks.fromfile(file, ksl)

            ta = array.array('f')
            ta.fromfile(file, ks[ksi])
            self.ta = np.array(ta)

            tmp = array.array('i')
            tmp.fromfile(file, 3)
            self.ml, self.mul, iseg = tmp

            if iseg > 0:
                tmp = array.array('f')
                tmp.fromfile(file, iseg * 6)

            tmp = array.array('f')
            tmp.fromfile(file, 5)
            self.totm, self.totn, self.avs, self.jo, self.tauo = tmp

            self.ids = array.array('c')
            self.ids.fromfile(file, 80)

            tmp = array.array('f')
            tmp.fromfile(file, 4)
            self.tcut = tmp[0]
            self.ttt = tmp[1:]

            ids = array.array('c')
            ids.fromfile(file, 80)

            self.ids = array.array('c')
            self.ids.fromfile(file, 80)

            self.igw = array.array('i')
            self.igw.fromfile(file, 1)

            tmp = array.array('i')
            tmp.fromfile(file, F_l)

            self.iw = array.array('i')
            self.iw.fromfile(file, 1)

            wave = array.array('f')
            wave.fromfile(file, self.iw[0])
            self.wave = np.array(wave)

            # SED Section
            self.F = array.array('i')
            self.F.fromfile(file, F_l)
            self.iw = self.F[F_i]  # Number of wavelength elements

            self.sed = np.zeros((self.iw, ks[ksi]), dtype=np.float32)
            G = array.array('f')
            G.fromfile(file, self.iw)
            self.sed[:, 0] = G
            ik = array.array('i')
            ik.fromfile(file, 1)

            self.h = np.empty((ik[0], ks[ksi]), 'f')
            H = array.array('f')
            H.fromfile(file, ik[0])
            self.h[:, 0] = H

            for i in range(1, ks[ksi]):  # Fill rest of array with SEDs
                F = array.array('i')
                F.fromfile(file, F_l)
                iw = F[F_i]

                G = array.array('f')
                G.fromfile(file, iw)
                self.sed[:, i] = G
                ik = array.array('i')
                ik.fromfile(file, 1)

                H = array.array('f')
                H.fromfile(file, ik[0])
                self.h[:, i] = H

            tmp = array.array('i')
            tmp.fromfile(file, F_l)

            self.bflx = array.array('f')
            self.bflx.fromfile(file, tmp[F_i])

            tmp = array.array('i')
            tmp.fromfile(file, F_l)

            strm = array.array('f')
            strm.fromfile(file, tmp[F_i])
            self.strm = np.array(strm)

            tmp = array.array('i')
            tmp.fromfile(file, F_l)

            self.evf = array.array('f')
            self.evf.fromfile(file, tmp[F_i])

            tmp = array.array('i')
            tmp.fromfile(file, F_l)

            self.evf = array.array('f')
            self.evf.fromfile(file, tmp[F_i])

            tmp = array.array('i')
            tmp.fromfile(file, F_l)

            self.snr = array.array('f')
            self.snr.fromfile(file, tmp[F_i])

            tmp = array.array('i')
            tmp.fromfile(file, F_l)

            self.pnr = array.array('f')
            self.pnr.fromfile(file, tmp[F_i])

            tmp = array.array('i')
            tmp.fromfile(file, F_l)

            self.sn = array.array('f')
            self.sn.fromfile(file, tmp[F_i])

            tmp = array.array('i')
            tmp.fromfile(file, F_l)

            self.bh = array.array('f')
            self.bh.fromfile(file, tmp[F_i])

            tmp = array.array('i')
            tmp.fromfile(file, F_l)

            self.wd = array.array('f')
            self.wd.fromfile(file, tmp[F_i])

            tmp = array.array('i')
            tmp.fromfile(file, F_l)

            rmtm = array.array('f')
            rmtm.fromfile(file, tmp[F_i])
            self.rmtm = np.array(rmtm)


class SSP(object):
    """ Base class for SSP objects
    """

    def __init__(self, path=None):
        self.SSPpath = path


class BC(SSP):
    """ Bruzual & Charlot 2003/07 SSP Models
    
    """

    def __init__(self, path='../ssp/bc03/salpeter/lr/'):
        """


        :type path: string
        :param path:
        :return:
        """

        self.SSPpath = path
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
        # Set up
        metallicities = np.zeros(len(self.files))

        for i, file in enumerate(self.files):
            ised_binary = Ised(file)

            self.ta_arr.append(ised_binary.ta)
            self.metal_arr.append(ised_binary.ids)
            self.iw_arr.append(ised_binary.iw)
            self.wave_arr.append(ised_binary.wave)
            self.sed_arr.append(ised_binary.sed)
            self.strm_arr.append(ised_binary.strm)
            self.rmtm_arr.append(ised_binary.rmtm)
            self.iseds.append(ised_binary)
            metal = str(ised_binary.ids)[12:-3].strip()
            metallicities[i] = float(re.split("Z=?", metal)[1])

        self.metallicities = metallicities / 0.02  # Normalise to solar metallicity
        self.ages = np.array(self.ta_arr[0]) * u.yr
        self.ages[0] = (10 * u.yr) # Make non-zero but insignificant

        self.wave_arr = np.array(self.wave_arr) * u.AA

        # Reshape to match later analysis and strip t=0 columns
        self.sed_arr = np.array(self.sed_arr).swapaxes(1, 2)[:, :]
        self.sed_arr *= u.Lsun / u.AA
        self.strm_arr = np.array(self.strm_arr)[:, :]
        self.rmtm_arr = np.array(self.rmtm_arr)[:, :]
        self.iseds = np.array(self.iseds)


class CSP:
    """ Class for building composite stellar populations from input SSPs
            
    """

    def __init__(self, ssp,
                 age=None, sfh=None, dust=None, metal_ind=None, f_esc=None,
                 sfh_law='exp', dust_model='calzetti', neb_cont=True, neb_met=True):
        """
        :type ssp: SSP object
        :type age:
        """

        self.beta = None
        self.SSP = ssp
        self.metallicities = self.SSP.metallicities  # Normalise to solar metallicity
        self.ages = self.SSP.ages

        self.wave = self.SSP.wave_arr[0]
        self.iw = self.SSP.iw_arr[0]
        self.sed_arr = self.SSP.sed_arr

        if hasattr(self.SSP, 'strm_arr') and hasattr(self.SSP, 'rmtm_arr'):
            self.strm_arr = np.array(self.SSP.strm_arr)
            self.rmtm_arr = np.array(self.SSP.rmtm_arr)

        self.iseds = np.array(self.SSP.iseds)

        # Find closest match for each tg value in ta - set tg to these values

        nebular = np.loadtxt('nebular_emission.dat', skiprows=1)
        self.neb_cont = nebular[:, 1]
        self.neb_hlines = nebular[:, 2]
        self.neb_metal = nebular[:, 3:]
        self.neb_wave = nebular[:, 0]

        if None not in (age, sfh, dust, metal_ind):
            if f_esc == None:
                self.build(age, sfh, dust, metal_ind, sfh_law=sfh_law, dust_model=dust_model,
                           neb_cont=neb_cont, neb_met=neb_met)

            else:
                self.build(age, sfh, dust, metal_ind, f_esc, sfh_law, dust_model, neb_cont, neb_met)

    @staticmethod
    def _sfh_exp(t, tau):
        sfh = np.exp(-1 * t / tau) / abs(tau)
        return sfh

    @staticmethod
    def _sfh_pow(t, alpha):
        sfh = np.power(t, alpha)
        return sfh

    @staticmethod
    def _sfh_del(t, tau):
        sfh = t / (tau ** 2) * np.exp(-t / tau)
        return sfh

    @staticmethod
    def _sfh_tru(t, tstop):
        sfh = np.ones_like(t)
        sfh[t > tstop * np.max(t)] = 0.

        sfh /= np.trapz(sfh, t)
        return sfh

    @staticmethod
    def dust_func(lam, ai, bi, ni, li):
        """
        Functional form for SMC, LMC and MW extinction curves of
        Pei et al. 1992

        """
        lam = np.array(lam) / 1e4
        ki = np.power((lam / li), ni) + np.power((li / lam), ni) + bi
        eta_i = ai / ki
        return eta_i

    def build(self, age, sfh, dust, metal, fesc=1.,
              sfh_law='exp', dust_model='calzetti',
              neb_cont=True, neb_met=True, timesteps = 500):
        """ Docs

        """

        self.tg = age.to(u.yr)
        if sfh_law == 'exp':
            self.tau = sfh.to(u.yr)
        elif sfh_law == 'del':
            self.tau = sfh.to(u.yr)
        else:
            self.tau = sfh

        self.tauv = dust
        self.mi = metal
        self.fesc = fesc
        self.sfh_law = sfh_law
        self.inc_cont = neb_cont
        self.inc_met = neb_met
        self.dust_model = dust_model

        mu = 0.3

        self.ta = self.ages

        # Set up nebular emission arrays -- WILL CHANGE
        if len(self.neb_wave) != len(self.wave):
            self.neb_cont = griddata(self.neb_wave, self.neb_cont, self.wave)
            self.neb_hlines = griddata(self.neb_wave, self.neb_hlines, self.wave)
            neb_metaln = np.zeros((len(self.wave), 3))
            for i in range(3):
                neb_metaln[:, i] = griddata(self.neb_wave, self.neb_metal[:, i], self.wave)
            self.neb_metal = neb_metaln
            self.neb_wave = self.wave

        self.neb_cont[self.wave <= 912. * u.AA] = 0.
        self.neb_hlines[self.wave <= 912. * u.AA] = 0.
        self.neb_metal[self.wave <= 912. * u.AA, :] = 0.

        self.Nly_arr = self.calc_lyman(self.wave, self.sed_arr)
        self.neb_sed_arr = np.ones_like(self.sed_arr.value)

        nebrange1 = (self.metallicities <= 0.02)
        nebrange2 = (self.metallicities > 0.02) * (self.metallicities <= 0.2)
        nebrange3 = (self.metallicities > 0.2)
        self.neb_sed_arr[nebrange1, :, :] *= ((self.neb_cont * self.inc_cont) +
                                              self.neb_hlines +
                                              (self.neb_metal[:, 0] * self.inc_met))

        self.neb_sed_arr[nebrange2, :, :] *= ((self.neb_cont * self.inc_cont) +
                                              self.neb_hlines +
                                              (self.neb_metal[:, 1] * self.inc_met))

        self.neb_sed_arr[nebrange3, :, :] *= ((self.neb_cont * self.inc_cont) +
                                              self.neb_hlines +
                                              (self.neb_metal[:, 2] * self.inc_met))

        self.neb_sed_arr *= c.c.to(u.AA / u.s) / (self.wave ** 2)  # Convert to Flambda
        self.neb_sed_arr *= (self.Nly_arr[:, :, None] * (1 - self.fesc))

        # SSP Interpolation Section
        self.ta_sfh = np.logspace(np.log10(self.ages / u.yr).min(), 
                                  np.log10(self.tg / u.yr), timesteps) * u.yr
        self.me_sfh = np.ones(len(self.ta_sfh)) * self.mi

        # Calculate Barycentric coordinates for all ages/metallicities in SFH.
        points = np.array(zip(np.log10(self.me_sfh), np.log10(self.ta_sfh / u.yr)))

        ti, mi = np.meshgrid(np.log10(self.ages / u.yr), np.log10(self.metallicities))
        self.grid = zip(mi.flatten(), ti.flatten())

        tri_grid = Delaunay(self.grid)
        ss = tri_grid.find_simplex(points)

        X = tri_grid.transform[ss, :2]
        # Offset of each target from the origin of its containing tetrahedron.
        Y = points - tri_grid.transform[ss, 2]

        b = np.einsum('ijk,ik->ij', X, Y)
        self.bc = np.c_[b, 1 - b.sum(axis=1)]
        self.simplices = tri_grid.simplices[ss]

        # Interpolate SED, stellar mass fraction and remnant fractions
        # for SFH age grid using calculated Barycentric coordinates (bc).

        self.sed_sfh = ((self.sed_arr + self.neb_sed_arr).reshape(len(self.metallicities) *
                                                                  len(self.ages), self.iw)[self.simplices] *
                        self.bc[:, :, None]).sum(1)

        self.neb_sed_sfh = np.array(self.neb_sed_arr.reshape(len(self.metallicities) *
                                                             len(self.ages), self.iw)[self.simplices]
                                    * self.bc[:, :, None]).sum(1)

        self.strm_sfh = np.array(self.strm_arr.reshape(len(self.metallicities) *
                                                       len(self.ages))[self.simplices]
                                 * self.bc).sum(1)

        self.rmtm_sfh = np.array(self.rmtm_arr.reshape(len(self.metallicities) *
                                                       len(self.ages))[self.simplices]
                                 * self.bc).sum(1)
        self.Nly_sfh = np.array(self.Nly_arr.reshape(len(self.metallicities) *
                                                     len(self.ages))[self.simplices]
                                * self.bc).sum(1)

        # Dust Section
        if self.dust_model == "charlot":
            self.Att = np.zeros((len(self.ta_sfh), len(self.wave)))
            tv = ((self.tauv / 1.0857) * np.ones(len(self.ta_sfh)))
            tv[self.ta_sfh > 1e7 * u.yr] = mu * self.tauv
            lam = np.array((5500 * u.AA / self.wave) ** 0.7)
            self.Att[:, :] = (np.exp(-1 * np.outer(tv, lam)))

        elif self.dust_model == "calzetti":
            self.Att = np.ones([len(self.ta_sfh), len(self.wave)])
            k = np.zeros_like(self.wave.value)

            w0 = [self.wave <= 1200 * u.AA]
            w1 = [self.wave < 6300 * u.AA]
            w2 = [self.wave >= 6300 * u.AA]
            w_u = self.wave.to(u.um).value

            x1 = np.argmin(np.abs(self.wave - 1200 * u.AA))
            x2 = np.argmin(np.abs(self.wave - 1250 * u.AA))

            k[w2] = 2.659 * (-1.857 + 1.040 / w_u[w2])
            k[w1] = 2.659 * (-2.156 + (1.509 / w_u[w1]) - (0.198 / w_u[w1] ** 2) + (0.011 / w_u[w1] ** 3))
            k[w0] = k[x1] + ((self.wave[w0] - 1200. * u.AA) * (k[x1] - k[x2]) / (self.wave[x1] - self.wave[x2]))

            k += 4.05
            k[k < 0.] = 0.

            tv = self.tauv * k / 4.05
            self.Att *= np.power(10, -0.4 * tv)


        elif self.dust_model == "calzetti2":
            self.Att = np.ones([len(self.ta_sfh), len(self.wave)])
            k = np.zeros_like(self.wave.value)

            w0 = [self.wave <= 1000 * u.AA]
            w1 = [(self.wave > 1000 * u.AA) * (self.wave < 6300 * u.AA)]
            w2 = [self.wave >= 6300 * u.AA]
            w_u = self.wave.to(u.um).value

            k[w2] = 2.659 * (-1.857 + 1.040 / w_u[w2])
            k[w1] = 2.659 * (-2.156 + (1.509 / w_u[w1]) - (0.198 / w_u[w1] ** 2) + (0.011 / w_u[w1] ** 3))

            p1 = self.dust_func(self.wave, 27, 4, 5.5, 0.08) + self.dust_func(self.wave, 185, 90, 2, 0.042)

            k[w0] = p1[w0] / (p1[w1][0] / k[w1][0])
            k += 4.05
            k[k < 0.] = 0.
            tv = self.tauv * k / 4.05
            self.Att *= np.power(10, -0.4 * tv)

        elif self.dust_model == "smc":
            ai = [185., 27., 0.005, 0.01, 0.012, 0.03]
            bi = [90., 5.5, -1.95, -1.95, -1.8, 0.]
            ni = [2., 4., 2., 2., 2., 2.]
            li = [0.042, 0.08, 0.22, 9.7, 18., 25.]

            eta = np.zeros_like(self.wave)
            for i in xrange(len(ai)):
                eta += self.dust_func(self.wave, ai[i], bi[i], ni[i], li[i])

            Rv = 2.93
            Ab = self.tauv * (1 + (1 / Rv))
            self.Att = np.ones([len(self.ta_sfh), len(self.wave)])
            self.Att *= np.power(10, -0.4 * Ab * eta)

        elif self.dust_model == "lmc":
            ai = [175., 19., 0.023, 0.005, 0.006, 0.02]
            bi = [90., 4.0, -1.95, -1.95, -1.8, 0.]
            ni = [2., 4.5, 2., 2., 2., 2.]
            li = [0.046, 0.08, 0.22, 9.7, 18., 25.]

            eta = np.zeros_like(self.wave)
            for i in xrange(len(ai)):
                eta += self.dust_func(self.wave, ai[i], bi[i], ni[i], li[i])

            Rv = 3.16
            Ab = self.tauv * (1 + (1 / Rv))

            self.Att = np.ones([len(self.ta_sfh), len(self.wave)])
            self.Att *= np.power(10, -0.4 * Ab * eta)
            # Offset added to renormalise from B to V band.

        elif self.dust_model == "mw":
            ai = [165., 14., 0.045, 0.002, 0.002, 0.012]
            bi = [90., 4., -1.95, -1.95, -1.8, 0.]
            ni = [2., 6.5, 2., 2., 2., 2.]
            li = [0.047, 0.08, 0.22, 9.7, 18., 25.]

            eta = np.zeros_like(self.wave)
            for i in xrange(len(ai)):
                eta += self.dust_func(self.wave, ai[i], bi[i], ni[i], li[i])

            Rv = 3.08
            Ab = self.tauv * (1 + (1 / Rv))

            self.Att = np.ones([len(self.ta_sfh), len(self.wave)])
            self.Att *= np.power(10, -0.4 * Ab * eta)
            # Offset added to renormalise from B to V band.

        # Star-formation history
        if self.sfh_law == 'exp':
            self.sfr_func = self._sfh_exp
        elif self.sfh_law == 'pow':
            self.sfr_func = self._sfh_pow
        elif self.sfh_law == 'del':
            self.sfr_func = self._sfh_del
        elif self.sfh_law == 'tru':
            self.sfr_func = self._sfh_tru

        sfr_hist = self.sfr_func(self.ta_sfh, self.tau)
        # Enforce integrated SFR = 1 Msol.
        norm = np.trapz(sfr_hist, self.ta_sfh)

        sfr_hist /= norm

        sfh_weights = self.sfr_func(self.tg - self.ta_sfh, self.tau) / norm

        self.SED = simps(sfh_weights[:, None] * self.sed_sfh * self.Att,
                         self.ta_sfh, axis=0) * self.sed_sfh.unit
        self.STR = simps(sfh_weights * self.strm_sfh, self.ta_sfh)

        self.SED = self.SED / self.STR
        self.SFR = sfr_hist[-1] / self.STR

        self.Nly = self.calc_lyman_f(self.wave, self.SED)

        if self.Nly > 0:
            self.Nly = np.log10(self.Nly.cgs * self.fesc * u.s)
        else:
            self.Nly = 0

        # self.beta = self.calc_beta(self.wave,self.SED)
        self.Ms = 1 * u.Msun

    def calc_beta(self, wave, sed):
        """
        wave = wavelength array
        SED = Rest-frame flux

        Returns UV slope index and the error on that fit
        """
        # new_wave = np.arange(1200,2650)
        #new_SED = griddata(wave,SED,new_wave)
        #wave = new_wave
        #SED = new_SED

        #window_lower = np.array([1268.,1309.,1342.,1407.,1562.,1677.,1760.,1866.,1930.,2400.])
        #window_upper = np.array([1284.,1316.,1371.,1515.,1583.,1740.,1833.,1890.,1950.,2580.])

        window_lower = np.array([1600, 2400])
        window_upper = np.array([1950, 2600])

        ww = np.zeros_like(wave, dtype=bool)
        for w in np.arange(len(window_lower)):
            ww[(wave >= window_lower[w]) * (wave < window_upper[w])] = True

        #window_mean = (window_lower+window_upper)/2 #midpoint for power-law fitting

        fluxes = np.zeros_like(window_lower)

        #for w, window_lower in enumerate(window_lower):
        #    wf = np.where((wave > window_lower) & (wave <= window_upper[w]))
        #    fluxes[w] = np.mean(SED[wf])

        #fluxes *= 2.997925e18/(window_mean**2)
        fluxerr = np.sqrt(fluxes)

        logx = np.log10(wave[ww])
        logy = np.log10(sed[ww])
        logyerr = 1.  #fluxerr/fluxes

        fitfunc = lambda p, x: p[0] + (x * p[1])
        errfunc = lambda p, x, y, err: (y - fitfunc(p, x)) / err

        pinit = [np.max(sed[ww]), -2.0]
        out = leastsq(errfunc, pinit, args=(logx, logy, logyerr))
        #out = leastsq(errfunc, pinit, args=(log,fluxes,fluxerr))

        pfinal = out[0]
        covar = out[1]
        #print pfinal
        index = pfinal[1]
        #indexerr = np.sqrt(covar[0])

        return (index)  #, indexerr)

    @staticmethod
    def calc_lyman(wave, seds):
        wly = 912. * u.AA
        const = (1e-8 / (u.AA / u.cm)) / c.h.cgs / c.c.cgs

        n = int(sum([wave < wly][0]))
        f = np.zeros((seds.shape[0], seds.shape[1], n + 1)) * seds.unit * u.AA
        w = np.zeros(n + 1) * u.AA

        for i in range(n + 1):
            if wave[i] < wly:
                w[i] = wave[i]
                f[:, :, i] = w[i] * seds[:, :, i]
            elif wave[i] == wly:
                w[i] = wave[i]
                f[:, :, i] = w[i] * seds[:, :, i]
            elif wave[i] > wly:
                w[i] = wly
                f[:, :, i] = w[i] * (seds[:, :, i - 1] + (
                    (w[i] - wave[i - 1]) * (seds[:, :, i] - seds[:, :, i - 1]) / (wave[i] - wave[i - 1])))
                # f[:,:,i] = f[:,:,i]

        nlyman = const * np.trapz(f, w, axis=2)
        # print np.log10(N_lyman)
        return nlyman

    @staticmethod
    def calc_lyman_f(x, s):
        """

        :param s:
        :type x: np.array
        """
        wly = 912. * u.AA
        const = (1e-8 / (u.AA / u.cm)) / c.h.cgs / c.c.cgs

        n = int(sum([x < wly][0]))
        f = np.zeros(n + 1) * s.unit * u.AA
        w = np.zeros(n + 1) * u.AA

        for i in range(n + 1):
            if x[i] < wly:
                w[i] = x[i]
                f[i] = w[i] * s[i]
            elif x[i] == wly:
                w[i] = x[i]
                f[i] = w[i] * s[i]
            elif x[i] > wly:
                w[i] = wly
                f[i] = w[i] * (s[i - 1] + ((w[i] - x[i - 1]) * (s[i] - s[i - 1]) / (x[i] - x[i - 1])))
                # f[i] = w[i]*f[i]

        nlyman = const * np.trapz(f, w)
        return nlyman

    def __add__(self, other):
        new = None
        if isinstance(other, CSP):
            new = copy.deepcopy(self)
            new.SED += other.SED
            new.SFR += other.SFR
            new.Ms += other.Ms
            if new.Nly == 0:
                new.Nly = other.Nly
            if other.Nly == 0:
                new.Nly = new.Nly
            if (new.Nly > 0.) and (other.Nly > 0.):
                new.Nly = np.log10(10 ** self.Nly + 10 ** other.Nly)
            new.beta = new.calc_beta(new.wave, new.SED)
        assert isinstance(new, CSP)
        return new

    def __iadd__(self, other):
        new = None
        if isinstance(other, CSP):
            new = copy.deepcopy(self)
            new.SED += other.SED
            new.SFR += other.SFR
            new.Ms += other.Ms
            if new.Nly == 0:
                new.Nly = other.Nly
            if other.Nly == 0:
                new.Nly = new.Nly
            if (new.Nly > 0.) and (other.Nly > 0.):
                new.Nly = np.log10(10 ** self.Nly + 10 ** other.Nly)
            new.beta = new.calc_beta(new.wave, new.SED)
        assert isinstance(new, CSP)
        return new

    def __mul__(self, other):
        new = copy.deepcopy(self)
        new.SED *= other
        new.SFR *= other
        new.Ms *= other
        if other == 0.:
            new.Nly = 0.
        else:
            new.Nly += np.log10(other)
            new.Nly = np.maximum(new.Nly, 0)
        return new

    def __imul__(self, other):
        new = copy.deepcopy(self)
        new.SED *= other
        new.SFR *= other
        new.Ms *= other
        if other == 0.:
            new.Nly = 0.
        else:
            new.Nly += np.log10(other)
            new.Nly = np.maximum(new.Nly, 0)
        return new

    def __div__(self, other):
        new = copy.deepcopy(self)
        new.SED /= other
        new.SFR /= other
        new.Ms /= other
        if other == 0.:
            new.Nly = 0.
        else:
            new.Nly -= np.log10(other)
            new.Nly = np.maximum(new.Nly, 0)
        return new

    def __idiv__(self, other):
        new = copy.deepcopy(self)
        new.SED /= other
        new.SFR /= other
        new.Ms /= other
        if other == 0.:
            new.Nly = 0.
        else:
            new.Nly -= np.log10(other)
            new.Nly = np.maximum(new.Nly, 0)
        return new

    def __rmul__(self, other):
        new = copy.deepcopy(self)
        new.SED *= other
        new.SFR *= other
        new.Ms *= other
        if other == 0.:
            new.Nly = 0.
        else:
            new.Nly += np.log10(other)
            new.Nly = np.maximum(new.Nly, 0)
        return new

    def addEmissionLine(self, wavelength, EqW):
        wbin = np.argmin(np.abs(self.wave - wavelength))
        binwidth = np.mean(np.diff(self.wave)[wbin - 1:wbin + 1])
        continuum = np.mean(self.SED[wbin:wbin + 1])
        lineluminosity = continuum * EqW
        self.Lalpha = lineluminosity
        self.SED[wbin] += lineluminosity / binwidth


class Filter(object):
    def __init__(self):
        self.wave = []
        self.freq = []
        self.response = []

        self.lambda_c = []
        self.nu_c = []


class FileFilter(Filter):
    def __init__(self, filepath, maxbins=1000):
        """

        :type filepath: string
        :type maxbins: int
        """
        super(FileFilter, self).__init__()
        self.path = filepath

        data = np.loadtxt(self.path)
        wf = data[:, 0]
        tp = data[:, 1]
        if len(data[:, 0]) < maxbins:  # Re-sample large filters for performance
            wfx = np.linspace(wf[0], wf[-1], maxbins)
            tpx = griddata(wf, tp, wfx)

            wf = wfx
            tp = tpx

        self.wave = wf * u.angstrom
        self.response = tp

        self.freq = (c.c / self.wave).to(u.Hz)

        nmax = np.argmax(self.response)
        halfmax_low = self.wave[:nmax][np.argmin(np.abs(self.response[nmax] - 2 * self.response[:nmax]))]
        halfmax_hi = self.wave[nmax:][np.argmin(np.abs(self.response[nmax] - 2 * self.response[nmax:]))]

        self.fwhm = halfmax_hi - halfmax_low

        self.lambda_c = (simps(self.wave * self.response, self.wave) /
                         simps(self.response, self.wave)) * u.angstrom

        self.nu_c = (simps(self.freq * self.response, self.freq) /
                     simps(self.response, self.freq)) * u.Hz


class TophatFilter(Filter):
    def __init__(self, centre, width, steps=200):
        """

        :type steps: int
        """
        super(TophatFilter, self).__init__()
        self.centre = centre * u.angstrom
        self.width = width * u.angstrom
        self.steps = steps

        upper, lower = self.centre + self.width, self.centre - self.width
        resp_upper, resp_lower = self.centre + (self.width * 0.5), self.centre - (self.width * 0.5)

        self.wave = np.linspace(lower, upper, steps)
        self.response = np.zeros_like(self.wave.value)

        tophat = (self.wave >= resp_lower) * (self.wave < resp_upper)
        self.response[tophat] = 1

        self.freq = (c.c / self.wave).to(u.Hz)

        self.lambda_c = (simps(self.wave * self.response, self.wave) /
                         simps(self.response, self.wave))

        self.nu_c = (simps(self.freq * self.response, self.freq) /
                     simps(self.response, self.freq))
        self.fwhm = self.width


class LoadEAZYFilters(object):
    def __init__(self, path):
        self.path = path

        self.filters = []
        self.filternames = []
        self.central_wlengths = []

        with open(self.path) as file:
            for f in xrange(1000):
                x = file.readline().split()
                if len(x) < 1:
                    break
                nwave, name, lambda_c = x[0], x[1], x[4]
                nwave = int(nwave)

                wavelength = []
                response = []
                # print nwave
                for i in xrange(nwave):
                    N, w, r = np.array(file.readline().split()).astype('float')
                    wavelength.append(w)
                    response.append(r)

                wavelength *= u.angstrom
                freq = (c.c / wavelength).to(u.Hz)

                lambda_c = (simps(wavelength * response, wavelength) /
                            simps(response, wavelength))

                nu_c = (simps(freq * response, freq) /
                        simps(response, freq))

                new_filter = Filter()
                new_filter.wave = np.array(wavelength) * u.angstrom
                new_filter.response = np.array(response)
                new_filter.freq = np.array(freq) * u.Hz
                new_filter.lambda_c = lambda_c * u.angstrom
                new_filter.nu_c = nu_c * u.Hz
                nmax = np.argmax(new_filter.response)
                halfmax_low = new_filter.wave[:nmax][
                    np.argmin(np.abs(new_filter.response[nmax] - 2 * new_filter.response[:nmax]))]
                halfmax_hi = new_filter.wave[nmax:][
                    np.argmin(np.abs(new_filter.response[nmax] - 2 * new_filter.response[nmax:]))]
                #print new_filter.wave[nmax],halfmax_low, halfmax_hi
                new_filter.fwhm = halfmax_hi - halfmax_low

                self.filters.append(new_filter)
                self.filternames.append(name)
                self.central_wlengths.append(lambda_c * u.angstrom)

        self.central_wlengths *= u.angstrom


class FilterSet:
    def __init__(self, path=None):
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

    def addFileFilter(self, path):
        self.filters.append(FileFilter(path))

    def addTophatFilter(self, centre, width, steps=200):
        self.filters.append(TophatFilter(centre, width, steps))

    def addEAZYFilter(self, EAZYset, N):
        if type(N) == int:
            self.filters.append(EAZYset.filters[N])
        elif type(N) == list:
            for x in N:
                self.filters.append(EAZYset.filters[x])


class Observe:
    """
    
    """

    def __init__(self, SED, Filters, redshift, v=1, force_age=True, madau=True, units=u.uJy):
        self.SED = SED
        self.F = Filters
        self.redshifts = np.array(redshift, ndmin=1)
        self.wave = self.SED.wave

        self.fluxes = np.zeros((len(self.redshifts), len(self.F.filters))) * units
        self.AB = np.zeros((len(self.redshifts), len(self.F.filters))) * u.mag
        self.wl = np.zeros(len(self.F.filters)) * u.AA
        self.fwhm = np.zeros(len(self.F.filters)) * u.AA

        self.dl = cosmo.luminosity_distance(self.redshifts).cgs
        self.dl[self.redshifts == 0] = 10 * c.pc

        self.dm = cosmo.distmod(self.redshifts).value
        self.dm[self.redshifts == 0] = 0.

        for i, z in enumerate(self.redshifts):
            self.lyman_abs = np.ones(len(self.wave))
            if madau:
                self.lyman_abs = np.clip(self.tau_madau(self.wave, z), 0., 1.)

            if (self.SED.tg.to(u.Gyr) > cosmo.age(z)) and force_age:
                print 'SSP age older than universe...stopping.'
            else:
                for j, filter in enumerate(self.F.filters):
                    flux, mag = 0, 0
                    if v == 1:
                        flux, mag = self.calcflux(filter, z, self.dl[i], units)
                    elif v == 2:
                        flux, mag = self.calcFlux2(filter, z, self.dm[i], units)

                    self.wl[j] = filter.lambda_c
                    self.fwhm[j] = filter.fwhm
                    self.fluxes[i, j] = flux
                    self.AB[i, j] = mag * u.mag

        self.wl *= u.angstrom
        self.fwhm *= u.angstrom
        self.fluxes = np.squeeze(self.fluxes)  # * units
        self.AB = np.squeeze(self.AB)  # * u.mag

    @staticmethod
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
        return lyl_absorption

    def calcflux(self, filt, z, dl, units):
        wf = filt.wave
        tp = filt.response

        # Find SED wavelength entries within filter range
        wff = np.array([wf[0] < self.wave[i] < wf[-1]
                        for i in range(len(self.wave))])
        wft = self.wave[wff]

        # Interpolate to find throughput values at new wavelength points
        tpt = griddata(wf, tp, wft)

        # Join arrays and sort w.r.t to wf
        # Also replace units stripped by concatenate
        wf = np.array(np.concatenate((wf, wft))) * u.AA
        tp = np.concatenate((tp, tpt))

        order = np.argsort(wf)
        wf = wf[order]
        tp = tp[order]

        # Interpolate redshifted SED and LyAbs at new wavelength points
        sed = griddata(self.wave * (1 + z), self.SED.SED, wf) * self.SED.SED.unit
        lyabs = griddata(self.wave , self.lyman_abs, wf)

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

        top = np.trapz(sed * lyabs * tp * wf / c.c.to(u.AA / u.s), wf)
        bottom = np.trapz(tp / wf, wf)

        area = (4 * np.pi * (dl ** 2))

        F_mean = top / bottom / (1 + z) / area

        # Set flux to appropriate units and calculate AB magnitude
        Flux = F_mean.to(u.erg / u.cm ** 2 / u.s / u.Hz)
        ABmag = -2.5 * np.log10(Flux.to(u.Jy) / (3631 * u.Jy))

        return Flux.to(units), ABmag

    def calcFlux2(self, filt, z, dm, units):
        wf = filt.wave.value
        tp = filt.response
        z1 = z + 1

        if len(wf) > 1000:  # Re-sample large filters for performance
            wfx = np.linspace(wf[0], wf[-1], 1000)
            tpx = griddata(wf, tp, wfx)

            wf = wfx
            tp = tpx


        # Find SED wavelength entries within filter range
        wff = np.array([wf[0] < self.wave.value[i] < wf[-1]
                        for i in range(len(self.wave.value))])
        wft = self.wave[wff]

        #Interpolate to find throughput values at new wavelength points
        tpt = griddata(wf, tp, wft)

        #Join arrays and sort w.r.t to wf
        wf = np.concatenate((wf, wft))
        tp = np.concatenate((tp, tpt))

        order = np.argsort(wf)
        wf = wf[order]
        tp = tp[order]
        sed = griddata(self.wave * (1 + z), self.SED.SED.value, wf)
        lyabs = griddata(self.wave, self.lyman_abs, wf)
        dwf = np.diff(wf)
        nwf = len(wf)

        tpwf = tp / wf
        f_mean2 = np.dot(dwf, (tpwf[:nwf - 1] + tpwf[1:]) / 2)
        tpwf = tp * wf  #Reassign tpwf as product

        wf1 = wf / z1

        WR = 0.
        for i in range(nwf):
            #Interpolation indices
            j = np.where(self.wave.value < wf1[i])[0][-1]

            a = (wf1[i] - self.wave.value[j]) / (self.wave.value[j + 1] - self.wave.value[j])
            tpa = (tpwf[i] * ((1 - a) * (self.SED.SED.value[j] * self.lyman_abs[j]) +
                              a * self.SED.SED.value[j + 1] * self.lyman_abs[j + 1]))
            if i != 0:
                WR += dwf[i - 1] * (tpb + tpa)

            tpb = tpa

        F_mean = WR / 2 / z1 / f_mean2 / 2.997925e18
        AB0 = 5 * np.log10(1.7684e8 * 1e-5)
        # dl = 10pc in Mpc
        # this factor is sqrt(4*pi*(3.0856e24)^2 Lsun)

        #Convert fluxes to AB magnitudes
        Mag = AB0 - 2.5 * np.log10(F_mean) - 48.6
        Mag += dm

        Flux = 10 ** ((23.9 - Mag) / 2.5) * u.uJy  #uJy

        return Flux.to(units), Mag


# must define cosmo before calling an Observe
