import numpy as np
import array
import copy
import re
import sys
from glob import glob
from scipy.interpolate import griddata
from scipy.optimize import leastsq
from scipy.spatial import Delaunay

from astropy import units as u
from astropy import constants as c
from astropy import cosmology as cos
from astropy.utils.console import ProgressBar

from ssp import Ised, SSP, BC
from dust import Charlot, Calzetti, Calzetti2, MW, LMC, SMC
from sfh import exponential

cosmo = cos.FlatLambdaCDM(H0=70, Om0=0.3)

f = open("error.log", "w")
original_stderr = sys.stderr
sys.stderr = f


class CSP:
    """ Class for building composite stellar populations from input SSPs
            
    """

    def __init__(self, ssp,
                 age=None, sfh=None, dust=None, metal_ind=None, f_esc=None,
                 sfh_law=exponential, dust_model=Calzetti, neb_cont=True, neb_met=True):
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

        nebular = np.loadtxt('data/nebular_emission.dat', skiprows=1)
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


    #@profile
    def build(self, age, sfh, dust, metal, fesc=1.,
              sfh_law=exponential, dust_model=Calzetti,
              neb_cont=True, neb_met=True, timesteps = 400, verbose=False):
        """ Docs

        """
        self.tau = u.Quantity(sfh, ndmin=1)
        self.tg = u.Quantity(age, ndmin=1).to(u.yr)

        self.tauv = np.array(dust, ndmin=1)
        self.mi = np.array(metal, ndmin=1)
        self.fesc = np.array(fesc, ndmin=1)
        self.sfh_law = sfh_law
        self.inc_cont = neb_cont
        self.inc_met = neb_met
        self.getAttenuation = dust_model
        self.sfr_func = sfh_law
        
        mu = 0.3

        self.ta = self.ages
        
        outshape = [len(self.mi), len(self.tg), 
                    len(self.tau), len(self.tauv), 
                    len(self.fesc), len(self.wave)]
        self.SED = np.zeros(outshape) * self.sed_arr.unit
        self.STR = np.zeros(self.SED.shape[: -1])
        self.SFR = np.zeros(self.SED.shape[: -1]) * u.solMass / u.yr

        # Set up grid for NDinterpolation
        ti, mi = np.meshgrid(np.log10(self.ages / u.yr), np.log10(self.metallicities))
        self.grid = zip(mi.flatten(), ti.flatten())

        tri_grid = Delaunay(self.grid)

        # Make cube of interpolated age and SFHs.
        # Uses array slicing and vectorisation to minimise
        # loops where possible.
        sfh_grid_shape = (len(self.mi), len(self.tg), 
                          timesteps)

        self.ta_sfh = np.ones(sfh_grid_shape) * u.yr
        self.me_sfh = np.ones(sfh_grid_shape)
                
        for met_idx, metal in enumerate(self.mi):
            self.me_sfh[met_idx] *= metal

        for age_idx, age in enumerate(self.tg):
            ta_range = np.logspace(np.log10(self.ages / u.yr).min(), 
                                   np.log10(age / u.yr), timesteps)
            self.ta_sfh[:, age_idx] *= ta_range[None, :]

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

        self.neb_sed_arr *= (c.c.to(u.AA / u.s) / (self.wave ** 2)).value  # Convert to Flambda
        self.neb_sed_arr =  self.neb_sed_arr * self.sed_arr.unit / self.Nly_arr.unit

        # Calculate Barycentric coordinates for all ages/metallicities in SFH.
        points = np.array(zip(np.log10(self.me_sfh.flatten()), np.log10(self.ta_sfh.flatten() / u.yr)))

        if verbose:
            print 'Interpolating SEDs at SFH timesteps'
        ss = tri_grid.find_simplex(points)

        X = tri_grid.transform[ss, :2]
        # Offset of each target from the origin of its containing tetrahedron.
        Y = points - tri_grid.transform[ss, 2]

        b = np.einsum('ijk,ik->ij', X, Y)
        self.bc = np.c_[b, 1 - b.sum(axis=1)]
        self.simplices = tri_grid.simplices[ss]

        # Interpolate SED, stellar mass fraction and remnant fractions
        # for SFH age grid using calculated Barycentric coordinates (bc).

        self.sed_sfh = (self.sed_arr.reshape(len(self.metallicities) *
                                             len(self.ages), self.iw)[self.simplices] *
                                             self.bc[:, :, None]).sum(1)
        self.sed_sfh = self.sed_sfh.reshape(np.append(sfh_grid_shape, self.iw))
        
        self.neb_sed_sfh = (self.neb_sed_arr.reshape(len(self.metallicities) *
                                                     len(self.ages), self.iw)[self.simplices]
                                                     * self.bc[:, :, None]).sum(1)
        self.neb_sed_sfh = self.neb_sed_sfh.reshape(np.append(sfh_grid_shape, self.iw))
        
        self.strm_sfh = np.array(self.strm_arr.reshape(len(self.metallicities) *
                                                       len(self.ages))[self.simplices]
                                                        * self.bc).sum(1)
        self.strm_sfh = self.strm_sfh.reshape(sfh_grid_shape)

        self.rmtm_sfh = np.array(self.rmtm_arr.reshape(len(self.metallicities) *
                                                       len(self.ages))[self.simplices]
                                                        * self.bc).sum(1)
        self.rmtm_sfh = self.rmtm_sfh.reshape(sfh_grid_shape)
        
        self.Nly_sfh = (self.Nly_arr.reshape(len(self.metallicities) *
                                            len(self.ages))[self.simplices]
                                            * self.bc).sum(1)
                                                     
        self.Nly_sfh = self.Nly_sfh.reshape(sfh_grid_shape) #* 
                        #(1 - self.fesc[None, None, :, None]))

        # Star-formation history

        if verbose:
            bar = ProgressBar(np.product(self.SED.shape[2:-1]))
            
        for idT, t in enumerate(self.tau):
            if type(t) == tuple:
                tau = t
            else:
                tau = tuple([t])

            self.sfh_weights = np.ones(sfh_grid_shape)

            self.sfr_hist = self.sfr_func(self.ta_sfh, *tau)
            # Enforce integrated SFR = 1 Msol.
            self.norm = np.trapz(self.sfr_hist, self.ta_sfh, axis=-1)[:, :, None]

            self.sfr_hist /= self.norm
            self.weights = self.sfr_func(self.tg[None, :, None] - self.ta_sfh, 
                                         *tau) / self.norm
            self.sfh_weights *= self.weights


                # Offset added to renormalise from B to V band.
            for idA, Av in enumerate(self.tauv):
                for idf, fesc in enumerate(self.fesc):
                    self.Att = self.getAttenuation(self.ta_sfh, self.wave, Av)
                    combined_sed = self.sed_sfh
                
                    # Absorbed LyC photons
                    combined_sed[:, :, :, self.wave <= 912 * u.AA] *= (1-fesc)
                    # Resulting nebular emission
                    combined_sed += ((1 - fesc) * self.Nly_sfh[:, :, :, None] * self.neb_sed_sfh)
                    combined_sed *= self.Att # Dust attenuated combined SED
                
                    # Integrate over star-formation history
                    self.SED[:, :, idT, idA, idf, :] = np.trapz(self.sfh_weights[:, :, :, None] * combined_sed,
                                                              self.ta_sfh[:, :, :, None].value,
                                                              axis=-2)
                    self.STR[:, :, idT, idA, idf] = np.trapz(self.sfh_weights * self.strm_sfh, self.ta_sfh)

                    self.SFR[:, :, idT, idA, idf] = self.sfr_hist[:, :, -1] * u.solMass
                    if verbose:
                        bar.update()

        self.SFR /= self.STR
        self.SED = self.SED / self.STR[:, :, :, :, :, None]

        self.Nly = self.calc_lyman_f(self.wave, self.SED).cgs


        # self.beta = self.calc_beta(self.wave,self.SED)
        self.Ms = np.ones_like(self.SFR) * u.Msun

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

        nlyman = const * np.trapz(f, w, axis=-1)
        # print np.log10(N_lyman)
        return nlyman
        
    @staticmethod
    def calc_lyman_f(wave, seds):
        wly = 912. * u.AA
        const = (1e-8 / (u.AA / u.cm)) / c.h.cgs / c.c.cgs

        n = int(sum([wave < wly][0]))
        S = np.array(seds.shape)
        S[-1] = n+1
        f = np.zeros(S) * seds.unit * u.AA
        w = np.zeros(n + 1) * u.AA

        for i in range(n + 1):
            if wave[i] < wly:
                w[i] = wave[i]
                f[:, :, :, :, :, i] = w[i] * seds[:, :, :, :, :, i]
            elif wave[i] == wly:
                w[i] = wave[i]
                f[:, :, :, :, :, i] = w[i] * seds[:, :, :, :, :, i]
            elif wave[i] > wly:
                w[i] = wly
                f[:, :, :, :, :, i] = w[i] * (seds[:, :, :, :, :, i - 1] + (
                    (w[i] - wave[i - 1]) * (seds[:, :, :, :, :, i] - seds[:, :, :, :, :, i - 1]) / (wave[i] - wave[i - 1])))
                # f[:,:,i] = f[:,:,i]

        nlyman = const * np.trapz(f, w, axis=-1)
        # print np.log10(N_lyman)
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

        self.lambda_c = (np.trapz(self.wave * self.response, self.wave) /
                         np.trapz(self.response, self.wave)) * u.angstrom

        self.nu_c = (np.trapz(self.freq * self.response, self.freq) /
                     np.trapz(self.response, self.freq)) * u.Hz


class TophatFilter(Filter):
    def __init__(self, centre, width, steps=200):
        """

        :type steps: int
        """
        super(TophatFilter, self).__init__()
        self.centre = centre.to(u.angstrom)
        self.width = width.to(u.angstrom)
        self.steps = steps

        upper, lower = self.centre + self.width, self.centre - self.width
        resp_upper, resp_lower = self.centre + (self.width * 0.5), self.centre - (self.width * 0.5)

        self.wave = np.linspace(lower, upper, steps)
        self.response = np.zeros_like(self.wave.value)

        tophat = (self.wave >= resp_lower) * (self.wave < resp_upper)
        self.response[tophat] = 1

        self.freq = (c.c / self.wave).to(u.Hz)

        self.lambda_c = (np.trapz(self.wave * self.response, self.wave) /
                         np.trapz(self.response, self.wave))

        self.nu_c = (np.trapz(self.freq * self.response, self.freq) /
                     np.trapz(self.response, self.freq))
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

                lambda_c = (np.trapz(wavelength * response, wavelength) /
                            np.trapz(response, wavelength))

                nu_c = (np.trapz(freq * response, freq) /
                        np.trapz(response, freq))

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
    """ Mock photometry class
    """

    def __init__(self, SED, Filters, redshift, v=1, force_age=True, madau=True, units=u.uJy, verbose=False):
        """

        :type self: object
        """
        self.F = Filters
        self.redshifts = np.array(redshift, ndmin=1)
        self.wave = SED.wave

        self.fluxes = np.zeros(np.append([len(self.redshifts), 
                                          len(self.F.filters)], SED.SED.shape[:-1])) * units
        self.AB = np.zeros_like(self.fluxes.value) * u.mag
        self.wl = np.zeros(len(self.F.filters)) * u.AA
        self.fwhm = np.zeros(len(self.F.filters)) * u.AA

        self.dl = cosmo.luminosity_distance(self.redshifts).cgs
        self.dl[self.redshifts == 0] = 10 * c.pc

        if verbose:
            bar = ProgressBar(len(self.redshifts))

        for i, z in enumerate(self.redshifts):
            self.lyman_abs = np.ones(len(self.wave))
            if madau:
                self.lyman_abs = np.clip(self.tau_madau(self.wave, z), 0., 1.)

            for j, filter in enumerate(self.F.filters):
                flux = self.calcflux(SED, filter, z, self.dl[i], units)

                self.wl[j] = filter.lambda_c
                self.fwhm[j] = filter.fwhm
                self.fluxes[i, j] = flux
            
                if not force_age:
                    agecut = (SED.tg.to(u.Gyr) > cosmo.age(z))
                    self.fluxes[i, j, :, agecut] = 0.    
            
            if verbose:
                assert isinstance(bar, object)
                bar.update()
        
        self.AB = (-2.5 * np.log10(self.fluxes.to(u.Jy) / (3631 * u.Jy))) * u.mag
        
        self.wl *= u.angstrom
        self.fwhm *= u.angstrom

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

    def calcflux(self, SED, filt, z, dl, units):
        """ Convolve synthetic SEDs with a given filter

            Arguments:
                SED (np.array): Grid of synthetic spectra
                filt (Filter object):
                z (float): Redshift at which models are being observed
                dl (astropy.quantity): luminosity distance
                     corresponding to redshift(z) in given cosmology.
                units (astropy.units) : Desired output flux units (in spectral flux density)

            Returns:
                Flux (astropy.quantity): with units as given by 'units'
                       (spectral flux density)
        """
        # Find SED wavelength entries within filter range
        wff = np.logical_and(filt.wave[0] < self.wave, self.wave < filt.wave[-1])
        wft = self.wave[wff]

        # Interpolate to find throughput values at new wavelength points
        tpt = griddata(filt.wave, filt.response, wft)

        # Join arrays and sort w.r.t to wf
        # Also replace units stripped by concatenate
        wf = np.array(np.concatenate((filt.wave, wft))) * u.AA
        tp = np.concatenate((filt.response, tpt))

        order = np.argsort(wf)
        wf = wf[order]
        tp = tp[order]

        # Interpolate redshifted SED and LyAbs at new wavelength points
        sed = griddata(self.wave * (1 + z), SED.SED.T, wf).T * SED.SED.unit
        lyabs = griddata(self.wave, self.lyman_abs, wf)

        # Calculate f_nu mean
        # Integrate SED through filter, as per BC03 Fortran
        # As: f_nu=int(dnu Fnu Rnu/h*nu)/int(dnu Rnu/h*nu)
        # ie: f_nu=int(dlm Flm Rlm lm / c)/int(dlm Rlm/lm)
        top = np.trapz(sed * lyabs[None, None, None, None, None, :] * tp * wf / c.c.to(u.AA / u.s), wf)
        bottom = np.trapz(tp / wf, wf)
        area = (4 * np.pi * (dl ** 2))
        Flux = top / bottom / (1 + z) / area

        return Flux.to(units)
