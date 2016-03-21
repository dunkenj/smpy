from __future__ import print_function, division
import numpy as np
import copy
import os
import sys
import h5py
import types
import glob

from scipy.interpolate import griddata
from scipy.spatial import Delaunay

from astropy import units as u
from astropy import constants as c
from astropy import cosmology as cos
from astropy.utils.console import ProgressBar

import data

from .dust import Calzetti
from .sfh import exponential
from .misc import tau_madau

cosmo = cos.FlatLambdaCDM(H0=70, Om0=0.3)

f = open("error.log", "w")
original_stderr = sys.stderr
sys.stderr = f

data_path = data.__path__[0]


class CSP:
    """ Class for building composite stellar populations from input SSPs
        
    Attributes
    ----------
    
    Methods
    -------
    
    """
    
    def __init__(self, ssp,
                 age=None, sfh=None, dust=None, metal_ind=None, f_esc=None,
                 sfh_law=exponential, dust_model=Calzetti, neb_cont=True,
                 neb_met=True):
        """
            
        Parameters
        ----------
            ssp :
            
            age : '~astropy.units.Quantity' array (with units of time)
                Desired stellar population age(s) since the onset of
                star-formation
            sfh :  array of float or '~astropy.units.Quantity',
                Star-formation history parameter/exponent
            dust : array of floats,
                Strength of dust extinction in Av
            metal_ind : float or float array,
                Metallicity or metallicity range of stellar population 
                relative to solar metallicity (Z/Z_sol), min and max
                allowed values set by range of input metallicities in 
                'ssp' class
            f_esc : float or array of float,
                Escape fraction(s) of nebular emission component
            sfh_law : '~smpy.sfh' or user defined function,
                Star-formation history parametrisation for composite
                stellar population.
            dust_model : '~smpy.dust' or user defined function,
                Dust extinction or attenuation law parametrisation
            neb_cont : boolean, default = True
                Include continuum emission component in nebular emission
                model
            neb_met : boolean, default = True
                Include metal line component in nebular emission model
        """
        
        self.beta = None
        self.SSP = ssp
        self.metallicities = self.SSP.metallicities
        # Normalise to solar metallicity
        self.ages = self.SSP.ages
        
        self.wave = self.SSP.wave_arr[0]
        self.iw = self.SSP.iw_arr[0]
        self.sed_arr = self.SSP.sed_arr
        
        if hasattr(self.SSP, 'strm_arr') and hasattr(self.SSP, 'rmtm_arr'):
            self.strm_arr = np.array(self.SSP.strm_arr)
            self.rmtm_arr = np.array(self.SSP.rmtm_arr)
        
        self.iseds = np.array(self.SSP.iseds)
        
        # Find closest match for each tg value in ta - set tg to these values
        
        nebular = np.loadtxt(data_path+'/nebular_emission.dat', skiprows=1)
        self.neb_cont = nebular[:, 1]
        self.neb_hlines = nebular[:, 2]
        self.neb_metal = nebular[:, 3:]
        self.neb_wave = nebular[:, 0]
        
        if None not in (age, sfh, dust, metal_ind):
            if f_esc == None:
                self.build(age, sfh, dust, metal_ind, sfh_law=sfh_law,
                           dust_model=dust_model,
                           neb_cont=neb_cont, neb_met=neb_met)
            
            else:
                self.build(age, sfh, dust, metal_ind,
                           f_esc, sfh_law, dust_model, neb_cont, neb_met)
    
    def build(self, age, sfh, dust, metal, fesc = 1.,
              sfh_law = exponential, dust_model = Calzetti,
              neb_dust_weight = 1.,
              neb_cont = True, neb_met = True,
              timesteps = 500, verbose = False):
        
        """ Build composite stellar population SED(s) from input SSP models
            
            Parameters
            ----------
                ssp :
                
                age : '~astropy.units.Quantity' array (with units of time)
                    Desired stellar population age(s) since the onset of
                    star-formation
                sfh :  array of float or '~astropy.units.Quantity',
                    Star-formation history parameter/exponent
                dust : array of floats,
                    Strength of dust extinction in Av
                metal_ind : float or float array,
                    Metallicity or metallicity range of stellar population
                    relative to solar metallicity (Z/Z_sol), min and max
                    allowed values set by range of input metallicities in 'ssp' 
                    class
                f_esc : float or array of float,
                    Escape fraction(s) of nebular emission component
                sfh_law : '~smpy.sfh' or user defined function,
                    Star-formation history parametrisation for composite
                    stellar population.
                dust_model : '~smpy.dust' or user defined function,
                    Dust extinction or attenuation law parametrisation
                neb_cont : boolean, default = True
                    Include continuum emission component in nebular emission
                    model
                neb_met : boolean, default = True
                    Include metal line component in nebular emission model
            
            Returns
            -------
        
        
        """
        try:
            self.tau = u.Quantity(sfh, ndmin=1)
        except:
            self.tau = sfh
        
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
        
        # Set up nebular emission arrays -- WILL CHANGE
        if len(self.neb_wave) != len(self.wave):
            self.neb_cont = griddata(self.neb_wave, 
                                     self.neb_cont,
                                     self.wave)
            
            self.neb_hlines = griddata(self.neb_wave, 
                                       self.neb_hlines,
                                       self.wave)
                                       
            neb_metaln = np.zeros((len(self.wave), 3))
            for i in range(3):
                neb_metaln[:, i] = griddata(self.neb_wave, 
                                            self.neb_metal[:, i], 
                                            self.wave)
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
                                              (self.neb_metal[:, 0] * 
                                              self.inc_met))
        
        self.neb_sed_arr[nebrange2, :, :] *= ((self.neb_cont * self.inc_cont) +
                                              self.neb_hlines +
                                              (self.neb_metal[:, 1] * 
                                              self.inc_met))
        
        self.neb_sed_arr[nebrange3, :, :] *= ((self.neb_cont * self.inc_cont) +
                                              self.neb_hlines +
                                              (self.neb_metal[:, 2] * 
                                              self.inc_met))
        
        self.neb_sed_arr *= (c.c.to(u.AA / u.s) / (self.wave ** 2)).value  
        # Convert to Flambda
        self.neb_sed_arr =  (self.neb_sed_arr * self.sed_arr.unit /
                             self.Nly_arr.unit)
        
        
        # Set up grid for ND-interpolation
        ti, mi = np.meshgrid(np.log10(self.ages / u.yr).value,
                             np.log10(self.metallicities))
        self.grid = zip(mi.flatten(), ti.flatten())
        
        tri_grid = Delaunay(self.grid)
        
        # Make cube of interpolated age and SFHs.
        # Uses array slicing and vectorisation to minimise
        # loops where possible.
        sfh_grid_shape = (len(self.mi),
                          timesteps)
        
        self.me_sfh = np.ones(sfh_grid_shape)
        
        for met_idx, metal in enumerate(self.mi):
            self.me_sfh[met_idx] *= metal
        
        if verbose:
            bar = ProgressBar(np.product(self.SED.shape[1:-1]))
        
        for idG, age in enumerate(self.tg):
            ta_range = np.logspace(np.log10(self.ages / u.yr).min(),
                                   np.log10(age / u.yr), timesteps)
            self.ta_sfh = np.ones(sfh_grid_shape) * u.yr
            self.ta_sfh *= ta_range[None, :]
        
        # Calculate Barycentric coordinates for all ages/metallicities in SFH.
            points = np.array(zip(np.log10(self.me_sfh.flatten()),
                              np.log10(self.ta_sfh.flatten() / u.yr)))
            
            #if verbose:
            #    print 'Interpolating SEDs at SFH timesteps'
            ss = tri_grid.find_simplex(points)
            
            X = tri_grid.transform[ss, :2]
            Y = points - tri_grid.transform[ss, 2]
            
            b = np.einsum('ijk,ik->ij', X, Y)
            self.bc = np.c_[b, 1 - b.sum(axis=1)]
            self.simplices = tri_grid.simplices[ss]
            
            #if verbose:
            #    print 'Reshaping grids'
            # Interpolate SED, stellar mass fraction and remnant fractions
            # for SFH age grid using calculated Barycentric coordinates (bc).
            #print self.sed_arr.shape
            #print (len(self.metallicities) * len(self.ages), self.iw)
            
            self.temp = self.sed_arr.reshape(len(self.metallicities) *
                                             len(self.ages),
                                             self.iw)[self.simplices]
                                             
            self.sed_sfh = (self.temp * self.bc[:, :, None]).sum(1)
            self.sed_sfh = self.sed_sfh.reshape(np.append(sfh_grid_shape,
                                                          self.iw))
            
            self.temp = self.neb_sed_arr.reshape(len(self.metallicities) *
                                                 len(self.ages), 
                                                 self.iw)[self.simplices]
                                                 
            self.neb_sed_sfh = (self.temp * self.bc[:, :, None]).sum(1)
            self.neb_sed_sfh = self.neb_sed_sfh.reshape(np.append(
                                                        sfh_grid_shape, 
                                                        self.iw))
            
            self.strm_sfh = np.array(self.strm_arr.reshape(
                                     len(self.metallicities) * 
                                     len(self.ages))[self.simplices] *
                                     self.bc).sum(1)
            self.strm_sfh = self.strm_sfh.reshape(sfh_grid_shape)
            
            self.rmtm_sfh = np.array(self.rmtm_arr.reshape(
                                     len(self.metallicities) *
                                     len(self.ages))[self.simplices] * 
                                     self.bc).sum(1)
                                     
            self.rmtm_sfh = self.rmtm_sfh.reshape(sfh_grid_shape)
            
            self.Nly_sfh = (self.Nly_arr.reshape(len(self.metallicities) *
                                                 len(self.ages))[self.simplices]
                                                 * self.bc).sum(1)
            self.Nly_sfh = self.Nly_sfh.reshape(sfh_grid_shape) #*
                            #(1 - self.fesc[None, None, :, None]))
            
            # Star-formation history
            
            
            
            for idT, t in enumerate(self.tau):
                if type(t) == tuple:
                    tau = t
                else:
                    tau = tuple([t])
                
                self.sfr_hist = self.sfr_func(self.ta_sfh, *tau)
                # Enforce integrated SFR = 1 Msol.
                self.norm = np.trapz(self.sfr_hist, 
                                     self.ta_sfh, axis=-1)[:, None]
                
                self.sfr_hist /= self.norm
                self.weights = np.abs(self.sfr_func(self.tg[None, idG, None] \
                                      - self.ta_sfh, *tau)) / self.norm
                self.sfh_weights = np.ones(sfh_grid_shape) * self.weights
                
                for idA, Av in enumerate(self.tauv):
                    for idf, fesc in enumerate(self.fesc):
                        self.Att = self.getAttenuation(self.ta_sfh, \
                                                       self.wave, Av)
                        
                        neb_att = self.getAttenuation(self.ta_sfh, \
                                                      self.wave, 
                                                      Av * neb_dust_weight)
                        
                        combined_sed = self.sed_sfh
                        
                        # Absorbed LyC photons
                        combined_sed[:, :, self.wave <= 912 * u.AA] *= fesc
                        
                        # Resulting nebular emission
                        combined_sed *= self.Att # Dust attenuated combined SED
                        combined_sed += (neb_att * ((1 - fesc) * 
                                        self.Nly_sfh[:, :, None] * 
                                        self.neb_sed_sfh))
                        
                        # Integrate over star-formation history
                        self.SED[:, idG, idT, idA, idf, :] = \
                            np.trapz(self.sfh_weights[:, :, None] * \
                            combined_sed, self.ta_sfh[:, :, None], axis=-2)
                        
                        self.STR[:, idG, idT, idA, idf] = \
                            np.trapz(self.sfh_weights * self.strm_sfh,
                                     self.ta_sfh)
                        
                        self.SFR[:, idG, idT, idA, idf] = \
                            self.sfr_hist[:, -1] * u.solMass
                            
                        if verbose:
                            bar.update()
        
        # Normalise by stellar mass fraction (stars + remnants in case of BC03)
        self.SFR /= self.STR
        self.SED = self.SED / self.STR[:, :, :, :, :, None]
        self.Ms = np.ones_like(self.SFR.value) * u.Msun
        
        
        self.Nly = self.calc_lyman(self.wave, self.SED).cgs
    
    
    @staticmethod
    def calc_lyman(wave, seds):
        """ Calculate total Lyman continuum photons for an SED
            
            Semi-pythonic version of the fortran routine in BC03 code,
            allows flexible input array shapes provided wavelength axis
            is the last one.
            
            Parameters
            ----------
            
            wave : array of floats,
                Wavelength array
            seds :
        
        """
        wly = 912. * u.AA
        const = (1e-8 / (u.AA / u.cm)) / c.h.cgs / c.c.cgs
        
        n = int(sum([wave < wly][0]))
        S = np.array(seds.shape)
        S[-1] = n+1
        f = np.zeros(S) * seds.unit * u.AA
        w = np.zeros(n + 1) * u.AA
        
        for i in range(n + 1):
            if wave[i] <= wly:
                w[i] = wave[i]
                np.rollaxis(f, -1)[i] = w[i] * np.rollaxis(seds, -1)[i]
            elif wave[i] > wly:
                w[i] = wly
                np.rollaxis(f, -1)[i] = w[i] * (np.rollaxis(seds,-1)[i - 1] + \
                 ((w[i] - wave[i - 1]) * (np.rollaxis(seds, -1)[i] - \
                 np.rollaxis(seds, -1)[i - 1]) / (wave[i] - wave[i - 1])))
        
        nlyman = const * np.trapz(f, w, axis=-1)
        return nlyman
    
    def __add__(self, other):
        new = None
        if isinstance(other, CSP):
            new = copy.deepcopy(self)
            new.SED += other.SED
            new.SFR += other.SFR
            new.Ms += other.Ms
            new.Nly = self.Nly + other.Nly
        
        assert isinstance(new, CSP)
        return new
    
    def __iadd__(self, other):
        new = None
        if isinstance(other, CSP):
            new = copy.deepcopy(self)
            new.SED += other.SED
            new.SFR += other.SFR
            new.Ms += other.Ms
            new.Nly = self.Nly + other.Nly
        
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
            new.Nly /= other
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
            new.Nly /= other
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
            new.Nly /= other
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
            new.Nly /= other
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
            new.Nly /= other
            new.Nly = np.maximum(new.Nly, 0)
        return new
    
    def addEmissionLine(self, wavelength, EqW):
        """ Add emission line to SED
        
        """
        wbin = np.argmin(np.abs(self.wave - wavelength))
        binwidth = np.mean(np.diff(self.wave)[wbin-1:wbin+1])
        continuum = self.SED[:, :, :, :, :, wbin:wbin+1].mean(-1)
        
        print(continuum.unit)
        lineluminosity = continuum * EqW
        
        print((lineluminosity / binwidth).unit)
        self.Lalpha = lineluminosity
        self.SED[:, :, :, :, :,wbin] += (lineluminosity / binwidth)

class Filter(object):
    def __init__(self):
        self.wave = []
        self.freq = []
        self.response = []
        
        self.lambda_c = []
        self.nu_c = []


class FileFilter(Filter):
    """
        
        Attributes
        ----------
        
        Methods
        -------
            
            __init__ :
    
    
    """
    def __init__(self, filepath, maxbins=500):
        """ Build filter class from ascii file
            
            Parameters
            ----------
            filepath : str
                Filename (including relative/absolute path) for filter file.
                Should be readable as ascii file and contain two columns,
                where column one is wavelength in angstroms and columns two
                is the filter response.
            maxbins : int
                Maximum number of wavelength bins for filter. Large (>1000)
                numbers of wavelength bins result in slower calculations for
                minimal gain in accuracy.
        
        """
        super(FileFilter, self).__init__()
        self.path = filepath
        
        data = np.loadtxt(self.path)
        wf = data[:, 0]
        tp = data[:, 1]
        if len(data[:, 0]) > maxbins:  # Re-sample large filters for performance
            wfx = np.linspace(wf[0], wf[-1], maxbins)
            tpx = griddata(wf, tp, wfx)
            
            wf = wfx
            tp = tpx
        
        self.wave = wf * u.angstrom
        self.response = tp
        
        self.freq = (c.c / self.wave).to(u.Hz)
        
        nmax = np.argmax(self.response)
        halfmax_low = (self.wave[:nmax][np.argmin(np.abs(self.response[nmax] -
                       2 * self.response[:nmax]))])
        halfmax_hi = (self.wave[nmax:][np.argmin(np.abs(self.response[nmax] -
                       2 * self.response[nmax:]))])
        
        self.fwhm = halfmax_hi - halfmax_low
        
        self.lambda_c = (np.trapz(self.wave * self.response, self.wave) /
                         np.trapz(self.response, self.wave))
        
        self.nu_c = (np.trapz(self.freq * self.response, self.freq) /
                     np.trapz(self.response, self.freq))


class TophatFilter(Filter):
    def __init__(self, centre, width, steps=200):
        """ Build filter class for any arbitrary tophat function
            
            Parameters
            ----------
            centre : int or float
                Central wavelength for tophat function in Angstroms
            width : int or float
                Width of tophat function in Angstroms
            steps : int or float
                Number of wavelength steps to use. Should be a relatively fine
                grid as SED interpolated to wavelength grid of the filter during
                convolution.
        
        """
        super(TophatFilter, self).__init__()
        self.centre = centre.to(u.angstrom)
        self.width = width.to(u.angstrom)
        self.steps = steps
        
        upper, lower = self.centre + self.width, self.centre - self.width
        resp_upper = self.centre + (self.width * 0.5)
        resp_lower = self.centre - (self.width * 0.5)
        
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
    """ Load EAZY filter file into easily accessible library
        
        Attributes
        ----------
        
        Methods
        -------
    
    """
    def __init__(self, path):
        """ Load in EAZY filters.res formatted files
            
            Parameters
            ----------
            path : str
                EAZY filter file to load
        
        """
        
        self.path = path
        
        self.filters = []
        self.filternames = []
        
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
                new_filter.wave = wavelength
                new_filter.response = np.array(response)
                new_filter.freq = np.array(freq)
                new_filter.lambda_c = lambda_c
                new_filter.nu_c = nu_c * u.Hz
                nmax = np.argmax(new_filter.response)
                halfmax_low = new_filter.wave[:nmax][
                    np.argmin(np.abs(new_filter.response[nmax] - 2 * 
                    new_filter.response[:nmax]))]
                halfmax_hi = new_filter.wave[nmax:][
                    np.argmin(np.abs(new_filter.response[nmax] - 2 * 
                    new_filter.response[nmax:]))]
                #print new_filter.wave[nmax],halfmax_low, halfmax_hi
                new_filter.fwhm = halfmax_hi - halfmax_low
                
                self.filters.append(new_filter)
                self.filternames.append(name)
                if f == 0:
                    self.central_wl = np.array(lambda_c)
                else:
                    self.central_wl = np.append(self.central_wl, lambda_c)
        
        #self.central_wl = self.central_wl*lambda_c.unit
    
    def search(self, searchterm, error=0.1):
        if type(searchterm) == str:
            matches = np.array([searchterm in name for 
                                name in self.filternames])
            sorted_id = np.arange(len(self.filternames))[matches]
        elif type(searchterm) == u.quantity.Quantity:
            sorted_id = np.where(np.abs(self.central_wl - searchterm) <
                                 error*searchterm)[0]
        
        print('{0:<5s} {1:<10s}' \
              'Angstrom {2}'.format('Index',
                                    'Lambda_c',
                                    'Name'))
        for i in sorted_id:
            print('{0:<5.0f} {1:<8.2f}' \
                  'Angstrom {2}'.format(i,
                                        self.central_wl[i].value,
                                        self.filternames[i]))
        return sorted_id


class FilterSet:
    def __init__(self, path=None):
        self.directory = path
        self.filters = []
        if type(self.directory) == str:
            try:
                self.files = glob.glob(self.directory)
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
        elif type(N) == list or type(N) == np.ndarray:
            for x in N:
                self.filters.append(EAZYset.filters[x])


class Observe(object):
    """ Mock photometry class
        
        Attributes
        ----------
        
        Methods
        -------
    
    """
    
    def __init__(self, SED, Filters, redshift, force_age=True,
                 madau=True, units=u.uJy, verbose=False):
        """
        
        Parameters
        ----------
        SED : '~smpy.CSP' object
            Built
        Filters : '~smpy.FilterSet' object
            Filter set through which to observe the set of models included
            in SED object
        redshift : float of numpy.array
            Redshift(s) at which models are to be observed
        force_age : boolean
            Require age of the stellar population to be younger than
            the age of the Universe at the desired redshift.
        madau : boolean
            Apply IGM absorption following Madau 1999 prescription
        units : '~astropy.units.Quantity'
            Desired output units, must be in spectral flux density equivalent
        verbose : boolean
            Add additional terminal outputs if true
        
        
        Attributes
        ----------
        fluxes : Apparent fluxes of CSP models observed through 'Filters' at
            the desired redshifts
        AB : Apparent AB magnitudes of CSP models observed through 'Filters' at
            the desired redshifts
        
        
        Examples
        --------
        
        >>> redshifts = np.linspace(0, 3, 10)
        >>> A = Observe(CSP, Filters, redshifts)
        
        """
        self.F = Filters
        self.redshifts = np.array(redshift, ndmin=1)
        self.wave = SED.wave
        self.Ms = SED.Ms
        self.SFR = SED.SFR
        
        self.fluxes = np.zeros(np.append([len(self.redshifts),
                                          len(self.F.filters)],
                                          SED.SED.shape[:-1])) * units
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
                self.lyman_abs = np.clip(tau_madau(self.wave, z), 0., 1.)
            
            for j, filter in enumerate(self.F.filters):
                self.wl[j] = filter.lambda_c
                self.fwhm[j] = filter.fwhm
                self.fluxes[i, j] = self.calcflux(SED, filter, z, 
                                                  self.dl[i], units)
                
                if not force_age:
                    # Set fluxes for ages older than universe to zero
                    agecut = (SED.tg.to(u.Gyr) > cosmo.age(z))
                    self.fluxes[i, j, :, agecut] = 0.
            
            if verbose:
                assert isinstance(bar, object)
                bar.update()
        
        # Convert spectral flux density to AB magnitudes
        self.AB = (-2.5 * np.log10(self.fluxes.to(u.Jy) / 
                   (3631 * u.Jy))) * u.mag
    
    def calcflux(self, SED, filt, z, dl, units):
        """ Convolve synthetic SEDs with a given filter
            
            Arguments
            ---------
                SED : numpy.array
                    Grid of synthetic spectra
                filt : '~smpy.Filter' class
                    Filter through which to convolve SED grid
                z : float
                    Redshift at which models are to be observed
                dl : '~astropy.units.Quantity'
                    Luminosity distance corresponding to redshift(z) in given
                    cosmology.
                units : '~astropy.units'
                    Desired output flux units (in spectral flux density)
            
            Returns
            -------
                Flux : '~astropy.units.Quantity'
                    Spectral flux density, with exact units as given by 'units'
        
        """
        # Find SED wavelength entries within filter range
        wff = np.logical_and(filt.wave[0] < self.wave, 
                             self.wave < filt.wave[-1])
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
        top = np.trapz(sed * lyabs[None, None, None, None, None, :] * tp * wf /
                       c.c.to(u.AA / u.s), wf)
        bottom = np.trapz(tp / wf, wf)
        area = (4 * np.pi * (dl ** 2))
        Flux = top / bottom / (1 + z) / area
        
        return Flux.to(units)
    
    def __getitem__(self, items):
        return self.fluxes[items]


class ObserveToFile(object):
    """ Mock photometry class
        
        Attributes
        ----------
        
        Methods
        -------
    
    """
    def __init__(self):
        self.test = 1.
    
    def build(self, SED, Filters, redshift, savepath, v=1,
              force_age=True, madau=True, units=u.uJy,
              verbose=False, clobber=True):
        """
            
            Parameters
            ----------
            SED : '~smpy.CSP' object
                Built
            Filters : '~smpy.FilterSet' object
                Filter set through which to observe the set of models included
                in SED object
            
            redshift :
            
            savepath : str,
                Filename for hdf5 save file
            
            v :
            
            force_age : boolean
                Require age of the stellar population to be younger than
                the age of the Universe at the desired redshift.
            madau : boolean
                Apply IGM absorption following Madau 1999 prescription
            units : '~astropy.units.Quantity'
                Desired output units, must be in spectral flux density
                equivalent
            verbose : boolean
                Add additional terminal outputs if true
        
        """
        self.F = Filters
        self.redshifts = np.array(redshift, ndmin=1)
        self.wave = SED.wave
        
        #self.fluxes = np.zeros(gridshape) * units
        #self.AB = np.zeros(gridshape) * u.mag
        self.wl = np.zeros(len(self.F.filters)) * u.AA
        self.fwhm = np.zeros(len(self.F.filters)) * u.AA
        
        self.dl = cosmo.luminosity_distance(self.redshifts).cgs
        self.dl[self.redshifts == 0] = 10 * c.pc
        
        
        assert type(savepath) is types.StringType, "File save path is not a string: %r" % savepath
        self.savepath = savepath
        if clobber:
            if os.path.isfile(self.savepath):
                os.remove(self.savepath)
        
        with h5py.File(savepath, 'w') as f:
            
            gridshape = np.append([len(self.redshifts), 
                                  len(self.F.filters)], 
                                  SED.SED.shape[:-1])
                                  
            self.fluxes = f.create_dataset("fluxes", gridshape, dtype='f')
            self.fluxes.attrs['unit'] = units.to_string()
            
            self.AB = f.create_dataset("mags", gridshape, dtype='f')
            
            for j, filt in enumerate(self.F.filters):
                self.wl[j] = filt.lambda_c
                self.fwhm[j] = filt.fwhm
                
                print('Filter {0}' \
                      '(Central WL = {1:.1f}):'.format(j+1,
                                                       filt.lambda_c))
                
                # Find SED wavelength entries within filter range
                wff = np.logical_and(filt.wave[0] < self.wave, 
                                     self.wave < filt.wave[-1])
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
                
                if verbose:
                    niters = len(self.redshifts)*len(SED.tg)
                    bar = ProgressBar(niters)
                
                for i, z in enumerate(self.redshifts):
                    self.lyman_abs = np.ones(len(self.wave))
                    if madau:
                        self.lyman_abs = np.clip(tau_madau(self.wave, z), 
                                                 0., 1.)
                    
                    for a, age in enumerate(SED.tg):
                        if np.logical_or(age < cosmo.age(z), force_age):
                            fluxes = self.calcflux(SED.SED[:,a], wf, tp,
                                                   z, self.dl[i], units)
                            self.fluxes[i, j, :, a] = fluxes.to(units)
                        else:
                            # Set fluxes for ages older
                            # than universe to zero
                            self.fluxes[i, j, :, a] = 0.
                        if verbose:
                            assert isinstance(bar, object)
                            bar.update()
                print('\n')
                
            for i, z in enumerate(self.redshifts):
                for j, filt in enumerate(self.F.filters):
                    self.AB[i,j] = (-2.5 * 
                                    np.log10((self.fluxes[i,j]*units).to(u.Jy) 
                                    / (3631 * u.Jy))).value
                    zeros = (self.fluxes[i,j] == 0.)
                    self.AB[i,j][zeros] = np.inf
            
            # Set up hdf5 dimensions based on SED ranges for convenient
            # slicing if not used for fitting.
            f.create_dataset('z', data = self.redshifts)
            f.create_dataset('ages', data = SED.tg.value)
            f['ages'].attrs['unit'] = SED.tg.unit.to_string()

            f.create_dataset('dust', data = SED.tauv)
            f.create_dataset('metallicities', data = SED.mi)
            f.create_dataset('fesc', data = SED.fesc)
            
            try:
                s = SED.taus.shape
                taus = f.create_dataset('sfh', data = SED.tau)
                taus.attrs['unit'] = SED.tau.unit.to_string()
            except:
                # For multi-parameter SFHs, use indices for scale
                taus = f.create_dataset('sfh', data = np.arange(len(SED.tau)))
                
            
            for dataset in ['fluxes', 'mags']:
                f[dataset].dims.create_scale(f['z'], 'z')
                f[dataset].dims.create_scale(f['sfh'], 'sfh')
                f[dataset].dims.create_scale(f['ages'], 'ages')
                f[dataset].dims.create_scale(f['dust'], 'dust')
                f[dataset].dims.create_scale(f['metallicities'], 'met')
                f[dataset].dims.create_scale(f['fesc'], 'fesc')
                
                f[dataset].dims[0].attach_scale(f['z'])
                f[dataset].dims[2].attach_scale(f['metallicities'])
                f[dataset].dims[3].attach_scale(f['ages'])
                f[dataset].dims[4].attach_scale(f['sfh'])
                f[dataset].dims[5].attach_scale(f['dust'])
                f[dataset].dims[6].attach_scale(f['fesc'])
            
            # Store all CSP attributes in case they are needed later
            for attribute in SED.__dict__.keys():
                if attribute == 'SED':
                    # SED excluded due to large size
                    continue
                try:
                    f.create_dataset(attribute, data=SED.__dict__[attribute])
                except:
                    # in cases of functions etc.
                    continue
        
        self.load(self.savepath)
    
    def load(self, loadpath):
        self.f = h5py.File(loadpath, 'r')
        
        self.fluxes = self.f['fluxes']
        self.AB = self.f['mags']
        
        self.Ms = self.f['Ms'].value
        self.SFR = self.f['SFR'].value
    
    
    def calcflux(self, SED, wf, tp, z, dl, units):
        """ Convolve synthetic SEDs with a given filter
            
            Arguments
            ---------
                SED : numpy.array
                    Grid of synthetic spectra
                filt : '~smpy.Filter' class
                    Filter through which to convolve SED grid
                z : float
                    Redshift at which models are to be observed
                dl : '~astropy.units.Quantity'
                    Luminosity distance corresponding to redshift(z) in 
                    given cosmology.
                units : '~astropy.units'
                    Desired output flux units (in spectral flux density)
            
            Returns
            -------
                Flux : '~astropy.units.Quantity'
                    Spectral flux density, with exact units as given by 'units'
        
        """
        
        # Interpolate redshifted SED and LyAbs at new wavelength points
        sed = griddata(self.wave * (1 + z), SED.T, wf).T * SED.unit
        lyabs = griddata(self.wave, self.lyman_abs, wf)
        
        # Calculate f_nu mean
        # Integrate SED through filter, as per BC03 Fortran
        # As: f_nu=int(dnu Fnu Rnu/h*nu)/int(dnu Rnu/h*nu)
        # ie: f_nu=int(dlm Flm Rlm lm / c)/int(dlm Rlm/lm)
        top = np.trapz(sed * lyabs[None, None, None, None, :] * tp * wf /
                       c.c.to(u.AA / u.s), wf)
        bottom = np.trapz(tp / wf, wf)
        area = (4 * np.pi * (dl ** 2))
        Flux = top / bottom / (1 + z) / area
        
        return Flux.to(units)
    
    def __getitem__(self, items):
        return self.fluxes[items] * u.Unit(self.f['fluxes'].attrs['unit'])
