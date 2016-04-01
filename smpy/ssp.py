import numpy as np
import array
import re
import os
from glob import glob

from astropy import units as u
from astropy.utils.console import ProgressBar

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

            Paramaters
            ----------
            
            filename : str, filename (including relative/full path) of
                Bruzual & Charlot model to be loaded into class.
    
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
            Paramaters
            ----------
            path : str,
                Path to desired subset of Bruzual & Charlot models, should include
                relevant wildcard to match multiple metallicities.
                
        """
        super(BC, self).__init__(path)

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
        self.ages[0] = (10000 * u.yr)  # Make non-zero but insignificant
        self.wave_arr = np.array(self.wave_arr) * u.AA
        sed_arr = np.array(self.sed_arr).swapaxes(1, 2)[:, :]
        self.sed_arr = sed_arr * (u.Lsun / u.AA)
        self.strm_arr = np.array(self.strm_arr)[:, :]
        self.rmtm_arr = np.array(self.rmtm_arr)[:, :]
        self.iseds = np.array(self.iseds)

class BPASS(SSP):
    """ BPASS
    
    www.bpass.org.uk
    Eldridge & Stanway, 2009, MNRAS, 400, 1019
    
    """

    def __init__(self, path='../ssp/bpass_v1.1/SEDS/', nsteps = 3000):
        """ 
            Paramaters
            ----------
            path : str,
                Path to desired subset of BPASS models, should include
                relevant wildcard to match multiple metallicities.
            
        """
        super(BPASS, self).__init__(path)
        head, tail = os.path.split(path)
        readme = open(head+'/'+'sed.bpass.readme.txt','r')
        ages = [1]
        for line in readme.readlines():
            if 'log(Age/yrs)=' in line:
                n = line.split(')')[2].strip('=').strip()
                ages.append(float(n))
            else:
                continue

        #self.ages[0] = 5.
        self.ages = np.power(10, ages)
        
        self.files = glob(self.SSPpath)
        self.files.sort()
        self.iseds = []
        self.ta_arr = []
        self.metal_arr = []
        self.iw_arr = []
        self.wave_arr = []
        self.sed_arr = []
        self.strm_arr = []
        self.rmtm_arr = []
        
        metallicities = np.zeros(len(self.files))
        with ProgressBar(len(self.files)) as bar:
            for i, file in enumerate(self.files):
                data = np.loadtxt(file)
                #self.wave_arr.append(data[:,0])
                #sed = np.copy(data)
                #sed[:,0] *= 0. # Make zero age column
    #

                wave_large = data[:,0] * u.AA
                wave_sht = np.unique(np.logspace(np.log10(wave_large.min()/u.AA),
                                                 np.log10(wave_large.max()/u.AA), nsteps).astype('int')) *u.AA

                sed = np.copy(data)
                sed[:,0] *= 0. # Make zero age column
                print('rebin')
                sed_sht = rebin_sed(wave_large, sed, wave_sht)
                sed_sht[np.isnan(sed_sht)] = 0.
                self.wave_arr.append(wave_sht)
                self.ta_arr.append(self.ages)
                self.metal_arr.append(file)
                self.iw_arr.append(sed_sht.shape[0])
                self.sed_arr.append(sed_sht)
                self.strm_arr.append(np.ones_like(self.ages))
                self.rmtm_arr.append(np.ones_like(self.ages))
                self.iseds.append(None)
                metallicities[i] = float(file.split('z')[-1])
                bar.update()
            
        self.metallicities = metallicities / 1000 / 0.02  # Normalise to solar metallicity
        self.ages = np.array(self.ta_arr[0]) * u.yr
        self.wave_arr = np.array(self.wave_arr) * u.AA
        self.sed_arr = np.array(self.sed_arr).swapaxes(1, 2)[:, :] / 1e6
        self.sed_arr *= u.Lsun / u.AA
        self.strm_arr = np.array(self.strm_arr)[:, :]
        self.rmtm_arr = np.array(self.rmtm_arr)[:, :]
        self.iseds = np.array(self.iseds)


class BPASS2(SSP):
    """ BPASS v2
    
    www.bpass.org.uk
    Eldridge & Stanway, 2009, MNRAS, 400, 1019
    
    """

    def __init__(self, path='../ssp/bpass_v1.1/SEDS/', nsteps = 3000):
        """ 
            Parameters
            ----------
            path : str,
                Path to desired subset of BPASS models, should include
                relevant wildcard to match multiple metallicities.
            nsteps : int
                Number of wavelength steps in which to interpolate 
            
        """
        super(BPASS2, self).__init__(path)

        #self.ages[0] = 5.
        self.ages = np.power(10, 6+(0.1*np.arange(41)))
        self.ages = np.insert(self.ages, 0, 1e5)
        
        self.files = glob(self.SSPpath)
        self.files.sort()
        self.iseds = []
        self.ta_arr = []
        self.metal_arr = []
        self.iw_arr = []
        self.wave_arr = []
        self.sed_arr = []
        self.strm_arr = []
        self.rmtm_arr = []
        
        metallicities = np.zeros(len(self.files))
        with ProgressBar(len(self.files)) as bar:
            for i, file in enumerate(self.files):
                data = np.loadtxt(file)
            
                wave_large = data[:,0] * u.AA
                wave_sht = np.unique(np.logspace(np.log10(wave_large.min()/u.AA),
                                                 np.log10(wave_large.max()/u.AA), nsteps).astype('int')) *u.AA

                sed = np.copy(data)
                sed[:,0] *= 0. # Make zero age column
                sed_sht = rebin_sed(wave_large, sed, wave_sht)
                self.wave_arr.append(wave_sht)
            

                self.ta_arr.append(self.ages)
                self.metal_arr.append(file)
                self.iw_arr.append(sed_sht.shape[0])
                self.sed_arr.append(sed_sht)
                self.strm_arr.append(np.ones_like(self.ages))
                self.rmtm_arr.append(np.ones_like(self.ages))
                self.iseds.append(None)
                metallicities[i] = float(os.path.splitext(file)[0].split('z')[-1])
                bar.update()
        self.metallicities = metallicities / 1000 / 0.02  # Normalise to solar metallicity
        self.ages = np.array(self.ta_arr[0]) * u.yr
        self.wave_arr = np.array(self.wave_arr) * u.AA
        self.sed_arr = np.array(self.sed_arr).swapaxes(1, 2)[:, :] / 1e6
        self.sed_arr *= u.Lsun / u.AA
        self.strm_arr = np.array(self.strm_arr)[:, :]
        self.rmtm_arr = np.array(self.rmtm_arr)[:, :]
        self.iseds = np.array(self.iseds)

def rebin_sed(wave_in, sed_in, wave_out):
    """ Rebin an SED to new wavelength grid with flux density conserved
    
        Parameters
        ----------
            wave_in : numpy.array or '~astropy.units.Quantity' array
                Input wavelength grid, with length N
            sed_in : numpy.array or '~astropy.units.Quantity' array
                Input grid of SEDs, with shape (N, x)
            wave_out : numpy.array or '~astropy.units.Quantity' array
                Desired wavelength grid for rebinned SED, with length M
    
        Returns
        -------
            sed_out : numpy.array or '~astropy.units.Quantity' array
                Rebinned SEDs, with shape (M, x)
    
    """

    edges = np.empty(len(wave_out) + 1, dtype=np.float64) * wave_in.unit
    edges[1:-1] = (wave_out[1:] + wave_out[:-1]) * 0.5

    # Compute the first and last by making them symmetric
    edges[0] = 2.0 * wave_out[0] - edges[1]
    edges[-1] = 2.0 * wave_out[-1] - edges[-2]

    indices = np.searchsorted(wave_in, edges)
    i_beg = indices[:-1]
    i_end = indices[1:]

    avflux = sed_in[:-1] #(sed_in[1:, :] + sed_in[:-1, :]) * 0.5
    deltaw = wave_in[1:] - wave_in[:-1]

    #print avflux.shape
    sed_out = np.empty(shape=(len(wave_out), sed_in.shape[1]), dtype=np.float64)

    for i in range(len(i_beg)):
        first = i_beg[i]
        last = i_end[i]
        cur_dw = deltaw[first:last]
        avg = np.sum(avflux[first:last, :] * cur_dw[:, None], 0) / cur_dw.sum()
        sed_out[i, :] = avg

    return sed_out
        
