import numpy as np
import array
import copy
import re
import sys
from glob import glob

from astropy import units as u
from astropy import constants as c

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
        self.ages[0] = (10 * u.yr)  # Make non-zero but insignificant
        self.wave_arr = np.array(self.wave_arr) * u.AA
        self.sed_arr = np.array(self.sed_arr).swapaxes(1, 2)[:, :]
        self.sed_arr *= u.Lsun / u.AA
        self.strm_arr = np.array(self.strm_arr)[:, :]
        self.rmtm_arr = np.array(self.rmtm_arr)[:, :]
        self.iseds = np.array(self.iseds)

