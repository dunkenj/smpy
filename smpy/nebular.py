import re
import numpy as np
import data
data_path = data.__path__[0]

from astropy.modeling import models
from astropy.table import Table
from astropy import units as u
from astropy import constants as c

metallicities = np.array([2e-2, 8e-3, 4e-3, 4e-4, 1e-5, 1e-7])
metallicities /= metallicities[0]

inoue = np.loadtxt('{0}/{1}'.format(data_path, 'LineRatio_nodust.txt'))

names = open('{0}/{1}'.format(data_path, 'LineList.txt'), 'r')
line_names = ['_'.join(re.split(' ', n.rstrip())[1:]) for n in
              names.readlines()]

line_ratios = inoue[:, 2:-2:2].T
line_lambdas = inoue[:,1] * u.AA

def inoue_lines_model(wave, fwhm=100*u.km/u.s, fwhm_lyman_m=5):
    """
    Evaluate Inoue 2011 Nebular line model for a given wavelength grid
    """
    # Gaussian StdDev in Angstrom from velocity FWHM
    line_sigmas = (fwhm/c.c.to(u.km/u.s) * line_lambdas.to(wave.unit) /
                   2*np.sqrt(2*np.log(2)))
    line_sigmas[0] *= fwhm_lyman_m
    line_models = []
    
    for ixr, ratios in enumerate(line_ratios):
        line_amplitudes = ratios / np.sqrt(2*np.pi*(line_sigmas**2))
        #Built the model set
        line_models.append(models.Gaussian1D(mean=line_lambdas.to(wave.unit), 
                                         amplitude=line_amplitudes, 
                                         stddev=line_sigmas))
    return line_models

def inoue_lines(wave, fwhm=100*u.km/u.s):
    """
    Evaluate Inoue 2011 Nebular line model for a given wavelength grid
    """
    # Gaussian StdDev in Angstrom from velocity FWHM
    line_models = inoue_lines_model(wave, fwhm)

    line_sed = np.zeros((len(wave), len(metallicities)))
    for ixr, models in enumerate(line_models):
        line_sed[:,ixr] = np.array(map(models, wave.value)).sum(1)
    return line_sed / wave.unit, metallicities