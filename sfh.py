import numpy as np
import array
import copy
import re
import sys
from glob import glob

from astropy import units as u
from astropy import constants as c

def exponential(t, tau):
    """ Exponential star formation history

    """
    sfh = np.exp(-1 * t / tau) / abs(tau)
    return sfh


def power(t, alpha):
    """ Power-law star formation history"""
    sfh = np.power(t, alpha)
    return sfh


def delayed(t, tau):
    """ 'Delated' star formation history"""
    sfh = t / (tau ** 2) * np.exp(-t / tau)
    return sfh


def truncated(t, tstop):
    """ Truncated star formation history

    Star-formation is continuous until some fraction, tstop, of the total
    time since onset of star-formation history np.max(t).

    """
    sfh = np.ones_like(t)
    cut = np.around(tstop*t.shape[1], 0) # Nearest whole timestep.
    sfh[:, cut:] = 0.

    sfh /= np.trapz(sfh, t)
    return sfh

