import numpy as np


def exponential(t, tau):
    """ Exponential star formation history

    """
    sfh = np.exp(-1 * t / tau) / np.abs(tau)
    return sfh


def power(t, alpha):
    """ Power-law star formation history
    """
    sfh = np.power(t, alpha)
    return sfh

def dblpower(t, sl):
    """ Power-law star formation history
    """
    alpha, beta, tau = sl
    sfh = (np.power(np.abs(t)/(tau*np.max(t)), alpha) + np.power(np.abs(t)/(tau*np.max(t)), beta))**-1
    return sfh

def delayed(t, tau):
    """ 'Delated' star formation history
    """
    sfh = t / (tau ** 2) * np.exp(-t / tau)
    return sfh


def truncated(t, tstop):
    """ Truncated star formation history

    Star-formation is continuous until some fraction, tstop, of the total
    time since onset of star-formation history np.max(t).

    """
    sfh = 1 / np.ones_like(t)
    cut = np.argmin(np.abs(t - tstop*np.max(t))) # Nearest whole timestep.
    np.rollaxis(sfh, -1)[cut:] = 0.
    return sfh


def truncated_exp(t, sl):
    """ Truncated exponential star formation history

    Star-formation is exponential, tstop, of the total
    time since onset of star-formation history np.max(t).

    """
    tau, tstop = sl
    sfh = np.exp(-1 * t / tau) / abs(tau)
    cut = np.argmin(np.abs(t - tstop*np.max(t))) # Nearest whole timestep.
    np.rollaxis(sfh, -1)[cut:] = 0.
    return sfh
