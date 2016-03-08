import numpy as np
import astropy.units as u

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
        tau = np.where(wave <= lyn_w[i] * (1 + z), 
                       tau + lyn_c[i] * (wave / lyn_w[i]) ** 3.46, 
                       tau)

    wly = wave / lymanlim
    wly3 = wly ** 3
    pe_abs = (0.25 * wly3 * ((1 + z) ** 0.46 - wly ** 0.46) +
              9.4 * wly ** 1.5 * ((1 + z) ** 0.18 - wly ** 0.18) -
              0.7 * wly3 * (wly ** (-1.32) - (1 + z) ** (-1.32)) -
              0.023 * ((1 + z) ** 1.68 - wly ** 1.68))
    tau = np.where(wave <= lymanlim * (1 + z), tau + pe_abs, tau)
    lyl_absorption = np.where(tau > 700., 0., np.exp(-tau))
    return np.clip(lyl_absorption, 0., 1.)