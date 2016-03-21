import numpy as np
from astropy import units as u
#from astropy import constants as c

""" 
Functional Forms:

"""

def dust_func(lam, ai, bi, ni, li):
    """
    Functional form for SMC, LMC and MW extinction curves of
    Pei et al. 1992

    """
    lam = np.array(lam) / 1e4
    ki = np.power((lam / li), ni) + np.power((li / lam), ni) + bi
    eta_i = ai / ki
    return eta_i

def drude_profile(lam, Eb, lambda_central=2175*u.AA, delta_lambda=350*u.AA):
    top = Eb * (lam*delta_lambda)**2
    bottom = (lam**2 - lambda_central**2)**2 + (lam*delta_lambda)**2
    return top/bottom


"""
Attenuation and Extinction Laws

"""

def Charlot(ta_grid, wave, Av, mu=0.3):
    Att = np.ones(np.append(ta_grid.shape, len(wave)))
    tv = ((Av / 1.0857) * np.ones_like(ta_grid))
    tv[ta_grid > 1e7 * u.yr] *= mu
    lam = ((5500 * u.AA / wave) ** 0.7)
    Att *= (np.exp(-1 * np.outer(tv, lam))).reshape(Att.shape)
    return Att

def Calzetti(ta_grid, wave, Av):
    Att = np.ones(np.append(ta_grid.shape, len(wave)))
    k = np.zeros_like(wave.value)

    w0 = [wave <= 1200 * u.AA]
    w1 = [wave < 6300 * u.AA]
    w2 = [wave >= 6300 * u.AA]
    w_u = wave.to(u.um).value

    x1 = np.argmin(np.abs(wave - 1200 * u.AA))
    x2 = np.argmin(np.abs(wave - 1250 * u.AA))

    k[w2] = 2.659 * (-1.857 + 1.040 / w_u[w2])
    k[w1] = 2.659 * (-2.156 + (1.509 / w_u[w1]) - (0.198 / w_u[w1] ** 2) + (0.011 / w_u[w1] ** 3))
    k[w0] = k[x1] + ((wave[w0] - 1200. * u.AA) * (k[x1] - k[x2]) / (wave[x1] - wave[x2]))

    k += 4.05
    k[k < 0.] = 0.

    tv = Av * k / 4.05
    Att *= np.power(10, -0.4 * tv)
    return Att


def Calzetti2(ta_grid, wave, Av):
    Att = np.ones(np.append(ta_grid.shape, len(wave)))
    k = np.zeros_like(wave.value)

    w0 = [wave <= 1000 * u.AA]
    w1 = [(wave > 1000 * u.AA) * (wave < 6300 * u.AA)]
    w2 = [wave >= 6300 * u.AA]
    w_u = wave.to(u.um).value

    k[w2] = 2.659 * (-1.857 + 1.040 / w_u[w2])
    k[w1] = 2.659 * (-2.156 + (1.509 / w_u[w1]) - (0.198 / w_u[w1] ** 2) + (0.011 / w_u[w1] ** 3))

    p1 = dust_func(wave, 27, 4, 5.5, 0.08) + dust_func(wave, 185, 90, 2, 0.042)

    k[w0] = p1[w0] / (p1[w1][0] / k[w1][0])
    k += 4.05
    k[k < 0.] = 0.
    tv = Av * k / 4.05
    Att *= np.power(10, -0.4 * tv)
    return Att


def SMC(ta_grid, wave, Av):
    Att = np.ones(np.append(ta_grid.shape, len(wave)))
    ai = [185., 27., 0.005, 0.01, 0.012, 0.03]
    bi = [90., 5.5, -1.95, -1.95, -1.8, 0.]
    ni = [2., 4., 2., 2., 2., 2.]
    li = [0.042, 0.08, 0.22, 9.7, 18., 25.]

    eta = np.zeros_like(wave)
    for i in xrange(len(ai)):
        eta += dust_func(wave, ai[i], bi[i], ni[i], li[i])

    Rv = 2.93
    Ab = Av * (1 + (1 / Rv))

    Att *= np.power(10, -0.4 * Ab * eta)
    return Att


def LMC(ta_grid, wave, Av):
    Att = np.ones(np.append(ta_grid.shape, len(wave)))
    ai = [175., 19., 0.023, 0.005, 0.006, 0.02]
    bi = [90., 4.0, -1.95, -1.95, -1.8, 0.]
    ni = [2., 4.5, 2., 2., 2., 2.]
    li = [0.046, 0.08, 0.22, 9.7, 18., 25.]

    eta = np.zeros_like(wave)
    for i in xrange(len(ai)):
        eta += dust_func(wave, ai[i], bi[i], ni[i], li[i])

    Rv = 3.16
    Ab = Av * (1 + (1 / Rv))

    Att *= np.power(10, -0.4 * Ab * eta)
    # Offset added to renormalise from B to V band.
    return Att


def MW(ta_grid, wave, Av):
    Att = np.ones(np.append(ta_grid.shape, len(wave)))    
    ai = [165., 14., 0.045, 0.002, 0.002, 0.012]
    bi = [90., 4., -1.95, -1.95, -1.8, 0.]
    ni = [2., 6.5, 2., 2., 2., 2.]
    li = [0.047, 0.08, 0.22, 9.7, 18., 25.]

    eta = np.zeros_like(wave)
    for i in xrange(len(ai)):
        eta += dust_func(wave, ai[i], bi[i], ni[i], li[i])

    Rv = 3.08
    Ab = Av * (1 + (1 / Rv))

    Att *= np.power(10, -0.4 * Ab * eta)
    return Att

def flexible(ta_grid, wave, Av, delta = 0., e_bump = 0., 
             lambda_central=2175*u.AA, delta_lambda=350*u.AA):
    Att = np.ones(np.append(ta_grid.shape, len(wave)))
    k = np.zeros_like(wave.value)

    w0 = [wave <= 1200 * u.AA]
    w1 = [wave < 6300 * u.AA]
    w2 = [wave >= 6300 * u.AA]
    w_u = wave.to(u.um).value

    x1 = np.argmin(np.abs(wave - 1200 * u.AA))
    x2 = np.argmin(np.abs(wave - 1250 * u.AA))

    k[w2] = 2.659 * (-1.857 + 1.040 / w_u[w2])
    k[w1] = 2.659 * (-2.156 + (1.509 / w_u[w1]) - (0.198 / w_u[w1] ** 2) + (0.011 / w_u[w1] ** 3))
    k[w0] = k[x1] + ((wave[w0] - 1200. * u.AA) * (k[x1] - k[x2]) / (wave[x1] - wave[x2]))

    k += 4.05
    k[k < 0.] = 0.

    bump = drude_profile(wave, e_bump, lambda_central, delta_lambda)

    tv = Av * (k + bump) * (wave/5500*u.AA)**delta / 4.05 
    Att *= np.power(10, -0.4 * tv)
    return Att
