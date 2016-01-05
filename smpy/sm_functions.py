import array
import numpy
import math

from scipy.integrate import trapz,simps
from scipy.interpolate import griddata
from scipy.optimize import leastsq, fsolve


def read_ised(filename):
    """
    This function reads data from Bruzual & Charlot binary format
    SSP files and returns the necessary data in an array The input files
    should be '.ised' files, either 2003 or 2007.
    
    'ks' in the binary files is slightly different between 03/07 files
    so the read length and index should be set appropriately, therefore 
    the function tries '03 format first and retries with the '07 format
    if the returned number of ages isn't as expected (e.g. 221 ages)
    """

    with open(filename,'rb') as f:
        check = array.array('i')
        check.fromfile(f,2)
        
    if check[1] == 221:
        ksl, ksi = 2, 1
    else:
        ksl, ksi = 3, 2
        
    with open(filename,'rb') as f:
        ks = array.array('i')
        ks.fromfile(f,ksl)

        ta = array.array('f')
        ta.fromfile(f,ks[ksi])
        ta = numpy.array(ta)

        tmp = array.array('i')
        tmp.fromfile(f,3)
        ml,mul,iseg = tmp

        if iseg > 0:
            tmp = array.array('f')
            tmp.fromfile(f,iseg*6)

        tmp = array.array('f')
        tmp.fromfile(f,5)
        totm, totn, avs, jo, tauo = tmp


        ids= array.array('c')
        ids.fromfile(f,80)

        tmp = array.array('f')
        tmp.fromfile(f,4)
        tcut = tmp[0]
        ttt = tmp[1:]

        ids = array.array('c')
        ids.fromfile(f,80)

        ids = array.array('c')
        ids.fromfile(f,80)

        igw = array.array('i')
        igw.fromfile(f,1)

        tmp = array.array('i')
        tmp.fromfile(f,3)

        iw = array.array('i')
        iw.fromfile(f,1)

        wave = array.array('f')
        wave.fromfile(f,iw[0])
        wave = numpy.array(wave)

        #SED Section
        F = array.array('i')
        F.fromfile(f,3)
        iw = F[2] #Number of wavelength elements

        sed = numpy.zeros((iw,ks[ksi]),dtype=numpy.float32)
        G = array.array('f')
        G.fromfile(f,iw)
        sed[:,0] = G
        ik = array.array('i')
        ik.fromfile(f,1)

        h = numpy.empty((ik[0],ks[ksi]),'f')
        H = array.array('f')
        H.fromfile(f,ik[0])
        h[:,0] = H

        for i in range(1,ks[ksi]): #Fill rest of array with SEDs
            F = array.array('i')
            F.fromfile(f,3)
            iw = F[2]

            G = array.array('f')
            G.fromfile(f,iw)
            sed[:,i] = G
            ik = array.array('i')
            ik.fromfile(f,1)

            H = array.array('f')
            H.fromfile(f,ik[0])
            h[:,i] = H

        tmp = array.array('i')
        tmp.fromfile(f,3)

        bflx = array.array('f')
        bflx.fromfile(f,tmp[2])

        tmp = array.array('i')
        tmp.fromfile(f,3)

        strm = array.array('f')
        strm.fromfile(f,tmp[2])
        strm = numpy.array(strm)

        tmp = array.array('i')
        tmp.fromfile(f,3)

        evf = array.array('f')
        evf.fromfile(f,tmp[2])

        tmp = array.array('i')
        tmp.fromfile(f,3)

        evf = array.array('f')
        evf.fromfile(f,tmp[2])

        tmp = array.array('i')
        tmp.fromfile(f,3)

        snr = array.array('f')
        snr.fromfile(f,tmp[2])

        tmp = array.array('i')
        tmp.fromfile(f,3)

        pnr = array.array('f')
        pnr.fromfile(f,tmp[2])

        tmp = array.array('i')
        tmp.fromfile(f,3)

        sn = array.array('f')
        sn.fromfile(f,tmp[2])

        tmp = array.array('i')
        tmp.fromfile(f,3)

        bh = array.array('f')
        bh.fromfile(f,tmp[2])

        tmp = array.array('i')
        tmp.fromfile(f,3)

        wd = array.array('f')
        wd.fromfile(f,tmp[2])

        tmp = array.array('i')
        tmp.fromfile(f,3)

        rmtm = array.array('f')
        rmtm.fromfile(f,tmp[2])
        rmtm = numpy.array(rmtm)

        data = [[ta,ids,iw,wave,sed,strm,rmtm],
                ['Ages(yr)', 'SSP ID', 
                 'No. wavelength entries',
                 'Wavelength array',
                 'SEDs', 'Stellar mass history',
                 'Remnant mass history']]

    return data

def read_ised2(filename):
    """
    This function reads data from Bruzual & Charlot binary format
    SSP files and returns the necessary data in an array

    The input files should be '.ised' files
    
    Original function, kept as backup in case the less verbose
    version is bugged in some way not yet found.
    """

    with open(filename,'rb') as f:
        check = array.array('i')
        check.fromfile(f,2)
 
    if check[1] == 221:
        with open(filename,'rb') as f:

            ks = array.array('i')
            ks.fromfile(f,2)

            ta = array.array('f')
            ta.fromfile(f,ks[1])
            ta = numpy.array(ta)

            tmp = array.array('i')
            tmp.fromfile(f,3)
            ml,mul,iseg = tmp

            if iseg > 0:
                tmp = array.array('f')
                tmp.fromfile(f,iseg*6)

            tmp = array.array('f')
            tmp.fromfile(f,5)
            totm, totn, avs, jo, tauo = tmp


            ids= array.array('c')
            ids.fromfile(f,80)
 
            tmp = array.array('f')
            tmp.fromfile(f,4)
            tcut = tmp[0]
            ttt = tmp[1:]

            ids = array.array('c')
            ids.fromfile(f,80)

            ids = array.array('c')
            ids.fromfile(f,80)

            igw = array.array('i')
            igw.fromfile(f,1)
 
            tmp = array.array('i')
            tmp.fromfile(f,3)

            iw = array.array('i')
            iw.fromfile(f,1)
 
            wave = array.array('f')
            wave.fromfile(f,iw[0])
            wave = numpy.array(wave)

            #SED Section
            F = array.array('i')
            F.fromfile(f,3)
            iw = F[2] #Number of wavelength elements

            sed = numpy.zeros((iw,ks[1]),dtype=numpy.float32)
            G = array.array('f')
            G.fromfile(f,iw)
            sed[:,0] = G
            ik = array.array('i')
            ik.fromfile(f,1)

            h = numpy.empty((ik[0],ks[1]),'f')
            H = array.array('f')
            H.fromfile(f,ik[0])
            h[:,0] = H

            for i in range(1,ks[1]): #Fill rest of array with SEDs
                F = array.array('i')
                F.fromfile(f,3)
                iw = F[2]

                G = array.array('f')
                G.fromfile(f,iw)
                sed[:,i] = G
                ik = array.array('i')
                ik.fromfile(f,1)

                H = array.array('f')
                H.fromfile(f,ik[0])
                h[:,i] = H

            tmp = array.array('i')
            tmp.fromfile(f,3)

            bflx = array.array('f')
            bflx.fromfile(f,tmp[2])

            tmp = array.array('i')
            tmp.fromfile(f,3)

            strm = array.array('f')
            strm.fromfile(f,tmp[2])
            strm = numpy.array(strm)

            tmp = array.array('i')
            tmp.fromfile(f,3)

            evf = array.array('f')
            evf.fromfile(f,tmp[2])

            tmp = array.array('i')
            tmp.fromfile(f,3)

            evf = array.array('f')
            evf.fromfile(f,tmp[2])

            tmp = array.array('i')
            tmp.fromfile(f,3)

            snr = array.array('f')
            snr.fromfile(f,tmp[2])

            tmp = array.array('i')
            tmp.fromfile(f,3)

            pnr = array.array('f')
            pnr.fromfile(f,tmp[2])

            tmp = array.array('i')
            tmp.fromfile(f,3)

            sn = array.array('f')
            sn.fromfile(f,tmp[2])

            tmp = array.array('i')
            tmp.fromfile(f,3)

            bh = array.array('f')
            bh.fromfile(f,tmp[2])

            tmp = array.array('i')
            tmp.fromfile(f,3)

            wd = array.array('f')
            wd.fromfile(f,tmp[2])

            tmp = array.array('i')
            tmp.fromfile(f,3)

            rmtm = array.array('f')
            rmtm.fromfile(f,tmp[2])
            rmtm = numpy.array(rmtm)

            data = [[ta,ids,iw,wave,sed,strm,rmtm],
                    ['Ages(yr)', 'SSP ID', 
                     'No. wavelength entries',
                     'Wavelength array',
                     'SEDs', 'Stellar mass history',
                     'Remnant mass history']]

    else:       
        with open(filename,'rb') as f:


            ks = array.array('i')
            ks.fromfile(f,3)

            ta = array.array('f')
            ta.fromfile(f,ks[2])
            ta = numpy.array(ta)


            tmp = array.array('i')
            tmp.fromfile(f,3)
            ml,mul,iseg = tmp

            if iseg > 0:
                tmp = array.array('f')
                tmp.fromfile(f,iseg*6)

            tmp = array.array('f')
            tmp.fromfile(f,5)
            totm, totn, avs, jo, tauo = tmp



            ids= array.array('c')
            ids.fromfile(f,80)


            tmp = array.array('f')
            tmp.fromfile(f,4)
            tcut = tmp[0]
            ttt = tmp[1:]

            ids = array.array('c')
            ids.fromfile(f,80)
 
            ids = array.array('c')
            ids.fromfile(f,80)
 
            igw = array.array('i')
            igw.fromfile(f,1)
 
            tmp = array.array('i')
            tmp.fromfile(f,5)
 

            iw = array.array('i')
            iw.fromfile(f,1)

            wave = array.array('f')
            wave.fromfile(f,iw[0])
            wave = numpy.array(wave)

            #SED Section
            F = array.array('i')
            F.fromfile(f,5)
            iw = F[4] #Number of wavelength elements


            sed = numpy.zeros((iw,ks[2]),dtype=numpy.float32)
            G = array.array('f')
            G.fromfile(f,iw)
            sed[:,0] = G
            ik = array.array('i')
            ik.fromfile(f,1)

            h = numpy.empty((ik[0],ks[2]),'f')
            H = array.array('f')
            H.fromfile(f,ik[0])
            h[:,0] = H

            for i in range(1,ks[2]): #Fill rest of array with SEDs
                F = array.array('i')
                F.fromfile(f,5)
                iw = F[4]

                G = array.array('f')
                G.fromfile(f,iw)
                sed[:,i] = G
                ik = array.array('i')
                ik.fromfile(f,1)

                H = array.array('f')
                H.fromfile(f,ik[0])
                h[:,i] = H

            tmp = array.array('i')
            tmp.fromfile(f,5)


            bflx = array.array('f')
            bflx.fromfile(f,tmp[4])
 
            tmp = array.array('i')
            tmp.fromfile(f,5)

            strm = array.array('f')
            strm.fromfile(f,tmp[4])
            strm = numpy.array(strm)

            tmp = array.array('i')
            tmp.fromfile(f,5)

            evf = array.array('f')
            evf.fromfile(f,tmp[4])

            tmp = array.array('i')
            tmp.fromfile(f,5)

            evf = array.array('f')
            evf.fromfile(f,tmp[4])

            tmp = array.array('i')
            tmp.fromfile(f,5)

            snr = array.array('f')
            snr.fromfile(f,tmp[4])

            tmp = array.array('i')
            tmp.fromfile(f,5)

            pnr = array.array('f')
            pnr.fromfile(f,tmp[4])

            tmp = array.array('i')
            tmp.fromfile(f,5)

            sn = array.array('f')
            sn.fromfile(f,tmp[4])

            tmp = array.array('i')
            tmp.fromfile(f,5)

            bh = array.array('f')
            bh.fromfile(f,tmp[4])

            tmp = array.array('i')
            tmp.fromfile(f,5)

            wd = array.array('f')
            wd.fromfile(f,tmp[4])

            tmp = array.array('i')
            tmp.fromfile(f,5)

            rmtm = array.array('f')
            rmtm.fromfile(f,tmp[4])
            rmtm = numpy.array(rmtm)

            data = [[ta,ids,iw,wave,sed,strm,rmtm],
                    ['Ages(yr)', 'SSP ID', 
                     'No. wavelength entries',
                     'Wavelength array',
                     'SEDs', 'Stellar mass history',
                     'Remnant mass history']]

    return data


def tl(z,OMEGA_M0=0.3,OMEGA_L=0.7,H0=70,N=1000): 
     """ Calculates the lookback time in Gyr to redshift z for the current set of cosmological 
      parameters. 
       
      @type z: float 
      @param z: redshift   
      @rtype: float 
      @return: lookback time in Gyr to redshift z 
       
      """ 
     OMEGA_K=1.0-OMEGA_M0-OMEGA_L 
     
     x=numpy.linspace((1.0/(1+z)),1.0,N)
     y=x/numpy.sqrt(OMEGA_M0*x+OMEGA_L*numpy.power(x, 4)+OMEGA_K*numpy.power(x, 2))  
     integralValue=simps(y,x)
     T0=(1.0/H0*integralValue*3.08568025e19)/3.15569e7/1e9 

     return numpy.nan_to_num(T0)

def t0(OMEGA_M0=0.3,OMEGA_L=0.7,H0=70,N=1000): 
     """ Calculates the lookback time in Gyr to redshift z for the current set of cosmological 
      parameters. 
       
      @type z: float 
      @param z: redshift   
      @rtype: float 
      @return: lookback time in Gyr to redshift z 
       
      """ 
     OMEGA_K=1.0-OMEGA_M0-OMEGA_L 
     
     x=numpy.linspace(0.0,1.0,N)[1:]
     y=x/numpy.sqrt(OMEGA_M0*x+OMEGA_L*numpy.power(x, 4)+OMEGA_K*numpy.power(x, 2))  
     integralValue=simps(y,x)
     T0=(1.0/H0*integralValue*3.08568025e19)/3.15569e7/1e9 

     return T0  

def dm(z,OMEGA_M0=0.3,OMEGA_L=0.7,H0=70,N=1000,C_LIGHT=299792.458):

    """Calculates the transverse comoving distance (proper motion distance) in Mpc at 
      redshift z. 
       
      @type z: float 
      @param z: redshift  
      @rtype: float 
      @return: transverse comoving distance (proper motion distance) in Mpc 
       
      """ 
       
    OMEGA_K=1.0-OMEGA_M0-OMEGA_L 
 
    x=numpy.linspace((1.0/(1+z)),1.0,N)
    y=1.0/numpy.sqrt(OMEGA_M0*x+OMEGA_L*numpy.power(x, 4)+OMEGA_K*numpy.power(x, 2)) 
                
    integralValue=simps(y,x)

    if OMEGA_K>0.0: 
        DM=C_LIGHT/H0*math.pow(abs(OMEGA_K),-0.5*math.sinh(math.sqrt(abs(OMEGA_K))*integralValue))
    elif OMEGA_K==0.0: 
        DM=C_LIGHT/H0*integralValue 
    elif OMEGA_K<0.0: 
        DM=C_LIGHT/H0*math.pow(abs(OMEGA_K),-0.5*math.sin(math.sqrt(abs(OMEGA_K))*integralValue))
       
    return DM

def dist(z):
    """
    Calculates the distance modulus for a given z using the proper distance dm(z)

    """
    return 5*numpy.log10(dm(z)*(1+z)*1e6/10)

def dV(z,z0=0.0001,OMEGA_M0=0.3,OMEGA_L=0.7,H0=70,N=1000,C_LIGHT=299792.458):
    """
    Calculates the comoving volume per solid angle at z, Mpc^3 per steradian

    """

    x = numpy.linspace(z0,z,N)
    y = numpy.array(map(dm,x))**2/numpy.array(map(Ez,x))
    
    integralValue = simps(y,x)
    dV = (C_LIGHT/H0)*integralValue
    return dV

def dVz(z,OMEGA_M0=0.3,OMEGA_L=0.7,H0=70,N=1000,C_LIGHT=299792.458):
    """
    Calculates the comoving volume per solid angle per z, Mpc^3 per steradian

    """

    return dm(z)**2/Ez(z)*(C_LIGHT/H0)

def Ez(z,OMEGA_M0=0.3,OMEGA_L=0.7,H0=70):
    Ezz = math.sqrt(OMEGA_M0*math.pow(1.0+z, 3)+(1.0-OMEGA_M0-OMEGA_L)*math.pow(1.0+z, 2)+OMEGA_L)
    return Ezz


def dec_a_func(wave_obs):
    return numpy.exp(-1*0.0036*(numpy.power(wave_obs/1216.,3.46)))

def dec_b_func(wave_obs):
    teff_beta=1.7e-3*(numpy.power(wave_obs/1026.,3.46))
    teff_gamma=1.2e-3*(numpy.power(wave_obs/973.,3.46))
    teff_delta=9.3e-4*(numpy.power(wave_obs/950.,3.46))
    #The above absorption lines dominate the average absorption over the range 
    #but higher order lines should be added in future
    
    teff_total=teff_beta+teff_gamma+teff_delta

    return numpy.exp(-1*teff_total)

def calc_lyman(x,s):
    wly = 912.
    const_c = 2.997925e10
    const_h = 6.6262e-27
    const = 1e-8/const_h/const_c

    n = int(sum([x < wly][0]))
    f = numpy.zeros(n+1)
    w = numpy.zeros(n+1)

    for i in range(n+1):
        if x[i] < wly:
            w[i] = x[i]
            f[i] = w[i]*s[i]
        elif x[i] == wly:
            w[i] = x[i]
            f[i] = w[i]*s[i]
        elif x[i] > wly:
            w[i] = wly
            f[i] = s[i-1] + ((w[i]-x[i-1])*(s[i]-s[i-1])/(x[i]-x[i-1]))
            f[i] = w[i]*f[i]

    nlyman = const*trapz(f,w)
    #print numpy.log10(N_lyman)
    return nlyman

def calc_beta(wave, SED):
    """
    wave = wavelength array
    SED = Rest-frame flux

    Returns UV slope index and the error on that fit
    """
    new_wave = numpy.arange(1200,2650)
    new_SED = griddata(wave,SED,new_wave)
    wave = new_wave
    SED = new_SED

    window_lower = numpy.array([1268.,1309.,1342.,1407.,1562.,1677.,1760.,1866.,1930.])#,2400.])
    window_upper = numpy.array([1284.,1316.,1371.,1515.,1583.,1740.,1833.,1890.,1950.])#,2580.])

    window_mean = (window_lower+window_upper)/2 #midpoint for power-law fitting

    fluxes = numpy.zeros_like(window_lower)

    for window in range(len(window_lower)):
        wf = numpy.arange(window_lower[window]-5,window_upper[window]+5)
        tp = numpy.zeros(len(wf))
        tp[(wf>=window_lower[window]) & (wf<window_upper[window])] = 1.0

        #Find SED wavelength entries within filter range
        wff = numpy.array([wf[0] < wave[i] < wf[-1] for i in range(len(wave))])
        wft = wave[wff]

        #Interpolate to find throughput values at new wavelength points
        tpt = griddata(wf,tp,wft)

        #Join arrays and sort w.r.t to wf
        wf = numpy.concatenate((wf,wft))
        tp = numpy.concatenate((tp,tpt))

        order = numpy.argsort(wf)
        wf = wf[order]
        tp = tp[order]
        
        dwf = numpy.diff(wf)
        nwf = len(wf)

        tpwf = tp/wf
        
        f_mean2 = numpy.dot(dwf,(tpwf[:nwf-1]+tpwf[1:])/2)
        tpwf = tp*wf

        WR = 0.

        for i in range(nwf):

            #Interpolation indices
            j = numpy.where(wave<wf[i])[0][-1]

            a = (wf[i] - wave[j])/(wave[j+1]-wave[j])
            tpa = tpwf[i]*((1-a)*(SED[j]) + a*SED[j+1])

            if i != 0:
                WR += dwf[i-1]*(tpb+tpa)

            tpb = tpa

        fluxes[window] = WR/2/f_mean2/2.997925e18

    fluxes *= 2.997925e18/(window_mean**2)
    fluxerr = numpy.sqrt(fluxes)

    logx = numpy.log10(window_mean/1e10)
    logy = numpy.log10(fluxes)
    logyerr = fluxerr/fluxes
    
    fitfunc = lambda p, x: p[0] + p[1]*x
    errfunc = lambda p, x, y, err: (y - fitfunc(p,x))/err
    
    pinit = [numpy.max(fluxes), -2.0]
    out = leastsq(errfunc, pinit, args=(logx,logy,logyerr))

    pfinal = out[0]
    covar = out[1]
    
    index = pfinal[1]
    #indexerr = numpy.sqrt(covar[0])
    
    return (index)#, indexerr)


def zmax_func(z,limitingMagnitude,absoluteMagnitude):
    return dist(z)+absoluteMagnitude-limitingMagnitude

def zmax_find(z0,observedRedshift,observedMagnitude,limitingMagnitude):
    absoluteMagnitude=observedMagnitude-dist(observedRedshift)
    findit = fsolve(zmax_func,z0,args=(limitingMagnitude,absoluteMagnitude))
    return findit[0]
