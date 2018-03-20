import numpy as np
import array
import os, sys
import re
import time
import multiprocessing
import h5py
import logging
from astropy.table import Table, Column
from astropy import units as u


import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-p","--params", type=str,
                    help = "Parameter file")
parser.add_argument("-q", "--quiet", help = "Suppress extra outputs",
                    action = "store_true")
args = parser.parse_args()
quiet = args.quiet

params_root = re.split(".py", args.params)[0]
if os.path.isfile(params_root+".pyc"):
    os.remove(params_root+".pyc")

import importlib
try:
    params = importlib.import_module(params_root)
    print('Successfully loaded "{0}" as params'.format(args.params))
    reload(params)
except:
    print('Failed to load "{0}" as params'.format(args.params))
    raise

if quiet:
    quietprint = lambda *a: None
else:
    def quietprint(*args):
        for arg in args:
            print arg,
        print

# Fitting function definition for later use by Processess

def galaxyFit(inputQueue, printQueue, printlock):
    for gal in iter(inputQueue.get, 'STOP'):
        j = np.argmin(np.abs(z-zobs[gal])) # Find closest model redshift


        flux_obs = obs[gal,:]
        flux_err = obs_err[gal,:]
        
        #flux_obs[fo <= 0.] = 0.       # Set negative fluxes to zero
        I = np.where(flux_err > 0.)[0] # Find bands with no observation
        
        if len(I) == 0:
            if include_rest:
                M_scaled = np.ones(len(fo)) * -99.
                restframe_output = ' '.join(M_scaled.astype('str'))
                output_string = '{0} {1} {2} {3} {4} {5} {6} {7}' \
                                ' {8} {9} {10} {11} {12} {13} {14} {15} {16}'.format(gal+1,ID[gal],zobs[gal],-99,-99,-99,-99,-99,-99, -99, -99,-99,len(I),-99,z[j],restframe_output,'\n')
            else:
                output_string = '{0} {1} {2} {3} {4} {5} {6} {7} {8} {9} {10} {11} {12} {13} {14}'.format(gal+1,ID[gal],zobs[gal],-99,-99,-99,-99,-99,-99,-99, -99,-99,len(I),-99,'\n')
            printQueue.put(output_string)
            continue
            
        flux_obs = flux_obs[I]                    # and exclude from fit
        flux_err = flux_err[I]
        flux_models = f[j,I,:]

        tot_err = np.sqrt(flux_err**2 + (0.1*flux_obs)**2)
        
        top = 0.
        bottom = 0.
    
        for i in range(len(flux_obs)):
            top += (flux_models[i,:]*flux_obs[i])/(tot_err[i]**2)
            bottom += (flux_models[i,:]**2)/(tot_err[i]**2)
    
        scale = top/bottom
        scale = np.reshape(scale, (n_metal, n_tg, n_tau, n_tauv, n_fesc))  

        chisq = 0.
        for i in range(len(flux_obs)):
            chisq += ((np.abs(scale*flux_models[i,:]-flux_obs[i])**2)/(flux_err[i])**2)

        chimin, minind = np.nanmin(chisq), np.nanargmin(chisq)
        
        if np.isinf(chimin) or np.isnan(minind):
            if include_rest:
                M_scaled = np.ones(len(flux_obs)) * -99.
                restframe_output = ' '.join(M_scaled.astype('str'))
                output_string = '{0} {1} {2} {3} {4} {5} {6} {7} {8} {9} {10} {11} {12} {13} {14} {15} {16}'.format(gal+1,ID[gal],zobs[gal],-99,-99,-99,-99,-99,-99, -99, -99,-99,len(I),-99,z[j],restframe_output,'\n')
            else:
                output_string = '{0} {1} {2} {3} {4} {5} {6} {7} {8} {9} {10} {11} {12} {13} {14}'.format(gal+1,ID[gal],zobs[gal],-99,-99,-99,-99,-99,-99,-99, -99,-99,len(I),-99,'\n')
            printQueue.put(output_string)
            continue


        #Find the coordinate of the model with the bestfit mass
        mi, tgi, ti, tvi, fi = np.unravel_index(minind, 
                                                   (n_metal, n_tg, 
                                                   n_tau, n_tauv, n_fesc)) 

        Bestfit_Mass = np.log10(scale[mi, tgi, ti, tvi, fi]*flux_corr)
        Bestfit_SFR = (scale[mi, tgi, ti, tvi, fi] * 
                       SFR[mi, tgi, ti, tvi, fi]*flux_corr)
        #Bestfit_Beta = beta[tgi,tvi,ti,mi]
        Bestfit_Beta = -99.

        #Scale the observed tot_mag band of the template to be the same as the observed tot_mag band of the galaxy
        #Convert the templates so they are no longer units of per stellar mass

        F_rest = f[0,:]*scale[mi, tgi, ti, tvi, fi]*flux_corr
        restframeMags = 23.9 - 2.5*np.log10(F_rest)
    
        #UV_rest = UV_flux[0]*scale[tgi,tvi,ti,mi]*flux_corr
        #restframeMUV = 23.9 - 2.5*np.log10(UV_rest)

        M_scaled = restframeMags[:, mi, tgi, ti, tvi, fi]
        #MUV_scaled = restframeMUV[tgi,tvi,ti,mi]
        MUV_scaled = -99.
        
        if np.isnan(Bestfit_Mass) or np.isinf(chimin):
            Bestfit_Mass = -99
            #M_scaled[:] = -99
            tgs = -99
            tvs = -99
            taus = -99
            mis = -99
            escape_fraction = -99

        else:
            tgs = tg[tgi]/1e9
            tvs = tv[tvi]
            taus = tau[ti]
            mis = metallicities[mi]
            escape_fraction = fesc[fi] 


        printlock.acquire()

        print '{0:6d} {1:8d} {2:>5.2f} {3:>7.2f} {4:>8.1f} {5:>8.3f} {6:>5.1f} {7:>8.2f} {8:>4.2f} {9:>5.2f}'.format(gal+1,ID[gal], zobs[gal],Bestfit_Mass,chimin,tgs,tvs,taus,mis,np.log10(Bestfit_SFR))

        if include_rest:
            restframe_output = ' '.join(M_scaled.astype('str'))
            output_string = '{0} {1} {2} {3} {4} {5} {6} {7} {8} {9} {10} {11} {12} {13} {14} {15} {16}'.format(gal+1,ID[gal],zobs[gal],Bestfit_Mass,chimin,tgs,tvs,taus,mis, MUV_scaled, minind,Bestfit_SFR,len(I),Bestfit_Beta,z[j],restframe_output,'\n')
        else:
            output_string = '{0} {1} {2} {3} {4} {5} {6} {7} {8} {9} {10} {11} {12} {13} {14}'.format(gal+1,ID[gal],zobs[gal],Bestfit_Mass,chimin,tgs,tvs,taus,mis, MUV_scaled, minind,Bestfit_SFR,len(I),Bestfit_Beta,'\n')

        printlock.release()
        printQueue.put(output_string)

def galaxyFit2(inputQueue, printQueue, printlock):
    for gal in iter(inputQueue.get, 'STOP'):
        
        
        output_string = '{0[0]} {0[1]} {0[2]} {0[3]} {0[4]} {0[5]} ' + \
                        '{0[6]} {0[7]} {0[8]} {0[9]} {0[10]} {0[11]} ' + \
                        '{0[12]} {0[13]} {0[14]}'
        
        j = np.argmin(np.abs(z-zobs[gal])) # Find closest model redshift

        log_mass_min, log_mass_max = 7, 13
        log_sfr_min, log_sfr_max = -3, 4

        flux_obs = obs[gal,:]
        flux_err = obs_err[gal,:]
        
        #flux_obs[fo <= 0.] = 0.       # Set negative fluxes to zero
        I = np.where(flux_err > 0.)[0] # Find bands with no observation
        
        if len(I) == 0:
            output_array = [gal+1, ID[gal], zobs[gal], z[j],
                            -99, -99, -99, -99, -99, -99, -99,
                            -99,-99,len(I),-99,'\n']
            output = output_string.format(output_array)
            
            if include_rest:
                M_scaled = np.ones(len(flux_obs)) * -99.
                restframe_output = ' '.join(M_scaled.astype('str'))
                output = output + restframe_output + ' \n'
                
            else:
                output = output + ' \n'
            printQueue.put(output_string)
            continue
            
        flux_obs = flux_obs[I]                    # and exclude from fit
        flux_err = flux_err[I]
        flux_models = f[j,I,:]

        tot_err = np.sqrt(flux_err**2 + (params.flux_err*flux_obs)**2)

        top = 0.
        bottom = 0.
    
        for i in range(len(flux_obs)):
            top += (flux_models[i,:]*flux_obs[i])/(tot_err[i]**2)
            bottom += (flux_models[i,:]**2)/(tot_err[i]**2)
    
        scale = top/bottom
        scale = np.reshape(scale, (n_metal, n_tg, n_tau, n_tauv, n_fesc))  

        chisq = 0.
        for i in range(len(flux_obs)):
            chisq += ((np.abs(scale*flux_models[i,:]-flux_obs[i])**2)/(tot_err[i])**2)

        chimin, minind = np.nanmin(chisq), np.nanargmin(chisq)
        likelihood = np.reshape(np.exp(-0.5*chisq), 
                                (n_metal, n_tg, n_tau, n_tauv, n_fesc))
        likelihood[np.isnan(likelihood)] = 0.
        likelihood = np.abs(likelihood/likelihood.sum())
        
        
        if np.isinf(chimin) or np.isnan(minind):
            output_array = [gal+1, ID[gal], zobs[gal], z[j],
                            -99, -99, -99, -99, -99, -99, -99,
                            -99,-99,len(I),-99,'\n']
            output = output_string.format(output_array)

        else:
            #Find the coordinate of the model with the bestfit mass
            mi, tgi, ti, tvi, fi = np.unravel_index(minind, 
                                                       (n_metal, n_tg, 
                                                       n_tau, n_tauv, n_fesc)) 

            Masses = np.abs(np.log10(scale*flux_corr))
            SFRs = np.log10(scale * SFR * flux_corr)
        
            mass_hist = np.histogram(Masses.flatten(),
                                     range = (log_mass_min, log_mass_max),
                                     bins = 120,
                                     weights = likelihood.flatten(),
                                     density = True)
                                 
            sfr_hist = np.histogram(SFRs.flatten(),
                                     range = (log_sfr_min, log_sfr_max),
                                     bins = 140,
                                     weights = likelihood.flatten(),
                                     density = True)
        
            Bestfit_Mass = np.abs(np.log10(scale[mi, tgi, ti, tvi, fi]*flux_corr))
            Bestfit_SFR = np.abs(np.log10(scale[mi, tgi, ti, tvi, fi] * 
                                   SFR[mi, tgi, ti, tvi, fi]*flux_corr))
        
            if np.isnan(Bestfit_Mass) or np.isinf(chimin):
                Bestfit_Mass = -99
                #M_scaled[:] = -99
                tgs = -99
                tvs = -99
                taus = -99
                mis = -99
                escape_fraction = -99

            else:
                tgs = tg[tgi]/1e9
                tvs = tv[tvi]
                taus = tau[ti]
                mis = metallicities[mi]
                escape_fraction = fesc[fi] 

            m16, m50, m84 = weighted_quantile(Masses.flatten(), 
                                              [0.16, 0.5, 0.84],
                                              sample_weight=likelihood.flatten(),
                                              values_sorted=False)
            s16, s50, s84 = weighted_quantile(SFRs.flatten(), 
                                              [0.16, 0.5, 0.84],
                                              sample_weight=likelihood.flatten(),
                                              values_sorted=False)

            printlock.acquire()

            MUV_scaled = -99. 
            Bestfit_Beta = -99.
        
            print_string = "{0[0]:6d} {0[1]:8d} {0[2]:>5.2f} " + \
                           "{0[3]:>7.2f} {0[4]:>8.1f} {0[5]:>8.3f} " + \
                           "{0[6]:>5.1f} {0[7]:>8.2f} {0[8]:>4.2f} " + \
                           "{0[9]:>5.2f}"
                       
            print_array = [gal+1, ID[gal], zobs[gal],
                           Bestfit_Mass, chimin, 
                           tgs, tvs, taus, mis, 
                           Bestfit_SFR]                  
            print(print_string.format(print_array))

            output_string = '{n} {id} {zobs} {ztemp} {mass_best} {sfr_best} '+ \
                            '{chi_best} {tg} {tvs} {taus} {mis} {fesc} '+ \
                            '{mass_med} {mass_l68} {mass_u68} ' + \
                            '{sfr_med} {sfr_l68} {sfr_u68} ' + \
                            '{nfilts} '
                        
            output_values = {'n': gal+1,
                             'id': ID[gal],
                             'zobs': zobs[gal], 'ztemp':z[j],
                             'mass_best': Bestfit_Mass,
                             'sfr_best': Bestfit_SFR,
                             'chi_best': chimin,
                             'tg': tgs, 'tvs': tvs, 'taus': taus, 
                             'mis': mis, 'fesc': escape_fraction,
                             'mass_med': m50, 'mass_l68': m16, 'mass_u68': m84,
                             'sfr_med': s50, 'sfr_l68': s16, 'sfr_u68': s84,
                             'nfilts': len(I)}

            output_array = [gal+1, ID[gal], zobs[gal], 
                            Bestfit_Mass, chimin, tgs, tvs, taus, mis, 
                            MUV_scaled, minind, Bestfit_SFR, len(I), -99., '\n']
            output = output_string.format(**output_values)
        
        if include_rest:
            if np.isinf(chimin) or np.isnan(minind):
                M_scaled = np.ones(len(flux_obs)) * -99.
                restframe_output = ' '.join(M_scaled.astype('str'))
                output = output + restframe_output + ' \n'
            
            else:
                F_rest = np.array(f[0, :, mi, tgi, ti, tvi, fi] * 
                                  scale[mi, tgi, ti, tvi, fi] * flux_corr)
                restframeMags = 23.9 - 2.5*np.log10(F_rest)
                restframe_output = ' '.join(restframeMags.astype('str'))
                output = output + restframe_output + ' \n'
        else:
            output = output + ' \n'
            
        printlock.release()
        printQueue.put([output, mass_hist, sfr_hist])

def galaxyFitPlus(inputQueue, printQueue, printlock):
    for gal in iter(inputQueue.get, 'STOP'):
        mass_range = 7, 13
        log_sfr_min, log_sfr_max = -3, 4
        
        j = np.argmin(np.abs(z-zobs[gal])) # Find closest model redshift

        fo = obs[gal,:]
        ferr = obs_err[gal,:]


        flux_obs[fo <= 0.] = 0.       # Set negative fluxes to zero
        #print fo
        I = (ferr > 0.)*(ferr < 1e6) # Find bands with no observation
        fo = flux_obs[I]                    # and exclude from fit
        ferr = flux_err[I]
        fm = f[I,j,:]
        #print flux_models[:,0,0,0,0]        

        top = 0.
        bottom = 0.
    
        for i in range(len(fo)):
            top += (flux_models[i,:]*flux_obs[i])/(flux_err[i]**2)
            bottom += (flux_models[i,:]**2)/(flux_err[i]**2)
    
        scale = top/bottom
        scale = np.reshape(scale, (n_metal, n_tg, n_tau, n_tauv, n_fesc))  

        chisq = 0.
        for i in range(len(fo)):
            chisq += ((np.abs(scale*flux_models[i,:]-flux_obs[i])**2)/(flux_err[i])**2)

        chimin, minind = np.nanmin(chisq), np.nanargmin(chisq)

        chisq -= (chisq.min() - 1)
        likelihood = np.exp(-0.5*chisq)
        likelihood /= likelihood.sum()
        
        if np.isinf(chimin) or np.isnan(minind) or len(fo) == 0:
            output_string = '{0} {1} {2} {3} {4} {5} {6} {7} {8} {9} \
            {10} {11} {12} {13} {14} {15} {16} {17} {18}'.format(gal+1,ID[gal],zobs[gal],
                                                                 -99,-99,-99,-99,-99,-99,
                                                                 -99, -99, -99, -99,-99,-99,-99,
                                                                 len(I),-99,'\n')
            
            massLikelihood = np.zeros(mass_bins+1)
            massLikelihood[0] = gal
            muvLikelihood = np.zeros(muv_bins+1)
            muvLikelihood[0] = gal
            betaLikelihood = np.zeros(beta_bins+1)
            betaLikelihood[0] = gal
            #tauLikelihood = np.zeros(n_tau)        
            #tauLikelihood = np.insert(tauLikelihood,0,gal)        
            printQueue.put([output_string,massLikelihood,muvLikelihood,betaLikelihood])
            continue

        #Find the coordinate of the model with the bestfit mass
        si,tgi,tvi,ti,mi = np.unravel_index(minind,(mass_bins,n_tg,n_tauv,n_tau,n_ssp))
        Bestfit_Mass = np.log10(mass_range[si]*flux_corr)
        Bestfit_SFR = (mass_range[si]*SFR[tgi,ti,mi]*flux_corr)
        Bestfit_Beta = beta[tgi,tvi,ti,mi]

        F_rest = f[:,0]*mass_range[likelihood.argmax(0)]*flux_corr
        restframeMags = 23.9 - 2.5*np.log10(F_rest)
    
        UV_rest = UV_flux[0]*mass_range[likelihood.argmax(0)]*flux_corr
        restframeMUV = 23.9 - 2.5*np.log10(UV_rest)

        Bestfit_restframeMags = restframeMags[:,tgi,tvi,ti,mi]
        Bestfit_restframeMUV = restframeMUV[tgi,tvi,ti,mi]

        if np.isnan(Bestfit_Mass) or np.isinf(chimin):
            Bestfit_Mass = -99
            #M_scaled[:] = -99
            tgs = -99
            tvs = -99
            taus = -99
            mis = -99

        else:
            tgs = tg[tgi]/1.e9
            tvs = tv[tvi]
            taus = tau[ti]/1.e9
            mis = mi
            
        """
        Likelihood array section:
        """
        mass_hist = np.histogram(np.log10(mass_))
        
        printlock.acquire()

        if calc_mode:
            print '{0:4d} {1:6d} {2:>6.2f} {3:>8.1f} {4:>6.2f}'.format(gal+1,ID[gal],Bestfit_Mass,chimin, np.log10(Mode_Mass), '/n')
        else:
            print '{0:6d} {1:8f} {2:>5.2f} {3:>7.2f} {4:>8.1f} {5:>8.3f} {6:>5.1f} {7:>8.2f} {8:>3d} {9:>5.2f}'.format(gal+1,int(ID[gal]),zobs[gal],Bestfit_Mass,chimin,tgs,tvs,taus,mis,np.log10(Bestfit_SFR))

        output_string = '{0} {1} {2} {3} {4} {5} {6} {7} {8} {9} {10} {11} {12} {13} {14} {15}'.format(gal+1,int(ID[gal]),zobs[gal],Bestfit_Mass,chimin,tgs,tvs,taus,mis,Bestfit_restframeMags[tot],Bestfit_restframeMUV,minind,Bestfit_SFR,len(I),Bestfit_Beta,'\n')

        printlock.release()
        printQueue.put([output_string, massLikelihoods, muvLikelihoods, betaLikelihoods])


def getObservations(inputpath):
    input_data = Table.read(inputpath,format=input_format)

    column_names = input_data.columns.keys()

    ID = input_data[ID_col]
    zobs = input_data[z_col]

    filter_names = []

    k,l = 0,0
    for ii in range(len(column_names)):
        if column_names[ii].lower().endswith(flux_col_end.lower()):
            if k == 0:
                fluxes = input_data[column_names[ii]]
            else:
                fluxes = np.column_stack((fluxes,input_data[column_names[ii]]))
            k+=1
            filter_names.append(column_names[ii])

        if column_names[ii].lower().endswith(fluxerr_col_end.lower()):
            if l == 0:
                fluxerrs = input_data[column_names[ii]]
            else:
                fluxerrs = np.column_stack((fluxerrs,input_data[column_names[ii]]))
            l+=1
    """        
    if filts_used != None:
        try:
            fluxes = fluxes[:,filts_used]
            fluxerrs = fluxerrs[:,filts_used]
        except:r
            print('Filter mismatch 1')
            # Array slicing fail
    """
    return ID, zobs, fluxes, fluxerrs, k, filter_names

class _function_wrapper(object):
    """
    This is a hack to make the likelihood function pickleable when ``args``
    or ``kwargs`` are also included.
    
    Stolen from emcee
    """
    def __init__(self, f, args, kwargs):
        self.f = f
        self.args = args
        self.kwargs = kwargs

    def __call__(self, x):
        try:
            return self.f(x, *self.args, **self.kwargs)
        except:
            import traceback
            print("emcee: Exception while calling your likelihood function:")
            print("  params:", x)
            print("  args:", self.args)
            print("  kwargs:", self.kwargs)
            print("  exception:")
            traceback.print_exc()
            raise

def weighted_quantile(values, quantiles, sample_weight=None, values_sorted=False, old_style=False):
    """ Very close to np.percentile, but supports weights.
    NOTE: quantiles should be in [0, 1]!
    :param values: np.array with data
    :param quantiles: array-like with many quantiles needed
    :param sample_weight: array-like of the same length as `array`
    :param values_sorted: bool, if True, then will avoid sorting of initial array
    :param old_style: if True, will correct output to be consistent with np.percentile.
    :return: np.array with computed quantiles.
    """
    values = np.array(values)
    quantiles = np.array(quantiles)
    if sample_weight is None:
        sample_weight = np.ones(len(values))
    sample_weight = np.array(sample_weight)
    assert np.all(quantiles >= 0) and np.all(quantiles <= 1), 'quantiles should be in [0, 1]'

    if not values_sorted:
        sorter = np.argsort(values)
        values = values[sorter]
        sample_weight = sample_weight[sorter]

    weighted_quantiles = np.cumsum(sample_weight) - 0.5 * sample_weight
    if old_style:
        # To be convenient with np.percentile
        weighted_quantiles -= weighted_quantiles[0]
        weighted_quantiles /= weighted_quantiles[-1]
    else:
        weighted_quantiles /= np.sum(sample_weight)
    return np.interp(quantiles, weighted_quantiles, values)

if __name__ == '__main__':
    
    logfile = open("error.log", "w")
    original_stderr = sys.stderr
    sys.stderr = logfile
    
    start = time.time()
    
    """
    SECTION 1

    """
    model_path = params.model_path

    input_catalog = params.input_catalog
    input_format = params.input_format
    z_col = params.z_col
    ID_col = params.ID_col
    flux_col_end = params.flux_col_end
    fluxerr_col_end = params.fluxerr_col_end

    ncpus = params.ncpus
    filts_used = params.filts_used
    include_rest = params.include_rest
    
    output_path = params.output_catalog_path
    output_format = params.output_format
    output_hdf_path = params.output_hdf_path
    
    calc_mode = params.fitting_mode
    flux_corr = params.flux_corr 

    
    ID, zobs, obs, obs_err, filters_found, filter_names = getObservations(input_catalog)
    
    """
    Section 2
    
    """


    print "Loading synthetic mags and mass array:"
    models = h5py.File(model_path, 'r')
    tg = models['ages'].value
    tv = models['dust'].value
    tau = models['sfh'].value
    metallicities = models['metallicities'].value
    fesc = models['fesc'].value

    Mshape = models['fluxes'].shape
    z = models['z']
    nfilts = Mshape[1]
    n_metal = Mshape[2]        
    n_tg = Mshape[3]
    n_tau = Mshape[4]
    n_tauv = Mshape[5]
    n_fesc = Mshape[6]
    
    #UV_flux = synmags['UV_flux']
    SFR = models['SFR']
    Ms = models['Ms']

    if (nfilts == filters_found) and (filts_used == None):
        f = models['fluxes']
    
    elif filts_used != None:
        try:        
            f = models['fluxes'][:,filts_used]
            obs = obs[:,filts_used]
            obs_err = obs_err[:,filts_used]
            filter_names = np.array(filter_names)[filts_used]
        except:
            print('Mis-match between model and observed filter numbers')
            raise
            # Slice fail

    print ("Done.")

    """
    SECTION 3
    """
    if os.path.isfile(output_path+".temp_output.txt"):
        os.remove(output_path+".temp_output.txt")
    temp_file = open(output_path+".temp_output.txt","w")
    
    """
    SECTION 4
    Chi-sq calculation

    """
    out_string = '{0:6s} {1:8s} {2:>5s} {3:>7s} {4:>8s} {5:>8s}' + \
                 '{6:>5s} {7:>8s} {8:>4s} {9:>5s}'

    print(out_string.format('N','ID','zobs','Best', 'chimin', 
                            'tg', 'tauv','tau','met', 'sfr'))

    loop_start = time.time()
    ncpus = np.clip(ncpus, 1, multiprocessing.cpu_count())

    inputQueue = multiprocessing.Queue()
    printQueue = multiprocessing.Queue()
    printlock = multiprocessing.Lock()
    
    if calc_mode == 'hist':
        output_hdf = h5py.File(output_hdf_path, 'w')
        output_hdf.create_dataset("mass_pdf", (len(ID), 120), dtype="f")
        output_hdf.create_dataset("sfr_pdf", (len(ID), 140), dtype="f")
        fitFunction = galaxyFit2
        
    else:
        fitFunction = galaxyFit
    
    for i in range( ncpus ):
        multiprocessing.Process(target = fitFunction,
                                args = (inputQueue, printQueue,
                                        printlock)).start()

    # Put elements in the send queue for processing
    for gal in range( len(ID) ):
        inputQueue.put( gal )

    if calc_mode == 'hist':
        for i, gal in enumerate(ID):
            printout, mass_hist, sfr_hist = printQueue.get()
            if i == 0:
                mass_centers = 0.5*(mass_hist[1][1:] + mass_hist[1][:-1])
                sfr_centers = 0.5*(sfr_hist[1][1:] + sfr_hist[1][:-1])
                
                output_hdf.create_dataset("mass_bins", data = mass_centers)
                output_hdf.create_dataset("sfr_bins", data = sfr_centers)            
            output_hdf["mass_pdf"][i] = mass_hist[0]
            output_hdf["sfr_pdf"][i] = sfr_hist[0]
            temp_file.write( printout )
        #tau_array.tofile(tau_file)
    else:
        for i, gal in enumerate(ID):
            printout = printQueue.get()
            temp_file.write( printout )
            #print len(mass_array), len(muv_array), len(beta_array)
       

    # Stop all the running processes
    for i in range( ncpus ):
        inputQueue.put( 'STOP' )

    # Close both send and receive queues
    inputQueue.close()
    printQueue.close()

    temp_file.close()    
    models.close()
    output_hdf.close()
    print("Fitting time taken: {0:.2f} {1}".format(time.time()-loop_start,
                                                   '\n'))
    
    """
    Section 3
    Reload, format and save output table
    """
    while temp_file.closed == False:
        pause(0.1)

    data = np.loadtxt(output_path+".temp_output.txt")
    try:
        rows, cols = data.shape
    except:
        cols = len(data)

    output = Table()

    names = ['N', 'ID', 'z', 'zmodel', 
             'Mass_best', 'SFR_best', 'chi_best', 
             'Age_best','Dust_best', 'SFH_best',
             'Metallicity_best', 'fesc_best',
             'Mass_median', 'Mass_l68', 'Mass_u68',
             'SFR_median', 'SFR_l68', 'SFR_u68',
             'Nfilts']
             
    units = [None, None, None, None,
             u.Msun, u.Msun/u.yr, None,
             u.Gyr, None, None, 
             None, None,
             u.Msun, u.Msun, u.Msun,
             u.Msun/u.yr, u.Msun/u.yr, u.Msun/u.yr,
             None]
             
    types = ['i4', 'i4', 'f4', 'f4',
             'f4', 'f4', 'f4',
             'f4', 'f4', 'f4', 
             'f4', 'f4', 
             'f4', 'f4', 'f4',
             'f4', 'f4', 'f4',
             'i4']
             
    if include_rest:
        for name in filter_names:
            names.append(name[:-len(flux_col_end)]+'_rest')
            units.append(u.mag)
            types.append('f4')
        
    for col in range(cols):
        column = Column( data[:,col], name = names[col], unit=units[col], dtype=types[col])
        output.add_column(column)

    table_format = params.output_format
    output.sort('ID')
    if os.path.isfile(output_path):
        os.remove(output_path)
    output.write(output_path,format=table_format, overwrite=True)
    print('Catalog saved')
    
    os.remove(temp_file.name)

    print
    print "Total time taken: "+str(time.time()-start)
    
    sys.stderr = original_stderr
    logfile.close()
    
