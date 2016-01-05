import numpy as np
import array
import os
import sys
import re
import time
import multiprocessing

#import atpy
import logging
from astropy.table import Table, Column
from astropy import units as u
from astropy.utils.console import human_time, ProgressBar


def galaxyFit(inputQueue, printQueue, printlock):
    """ Fitting function definition for later use by Processes
    """
    for gal in iter(inputQueue.get, 'STOP'):
        j = np.argmin(np.abs(z-zobs[gal]))  # Find closest model redshift

        fo = obs[gal,:]
        ferr = obs_err[gal,:]

        #mtoterr = np.log10(np.e)*fo[params.tot]/ferr[params.tot]

        fo[fo <= 0.] = 0.       # Set negative fluxes to zero
        #print fo
        I = np.where(ferr > 0.)[0] # Find bands with no observation
        
        if len(I) == 0:
            if params.include_rest:
                M_scaled = np.ones(len(fo)) * -99.
                restframe_output = ' '.join(M_scaled.astype('str'))
                output_string = '{0} {1} {2} {3} {4} {5} {6} {7} {8} {9} {10} {11} {12} {13} {14} {15} {16}'.format(gal+1,ID[gal],zobs[gal],-99,-99,-99,-99,-99,-99, -99, -99,-99,len(I),-99,z[j],restframe_output,'\n')
            else:
                output_string = '{0} {1} {2} {3} {4} {5} {6} {7} {8} {9} {10} {11} {12} {13} {14}'.format(gal+1,ID[gal],zobs[gal],-99,-99,-99,-99,-99,-99,-99, -99,-99,len(I),-99,'\n')
            printQueue.put(output_string)
            continue
            
        fo = fo[I]                    # and exclude from fit
        ferr = ferr[I]
        fm = f[I, j, :]

        top = 0.
        bottom = 0.
    
        for i in range(len(fo)):
            top += np.abs(fm[i, :] * fo[i])/(ferr[i]**2)
            bottom += (fm[i, :]**2)/(ferr[i]**2)
    
        scale = top/bottom
        scale = np.reshape(scale,(n_tg,n_tauv,n_tau,n_metal))  

        chisq = 0.
        for i in range(len(fo)):
            chisq += (((scale*fm[i,:]-fo[i])**2)/(ferr[i])**2)

        chimin,minind = np.nanmin(chisq), np.nanargmin(chisq)
        if np.isinf(chimin) or np.isnan(minind):
            if params.include_rest:
                M_scaled = np.ones(len(fo)) * -99.
                restframe_output = ' '.join(M_scaled.astype('str'))
                output_string = '{0} {1} {2} {3} {4} {5} {6} {7} {8} {9} {10} {11} {12} {13} {14} {15} {16}'.format(gal+1,ID[gal],zobs[gal],-99,-99,-99,-99,-99,-99, -99, -99,-99,len(I),-99,z[j],restframe_output,'\n')
            else:
                output_string = '{0} {1} {2} {3} {4} {5} {6} {7} {8} {9} {10} {11} {12} {13} {14}'.format(gal+1,ID[gal],zobs[gal],-99,-99,-99,-99,-99,-99,-99, -99,-99,len(I),-99,'\n')
            printQueue.put(output_string)
            continue


        #Find the coordinate of the model with the bestfit mass
        tgi,tvi,ti,mi = np.unravel_index(minind,(n_tg,n_tauv,n_tau,n_metal))
        Bestfit_Mass = np.log10(scale[tgi,tvi,ti,mi]*params.flux_corr)
        Bestfit_SFR = (scale[tgi,tvi,ti,mi]*SFR[tgi,ti,mi]*params.flux_corr)
        Bestfit_Beta = beta[tgi,tvi,ti,mi]

        M_chisq_plus1 = np.log10(scale.flatten()[chisq.flatten() < chimin+1])
        Bestfit_Min, Bestfit_Max = np.min(M_chisq_plus1), np.max(M_chisq_plus1)


        #Scale the observed tot_mag band of the template to be the same as the observed tot_mag band of the galaxy
        #Convert the templates so they are no longer units of per stellar mass

        F_rest = f[:, 0] * scale[tgi, tvi, ti, mi] * params.flux_corr
        restframeMags = 23.9 - 2.5*np.log10(F_rest)
    
        UV_rest = UV_flux[0] * scale[tgi, tvi, ti, mi] * params.flux_corr
        restframeMUV = 23.9 - 2.5*np.log10(UV_rest)

        M_scaled = restframeMags[:, tgi, tvi, ti, mi]
        MUV_scaled = restframeMUV[tgi, tvi, ti, mi]
        
        if np.isnan(Bestfit_Mass) or np.isinf(chimin):
            Bestfit_Mass = -99
            #M_scaled[:] = -99
            tgs = -99
            tvs = -99
            taus = -99
            mis = -99

        else:
            tgs = tg[tgi]/1e9
            tvs = tv[tvi]
            taus = tau[ti]/1e9
            mis = mi

        if calc_mode:
            a = len(np.ravel(chisq))
            b = sum(np.isnan(np.ravel(chisq)))
            num_chosen = round(((a-b)/100.)*params.mode_mass_percentage)

            chisq_order = np.argsort(chisq)
            chisq_sorted = chisq[chisq_order]

            Mass_sorted = scale.flatten()[chisq_order]
            Mass_sorted[num_chosen:] = np.nan

            Mass_sorted = Mass_sorted*params.flux_corr

            Binned_Mass  = np.arange(np.min(np.log10(Mass_sorted)),
                                                                    np.max(np.log10(Mass_sorted)),
                                                                    0.05)

            No_in_bin, bin_edges = np.histogram(np.log10(Mass_sorted),Binned_Mass)
            mode_bin = np.argmax(No_in_bin)
            Mode_Mass = np.power(10,bin_edges[mode_bin]+0.025)
            print np.log10(Mode_Mass)
            if np.isnan(Mode_Mass):
                Mode_Mass = -99

        printlock.acquire()

        if calc_mode:
            print '{0:4d} {1:6d} {2:>6.2f} {3:>8.1f} {4:>6.2f}'.format(gal+1,ID[gal],Bestfit_Mass,chimin, np.log10(Mode_Mass), '/n')
        else:
            print '{0:6d} {1:8f} {2:>5.2f} {3:>7.2f} {4:>8.1f} {5:>8.3f} {6:>5.1f} {7:>8.2f} {8:>3d} {9:>5.2f}'.format(gal+1,ID[gal],zobs[gal],Bestfit_Mass,chimin,tgs,tvs,taus,mis,np.log10(Bestfit_SFR))

        if params.include_rest:
            restframe_output = ' '.join(M_scaled.astype('str'))
            output_string = '{0} {1} {2} {3} {4} {5} {6} {7} {8} {9} {10} {11} {12} {13} {14} {15} {16}'.format(gal+1,ID[gal],zobs[gal],Bestfit_Mass,chimin,tgs,tvs,taus,mis, MUV_scaled, minind,Bestfit_SFR,len(I),Bestfit_Beta,z[j],restframe_output,'\n')
        else:
            output_string = '{0} {1} {2} {3} {4} {5} {6} {7} {8} {9} {10} {11} {12} {13} {14}'.format(gal+1,ID[gal],zobs[gal],Bestfit_Mass,chimin,tgs,tvs,taus,mis, MUV_scaled, minind,Bestfit_SFR,len(I),Bestfit_Beta,'\n')

        printlock.release()
        printQueue.put(output_string)

def galaxyFitPlus(inputQueue, printQueue, printlock):
    for gal in iter(inputQueue.get, 'STOP'):
        mass_range = np.logspace(params.mass_min,params.mass_max,params.mass_bins)
        muvBins = np.linspace(params.muv_max,params.muv_min,params.muv_bins)
        betaBins = np.linspace(params.beta_min,params.beta_max,params.beta_bins)

        j = np.argmin(np.abs(z-zobs[gal])) # Find closest model redshift

        fo = obs[gal,:]
        ferr = obs_err[gal,:]


        fo[fo <= 0.] = 0.       # Set negative fluxes to zero
        #print fo
        I = (ferr > 0.)*(ferr < 1e6) # Find bands with no observation
        fo = fo[I]                    # and exclude from fit
        ferr = ferr[I]
        fm = f[I,j,:]
        #print fm[:,0,0,0,0]        

        chisq = np.zeros((params.mass_bins,n_tg,n_tauv,n_tau,n_metal))

        for m, mass in enumerate(mass_range):
            for i in range(len(fo)):
                chisq[m] += (((mass*fm[i,:]-fo[i])**2)/(ferr[i])**2)

        chimin,minind = np.nanmin(chisq), np.nanargmin(chisq)

        likelihood = np.exp(-0.5*chisq)
        likelihood /= likelihood.sum()
        
        if np.isinf(chimin) or np.isnan(minind) or len(fo) == 0:
            output_string = '{0} {1} {2} {3} {4} {5} {6} {7} {8} {9} \
            {10} {11} {12} {13} {14} {15} {16} {17} {18}'.format(gal+1,ID[gal],zobs[gal],
                                                                 -99,-99,-99,-99,-99,-99,
                                                                 -99, -99, -99, -99,-99,-99,-99,
                                                                 len(I),-99,'\n')
            
            massLikelihood = np.zeros(params.mass_bins+1)
            massLikelihood[0] = gal
            muvLikelihood = np.zeros(params.muv_bins+1)
            muvLikelihood[0] = gal
            betaLikelihood = np.zeros(params.beta_bins+1)
            betaLikelihood[0] = gal
            #tauLikelihood = np.zeros(n_tau)        
            #tauLikelihood = np.insert(tauLikelihood,0,gal)        
            printQueue.put([output_string,massLikelihood,muvLikelihood,betaLikelihood])
            continue

        #Find the coordinate of the model with the bestfit mass
        si,tgi,tvi,ti,mi = np.unravel_index(minind,(params.mass_bins,n_tg,n_tauv,n_tau,n_metal))
        Bestfit_Mass = np.log10(mass_range[si]*params.flux_corr)
        Bestfit_SFR = (mass_range[si]*SFR[tgi,ti,mi]*params.flux_corr)
        Bestfit_Beta = beta[tgi,tvi,ti,mi]

        F_rest = f[:,0]*mass_range[likelihood.argmax(0)]*params.flux_corr
        restframeMags = 23.9 - 2.5*np.log10(F_rest)
    
        UV_rest = UV_flux[0]*mass_range[likelihood.argmax(0)]*params.flux_corr
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
            tgs = tg[tgi]/1e9
            tvs = tv[tvi]
            taus = tau[ti]/1e9
            mis = mi
            
        """
        Likelihood array section:
        """
        
        muvSize = np.diff(muvBins)[0]
        betaSize = np.diff(betaBins)[0]
        
        massLikelihoods = likelihood.sum(axis=(1,2,3,4))            
        massLikelihoods /= np.trapz(massLikelihoods,np.log10(mass_range))
        massLikelihoods = np.insert(massLikelihoods,0,gal)
            
        muvLikelihoods = array.array('d')
        for muv in muvBins:
            window = ((restframeMUV.flatten() > muv - 0.5*muvSize) * 
                        (restframeMUV.flatten() <= muv + 0.5*muvSize))
            muvLikelihoods.append(np.sum(likelihood.sum(0).flatten()[window]))
        
        muvLikelihoods /= np.trapz(muvLikelihoods,muvBins)
        muvLikelihoods = np.insert(muvLikelihoods,0,gal)
        
        betaLikelihoods = array.array('d')
        for fbeta in betaBins:
            window = ((beta.flatten() > fbeta - 0.5*betaSize)*
                        (beta.flatten() <= fbeta + 0.5*betaSize))
            betaLikelihoods.append(np.sum(likelihood.sum(0).flatten()[window]))
            
        betaLikelihoods /= np.trapz(betaLikelihoods,betaBins)        
        betaLikelihoods = np.insert(betaLikelihoods,0,gal)
        
        """
        #betaQueue.put(betaLikelihoods)
        likelihood_shaped[np.invert(np.isfinite(likelihood_shaped))] = 0.
        tauLikelihood = np.sum(likelihood_shaped,axis=(0,1,3))
        tauLikelihood /= sum(tauLikelihood)
        print tauLikelihood
        tauLikelihood = np.insert(tauLikelihood,0,gal)
        """
        
        printlock.acquire()

        if calc_mode:
            print '{0:4d} {1:6d} {2:>6.2f} {3:>8.1f} {4:>6.2f}'.format(gal+1,ID[gal],Bestfit_Mass,chimin, np.log10(Mode_Mass), '/n')
        else:
            print '{0:6d} {1:8f} {2:>5.2f} {3:>7.2f} {4:>8.1f} {5:>8.3f} {6:>5.1f} {7:>8.2f} {8:>3d} {9:>5.2f}'.format(gal+1,int(ID[gal]),zobs[gal],Bestfit_Mass,chimin,tgs,tvs,taus,mis,np.log10(Bestfit_SFR))

        output_string = '{0} {1} {2} {3} {4} {5} {6} {7} {8} {9} {10} {11} {12} {13} {14} {15}'.format(gal+1,int(ID[gal]),zobs[gal],Bestfit_Mass,chimin,tgs,tvs,taus,mis,Bestfit_restframeMags[params.tot],Bestfit_restframeMUV,minind,Bestfit_SFR,len(I),Bestfit_Beta,'\n')

        printlock.release()
        printQueue.put([output_string, massLikelihoods, muvLikelihoods, betaLikelihoods])


class Observations:
    """ Class to parse and wrap a catalog of observed photometry"""
    def __init__(self, inputpath, cat_format='ascii', id_col='ID', z_col='z',
                 flux_suffix='_flux', fluxerr_suffix='_fluxerr', filts_used=None):

        input_data = Table.read(inputpath, format=cat_format)

        column_names = input_data.colnames

        self.ID = input_data[id_col]
        self.zobs = input_data[z_col]

        self.filter_names = []

        self.fluxes = None
        self.fluxerrs = None

        k,l = 0,0
        for ii in range(len(column_names)):
            if column_names[ii].lower().endswith(flux_suffix.lower()):
                if k == 0:
                    fluxes = input_data[column_names[ii]]
                else:
                    fluxes = np.column_stack((fluxes, input_data[column_names[ii]]))
                k += 1
                filter_names.append(column_names[ii])

            if column_names[ii].lower().endswith(fluxerr_suffix.lower()):
                if l == 0:
                    fluxerrs = input_data[column_names[ii]]
                else:
                    fluxerrs = np.column_stack((fluxerrs, input_data[column_names[ii]]))
                l += 1


        if filts_used != None:
            fluxes = fluxes[:,filts_used]
            fluxerrs = fluxerrs[:,filts_used]
            filter_names = filer_names[filts_used]
    

class Fitting(object):
    """ Base class for fitting model SEDS to observed photometry 
    """
    def __init__(self, observations, models, ncpus=1):
        self.obs = observations
        self.models = models
        self.tg = self.models.f['tg'].value
        self.tv = self.models.f['Av'].value
        self.tau = self.models.f['SFH'].value
        self.SSP = self.models.f['metallicities'].value
        self.fesc = self.models.f['fesc'].value

        Mshape = models.fluxes.shape
        self.z = models.redshifts
        self.nfilts = Mshape[1]
        self.n_metal = Mshape[2]        
        self.n_tg = Mshape[3]
        self.n_tau = Mshape[4]
        self.n_tauv = Mshape[5]
        self.n_fesc = Mshape[6]
        
        #UV_flux = synmags['UV_flux']
        self.SFR = self.models.SFR
        self.Ms = self.models.Ms

        # beta = obs.beta

        if self.obs.filts_used != None:
            self.f = self.models.fluxes[self.obs.filts_used]
        else:
            self.f = self.models.fluxes

        self.ncpus = np.clip(ncpus, 1, multiprocessing.cpu_count())



class Simple(Fitting):
    def __init__(self, observations, models, output_path, ncpus=1):
        super(Simple, self).__init__(observations, models, ncpus)


        if os.path.isfile(params.output_name+".temp_output.txt"):
            os.remove(params.output_name+".temp_output.txt")
        temp_file = open(params.output_name+".temp_output.txt","w")

        if calc_mode:
            print('{0:>4s} {1:>8s} {2:>6s} {3:>8s} {4:>6s}'.format('N','ID','Best', 'chimin', 'Mode'))

        else:
            print('{0:>4s} {1:>8s} {2:>5s} {3:>6s} {4:>8s} {5:>10s} {6:>5s} {7:>12s} {8:>4s}'.format('N','ID','zobs','Best', 'chimin', 'tg', 'tauv','tau','SSP'))

        loop_start = time.time()

        inputQueue = multiprocessing.Queue()
        printQueue = multiprocessing.Queue()

        printlock = multiprocessing.Lock()

        fitFunction = galaxyFit

        for i in range( ncpus ):
            multiprocessing.Process(target=fitFunction,
                                    args=(inputQueue, printQueue, printlock)).start()

        # Put elements in the send queue for processing
        for gal in range(len(ID)):
            inputQueue.put(gal)


        for gal in range( len(ID) ):
            printout = printQueue.get()
            temp_file.write( printout )
            #print len(mass_array), len(muv_array), len(beta_array)

        # Stop all the running processes
        for i in range(ncpus):
            inputQueue.put('STOP')

        # Close both send and receive queues
        inputQueue.close()
        printQueue.close()

        temp_file.close()

        while not temp_file.closed:
            pause(0.1)

        data = np.loadtxt(params.output_name+".temp_output.txt")
        try:
            rows, cols = data.shape
        except:
            cols = len(data)

        self.results = Table(meta=parameters)

        names = ['N','ID','z','Bestfit_Mass','Bestfit_chi2','Age','Dust_Tau','SFH_Tau','SSP_Number','TotCol_rest', 'M1500','temp_index','SFR','nfilts','Beta','z_model']
        units = [None,None,None, u.Msun ,None, u.Gyr, None, u.Gyr, None, u.mag, u.mag,None,u.Msun/u.yr,None,None,None]
        types = ['i4','i4','f4','f4','f4','f4','f4','f4','i4', 'f4', 'f4','f4','f4','i4','f4','f4']
        if params.include_rest:
            for name in filter_names:
                names.append(name[:-len(params.flux_col_end)]+'_rest')
                units.append(u.mag)
                types.append('f4')

        for col in range(cols):
            column = Column(data[:, col], name=names[col], unit=units[col], dtype=types[col])
            self.results.add_column(column)

        self.results.sort('ID')
        if os.path.isfile(params.output_name):
            os.remove(params.output_name)
        self.results.write(params.output_name, format=params.table_format)
        print('Catalog saved')

        os.remove(temp_file.name)

        print("Total time taken: {0}".format(human_time(time.time()-start)))

        sys.stderr = original_stderr
        logfile.close()


if __name__ == '__main__':
    
    logfile = open("error.log", "w")
    original_stderr = sys.stderr
    sys.stderr = logfile
    
    start = time.time()


    """
    SECTION 1C
    Setting up output table
    """
    if os.path.isfile(params.output_name+".temp_output.txt"):
        os.remove(params.output_name+".temp_output.txt")
    temp_file = open(params.output_name+".temp_output.txt","w")
    
    mass_file = open(params.output_name+".temp_masses.prob", "wb")
    muv_file = open(params.output_name+".temp_muv.prob", "wb")
    beta_file = open(params.output_name+".temp_betas.prob", "wb")
    tau_file = open(params.output_name+".temp_taus.prob","wb")
    

    """
    SECTION 2
    Chi-sq calculation

    """
    if calc_mode:
        print '{0:>4s} {1:>8s} {2:>6s} {3:>8s} {4:>6s}'.format('N','ID','Best', 'chimin', 'Mode')

    else:
        print '{0:>4s} {1:>8s} {2:>5s} {3:>6s} {4:>8s} {5:>10s} {6:>5s} {7:>12s} {8:>4s}'.format('N','ID','zobs','Best', 'chimin', 'tg', 'tauv','tau','SSP')

    loop_start = time.time()
    ncpus = np.clip(params.ncpus,1,multiprocessing.cpu_count())

    inputQueue = multiprocessing.Queue()
    printQueue = multiprocessing.Queue()
    #massQueue = multiprocessing.Queue()
    #muvQueue = multiprocessing.Queue()
    #betaQueue = multiprocessing.Queue()
    
    printlock = multiprocessing.Lock()
    
    if params.calc_pdf:
        fitFunction = galaxyFitPlus
    else:
        fitFunction = galaxyFit
    
    for i in range( ncpus ):
        #multiprocessing.Process( target = galaxyFit, args = ( inputQueue , printQueue, printlock ) ).start()
        multiprocessing.Process( target = fitFunction,
                                args = (inputQueue, printQueue, printlock ) ).start()

    # Put elements in the send queue for processing
    for gal in range( len(ID) ):
        inputQueue.put( gal )

    if params.calc_pdf:
        for gal in range( len(ID) ):
            printout, mass_array, muv_array, beta_array = printQueue.get()
            temp_file.write( printout )
            #print len(mass_array), len(muv_array), len(beta_array)
            mass_array.tofile(mass_file)
            muv_array.tofile(muv_file)
            beta_array.tofile(beta_file)
        #tau_array.tofile(tau_file)
    else:
        for gal in range( len(ID) ):
            printout = printQueue.get()
            temp_file.write( printout )
            #print len(mass_array), len(muv_array), len(beta_array)
       

    # Stop all the running processes
    for i in range( ncpus ):
        inputQueue.put( 'STOP' )

    # Close both send and receive queues
    inputQueue.close()
    printQueue.close()
    #massQueue.close()
    #muvQueue.close()
    #betaQueue.close()


    temp_file.close()
    mass_file.close()
    muv_file.close()
    beta_file.close()
    tau_file.close()
    #print results
    print "Fitting time taken: "+str(time.time()-loop_start)
    print


    """
    Section 3
    Reload, format and save output table
    """
    while temp_file.closed == False:
        pause(0.1)

    data = np.loadtxt(params.output_name+".temp_output.txt")
    try:
        rows, cols = data.shape
    except:
        cols = len(data)

    output = Table(meta=parameters)

    names = ['N','ID','z','Bestfit_Mass','Bestfit_chi2','Age','Dust_Tau','SFH_Tau','SSP_Number','TotCol_rest', 'M1500','temp_index','SFR','nfilts','Beta','z_model']
    units = [None,None,None, u.Msun ,None, u.Gyr, None, u.Gyr, None, u.mag, u.mag,None,u.Msun/u.yr,None,None,None]
    types = ['i4','i4','f4','f4','f4','f4','f4','f4','i4', 'f4', 'f4','f4','f4','i4','f4','f4']
    if params.include_rest:
        for name in filter_names:
            names.append(name[:-len(params.flux_col_end)]+'_rest')
            units.append(u.mag)
            types.append('f4')
        
    for col in range(cols):
        column = Column( data[:,col], name = names[col], unit=units[col], dtype=types[col])
        output.add_column(column)

    output.sort('ID')
    if os.path.isfile(params.output_name):
        os.remove(params.output_name)
    output.write(params.output_name,format=params.table_format)
    print('Catalog saved')
    
    if params.calc_pdf:
        print
        print('Opening, re-ordering and saving binaries:')
    
        print('{0:<20s}'.format('Masses')),
        with open(mass_file.name) as mass_binary:
            mass_likelihood = array.array('d')
            mass_likelihood.fromfile(mass_binary,len(ID)*(params.mass_bins + 1))
            mass_likelihood = np.reshape(mass_likelihood,(len(ID),(params.mass_bins + 1)))
            mass_likelihood = mass_likelihood[np.argsort(mass_likelihood[:,0]),1:]
    
            output_binary = open(params.output_name+".masses.prob","wb")
            mass_params = array.array('i',[len(ID),params.mass_bins])
            mass_params.tofile(output_binary)
        
            mass_bins = np.linspace(params.mass_min,params.mass_max,params.mass_bins)
            mass_bins.tofile(output_binary)
        
            for gal in range(len(ID)):
                mass_likelihood[gal,:].tofile(output_binary)
        
            output_binary.close()
            print('Done')
    
        print('{0:<20s}'.format('MUVs')),
        with open(muv_file.name) as muv_binary:
            muv_likelihood = array.array('d')
            muv_likelihood.fromfile(muv_binary,len(ID)*(params.muv_bins + 1))
            muv_likelihood = np.reshape(muv_likelihood,(len(ID),(params.muv_bins + 1)))
            muv_likelihood = muv_likelihood[np.argsort(muv_likelihood[:,0]),1:]
    
            output_binary = open(params.output_name+".muv.prob","wb")
            muv_params = array.array('i',[len(ID),params.muv_bins])
            muv_params.tofile(output_binary)
        
            muv_bins = np.linspace(params.muv_max,params.muv_min,params.muv_bins)
            muv_bins.tofile(output_binary)
        
            for gal in range(len(ID)):
                muv_likelihood[gal,:].tofile(output_binary)
        
            output_binary.close()
            print('Done')
        
        print('{0:<20s}'.format('Betas')),
        with open(beta_file.name) as beta_binary:
            beta_likelihood = array.array('d')
            beta_likelihood.fromfile(beta_binary,len(ID)*(params.beta_bins + 1))
            beta_likelihood = np.reshape(beta_likelihood,(len(ID),(params.beta_bins + 1)))
            beta_likelihood = beta_likelihood[np.argsort(beta_likelihood[:,0]),1:]
    
            output_binary = open(params.output_name+".betas.prob","wb")
            beta_params = array.array('i',[len(ID),params.beta_bins])
            beta_params.tofile(output_binary)
        
            beta_bins = np.linspace(params.beta_min,params.beta_max,params.beta_bins)
            beta_bins.tofile(output_binary)
        
            for gal in range(len(ID)):
                beta_likelihood[gal,:].tofile(output_binary)
        
            output_binary.close()
            print('Done')


    os.remove(mass_file.name)
    os.remove(muv_file.name)
    os.remove(beta_file.name)
    os.remove(tau_file.name)
    os.remove(temp_file.name)

    print
    print "Total time taken: "+str(time.time()-start)
    
    sys.stderr = original_stderr
    logfile.close()
    
