import numpy
import array
import os, sys
import re
import time
import multiprocessing
#import atpy
import logging
from astropy.table import Table, Column
from astropy import units as u

from sm_functions import dist

"""
Command line arguments and parameter import section
"""
version = sys.version_info[1]
if version == 7:
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-p","--params", type=str, default="sm_params",
                        help = "Parameter file, default = sm_params")
    parser.add_argument("-q", "--quiet", help = "Suppress extra outputs",
                        action = "store_true")
    args = parser.parse_args()
    quiet = args.quiet

    params_root = re.split(".py",args.params)[0]
    if os.path.isfile(params_root+".pyc"):
        os.remove(params_root+".pyc")

    import importlib
    try:
        params = importlib.import_module(params_root)
        print 'Loaded '+args.params+' as params'
    except:
        print 'Failed to load "'+args.params+'" as params, loading default instead'
        import sm_params as params

if version == 6:
    import optparse
    parser = optparse.OptionParser()
    parser.add_option("-p","--params", type=str, default="sm_params",
                                            dest="params",help = "Parameter file, default = sm_params")
    parser.add_option("-q", "--quiet", help = "Suppress extra outputs",
                                            dest="quiet", default=False,
                                            action = "store_true")
    args, dump = parser.parse_args()
    quiet = args.quiet

    params_root = re.split(".py",args.params)[0]
    if os.path.isfile(params_root+".pyc"):
        os.remove(params_root+".pyc")
    import imp
    try:
        fp, pathname, description = imp.find_module(params_root)
        params = imp.load_module(params_root, fp, pathname, description)
    except:
        print 'Failed to load "'+args.params+'" as params, loading default instead'
        import sm_params as params

elif version not in [6, 7]:
    print 'Import option only coded for python versions 2.6 and 2.7... \n',
    print 'Loading default instead'
    import sm_params as params
    quiet = False

if quiet:
    print "Shhhh"

if quiet:
    quietprint = lambda *a: None
else:
    def quietprint(*args):
        for arg in args:
            print arg,
        print


"""
Fitting function definition for later use by Processes

"""

def galaxyFit(inputQueue, printQueue, printlock):
    for gal in iter(inputQueue.get, 'STOP'):
        j = numpy.argmin(numpy.abs(z-zobs[gal])) # Find closest model redshift


        fo = obs[gal,:]
        ferr = obs_err[gal,:]

        #mtoterr = numpy.log10(numpy.e)*fo[params.tot]/ferr[params.tot]

        fo[fo <= 0.] = 0.       # Set negative fluxes to zero
        #print fo
        I = numpy.where(ferr > 0.)[0] # Find bands with no observation
        
        if len(I) == 0:
            if params.include_rest:
                M_scaled = numpy.ones(len(fo)) * -99.
                restframe_output = ' '.join(M_scaled.astype('str'))
                output_string = '{0} {1} {2} {3} {4} {5} {6} {7} {8} {9} {10} {11} {12} {13} {14} {15} {16}'.format(gal+1,ID[gal],zobs[gal],-99,-99,-99,-99,-99,-99, -99, -99,-99,len(I),-99,z[j],restframe_output,'\n')
            else:
                output_string = '{0} {1} {2} {3} {4} {5} {6} {7} {8} {9} {10} {11} {12} {13} {14}'.format(gal+1,ID[gal],zobs[gal],-99,-99,-99,-99,-99,-99,-99, -99,-99,len(I),-99,'\n')
            printQueue.put(output_string)
            continue
            
        fo = fo[I]                    # and exclude from fit
        ferr = ferr[I]
        fm = f[I,j,:]

            #print fm[:,0,0,0,0]

        top = 0.
        bottom = 0.
    
        for i in range(len(fo)):
            top += (fm[i,:]*fo[i])/(ferr[i]**2)
            bottom += (fm[i,:]**2)/(ferr[i]**2)
    
        scale = top/bottom
        scale = numpy.reshape(scale,(n_tg,n_tauv,n_tau,n_ssp))  

        chisq = 0.
        for i in range(len(fo)):
            chisq += (((scale*fm[i,:]-fo[i])**2)/(ferr[i])**2)

        chimin,minind = numpy.nanmin(chisq), numpy.nanargmin(chisq)
        if numpy.isinf(chimin) or numpy.isnan(minind):
            if params.include_rest:
                M_scaled = numpy.ones(len(fo)) * -99.
                restframe_output = ' '.join(M_scaled.astype('str'))
                output_string = '{0} {1} {2} {3} {4} {5} {6} {7} {8} {9} {10} {11} {12} {13} {14} {15} {16}'.format(gal+1,ID[gal],zobs[gal],-99,-99,-99,-99,-99,-99, -99, -99,-99,len(I),-99,z[j],restframe_output,'\n')
            else:
                output_string = '{0} {1} {2} {3} {4} {5} {6} {7} {8} {9} {10} {11} {12} {13} {14}'.format(gal+1,ID[gal],zobs[gal],-99,-99,-99,-99,-99,-99,-99, -99,-99,len(I),-99,'\n')
            printQueue.put(output_string)
            continue


        #Find the coordinate of the model with the bestfit mass
        tgi,tvi,ti,mi = numpy.unravel_index(minind,(n_tg,n_tauv,n_tau,n_ssp))
        Bestfit_Mass = numpy.log10(scale[tgi,tvi,ti,mi]*params.flux_corr)
        Bestfit_SFR = (scale[tgi,tvi,ti,mi]*SFR[tgi,ti,mi]*params.flux_corr)
        Bestfit_Beta = beta[tgi,tvi,ti,mi]

        M_chisq_plus1 = numpy.log10(scale.flatten()[chisq.flatten() < chimin+1])
        Bestfit_Min, Bestfit_Max = numpy.min(M_chisq_plus1), numpy.max(M_chisq_plus1)


        #Scale the observed tot_mag band of the template to be the same as the observed tot_mag band of the galaxy
        #Convert the templates so they are no longer units of per stellar mass

        F_rest = f[:,0]*scale[tgi,tvi,ti,mi]*params.flux_corr
        restframeMags = 23.9 - 2.5*numpy.log10(F_rest)
    
        UV_rest = UV_flux[0]*scale[tgi,tvi,ti,mi]*params.flux_corr
        restframeMUV = 23.9 - 2.5*numpy.log10(UV_rest)

        M_scaled = restframeMags[:,tgi,tvi,ti,mi]
        MUV_scaled = restframeMUV[tgi,tvi,ti,mi]
        
        if numpy.isnan(Bestfit_Mass) or numpy.isinf(chimin):
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
            a = len(numpy.ravel(chisq))
            b = sum(numpy.isnan(numpy.ravel(chisq)))
            num_chosen = round(((a-b)/100.)*params.mode_mass_percentage)

            chisq_order = numpy.argsort(chisq)
            chisq_sorted = chisq[chisq_order]

            Mass_sorted = scale.flatten()[chisq_order]
            Mass_sorted[num_chosen:] = numpy.nan

            Mass_sorted = Mass_sorted*params.flux_corr

            Binned_Mass  = numpy.arange(numpy.min(numpy.log10(Mass_sorted)),
                                                                    numpy.max(numpy.log10(Mass_sorted)),
                                                                    0.05)

            No_in_bin, bin_edges = numpy.histogram(numpy.log10(Mass_sorted),Binned_Mass)
            mode_bin = numpy.argmax(No_in_bin)
            Mode_Mass = numpy.power(10,bin_edges[mode_bin]+0.025)
            print numpy.log10(Mode_Mass)
            if numpy.isnan(Mode_Mass):
                Mode_Mass = -99

        printlock.acquire()

        if calc_mode:
            print '{0:4d} {1:6d} {2:>6.2f} {3:>8.1f} {4:>6.2f}'.format(gal+1,ID[gal],Bestfit_Mass,chimin, numpy.log10(Mode_Mass), '/n')
        else:
            print '{0:6d} {1:8f} {2:>5.2f} {3:>7.2f} {4:>8.1f} {5:>8.3f} {6:>5.1f} {7:>8.2f} {8:>3d} {9:>5.2f}'.format(gal+1,ID[gal],zobs[gal],Bestfit_Mass,chimin,tgs,tvs,taus,mis,numpy.log10(Bestfit_SFR))

        if params.include_rest:
            restframe_output = ' '.join(M_scaled.astype('str'))
            output_string = '{0} {1} {2} {3} {4} {5} {6} {7} {8} {9} {10} {11} {12} {13} {14} {15} {16}'.format(gal+1,ID[gal],zobs[gal],Bestfit_Mass,chimin,tgs,tvs,taus,mis, MUV_scaled, minind,Bestfit_SFR,len(I),Bestfit_Beta,z[j],restframe_output,'\n')
        else:
            output_string = '{0} {1} {2} {3} {4} {5} {6} {7} {8} {9} {10} {11} {12} {13} {14}'.format(gal+1,ID[gal],zobs[gal],Bestfit_Mass,chimin,tgs,tvs,taus,mis, MUV_scaled, minind,Bestfit_SFR,len(I),Bestfit_Beta,'\n')

        printlock.release()
        printQueue.put(output_string)

def galaxyFitPlus(inputQueue, printQueue, printlock):
    for gal in iter(inputQueue.get, 'STOP'):
        mass_range = numpy.logspace(params.mass_min,params.mass_max,params.mass_bins)
        muvBins = numpy.linspace(params.muv_max,params.muv_min,params.muv_bins)
        betaBins = numpy.linspace(params.beta_min,params.beta_max,params.beta_bins)

        j = numpy.argmin(numpy.abs(z-zobs[gal])) # Find closest model redshift

        fo = obs[gal,:]
        ferr = obs_err[gal,:]


        fo[fo <= 0.] = 0.       # Set negative fluxes to zero
        #print fo
        I = (ferr > 0.)*(ferr < 1e6) # Find bands with no observation
        fo = fo[I]                    # and exclude from fit
        ferr = ferr[I]
        fm = f[I,j,:]
        #print fm[:,0,0,0,0]        

        chisq = numpy.zeros((params.mass_bins,n_tg,n_tauv,n_tau,n_ssp))

        for m, mass in enumerate(mass_range):
            for i in range(len(fo)):
                chisq[m] += (((mass*fm[i,:]-fo[i])**2)/(ferr[i])**2)

        chimin,minind = numpy.nanmin(chisq), numpy.nanargmin(chisq)

        likelihood = numpy.exp(-0.5*chisq)
        likelihood /= likelihood.sum()
        
        if numpy.isinf(chimin) or numpy.isnan(minind) or len(fo) == 0:
            output_string = '{0} {1} {2} {3} {4} {5} {6} {7} {8} {9} \
            {10} {11} {12} {13} {14} {15} {16} {17} {18}'.format(gal+1,ID[gal],zobs[gal],
                                                                 -99,-99,-99,-99,-99,-99,
                                                                 -99, -99, -99, -99,-99,-99,-99,
                                                                 len(I),-99,'\n')
            
            massLikelihood = numpy.zeros(params.mass_bins+1)
            massLikelihood[0] = gal
            muvLikelihood = numpy.zeros(params.muv_bins+1)
            muvLikelihood[0] = gal
            betaLikelihood = numpy.zeros(params.beta_bins+1)
            betaLikelihood[0] = gal
            #tauLikelihood = numpy.zeros(n_tau)        
            #tauLikelihood = numpy.insert(tauLikelihood,0,gal)        
            printQueue.put([output_string,massLikelihood,muvLikelihood,betaLikelihood])
            continue

        #Find the coordinate of the model with the bestfit mass
        si,tgi,tvi,ti,mi = numpy.unravel_index(minind,(params.mass_bins,n_tg,n_tauv,n_tau,n_ssp))
        Bestfit_Mass = numpy.log10(mass_range[si]*params.flux_corr)
        Bestfit_SFR = (mass_range[si]*SFR[tgi,ti,mi]*params.flux_corr)
        Bestfit_Beta = beta[tgi,tvi,ti,mi]

        F_rest = f[:,0]*mass_range[likelihood.argmax(0)]*params.flux_corr
        restframeMags = 23.9 - 2.5*numpy.log10(F_rest)
    
        UV_rest = UV_flux[0]*mass_range[likelihood.argmax(0)]*params.flux_corr
        restframeMUV = 23.9 - 2.5*numpy.log10(UV_rest)

        Bestfit_restframeMags = restframeMags[:,tgi,tvi,ti,mi]
        Bestfit_restframeMUV = restframeMUV[tgi,tvi,ti,mi]

        if numpy.isnan(Bestfit_Mass) or numpy.isinf(chimin):
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
        
        muvSize = numpy.diff(muvBins)[0]
        betaSize = numpy.diff(betaBins)[0]
        
        massLikelihoods = likelihood.sum(axis=(1,2,3,4))            
        massLikelihoods /= numpy.trapz(massLikelihoods,numpy.log10(mass_range))
        massLikelihoods = numpy.insert(massLikelihoods,0,gal)
            
        muvLikelihoods = array.array('d')
        for muv in muvBins:
            window = ((restframeMUV.flatten() > muv - 0.5*muvSize) * 
                        (restframeMUV.flatten() <= muv + 0.5*muvSize))
            muvLikelihoods.append(numpy.sum(likelihood.sum(0).flatten()[window]))
        
        muvLikelihoods /= numpy.trapz(muvLikelihoods,muvBins)
        muvLikelihoods = numpy.insert(muvLikelihoods,0,gal)
        
        betaLikelihoods = array.array('d')
        for fbeta in betaBins:
            window = ((beta.flatten() > fbeta - 0.5*betaSize)*
                        (beta.flatten() <= fbeta + 0.5*betaSize))
            betaLikelihoods.append(numpy.sum(likelihood.sum(0).flatten()[window]))
            
        betaLikelihoods /= numpy.trapz(betaLikelihoods,betaBins)        
        betaLikelihoods = numpy.insert(betaLikelihoods,0,gal)
        
        """
        #betaQueue.put(betaLikelihoods)
        likelihood_shaped[numpy.invert(numpy.isfinite(likelihood_shaped))] = 0.
        tauLikelihood = numpy.sum(likelihood_shaped,axis=(0,1,3))
        tauLikelihood /= sum(tauLikelihood)
        print tauLikelihood
        tauLikelihood = numpy.insert(tauLikelihood,0,gal)
        """
        
        printlock.acquire()

        if calc_mode:
            print '{0:4d} {1:6d} {2:>6.2f} {3:>8.1f} {4:>6.2f}'.format(gal+1,ID[gal],Bestfit_Mass,chimin, numpy.log10(Mode_Mass), '/n')
        else:
            print '{0:6d} {1:8f} {2:>5.2f} {3:>7.2f} {4:>8.1f} {5:>8.3f} {6:>5.1f} {7:>8.2f} {8:>3d} {9:>5.2f}'.format(gal+1,int(ID[gal]),zobs[gal],Bestfit_Mass,chimin,tgs,tvs,taus,mis,numpy.log10(Bestfit_SFR))

        output_string = '{0} {1} {2} {3} {4} {5} {6} {7} {8} {9} {10} {11} {12} {13} {14} {15}'.format(gal+1,int(ID[gal]),zobs[gal],Bestfit_Mass,chimin,tgs,tvs,taus,mis,Bestfit_restframeMags[params.tot],Bestfit_restframeMUV,minind,Bestfit_SFR,len(I),Bestfit_Beta,'\n')

        printlock.release()
        printQueue.put([output_string, massLikelihoods, muvLikelihoods, betaLikelihoods])


def getObservations(inputpath):
    input_data = Table.read(inputpath,format=params.input_format)

    column_names = input_data.columns.keys()

    ID = input_data[params.ID_col]
    zobs = input_data[params.z_col]

    filts_used = params.filts_used
    filter_names = []

    k,l = 0,0
    for ii in range(len(column_names)):
        if column_names[ii].lower().endswith(params.flux_col_end.lower()):
            if k == 0:
                fluxes = input_data[column_names[ii]]
            else:
                fluxes = numpy.column_stack((fluxes,input_data[column_names[ii]]))
            k+=1
            filter_names.append(column_names[ii])

        if column_names[ii].lower().endswith(params.fluxerr_col_end.lower()):
            if l == 0:
                fluxerrs = input_data[column_names[ii]]
            else:
                fluxerrs = numpy.column_stack((fluxerrs,input_data[column_names[ii]]))
            l+=1

    fluxes = fluxes[:,filts_used]
    fluxerrs = fluxerrs[:,filts_used]

    return ID, zobs, fluxes, fluxerrs, k, filter_names


if __name__ == '__main__':
    
    logfile = open("error.log", "w")
    original_stderr = sys.stderr
    sys.stderr = logfile
    
    start = time.time()
    """
    SECTION 1A

    Loads in observed magnitudes/fluxes from catalog using ATPY
    From the catalog column names, it then picks out all the columns
    which it believes are magnitude/flux columns and their corresponding
    errors.

    """

    ID, zobs, obs, obs_err, filters_found, filter_names = getObservations(params.input_catalog)

    """
    Section 1C

    Loads the various numpy arrays from the synthesised magnitude binary

    """

    #Load in synthesised data
    input_binary = params.synmag_output

    print "Loading synthetic mags and mass array:"
    parameters = numpy.load(params.ssp_output+'.par.npy').item()
    tg = parameters['tg']
    tv = parameters['Av']
    tau = parameters['SFH']
    SSP = parameters['SSPs']

    with numpy.load(input_binary+'.main.npz') as synmags:
        Mshape = synmags['Mshape']
        z = synmags['z']
        
        nfilts = Mshape[0]
        n_tg = Mshape[2]
        n_tauv = Mshape[3]
        n_tau = Mshape[4]
        n_ssp = Mshape[5]
        
        UV_flux = synmags['UV_flux']
        SFR = synmags['SFR']

    with numpy.load(params.ssp_output) as beta_data:
        beta = beta_data['beta']

    calc_mode = params.calc_mode

    if nfilts != filters_found:
        f = numpy.load(input_binary+'.fluxes.npy')[params.filts_used]
    else:
        f = numpy.load(input_binary+'.fluxes.npy')


    print "Done."


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
    ncpus = numpy.clip(params.ncpus,1,multiprocessing.cpu_count())

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

    data = numpy.loadtxt(params.output_name+".temp_output.txt")
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
            mass_likelihood = numpy.reshape(mass_likelihood,(len(ID),(params.mass_bins + 1)))
            mass_likelihood = mass_likelihood[numpy.argsort(mass_likelihood[:,0]),1:]
    
            output_binary = open(params.output_name+".masses.prob","wb")
            mass_params = array.array('i',[len(ID),params.mass_bins])
            mass_params.tofile(output_binary)
        
            mass_bins = numpy.linspace(params.mass_min,params.mass_max,params.mass_bins)
            mass_bins.tofile(output_binary)
        
            for gal in range(len(ID)):
                mass_likelihood[gal,:].tofile(output_binary)
        
            output_binary.close()
            print('Done')
    
        print('{0:<20s}'.format('MUVs')),
        with open(muv_file.name) as muv_binary:
            muv_likelihood = array.array('d')
            muv_likelihood.fromfile(muv_binary,len(ID)*(params.muv_bins + 1))
            muv_likelihood = numpy.reshape(muv_likelihood,(len(ID),(params.muv_bins + 1)))
            muv_likelihood = muv_likelihood[numpy.argsort(muv_likelihood[:,0]),1:]
    
            output_binary = open(params.output_name+".muv.prob","wb")
            muv_params = array.array('i',[len(ID),params.muv_bins])
            muv_params.tofile(output_binary)
        
            muv_bins = numpy.linspace(params.muv_max,params.muv_min,params.muv_bins)
            muv_bins.tofile(output_binary)
        
            for gal in range(len(ID)):
                muv_likelihood[gal,:].tofile(output_binary)
        
            output_binary.close()
            print('Done')
        
        print('{0:<20s}'.format('Betas')),
        with open(beta_file.name) as beta_binary:
            beta_likelihood = array.array('d')
            beta_likelihood.fromfile(beta_binary,len(ID)*(params.beta_bins + 1))
            beta_likelihood = numpy.reshape(beta_likelihood,(len(ID),(params.beta_bins + 1)))
            beta_likelihood = beta_likelihood[numpy.argsort(beta_likelihood[:,0]),1:]
    
            output_binary = open(params.output_name+".betas.prob","wb")
            beta_params = array.array('i',[len(ID),params.beta_bins])
            beta_params.tofile(output_binary)
        
            beta_bins = numpy.linspace(params.beta_min,params.beta_max,params.beta_bins)
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
    
