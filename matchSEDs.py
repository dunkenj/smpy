import numpy
import array
import os, sys
import re
import time
import multiprocessing
import atpy

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

def galaxyFit(send_q, recv_q, printlock):
    for gal in iter(send_q.get, 'STOP'):
        j = numpy.argmin(numpy.abs(z-zobs[gal])) # Find closest model redshift


        fo = obs[gal,:]
        ferr = obs_err[gal,:]

        #mtoterr = numpy.log10(numpy.e)*fo[params.tot]/ferr[params.tot]

        fo[fo <= 0.] = 0.       # Set negative fluxes to zero
        #print fo
        I = numpy.where(ferr > 0.)[0] # Find bands with no observation
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
        if numpy.isinf(chimin) or numpy.isnan(minind) or len(fo) == 0:
            output_string = '{0} {1} {2} {3} {4} {5} {6} {7} {8} {9} \
            {10} {11} {12} {13} {14} {15} {16} {17} {18}'.format(gal+1,ID[gal],zobs[gal],
                                                                 -99,-99,-99,-99,-99,-99,
                                                                 -99, -99, -99, -99,-99,-99,-99,
                                                                 len(I),-99,'\n')
            printQueue.put(output_string)
            continue


        #Find the coordinate of the model with the bestfit mass
        tgi,tvi,ti,mi = numpy.unravel_index(minind,(n_tg,n_tauv,n_tau,n_ssp))
        Bestfit_Mass = numpy.log10(scale[tgi,tvi,ti,mi]*params.flux_corr)
        Bestfit_BMass = numpy.log10(Mass2[tgi,tvi,ti,mi]*params.flux_corr)
        Bestfit_SFR = (scale[tgi,tvi,ti,mi]*SFR[tgi,ti,mi]*params.flux_corr)
        Bestfit_Beta = beta[tgi,tvi,ti,mi]

        M_chisq_plus1 = numpy.log10(scale.flatten()[chisq < chimin+1])
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

        output_string = '{0} {1} {2} {3} {4} {5} {6} {7} {8} {9} {10} {11} {12} {13} {14} {15} {16} {17} {18}'.format(gal+1,ID[gal],zobs[gal],Bestfit_Mass,chimin,tgs,tvs,taus,mis, M_scaled[params.tot], MUV_scaled, Bestfit_Min, Bestfit_Max,minind,Bestfit_BMass,Bestfit_SFR,len(I),Bestfit_Beta,'\n')

        printlock.release()
        recv_q.put(output_string)

def getObservations(inputpath):
    input_data = atpy.Table(inputpath,type=params.input_format)

    column_names = input_data.columns.keys

    ID = input_data[params.ID_col]
    zobs = input_data[params.z_col]

    filts_used = params.filts_used

    k,l = 0,0
    for ii in range(len(column_names)):
        if column_names[ii].lower().endswith(params.flux_col_end):
            if k == 0:
                fluxes = input_data[column_names[ii]]
            else:
                fluxes = numpy.column_stack((fluxes,input_data[column_names[ii]]))
            k+=1

        if column_names[ii].lower().endswith(params.fluxerr_col_end):
            if l == 0:
                fluxerrs = input_data[column_names[ii]]
            else:
                fluxerrs = numpy.column_stack((fluxerrs,input_data[column_names[ii]]))
            l+=1

    tot_mag = -2.5*numpy.log10(fluxes[:,params.tot]) + 23.9

    fluxes = fluxes[:,filts_used]
    fluxerrs = fluxerrs[:,filts_used]

    return ID, zobs, tot_mag, fluxes, fluxerrs, k


if __name__ == '__main__':
    """
    SECTION 1A

    Loads in observed magnitudes/fluxes from catalog using ATPY
    From the catalog column names, it then picks out all the columns
    which it believes are magnitude/flux columns and their corresponding
    errors.

    """

    ID, zobs, tot_mag, obs, obs_err, filters_found = getObservations(params.input_catalog)

    """
    Section 1C

    Loads the various numpy arrays from the synthesised magnitude binary

    """

    #Load in synthesised data
    input_binary = params.synmag_output

    print "Loading synthetic mags and mass array:"
    with numpy.load(input_binary+'.main.npz') as synmags:
        parameters = synmags['parameters']
        Mshape = synmags['Mshape']
        z = synmags['z']
        tg, tv, tau, SSP =  parameters[1:5]

        nfilts = Mshape[0]
        n_tg = Mshape[2]
        n_tauv = Mshape[3]
        n_tau = Mshape[4]
        n_ssp = Mshape[5]

        MS = synmags['MS']
        MB = synmags['MB']
        MUV = synmags['MUV']
        SFR = synmags['SFR']

    with numpy.load(params.ssp_output) as beta_data:
        beta = beta_data['beta']

    calc_mode = params.calc_mode

    if nfilts != filters_found:
        f = numpy.load(input_binary+'.fluxes.npy')[params.filts_used]
        f_tot = numpy.load(input_binary+'.fluxes.npy')[params.tot]
    else:
        f = numpy.load(input_binary+'.fluxes.npy')


    print "Done."


    """
    SECTION 1C
    Setting up output table
    """
    if os.path.isfile("temp_output.txt"):
        os.remove("temp_output.txt")
    temp_file = open("temp_output.txt","w")

    """
    SECTION 2
    Chi-sq calculation

    """
    if calc_mode:
        print '{0:>4s} {1:>6s} {2:>6s} {3:>8s} {4:>6s}'.format('N','ID','Best', 'chimin', 'Mode')

    else:
        print '{0:>4s} {1:>6s} {2:>5s} {3:>6s} {4:>8s} {5:>10s} {6:>5s} {7:>12s} {8:>4s}'.format('N','ID','zobs','Best', 'chimin', 'tg', 'tauv','tau','SSP')

    start = time.time()
    ncpus = numpy.clip(params.ncpus,1,multiprocessing.cpu_count())

    send_q = multiprocessing.Queue()
    recv_q = multiprocessing.Queue()
    printlock = multiprocessing.Lock()
    for i in range( ncpus ):
        multiprocessing.Process( target = galaxyFit, args = ( send_q , recv_q, printlock ) ).start()


    # Put elements in the send queue for processing
    for gal in range( len(ID) ):
        send_q.put( gal )

    for gal in range( len(ID) ):
        temp_file.write( recv_q.get() )

    # Stop all the running processes
    for i in range( ncpus ):
        send_q.put( 'STOP' )

    # Close both send and receive queues
    send_q.close()
    recv_q.close()


    temp_file.close()
    #print results
    print "Time taken: "+str(time.time()-start)


    """
    Section 3
    Reload, format and save output table
    """
    while temp_file.closed == False:
        pause(0.1)

    data = numpy.loadtxt("temp_output.txt")
    try:
        rows, cols = data.shape
    except:
        cols = len(data)

    output = atpy.Table(name='Results')

    names = ['N','ID','z','Bestfit_Mass','Bestfit_chi2','Age','Dust_Tau','SFH_Tau','SSP_Number','H_rest', 'M1500','Bestfit_Min', 'Bestfit_Max','temp_index','Full_Mass','SFR','nfilts','Beta','MeanMass','MeanBeta','MeanM1500']
    units = [None,None,None,'log(Ms)',None,'Gyr',None,'Gyr',None, 'AB_mags', 'AB_mags','log(Ms)','log(Ms)',None,'log(Ms)','Ms/yr',None,None,'log(Ms)',None,'AB_mags']
    types = ['i4','i4','f4','f4','f4','f4','f4','f4','i4', 'f4', 'f4','f4','f4','i4','f4','f4','i4','f4','f4','f4','f4']
    for col in range(cols):
        output.add_column(names[col], data[:,col], unit=units[col], dtype=types[col])

    output.sort('ID')
    if os.path.isfile(params.output_name):
        os.remove(params.output_name)
    output.write(params.output_name,type=params.table_format)
