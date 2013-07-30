import numpy
import array
import os

import atpy

from sm_functions import dist
import sm_params as params

"""
SECTION 1A

Loads in observed magnitudes/fluxes from catalog using ATPY
From the catalog column names, it then picks out all the columns
which it believes are magnitude/flux columns and their corresponding
errors.

"""

#input_catalog = '/data/candels/catalogs/gs_all_tf_h_120607a_multi_highz_full.fits'

input_catalog = params.input_catalog
if input_catalog[-3:] == 'cat': 
    input_data = atpy.Table(input_catalog,type='ascii')
else:
    input_data = atpy.Table(input_catalog)

column_names = input_data.columns.keys

ID = input_data[params.ID_col]
zobs = input_data[params.z_col]

filts_used = params.filts_used

if params.fit_mode == "colours":
    i,j = 0,0
    for ii in range(len(column_names)):   
            if column_names[ii].lower().endswith(params.mag_col_end):
                    if i == 0:
                            mags = input_data[column_names[ii]]
                    else:
                            mags = numpy.column_stack((mags,input_data[column_names[ii]]))
                    i+=1
    
            if column_names[ii].lower().endswith(params.magerr_col_end):
                    if j == 0:
                            PSI = input_data[column_names[ii]]
                    else:
                            PSI = numpy.column_stack((PSI,input_data[column_names[ii]]))
                    j+=1
            
    tot_mag = mags[:,params.tot]

if params.fit_mode == "flux":
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
    i = k

	



"""
Section 1C

Loads the various numpy arrays from the synthesised magnitude binary

"""

#Load in synthesised data
input_binary = params.synmag_output

print "Loading synthetic mags and mass array:"
synmags = numpy.load(input_binary+'.main.npz')
parameters = synmags['parameters']
Mshape = synmags['Mshape']
z = synmags['z']
tg, tv, tau, SSP =  parameters[1:5]

nfilts = Mshape[0]
n_tg = Mshape[2]
n_tauv = Mshape[3]
n_tau = Mshape[4]
n_ssp = Mshape[5]

calc_mode = params.calc_mode

MS = synmags['MS']
MB = synmags['MB']
MUV = synmags['MUV']
SFR = synmags['SFR']
synmags.close()

beta_data = numpy.load(params.ssp_output)
beta = beta_data['beta']

if nfilts != i:
    M = numpy.load(input_binary+'.mags.npy')[filts_used]
else:
    M = numpy.load(input_binary+'.mags.npy')

if nfilts != i:
    f = numpy.load(input_binary+'.fluxes.npy')[filts_used]
    f_tot = numpy.load(input_binary+'.fluxes.npy')[params.tot]
    f /= f_tot
else:
    f = numpy.load(input_binary+'.fluxes.npy')
    f /= f[params.tot]
# Normalise template flux values to tot_mag filter fluxes


print "Done."


"""
SECTION 1C
Setting up output table
"""
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

for gal in range(len(ID)):
    j = numpy.argmin(numpy.abs(z-zobs[gal])) # Find closest model redshift

    # Turn MS into full mass
    Mass = MS[j,:]*numpy.power(10,-1*tot_mag[gal]/2.5)
    Mass2 = MB[j,:]*numpy.power(10,-1*tot_mag[gal]/2.5)
    cSFR = SFR[j,:]*numpy.power(10,-1*tot_mag[gal]/2.5)

    if params.fit_mode == "colours":
        psi = PSI[gal,:]
        mo = mags[gal,:]

        I = numpy.where(psi > -90)[0]

        I1 = I[:-1]
        I2 = I[1:]

        Co = mo[I1] - mo[I2] # Observed Colours
        psiqu = psi[I1]**2 + psi[I2]**2 # Errors
       
        Cm = M[I1,j,:]-M[I2,j,:] # Model colours

        # Calculate chi-sq
        chisq = 0.
        for i in range(len(Co)):
            chisq += (((Cm[i,:]-Co[i])**2)/psiqu[i])

        chimin,minind = numpy.nanmin(chisq), numpy.nanargmin(chisq)

        if numpy.isinf(chimin) or numpy.isnan(minind) or len(Co) == 0:
            continue

    elif params.fit_mode == "flux":
        fo = fluxes[gal,:]
        ferr = fluxerrs[gal,:]

        ferr /= fluxes[gal,params.tot] # Normalise fluxes and errors
        fo /= fluxes[gal,params.tot]   # to Total Mag column

        fo[fo <= 0.] = 0.       # Set negative fluxes to zero
        #print fo
        I = numpy.where(ferr > 0.)[0] # Find bands with no observation
        fo = fo[I]                    # and exclude from fit
        ferr = ferr[I]
        fm = f[I,j,:]
        #print fm[:,0,0,0,0]
        chisq = 0.
        for i in range(len(fo)):
            chisq += (((fm[i,:]-fo[i])**2)/(ferr[i])**2)

        chimin,minind = numpy.nanmin(chisq), numpy.nanargmin(chisq)
        if len(fo) == 0:
            continue
       
    isreal = numpy.isfinite(chisq)
    likelihood = numpy.exp(-0.5*chisq[isreal])
    wmean = numpy.sum(likelihood*Mass[isreal])/numpy.sum(likelihood)
    wmean = numpy.log10(wmean)
    wbeta = numpy.sum(likelihood*beta[isreal])/numpy.sum(likelihood)

    #Find the coordinate of the model with the bestfit mass
    zi,tgi,tvi,ti,mi = numpy.unravel_index(minind,(1,n_tg,n_tauv,n_tau,n_ssp))
    Bestfit_Mass=numpy.log10(Mass[tgi,tvi,ti,mi]*params.flux_corr)
    Bestfit_BMass=numpy.log10(Mass2[tgi,tvi,ti,mi]*params.flux_corr)
    Bestfit_SFR = (cSFR[tgi,tvi,ti,mi]*params.flux_corr)
    Bestfit_Beta = beta[tgi,tvi,ti,mi]

    M_chisq_plus1 = numpy.log10(Mass[chisq < chimin+1])
    Bestfit_Min, Bestfit_Max = numpy.min(M_chisq_plus1), numpy.max(M_chisq_plus1)
    
    #Scale the observed tot_mag band of the template to be the same as the observed tot_mag band of the galaxy
    #Convert the templates so they are no longer units of per stellar mass

    M_no_per_sm = M[params.tot,j,tgi,tvi,ti,mi]-(2.5*Bestfit_Mass)-(2.5*numpy.log10(1+zobs[gal]))
    Temp_flux = numpy.power(10,M_no_per_sm/-2.5)
    gal_totabs = tot_mag[gal] - dist(zobs[gal])
    Gal_flux = numpy.power(10,(gal_totabs/-2.5))
    flux_scale = Gal_flux/Temp_flux

    M_no_per_sm = (numpy.sum(M[params.tot,j][isreal]*likelihood)/numpy.sum(likelihood)) -(2.5*Bestfit_Mass)-(2.5*numpy.log10(1+zobs[gal]))
    wTemp_flux = numpy.power(10,M_no_per_sm/-2.5)
    wgal_totabs = tot_mag[gal] - dist(zobs[gal])
    wGal_flux = numpy.power(10,(wgal_totabs/-2.5))
    wflux_scale = wGal_flux/wTemp_flux


    M_scaled = M[:,0,tgi,tvi,ti,mi] - (2.5*Bestfit_Mass) - (2.5*numpy.log10(1+zobs[gal]))-2.5*numpy.log10(flux_scale)
    MUV_scaled = MUV[0,tgi,tvi,ti,mi] - (2.5*Bestfit_Mass) - (2.5*numpy.log10(1+zobs[gal]))-2.5*numpy.log10(flux_scale)

    wMUV_scaled = MUV[0,:,:,:,:] - (2.5*wmean) - (2.5*numpy.log10(1+zobs[gal])) - 2.5*numpy.log10(wflux_scale)
    wMUV_rest = numpy.sum(likelihood*wMUV_scaled[isreal])/numpy.sum(likelihood)

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
        
        Mass_sorted = Mass[chisq_order]
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

    if calc_mode:
        print '{0:4d} {1:6d} {2:>6.2f} {3:>8.1f} {4:>6.2f}'.format(gal+1,ID[gal],Bestfit_Mass,chimin, numpy.log10(Mode_Mass), '/n')
    else:
        print '{0:6d} {1:8f} {2:>5.2f} {3:>7.2f} {4:>8.1f} {5:>8.3f} {6:>5.1f} {7:>8.2f} {8:>3d} {9:>5.2f}'.format(gal+1,ID[gal],zobs[gal],Bestfit_Mass,chimin,tgs,tvs,taus,mis,wmean)


    temp_file.write('{0} {1} {2} {3} {4} {5} {6} {7} {8} {9} {10} {11} {12} {13} {14} {15} {16} {17} {18} {19} {20} {21}'.format(gal+1,ID[gal],zobs[gal],Bestfit_Mass,chimin,tgs,tvs,taus,mis, M_scaled[params.tot], MUV_scaled, Bestfit_Min, Bestfit_Max,minind,Bestfit_BMass,Bestfit_SFR,len(I),Bestfit_Beta, wmean, wbeta, wMUV_rest, '\n'))

                   
temp_file.close()


"""
Section 3
Reload, format and save output table
"""

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

#for col in range(cols):
#    output.add_column(names[col], [data[col],0], unit=units[col], dtype=types[col])

if os.path.isfile(params.output_name):
    os.remove(params.output_name)
output.write(params.output_name,type=params.table_format)
