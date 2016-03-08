
SMpy Observe
============

User's Guide
============

The SMpy observation tools are intended to be flexible and easy to use,
allowing for a wide range of SED and observation planning uses. But...if
there are any problems, bugs or suggestions, please email ppxkd at
nottingham.

The primary python pre-requisites are as follows: \* numpy \* matplotlib
(for plotting) \* astropy v4.0+

The required Bruzual & Charlot models are included in the SSP folder
along with an example set of file filters and an EAZY filter set, the
usage of which will be explained later.

.. code:: python

    from mkCSPs import *
    from IPython.display import Math
    %pylab inline
     
    from pylab import *

.. parsed-literal::

    Populating the interactive namespace from numpy and matplotlib


Populating the interactive namespace from numpy and matplotlib

Building composite stellar populations
--------------------------------------

Stellar populations are built using the CSP class. A CSP object can be
initiated empty:

.. code:: python

    Galaxy1 = CSP()
And the desired SED can be built using the build function. ####
build(age, sfh, dust, metal, fesc=1.0, sfh\_law ='exp'): \* age - Time
since the onset of star-formation, in Gyr \* sfh - Star-formation
history parameter. \* For sfh\_law = 'exp', sfh is the exponentially
declining(+ve)/increasing(-ve) timescale in Gyr. \* For sfh\_law =
'pow', sfh is the power law exponent. \* dust - Dust extinction A\_V \*
metal - Bruzual & Charlot metalliticity index, 0 -> 5 \* fesc - Nebular
emission escape fraction \* 'sfh\_law' - Desired star-formation history
parametrisation

.. code:: python

    Galaxy1.build(0.5, 0.05, 1., 4) # A solar metallicity, 0.5 Gyr old stellar population with a short burst 
                                    # and 1 Av of dust extinction
Or, SEDs can be initiated directly:

.. code:: python

    Galaxy2 = CSP(0.1,1.5,0.,2, fesc=0.,sfh_law='pow') # A dust-free, 100 Myr old, t^1.5 SFH, 
                                                       # 0.2 Z_sol and strong nebular emission
Once a CSP object has been built, you can easily access the calculated
SED properties by:

.. code:: python

    print Galaxy1
    print ''
    print Galaxy2

.. parsed-literal::

               Age: 0.509 Gyr
           SFH Tau: 0.05 Gyr
          Dust Tau: 1 Av
               SFR: 1.17e-12 Ms/yr
      Stellar Mass: 1 Msol
              Beta: 4.4 
    
               Age: 0.102 Gyr
           SFH Tau: 1.5e-09 Gyr
          Dust Tau: 0 Av
               SFR: 2.94e-08 Ms/yr
      Stellar Mass: 1 Msol
              Beta: -2.6 


To access the wavelength and SED arrays, simply use:

.. code:: python

    print Galaxy1.SED
    print Galaxy1.wave

.. parsed-literal::

    [  0.00000000e+00   0.00000000e+00   0.00000000e+00 ...,   5.00384043e-12
       2.69361049e-12   1.57482822e-12]
    [  9.10000000e+01   9.40000000e+01   9.60000000e+01 ...,   1.20000000e+06
       1.40000000e+06   1.60000000e+06]


Manipulating and combining CSPs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Any

.. code:: python

    Galaxy1 = Galaxy1 * 1e11
    Galaxy2 = Galaxy2 * 1e9
    
    print log10(Galaxy1.Ms), Galaxy1.SFR
    print log10(Galaxy2.Ms), Galaxy2.SFR
    
    Galaxy3 = Galaxy1 + Galaxy2
    print log10(Galaxy3.Ms), Galaxy3.SFR
    
    Galaxy3 *= 1.5
    print log10(Galaxy3.Ms), Galaxy2.SFR

.. parsed-literal::

    11.0 0.117306215484
    9.0 29.4069594418
    11.0043213738 29.5242656573
    11.1804126328 29.4069594418


The CSP objects can be combined and normalised in more complicated ways:

.. code:: python

    Galaxy4 = ((Galaxy1 / Galaxy1.Ms) + (Galaxy2 / Galaxy2.Ms)) * 1e9
    print Galaxy4

.. parsed-literal::

               Age: 0.509 Gyr
           SFH Tau: 0.05 Gyr
          Dust Tau: 1 Av
               SFR: 29.4 Ms/yr
      Stellar Mass: 2e+09 Msol
              Beta: -2.6 


Plotting SEDs
^^^^^^^^^^^^^

.. code:: python

    loglog(Galaxy1.wave,Galaxy1.SED,label='Galaxy 1')
    loglog(Galaxy2.wave,Galaxy2.SED,label='Galaxy 2')
    loglog(Galaxy3.wave,Galaxy3.SED,label='Galaxy 1+2')
    
    ylim([1e3,1e9])
    xlim([500,5e4])
    xlabel(r'Wavelength [$\AA$]')
    ylabel(r'Flux [$L_{\odot} \AA^{-1}$]')
    Leg = legend(loc='lower right')


.. image:: output_22_0.png


--------------

Creating filter sets
--------------------

In order to do the synthetic photometry on the SEDs we have created, we
need to create a filter set. A *FilterSet* object is a collection of
individual *Filter* objects which can be incorporated or created in a
range of ways.

An empty *FilterSet* can be initiated as...

.. code:: python

    FiltSet1 = FilterSet()
One option for creating a filter is to generate a simple tophat filter
at a desired wavelength and width:

.. code:: python

    FiltSet1.addTophatFilter(1500,150) # 150 Angstrom wide tophat filter centered at 1500A
Filters can also be included from a library of commonly used filters
based on the format of EAZY (the latest EAZY filter compilation is
included in the github). The EAZY filter library must first be loaded
into its own class:

.. code:: python

    EAZYfilters = LoadEAZYFilters('FILTER.RES.V8')
The filters included in the set can be accessed through the
*.filternames* attribute, along with the corresponding central
wavelengths. E.g.

.. code:: python

    print EAZYfilters.filternames[125:127]
    print EAZYfilters.central_wlengths[125:127].to(U.micron)

.. parsed-literal::

    ['UKIDSS/J.txt', 'UKIDSS/K.txt']
    [ 1.25106196  2.20848924] micron


To add the desired EAZY filters to the filter set:

.. code:: python

    FiltSet1.addEAZYFilter(EAZYfilters,[125,126]) # where [125, 126] is any length list with indices 
                                                  # in the desired order
    print '\n'

.. parsed-literal::

    
    


Filter sets can also be generated directly from a set of filter files:

.. code:: python

    filt_dir = 'filters/filt_0*.txt'
    FiltSet2 = FilterSet(filt_dir)
Or, to add a filter file to an existing *FilterSet*:

.. code:: python

    FiltSet1.addFileFilter('filters/filt_01_bessell_U.txt')
\_ \_

Any of the filters added to a filterset can be accessed through the
*.filters* attribute. For example, plotting the filter responses:

.. code:: python

    for filt in FiltSet1.filters:
        semilogx(filt.wave,filt.response)
    for filt in FiltSet2.filters:
        semilogx(filt.wave,filt.response,'--')
        
    xlabel(r'Wavelength [$\AA$]')
    ylabel('Response')

.. parsed-literal::

    
    



.. image:: output_41_1.png


Making synthetic observations
-----------------------------

The final class is the *Observe* class which takes a *CSP* SED and
convolves it with a *FilterSet* at any desired redshift. The result is
an object which contains the resulting fluxes (in Jy) and magnitudes (in
AB).

.. code:: python

    Obs1 = Observe(Galaxy1, FiltSet1, 0.3) # Galaxy1 through FiltSet1 at z = 0.3
    
    print Obs1.AB
    print Obs1.fluxes

.. parsed-literal::

    [ 26.63315295  16.72814137  16.32333272  20.5608554 ] mag
    [  8.06749711e-08   7.39168505e-04   1.07316602e-03   2.16599694e-05] Jy
    
    


.. code:: python

    Obs2 = Observe(Galaxy2, FiltSet1, 1.5)
    
    print Obs2.AB # Note how the 1500A has no flux because the 
                  # Lyman break has redshift outside the filter
    print '\n'

.. parsed-literal::

    [         inf  21.61780511  22.58267031  22.03059338] mag
    
    


By default, the code will not 'observe' a galaxy which is older than the
age of the universe at the desired redshift. However, this can be
overridden with the keyword *force\_age = False*.

.. code:: python

    Galaxy5 = CSP(10.,0.5,0.,4) 
    Galaxy5 *= 5e10
    Obs3 = Observe(Galaxy5, FiltSet1, 2.) 
    Obs3 = Observe(Galaxy5, FiltSet1, 2., force_age=False)
    
    print Obs3.AB[2]

.. parsed-literal::

    SSP age older than universe...stopping.
    23.452067506 mag

