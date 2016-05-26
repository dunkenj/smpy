import numpy as np
import array
import os, sys
import re
import time
import multiprocessing
import h5py
#import atpy
import logging
from astropy.table import Table, Column
from astropy import units as u

def fit_function(flux, models):
    return np.argmin(np.abs(models - flux))

class Fitter(object):
    def __init__(self, models, fluxes, nthreads=2):
        self.models = models
        self.fluxes = fluxes
        self.nthreads = nthreads
    
        self.args = [self.models]
        self.kwargs = {}
    
    def fit(self):
        self._fit = _function_wrapper(fit_function, self.args, self.kwargs)
        
        self.pool = multiprocessing.Pool(self.nthreads)
        self.results = self.pool.map(self._fit, self.fluxes)
        print('Done')
        
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
            
if __name__ == '__main__':
    models = np.arange(10)
    fluxes = np.random.random(100) * 10.
    
    fit = Fitter(models, fluxes)
    fit.fit()
    
