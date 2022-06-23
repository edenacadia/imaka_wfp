# Eden McEwen
# test_lengths_script.py
# June 2021

# imports
import numpy as np
import pandas as pd
import numpy.ma as ma
import importlib
import matplotlib
from astropy.stats import sigma_clipped_stats
from importlib import reload
import time
from astropy.io import fits

# ignore warnings
import warnings
warnings.filterwarnings("ignore")

#personal
from pipeline.est_pipeline import *
from pipeline.code.file_reader import *
from pipeline.code.corr_plots import *
from pipeline.code.cluster import *

import pipeline.code.Estimator as es
import pipeline.code.Estimator_R as er
import pipeline.code.data_table as d_t
import pipeline.code.graph_code as gc

import pipeline.code.Correlator as Cor


out_dir = "/home/emcewen/test_lengths/"
data_path = "/home/imaka/data/"

# aocb file:
date = '20210502'
aocb = 'aocb0012o'

########### Calculating the long circular buffers

def calc_par(date, aocb, tmax):
    # file things
    name = f'{date}_{aocb}'
    file_path = f'{data_path}{date}/ao/{aocb}.fits'
    ## a parallel correlation
    curr_data = Cor.Correlator(name, file_path, out_dir, tmax=tmax, s_sub=True, tt_sub=True)

    t0 = time.time()
    print("Starting parallel ACor")
    curr_data.acor_gen_par()
    t1 = time.time()
    print("... Finished in ", str(t1-t0))

    print("Starting XCor")
    curr_data.ccor_gen_par()
    t2 = time.time()
    print("... Finished in ", str(t2-t1))

    print("Writing Fits File")
    curr_data.fits_write()
    t3 = time.time()
    print("... Finished in ", str(t3-t2))
    
    
def date_iter_short(tmax):
    
    dates = [ '20180531',  '20180531',  '20180601',  '20180601', '20181219', '20181221', '20181221', '20181222', '20181223']
    aocbs = ['aocb0081o', 'aocb0033o', 'aocb0040o', 'aocb0130o','aocb0086o','aocb0156o','aocb0197o','aocb0133o','aocb0151o']
    
    for ex, i in enumerate(dates):
        calc_par(i, aocbs[ex], tmax)
        
def run12_iter_short(tmax):
    dates = ['20210429', '20210430', '20210502' ]
    aocbs = ['aocb0049o','aocb0031o','aocb0012o']
    
    for ex, i in enumerate(dates):
        calc_par(i, aocbs[ex], tmax)
        
################# Testing trange feature
# assumes that previous function has already been run

def calc_par_trange(date, aocb, tmax, tr1, tr2):
    # file things
    name = f'{date}_{aocb}'
    file_path = f'{data_path}{date}/ao/{aocb}.fits'
    out_dir_t = out_dir + "trange/"
    ## a parallel correlation
    curr_data = Cor.Correlator(name, file_path, out_dir_t, tmax=tmax, s_sub=True, tt_sub=True)
    curr_data.set_trange([tr1, tr2])
    
    print(f"{tr1} to {tr2}")
    t0 = time.time()
    print("Starting parallel ACor")
    curr_data.acor_gen_par()
    t1 = time.time()
    print("... Finished in ", str(t1-t0))

    print("Starting XCor")
    curr_data.ccor_gen_par()
    t2 = time.time()
    print("... Finished in ", str(t2-t1))

    print("Writing Fits File")
    curr_data.fits_write()
    t3 = time.time()
    print("... Finished in ", str(t3-t2))

        
def run12_trange(tmax):
    # files
    dates = ['20210429', '20210430', '20210502' ]
    aocbs = ['aocb0049o','aocb0031o','aocb0012o']
    
    # iter on files
    for i, d in enumerate(dates):
        date = d
        aocb = aocbs[i]
        # find the cb length
        file_path = f'{data_path}{date}/ao/{aocb}.fits'
        hdu_list = fits.open(file_path)
        tr_max = hdu_list[0].header["NAXIS2"]
        hdu_list.close()
        # send into the recursion funct.
        run12_trange_help(date, aocb, tmax, 0, tr_max)
    
    
def run12_trange_help(date, aocb, tmax, start, stop):
    if stop-start > 2*tmax:
        mid_t = stop//2
        calc_par_trange(date, aocb, tmax, start, mid_t)
        calc_par_trange(date, aocb, tmax, mid_t, stop)
        run12_trange_help(date, aocb, tmax, start, mid_t)
        run12_trange_help(date, aocb, tmax, mid_t, stop)
        return True
    return False


######### Long subs, same aocb length
### 2000 length subtractions
### redoing to get the right order of pre-subtractions
### keeping all aocbs in the same length (short)


def same_length_iter():
    
    dates = [ '20180531',  '20180531',  '20180601',  '20180601', '20181219', '20181221', '20181221', '20181222', '20181223', '20210429', '20210430', '20210502']
    aocbs = ['aocb0081o', 'aocb0033o', 'aocb0040o', 'aocb0130o','aocb0086o','aocb0156o','aocb0197o','aocb0133o','aocb0151o', 'aocb0049o','aocb0031o','aocb0012o']
    
    tmax = 2000
    tr1 = 0
    tr2 = 5000
    
    for i, d in enumerate(dates):
        date = d
        aocb = aocbs[i]
        calc_par_trange(date, aocb, tmax, tr1, tr2)