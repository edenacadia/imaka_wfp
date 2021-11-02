# Eden McEwen
# test_timer.py
# October 2021

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

#import pipeline.code.Estimator as es
#import pipeline.code.Estimator_R as er
#import pipeline.code.data_table as d_t
#import pipeline.code.graph_code as gc

import pipeline.code.Correlator as Cor



# Out directory specific for this test
out_dir = "/home/emcewen/test_timer/test2/"
data_path = "/home/imaka/data/"
target_path = "home/imaka_wfp/inputs/targets/"


####### This function takes in a range and only calculates the correlation in that range

def calc_par_trange(name, tmax, tr2, tr1=0):
    # file things
    aocb = name[9:]
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
    
    
########### Timing test   
if __name__ == '__main__':
    # finds the first ten aocbs in that date
    date = sys.argv[1]
    
    names, d_files, o_dirs, t_files = read_d([date], data_path, out_dir, target_path)
    
    # for an itterations on tmax length
    range_lst = [100, 1000, 2000, 5000, 10000, 15000, 20000, 25000]
    tmax_lst = [10, 50, 100, 200, 500, 1000]
    
    #test one: range, tmax, file
    #test two: tmax, range, file
    
    for tmx in tmax_lst:
        for r in range_lst:
            for n in names[:10]:
                print("STARTING:", n, r, tmx)
                #print(n[9:])
                calc_par_trange(n, tmx, r, tr1=0)