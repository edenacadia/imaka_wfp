### test_sim.py
# for use in test_simulation.ipynb
# Eden McEwen
# November 2021

# February 2022
# edited to include the most recent simulations of correlations

# imports
import numpy as np
import pandas as pd
import numpy.ma as ma
import multiprocessing as mp
import importlib
import matplotlib
from astropy.stats import sigma_clipped_stats
from importlib import reload
import time
from astropy.io import fits
import re
import itertools

# ignore warnings
import warnings
warnings.filterwarnings("ignore")

#personal
#from pipeline.est_pipeline import *
from pipeline.code.file_reader import *
from pipeline.code.corr_plots import *
from pipeline.code.cluster import *

import pipeline.code.Estimator as es
import pipeline.code.Estimator_R as er
import pipeline.code.data_table as d_t
import pipeline.code.graph_code as gc

import pipeline.code.Correlator as Cor
from importlib import reload

import t_const_code as tc

# Out directory specific for this test
out_dir = "/home/emcewen/test_simul/"
data_path = "/home/imaka/data/"
target_path = "home/imaka_wfp/inputs/targets/"

# for test injections:
out_dir = "/home/emcewen/out/injsim/" 
dates = ["20180531", "20180601", "20181219", "20181221", "20181222", "20181223", "20210429", "20210430", "20210501", "20210502"]
files = ["0069", "0072", "0127", "0323", "0133", "0181", "0013", "0023", "0025", "0076"]

file_opdir_idx = [0,0,0,2,2,2,3,3,2,2]

r0s = ["0.10", "0.15", "0.20", "0.25" ]
vspds = ["001", "003", "006", "012", "018", "024"]
vdirs = ["0", "1", "2", "3"]

def inj_file(date_i, r0_j, vspd_k, vdir_m):
    fmt = "{out_dir}fits/{date}_aocb{file}o+injsim.r0_{r0}.vspd_{vspd}.0.vdir_{vdir}_tmax1000_tts.fits"
    return fmt.format(out_dir=out_dir, date=dates[date_i], 
                      file=files[date_i], r0=r0s[r0_j], 
                      vspd=vspds[vspd_k], vdir=vdirs[vdir_m]) 

def inj_name(date_i, r0_j, vspd_k, vdir_m):
    fmt_name = "{date}_aocb{file}o+injsim.r0_{r0}.vspd_{vspd}.0.vdir_{vdir}"
    return fmt_name.format(date=dates[date_i], 
                      file=files[date_i], r0=r0s[r0_j], 
                      vspd=vspds[vspd_k], vdir=vdirs[vdir_m]) 

######### ESTIMATING #############


def iter_est_par(out_file="out.txt", save=True, cpu_count=8):
    """
    Iterate through a list of parameters, from begining 
    uses parallel on each file
    currently not customizable
    """
    
    date_idxs = [0, 4, 6] # pre-selected
    r0s_idxs = [1,2]
    vspd_idxs = [1,3,5]
    
    # generate a file list
    file_idxs = list(itertools.product(date_idxs, r0s_idxs, vspd_idxs))
    
    cor_files = [inj_file(i, j, k, file_opdir_idx[i]) for i,j,k in file_idxs]
    cor_names = [inj_name(i, j, k, file_opdir_idx[i]) for i,j,k in file_idxs]
    N_files = len(cor_files)
    
    # Start the pool
    pool = mp.Pool(cpu_count)
    results_async = []
    print(f'Estimate_R in parallel with {cpu_count} cores.')
    
    for ii in range(N_files):
        cor_file = cor_files[ii]
        cor_name = cor_names[ii]
        
        results = pool.apply_async(file_run, (cor_file, cor_name))
        results_async.append(results)

    pool.close()
    pool.join()

    df_sum = pd.DataFrame()
    for ii in range(N_files):
        try:
            results = results_async[ii].get()
            #print(results)
            # add file parameters
            n = results.shape[0]
            date_ii = dates[file_idxs[ii][0]]
            r0_ii = r0s[file_idxs[ii][1]]
            vspd_ii = vspds[file_idxs[ii][2]]
            vdir_ii = vdirs[file_opdir_idx[file_idxs[ii][0]]]
            results = results.assign(date=[date_ii]*n, r0=[r0_ii]*n, vspd=[vspd_ii]*n, vdir=[vdir_ii]*n)
            # add to larger table
            df_sum = pd.concat([df_sum, results])
        except exception as e:
            print("Error with file: ", cor_files[ii], e)
    if save:
        df_sum.to_csv(out_file, sep='\t')
        print(f"Saved to: {out_file}")
    
    return df_sum


def file_run(cor_file, name):
    """
    Runs multiple parameters of the estimator on a single file
    """
    pid = mp.current_process().pid
    df_file = pd.DataFrame()
    ## CHANGE: choice of estimator iterations
    wind_sub = [5, 10, 25, 50, 100, 200]
    d_clip = [2,3,4,5,6]
    est_params = file_idxs = list(itertools.product(wind_sub, d_clip))
    ## iteration to estimator
    for w, d, in est_params:
        t0 = time.time()
        print(f"p{pid}: {name[:18]} => START, sub_len:{w} detection_clip:{d}")
        
        df_tmp = single_est(cor_file, name, sub_len=w, detect_clp=d, 
                            bg_sub=True, move_sub=True, est_len=250)
        df_file = pd.concat([df_file, df_tmp])
        
        t1 = time.time()
        print(f"p{pid}: {name[:18]} ==> END, %s s"% str(t1-t0))
    return df_file
    

def single_est(cor_file, name, **kwargs):
    """
    Takes in a file and various kwargs
    """
    df_test = pd.DataFrame()

    er_pipe = er.Estimate_simple(cor_file)
    er_pipe.update_wfs([True, True, True, False, False])
    er_pipe.update_params(**kwargs)
    
    # Estimatior call 
    table_tmp = er_pipe.return_table(sdv_cnd = 2, n_filter = 10)
    
    #add this table to a list of tables 
    n = table_tmp.shape[0]
    table_tmp = table_tmp.assign(est_len=[er_pipe.est_len]*n, bg_sub=[er_pipe.bg_sub]*n, sub_len=[er_pipe.sub_len]*n)       
    df_test = pd.concat([df_test, table_tmp])
    
    return df_test
    

######### CORRELATING ############

def calc_par(name, file, tmax):
    # file things
    ## a parallel correlation
    curr_data = Cor.Correlator(name, file, out_dir, tmax=tmax, s_sub=False, tt_sub=False)
    curr_data.set_trange([0, 1024])
    
    print(f"FILE: {name}")
    t0 = time.time()
    print("Starting parallel ACor")
    curr_data.acor_gen_par()
    t1 = time.time()
    print("... Finished in ", str(t1-t0))

    print("Starting parallel XCor")
    curr_data.ccor_gen_par()
    t2 = time.time()
    print("... Finished in ", str(t2-t1))

    print("Writing Fits File")
    curr_data.fits_write()
    t3 = time.time()
    print("... Finished in ", str(t3-t2))
    return curr_data



if __name__ == '__main__':
    
    name1 = sys.argv[1]
    file1 = sys.argv[2]
    
    print("hello")
    
    
# Old functions:
# pre multi simulation run (-2.28.22)


def run1():
    file1 = '/data/emcewen/sim/20211101/simul.20211101_203350.aocb.fits'
    name1 = 'simul.20211101_203350'
    cor1 = calc_par(name1, file1, 1000)
    
    file2 = '/data/emcewen/sim/20211102/simul.20211102_105221.aocb.fits'
    name2 = 'simul.20211102_105221'
    cor2 = calc_par(name2, file2, 1000)

def run2():
    file1 = '/data/emcewen/sim/20211118/20211118_143311/simul.20211118_143311.aocb.fits'
    name1 = 'simul.20211118_143311'
    cor1 = calc_par(name1, file1, 500)
    
    file2 = '/data/emcewen/sim/20211118/20211118_144811/simul.20211118_144811.aocb.fits'
    name2 = 'simul.20211118_144811'
    cor2 = calc_par(name2, file2, 500)