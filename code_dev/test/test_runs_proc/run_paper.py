# run_paper.py
## Eden McEwen
## organizational Code for the paper

#### inputs ####
import os
import fnmatch
import numpy as np

import pipeline.code.Correlator as cor
import pipeline.code.Estimator_R as er


def run7():
    dates = ["20180531", "20180601"]
    active_wfs = [True, True, True, False, False]
    return dates, active_wfs
    
def run9():
    dates = ["20181219", "20181221", "20181222", "20181223"]
    active_wfs = [True, True, True, False, True]
    return dates, active_wfs
    
def run12():
    dates = ["20210429", "20210430", "20210501", "20210502"]
    active_wfs = [True, True, True, False, False]
    return dates, active_wfs

def run13():
    dates = ["20210827", "20210828", "20210829", "20210830"]
    active_wfs = [True, True, True, False, False]
    return dates, active_wfs

def run15():
    dates = ["20220120", "20220121", "20220122", "20220123", "20220124", "20220125"]
    active_wfs = [True, True, True, True, True]
    return dates, active_wfs

    
def fits_file_pull(DATE, suff='*stt.fits'):
    # pulling all fits given a date
    dir_pre = f"/home/emcewen/out/{DATE}/fits/"
    fits_in = []
    #check this directory is valid
    if os.path.isdir(dir_pre):
        files = os.listdir(dir_pre)   # list all files
        # collect fits with right suffix
        fits_in = [dir_pre + fn for fn in files if fnmatch.fnmatch(fn, suff)]
    fits_in.sort()
    return fits_in


################## CALLABLE FUNCTIONS ##################
    
def plots_from_corr(DATE, active_wfs, suff='*tts.fits', avg_sub=False, sub_len=200):
    '''
    Takes a date, and generate plots given an ending
    inputs: DATE - a string, active_wfs - boolean string
    outputs: none (just saved files)
    '''
    fits_in = fits_file_pull(DATE, suff=suff)
    for p_file in fits_in:
        curr_data = cor.Correlator("", "", "", f_file=p_file)
        curr_data.active_wfs = active_wfs
        curr_data.n_wfs = sum(active_wfs)
        print("======== %s ========" % curr_data.name)
        plot_all(curr_data, avg_sub=True, sub_len=5)
        plot_all(curr_data, avg_sub=True, sub_len=200)
        print("===> complete")

def plots_from_corr_files(fits_in, active_wfs, suff='*tts.fits', avg_sub=False, sub_len=200):
    '''
    Takes a date, and generate plots given an ending
    inputs: DATE - a string, active_wfs - boolean string
    outputs: none (just saved files)
    '''
    for p_file in fits_in:
        curr_data = cor.Correlator("", "", "", f_file=p_file)
        curr_data.active_wfs = active_wfs
        curr_data.n_wfs = sum(active_wfs)
        print("======== %s ========" % curr_data.name)
        plot_all(curr_data, avg_sub=True, sub_len=5)
        plot_all(curr_data, avg_sub=True, sub_len=200)
        print("===> complete")
        

def plot_all(data, avg_sub=False, sub_len=200):
    g_out = data.acor_graph(t_list=[0,5,10,20,30], avg_sub=avg_sub, avg_len=sub_len)
    print(g_out)
    g_out = data.acor_animate_avg(dt_max=40, avg_sub=avg_sub, avg_len=sub_len)
    print(g_out)
    g_out = data.ccor_graph_all(avg_sub=avg_sub, avg_len=sub_len)
    print(g_out)
    g_out = data.cor_animate_all(dt_max=40, avg_sub=avg_sub, avg_len=sub_len) 
    print(g_out)
    
def est_plots_from_corr(DATE, active_wfs):
    '''
    Given that a aocb has already been correlated, this generates plots of interest
    '''
    
    #TODO: setting active wavefronts
    
    # pulling all fits given a date
    dir_pre = f"/home/emcewen/out/{DATE}/fits/"
    #check this directory is valid
    if os.path.isdir(dir_pre):
        # list all files
        files = os.listdir(dir_pre)
        # collect fits with right suffix
        fits_in = [dir_pre + fn for fn in files if fnmatch.fnmatch(fn, '*stt.fits')]
    fits_in.sort()
    
    # estimating part
    df_ms = pd.DataFrame()
    for fts in fits_in[:]:
        print(fts)
        er_pipe = er.Estimate_simple(fts)
        table = er_pipe.return_table(sdv_cnd = 2, n_filter = 10)
        df_ms = pd.concat([df_ms, table])
    df_ms["method"] = "meanshift"
    
    # TODO: save instead of printing
    date_name = f"{DATE} mean shift"
    er.plot_clstr_table(df_ms[df_ms["dir_std"] < 50 ], min_count=1, name=date_name).show
    
    # TODO: save instead of printing
    plt.subplot(2,1,1)
    plt.hist(df_ms["spd"], bins = 25)
    plt.title(f"{DATE} spd")

    plt.subplot(2,1,2)
    plt.hist(df_ms["dir"], bins = 25)
    plt.title(f"{DATE} dir")
