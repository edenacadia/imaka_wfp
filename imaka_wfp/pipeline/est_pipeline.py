## estr_pipeline.py
## Eden McEwen
## Created: April 5th 2020
# a pipeline structure for the radial estimator code. 

import os
import sys
import fnmatch

import time
import numpy as np
import pandas as pd

# Personal code
from pipeline.code.file_reader import *
import pipeline.code.Estimator as est
import pipeline.code.Estimator_R as er

# Filepath defaults
# updated with fn set_global_paths
data_path = "../"
out_path = "../"
target_path = "../"
pipe_f = "pipe_iteration"
#parallel = False

########## Estimator Cluster
# preserving old file names
def df_iterator(df, xcor=False, sub_len=5, sig=3, thresh=3, c_size=4):
    return est.df_iterator(df, xcor=False, sub_len=5, sig=3, thresh=3, c_size=4)

def df_iterator_mult(df, sl_list, sig_list, thresh_list, c_list):
    return est.df_iterator_mult(df, sl_list, sig_list, thresh_list, c_list)


########### Estimator Radial

def df_iterator_estr(df, **kwargs):
    # TODO: set up **kwargs
    er.est_file(cor_f, method="meanshift",  text = True, plot = False)
    return True

########################### PIPELINE #########################
    
def pipe_run(pipe_funct, names, d_files, o_dirs, t_files):
    """ 
    Applying Pipeline code
        input: lists of names, data paths, output dirs, target paths
        generates: whatever cor_pipeline function called
        output: none
    """
    logging.info("============ STARTING %s ============")
    start = time.time()
    
    for i, f in enumerate(d_files):
        logging.info("==|==|==|==|==|== Starting file %s of %s ==|==|==|==|==|=="%(i + 1, len(d_files)))
        name = names[i]
        data_f = f
        out_d = o_dirs[i]
        target_f = t_files[i]
        
        #this is what change based on the type of pipeline you're trying to run
        pipe_funct(i, name, data_f, out_d, target_f)
        logging.info("==|==|==|==|==|==> TIME SO FAR: %s"% (time.time() - start))
        
    logging.info("============ RUN FINISHED IN: %s ============"%(time.time() - start))
    
########################### PIPELINE ITERS #########################    
    
def pipe_iteration(i, name, data_f, out_d, target_f):
    """ 
    Applies the file to data_proc_all
    Generates corr fits and corr graphs
    """
    thread_log_handler = start_thread_logging()
    try:
        logging.info("File %s: starting %s"% (i, name))
        logging.info("   %s name: %s"% (i, name))
        logging.info("   %s d_file: %s"% (i, data_f))
        logging.info("   %s o_dirs: %s"% (i, out_d))
        logging.info("   %s t_file: %s"% (i, target_f))
        data_proc_all(name, data_f, out_d)
        data_proc_all(name, data_f, out_d, sub_len=5)
        data_proc_all(name, data_f, out_d, s_sub=True)
        data_proc_all(name, data_f, out_d, s_sub=True, sub_len=5)
        data_proc_all(name, data_f, out_d, s_sub=True, tt_sub=True)
        data_proc_all(name, data_f, out_d, s_sub=True, tt_sub=True, sub_len=5)
    except Exception as e:
        logging.error("Iteration %s error: %s"%(i, e))
    logging.info("File %s: ending", i)
    stop_thread_logging(thread_log_handler)

######################

def set_global_paths(conf_dict):
    #set data_path if available
    if conf_dict["data_path"]:
        global data_path
        data_path = conf_dict["data_path"]
    #set out_path if available
    if conf_dict["out_path"]:
        global out_path
        out_path = conf_dict["out_path"]
    #set target_path if available
    if conf_dict["target_path"]:
        global target_path
        target_path = conf_dict["target_path"]

def start_run(conf_f, run_f):
    """Sets up the pipeline run from a configuration file and a list of date names
    --------
    conf_f : string, file.txt
        Contains file paths 
    run_f : string, file.txt
        Contains dates selected to be run in pipeline
    Outputs
    ----------
    darkname : string
        The name of the output dark FITS file.
    """
    
    ##### Pipeline inputs
    names=[]   # names per each file
    d_files=[] # aocb files paths
    o_dirs=[]  # out dirs per file
    t_files=[] # target dirs per file
    
    ##### CONF FILE: sets  #####
    # Read in configuration file 
    conf_d = read_file_dic(conf_f)
    print(conf_d)
    
    # Set the paths
    set_global_paths(conf_d)
    
    # Chose pipe_fn
    dispatcher = {'all':pipe_iteration, 'data':pipe_data_iteration, 'plots':pipe_plot_iteration}
    pipe_fn = pipe_iteration
    try:
        pipe_fn = dispatcher[conf_d["pipe_fn"]]
    except KeyError:
        logging.error("No valid function input: using all")
        pipe_fn = pipe_iteration
    
    ##### RUN FILE: made into list entries #####
    # look at infile, error if infile 
    entries = read_file(run_f)
    
    #search for all needed files:
    names, d_files, o_dirs, t_files = read_d(entries, data_path, out_path, target_path)
        
    ##### START PIPE
    # num_files = len(d_files)
    # logging.info("Number of files: %s", str(num_files))
    #logging.info('Name of files: %s', names)
    
    #### Applying Pipeline code ###
    pipe_run(pipe_fn, names, d_files, o_dirs, t_files)
    
    # BUG: waiting on making this parallel
    #if parallel:
        # parallel pipeline init
    #    pipe_run_parallel(pipe_fn, names, d_files, o_dirs, t_files)
    #else:
        # reg pipeline init
    #    pipe_run(pipe_fn, names, d_files, o_dirs, t_files)


###################################################################
########################### MAIN FUNCTION #########################
###################################################################


if __name__ == '__main__':
    """
    arg1: Conf file (txt file with paths, function, parallel option)
    arg2: text file (with it's extension)
    """
    conf_file = sys.argv[1]
    run_file = sys.argv[2]
    
    ## Starts logger
    #dir_in = run_file.replace(".txt", "_log/")
    #if not os.path.exists(dir_in):
    #    os.makedirs(dir_in)
    #config_root_logger(dir_in)
    
    ## Runs Pipeline
    start_run(conf_file, run_file)
    

    
    
############################## Class

class Est_pipe(object):
    
    def __init__(self, files, *args):
        self.out_fits = files
        self.est_obj = []
        self.df_total = pd.DataFrame()
    
    def run_file(self, file):
        # some sort of inputing on files
        er_pipe = er.Estimate_simple(file)
        self.est_obj.append(er_pipe)
        return er_pipe.return_table()
    
    def run_files(self):
        df_total = pd.DataFrame()
        for file in self.out_fits:
            table = self.run_file(file)
            df_total = pd.concat([df_total, table])
        self.df_total = df_total
        
    def save_table(self, ):
        if name == "":
            name = "0.txt"
        self.df_total.to_csv(r'c:\data\pandas.txt', index=None, sep=' ', mode='a')
        