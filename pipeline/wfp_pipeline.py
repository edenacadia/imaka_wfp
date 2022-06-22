## Eden McEwen
## Created: November 2021
## UH REU
# This file takes in a text file of dates, finds them on ehu, and then returns correlation outputs
## TODO: redid conf file to hold a bunch of variables, need to figure out the way to refactor these into 

import os
import sys
import fnmatch
import logging
import threading
import time
import concurrent.futures
from threading import Lock

lock = Lock()

sys.path.append('/home/emcewen/code_dev')
from pipeline.cor_pipeline import *
from pipeline.code.log_ex import *
from pipeline.code.file_reader import *

# Self-witten code
from pipeline.code.Correlator import *

data_path = "../"
out_path = "../"
target_path = "../"
pipe_f = "pipe_iteration"

# A class to store data for pipeline operations
class PipeData(object):
    ### Correlating options ###
    # Many of these are "false" for settign reasons
    tmax = 1000
    s_sub = False
    tt_sub = False
    ### Change len of AOCB ###
    trange = 0
    ### Processing options ###
    parallel_f = False # files in parallel
    parallel_c = False # wfs in parallel
    inj_sim = False #injected layer
    ### Potting options ###
    plotting = False
    sub_len = 0
    avg_sub = False
    med_sub = False
    mov_sub = False
    sig_clip= False
    ### datapaths ###
    data_path = "../"
    out_path = "../"
    target_path = "../"
    
    def __init__(self, conf_f="", run_f=""):
        ## initilize
        self.conf_f = ""
        self.run_f = ""
        #### Pipeline inputs (generated)
        self.names=[]   # names per each aocb
        self.d_files=[] # aocb files paths
        self.o_dirs=[]  # out dirs per file
        self.t_files=[] # target dirs per file
        self.t_ranges =[] # pairs of tranges
        ### Generating functions
        self.gen_conf(conf_f) # uses conf file to change variables
        self.gen_run(run_f) # uses run file to generate inputs
        if self.trange > 0:
            self.gen_trange(self.trange)
    
    def gen_conf(self, conf_new):
        ##### CONF FILE: sets variables #####
        # Read in configuration file
        if conf_new == "" or not os.path.isfile(conf_new):
            ## throw error and exit
            print("Given CONF file DNE: ", conf_new)
            return
        self.conf_f = conf_new
        conf_dict = read_file_dic(conf_new)
        self.set_vars(conf_dict) # Set the variables for pipeline
        #TODO: set the file type, d or e options
        
    def gen_run(self, run_new):
        ##### RUN FILE: made into list entries #####
        # look at infile, error if infile
        if self.inj_sim:
            # in the injected case, we don't need a runfile
            self.gen_files_inj()
            return
        if run_new == "" or not os.path.isfile(run_new):
            ## throw error and exit
            print("Given RUN file DNE: ", run_new)
            return
        self.run_f = run_new
        entries = read_file(run_new)
        # search for all needed files:
        self.gen_files(entries)
        
    def gen_trange(self, trange, num = -1):
        # set trange variable to a list of pairs
        self.trange = trange
        if trange == 0: # no range
            self.t_ranges = []
            return
        hdulist = fits.open(self.d_files[0])
        WFS_data = hdulist[3] #the slopes
        WFS_shape = WFS_data.header['NAXIS3']
        hdulist.close()
        # first, floor div by trange
        nums = WFS_shape // trange
        # if zero, return
        if nums == 0:
            return
        # cap at the specified nums
        if num != -1 and num > 0 and num <= nums:
            nums = num
        #else, iter through that amount, adding to lists
        t1_lst, t2_lst = [], []
        for i in range(nums):
            t1_lst.append(i*trange)
            t2_lst.append((i+1)*trange)
        #zip lists and set to the class variable
        self.t_ranges = [[t1_lst[i], t2_lst[i]] for i in range(nums)]
        
    #####################
    ### Manual inits: ###
    #####################
    
    def init_files(self, conf_f, run_f, trange=1000):
        ## Generating functions
        self.gen_conf(conf_f)
        self.gen_run(run_f)
        self.gen_trange(trange)
    
    def init_manual(self, date, data_path, out_dir, target_path):
        self.out_path = out_dir
        self.data_path = data_path
        self.target_path = target_path
        self.gen_files([date])   
    
    def init_date(self, date):
        self.gen_files([date]) 
        
    #####################
    ### generate      ###
    #####################
    
    def gen_files(self, entries):
        # make iterables for the pipeline
        # uses the date entries from the runfile read
        self.names, self.d_files, self.o_dirs, self.t_files = read_d(entries, self.data_path, self.out_path, self.target_path)
        
    def gen_files_inj(self):
        # make iterables for the pipeline
        # implified imput where all fits in datapath are meant to be read in
        # no targets given
        self.names, self.d_files, self.o_dirs, self.t_files = read_d_inj(self.data_path, self.out_path)
        
    def set_vars(self, conf_dict):
        #set data_path if available
        if conf_dict.get("data_path"): self.data_path = conf_dict["data_path"]
        #set out_path if available
        if conf_dict.get("out_path"): self.out_path = conf_dict["out_path"] 
        #set target_path if available
        if conf_dict.get("target_path"): self.target_path = conf_dict["target_path"] 
        #### Set Parallel vars ####
        if conf_dict.get("parallel"): self.parallel_c = eval(conf_dict["parallel"])
        if conf_dict.get("parallel_files"): self.parallel_f = eval(conf_dict["parallel_files"]) 
        #### Set Corr vars ####
        if conf_dict.get("tmax"): self.tmax = int(conf_dict["tmax"])
        if conf_dict.get("s_sub"): self.s_sub = eval(conf_dict["s_sub"])
        if conf_dict.get("tt_sub"): self.tt_sub = eval(conf_dict["tt_sub"])
        if conf_dict.get("trange"): self.trange = int(conf_dict["trange"])
        #### Set PLOT vars ####
        if conf_dict.get("plotting"): self.plotting = eval(conf_dict["plotting"])
        if conf_dict.get("sub_len"): self.sub_len = int(conf_dict["sub_len"])
        if conf_dict.get("avg_sub"): self.avg_sub = eval(conf_dict["avg_sub"])
        if conf_dict.get("med_sub"): self.med_sub = eval(conf_dict["med_sub"])
        if conf_dict.get("mov_sub"): self.mov_sub = eval(conf_dict["mov_sub"])
        if conf_dict.get("sig_clip"): self.sig_clip = eval(conf_dict["sig_clip"])
        #check to see if we can ignore runfile
        if conf_dict.get("inj_sim"): self.inj_sim = eval(conf_dict["inj_sim"])
    
    ###########################
    ###### PIPELINE INITS #####  
    ###########################
    
    def start_run(self, iter_lim=0):
        if self.parallel_f:
            self.pipe_run_par(iter_lim = iter_lim)
        else:
            self.pipe_run(iter_lim = iter_lim)
    
    def pipe_run(self, iter_lim = 0):
        """ 
        Applying Pipeline code
            input: lists of names, data paths, output dirs, target paths
            generates: whatever cor_pipeline function called
            output: none
        """
        print("============ STARTING RUN ============")
        start = time.time()
        # Limit the number of files for testing:
        d_files = self.d_files
        if iter_lim != 0:
            d_files = d_files[:iter_lim]
        for i, f in enumerate(d_files):
            print("==|==|==|==|==|== Starting file %s of %s ==|==|==|==|==|=="%(i + 1, len(d_files)))
            name = self.names[i]
            data_f = f
            out_d = self.o_dirs[i]
            target_f = self.t_files[i]
            self.log_iter(i, name, data_f, out_d, target_f)
            if not self.trange:
                self.pipe_iter(i, name, data_f, out_d, target_f, [0,0])
            else:
                # iter over all trange lengths
                for t1, t2 in self.t_ranges:
                    self.pipe_iter(i, name, data_f, out_d, target_f, [t1, t2])
            print("==|==|==|==|==|==> TIME SO FAR: %s"% (time.time() - start))
        print("============ RUN FINISHED IN: %s ============"%(time.time() - start))
    
    def pipe_run_par(self, iter_lim = 0):
        """ 
        Applying Pipeline code in parallel threads
        Uses executor to map 5 workers to each call of pipe iteration
            input: lists of names, data paths, output dirs, target paths
            generates: whatever cor_pipeline function called
            output: none
        """
        # TODO: make this run happily with class function
        print("============ STARTING THREAD RUN ============")
        print("BEGIN THREADS")
        d_files = self.d_files
        if iter_lim != 0:
            num_files = iter_lim
        else:
            num_files = len(d_files)
            
        names = self.names[:num_files]
        o_dirs = self.o_dirs[:num_files]
        t_files = self.t_files[:num_files]
        self_list = [self for i in range(num_files)]
        threaded_start = time.time()
        
        if not self.trange:
            tr_list = [[0,0] for i in range(num_files)]
            with concurrent.futures.ThreadPoolExecutor(max_workers=5) as executor:
                executor.map(self.pipe_iter, range(num_files), 
                             names, d_files, o_dirs, t_files, tr_list)
        else:
            # iter over all trange lengths
            for t1, t2 in self.t_ranges:
                with concurrent.futures.ThreadPoolExecutor(max_workers=5) as executor:
                    tr_list = [[t1,t2] for i in range(num_files)]
                    executor.map(self.pipe_iter, range(num_files), 
                                 names, d_files, o_dirs, t_files, tr_list)
        
        print("END THREADS")
        time_diff = time.time() - threaded_start
        print("THREAD Time: %s"%str(time_diff))
        
    ### PRINT FUNCTIONS ###

    def log_iter(self, i, name, data_f, out_d, target_f):
        """Standard logging format"""
        print("File %s: starting %s"% (i, name))
        print("   %s name: %s"% (i, name))
        print("   %s d_file: %s"% (i, data_f))
        print("   %s o_dirs: %s"% (i, out_d))
        print("   %s t_file: %s"% (i, target_f)) 
    
    def iter_print(self, name, t_range):
        # create an object
        #lock.acquire() # will block if lock is already held
        print("======== %s ========" % name)
        print("== tmax    = %s " % self.tmax)
        print("== sub_len = %s " % self.sub_len)
        print("== s_sub   = %s " % self.s_sub)
        print("== tt_sub  = %s " % self.tt_sub)
        if t_range[1] != 0:
            print("== t_range = %s " % t_range)
        #lock.release()
    
    ### MAIN PIPE FUNCTION ###
    
    def pipe_iter(self, i, name, data_f, out_d, target_f, t_range):
        """ 
        Correlates given class variables
        """
        try:
            self.iter_print(name, t_range)
            if not self.sub_len: sub_len=self.tmax
            # Create corr object
            curr_data = Correlator(name, data_f, out_d, 
                                   tmax=self.tmax, s_sub=self.s_sub, tt_sub= self.tt_sub)
            if target_f and not curr_data.target_file: 
                curr_data.set_target(target_f)
            if t_range[1] != 0:
                curr_data.set_trange(t_range)
            
            # Start correlations
            if self.parallel_c:
                # generate acor
                logging.info("====> Starting parallel acor")
                t0 = time.time()
                curr_data.acor_gen_par()
                t1 = time.time()
                logging.info("== finished in %s s"% str(t1-t0))
                # generate ccor
                logging.info("====> Starting parallel xcor")
                t2 = time.time()
                curr_data.ccor_gen_par()
                t3 = time.time()
                logging.info("== finished in %s s"% str(t3-t2))   
            else:
                # generate acor
                print("====> Starting acor")
                t0 = time.time()
                curr_data.acor_gen()
                t1 = time.time()
                print("== finished in %s s"% str(t1-t0))
                # generate ccor
                print("====> Starting ccor")
                t2 = time.time()
                curr_data.ccor_gen()
                t3 = time.time()
                print("== finished in %s s"% str(t3-t2))
                
            # create fits
            logging.info("====> writing fits")
            out = curr_data.fits_write()
            logging.info("== %s"% out)
            logging.info("====> Graphing")  
                
            if self.plotting:
                logging.info("====> Plotting")
                # Graph acor
                self.plot_all(curr_data, avg_sub = self.avg_sub, avg_len = self.sub_len)
                logging.info("===> complete")
            
        except Exception as e:
            print("Iteration %s error: %s"%(i, e))
        print("File %s: ending"%i)
    
    def plot_all(data, avg_sub=False, sub_len=200):
        g_out = data.acor_graph(t_list=[0,5,10,20,30], avg_sub=avg_sub, avg_len=sub_len)
        print(g_out)
        g_out = data.acor_animate_avg(dt_max=40, avg_sub=avg_sub, avg_len=sub_len)
        print(g_out)
        g_out = data.ccor_graph_all(avg_sub=avg_sub, avg_len=sub_len)
        print(g_out)
        g_out = data.cor_animate_all(dt_max=40, avg_sub=avg_sub, avg_len=sub_len) 
        print(g_out)
       
    
########################### PIPE ITERS #########################
########################### HELPER FUNCTS #########################

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

def plots_from_corr(DATE, active_wfs, suff='*stt.fits', avg_sub=False, sub_len=200):
    '''
    Takes a date, and generate plots given an ending
    inputs: DATE - a string, active_wfs - boolean string
    outputs: none (just saved files)
    '''
    fits_in = fits_file_pull(DATE, suff=suff)
    for p_file in fits_in:
        curr_data = cor.Correlator("", "", "", f_file=p_file)
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
    
    # init the pipe class
    pipe = PipeData(conf_file, run_file)
    
    ### Edits to variables, should be done in CONF file    
    #pipe.parallel_c = True
    #pipe.trange = True
    #pipe.tmax = 200
    #pipe.sub_len = 200
    #pipe.parallel_f = False
    #pipe.gen_trange(5000)
    
    ### Start the pipeline run ###
    pipe.start_run()
    #pipe.pipe_run()
    #pipe.pipe_run_par()
    # run an iteration all the way through
    # pipe.pipe_run(iter_lim = 1)
      
    ### Starts logger
    #logging.getLogger()
    #dir_in = run_file.replace(".txt", "_log/")
    #if not os.path.exists(dir_in):
    #    os.makedirs(dir_in)
    #config_root_logger()
    
    ## Runs Pipeline
    #start_run(conf_file, run_file)
    

    
    

    
    