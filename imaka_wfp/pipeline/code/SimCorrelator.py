## Eden McEwen
## December 25, 2021
# SimCorrelator
# This code takes in an imaka formated aocb file
# Will return auto correlations (acor) or cross correlations (ccor)
# as of feb 2022 this isn't really used, taking in files pre=simulated

import math
import imageio
import time
import os
import sys
import itertools

import pandas as pd
import numpy as np
from scipy.io import readsav
import multiprocessing as mp
import functools as fct

import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib as mpl
import matplotlib.animation as animation

from astropy.table import Table
from astropy.stats import sigma_clip
from astropy.io import fits

# Self-witten code
from pipeline.code.corr_code import *
from pipeline.code.corr_plots import *
from pipeline.code.CorrClock import CorrClock
from pipeline.code.Correlator import Correlator


mask_8_8_center = [[0,0,1,1,1,1,0,0],
           [0,1,1,1,1,1,1,0],
           [1,1,1,1,1,1,1,1],
           [1,1,1,0,0,1,1,1],
           [1,1,1,0,0,1,1,1],
           [1,1,1,1,1,1,1,1],
           [0,1,1,1,1,1,1,0],
           [0,0,1,1,1,1,0,0]]

mask_8_8_agg =  [[0,0,0,1,1,0,0,0],
           [0,0,1,1,1,1,0,0],
           [0,1,1,0,0,1,1,0],
           [1,1,0,0,0,0,1,1],
           [1,1,0,0,0,0,1,1],
           [0,1,1,0,0,1,1,0],
           [0,0,1,1,1,1,0,0],
           [0,0,0,1,1,0,0,0]]

############################################################## 
####################### Correlator Simulation Class #####################
##############################################################
        
class SimCorrelator(Correlator):
    def __init__(self, name, sim_f, data_f, out_d, sim_spd = 0, sim_dir = 0, f_file=None, tmax=200, s_sub=False, tt_sub=False):
        self.sim_file = sim_f
        self.sim_spd = sim_spd
        self.sim_dir = sim_dir
        self.sim_r0 = sim_r0
        # doing correlation
        super(SimCorrelator, self).__init__(name, data_f, out_d, f_file=f_file,
                                           tmax=tmax, s_sub = s_sub, tt_sub=tt_sub)
            
    # overwrite: x and y slope generations
    # file naming (include a sim regerence)
    # fits file generation (include sim file and strength)
    
    def slopes_gen(self):
        """
        If data file exists, pulls x and y slopes from datafile
        If fits is formatted differently, you'll want to change this file
            args: N/A
            return: T/F (whether or not slopes were updated, and thus data valid)
        """
        try:
            ### DATA pull
            if self.data_check(self.data_file):
                self.x_slopes, self.y_slopes, self.n_wfs, self.date = self.slopes_pull(self.data_file)
                if self.trange is not None:
                    self.x_slopes, self.y_slopes = self.slopes_trange_apply(self.x_slopes , self.y_slopes)
            else:
                print("No valid data")
                return False
            ### SIM pull
            if self.data_check(self.sim_file):
                x_sim, y_sim, _, _ = self.slopes_pull(self.data_file)
                # TODO: add ro scaling here?
                # TODO: add dir shifting here
                # TODO: add a check for lengths of the files?
                np.rot90(x_sim, k=1, axes=(0, 1))
                
                if self.trange is not None:
                    x_sim, y_sim = self.slopes_trange_apply(x_sim, y_sim)
                if self.sim_dir is not 0:
                    x_sim = np.rot90(x_sim, k=self.sim_dir, axes=(2, 3))
                    y_sim = np.rot90(y_sim, k=self.sim_dir, axes=(2, 3))
                self.x_slopes += x_sim 
                self.y_slopes += y_sim 
                # BUG: what if not the right number of wavefronts?
            if self.tt_sub:
                self.x_slopes, self.y_slopes = self.get_tiptilt_sub()
            if self.s_sub:
                self.x_slopes, self.y_slopes = self.get_statics_sub()
            return True
        except Exception as e: 
            print(e)
            return False
        
    def slopes_pull(self, file):
        hdulist = fits.open(file)
        WFS_data = hdulist[3] #the slopes
        WFS_list = np.arange(WFS_data.header['NAXIS2'])
        half_slopes = WFS_data.header['NAXIS1']//2
        WFS_shape = (WFS_data.header['NAXIS3'],8,8)
        x_wfs_slopes = np.array([np.array(WFS_data.data[:,wfs, :half_slopes]).reshape(WFS_shape) for wfs in WFS_list])
        y_wfs_slopes = np.array([np.array(WFS_data.data[:,wfs, half_slopes:]).reshape(WFS_shape) for wfs in WFS_list])
        n_wfs = hdulist[0].header['NWFS'] 
        date = hdulist[0].header['OBSDATE']
        hdulist.close()
        return x_wfs_slopes, y_wfs_slopes, n_wfs, date
    
    def slopes_trange_apply(self, x_slopes, y_slopes):
        x_slopes = x_slopes[:, self.trange[0]:self.trange[1], :, :]
        y_slopes = y_slopes[:, self.trange[0]:self.trange[1], :, :]
        return x_slopes, y_slopes
    
    ################ Naming helpers #############
    
    def fits_name_gen(self):
        """
        Gives back the fits file name associated with parameters
            args: none
            returns: string (fits file)
        """
        out_dir = self.out_dir + "fits/"
        out_file = f"sim{self.name}_spd{self.sim_spd}_dir{self.sim_dir}_tmax_{self.tmax}"
        # generating file structure more finely
        if self.s_sub and self.tt_sub:
            out_file = out_file + "_tts"
        elif self.tt_sub:
            out_file = out_file + "_tt"
        elif self.s_sub:
            out_file = out_file + "_s"
        else:
            out_file = out_file + "_raw"
        if self.trange is not None:
            out_file = out_file + f"_f{self.trange[0]}_{self.trange[1]}"
        out_file = out_file + ".fits"
        #checking to make sure this location exists
        path = self.path_check(out_dir, loud=False)
        if path: 
            fits_file = out_dir + out_file
            logging.debug("fits name: %s"% fits_file)
            return out_dir + out_file
        else:
            logging.warning("ERROR: fits path not created: %s"% out_dir)
            return out_file
    
    
    def plot_file_gen(self, plot_dir, plot_type, file_type, med_sub=False, avg_sub=False, avg_len=0):
        """
        Generates a file name for a plot, renumbers if name already taken
            args: plot_type (varies by plot func), file_type ("png" or "gif"), Graphing specifications
            retuns: out_file (descriptive string)
        """
        out_dir = f"{self.out_dir}plots/{plot_dir}/"
        #special naming bc sim correlator
        out_file = f"sim{self.name}_spd{self.sim_spd}_dir{self.sim_dir}"
        
        # generating file structure more finely
        if self.s_sub and self.tt_sub:
            out_file = out_file + "_tts"
        elif self.tt_sub:
            out_file = out_file + "_tt"
        elif self.s_sub:
            out_file = out_file + "_s"
        else:
            out_file = out_file + "_raw"
        out_file = out_file + plot_type
        #check if path exists
        path = self.path_check(out_dir, loud=False)
        # Renumber if repeats 
        out_base = out_dir + out_file
        i = 0
        while os.path.exists(out_base + "_%s.%s" % (i, file_type)):
            i += 1
        out = out_base + "_%s.%s" % (i, file_type)
        return out