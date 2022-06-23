### Def Estimate R
## Eden McEwen
## February 2, 2020

import math
import time
import os
import sys
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
import matplotlib
import statistics as stat
from scipy import signal
from astropy.io import fits
import pandas as pd
from astropy.stats import sigma_clipped_stats
from celluloid import Camera

# Windplot for detections
from windrose import WindroseAxes

# Using Clustering
import hdbscan
import seaborn as sns
from numpy import unique
from numpy import where
from sklearn.datasets import make_classification
from sklearn.cluster import MeanShift
from sklearn.cluster import DBSCAN
from sklearn.cluster import AffinityPropagation

# Self-witten code
from pipeline.code.tpoint import TPoint
from pipeline.code.corr_code import *
from pipeline.code.Correlator import *
from pipeline.code.est_plots import *
from pipeline.code.cluster import *

#generating a mask for corr data
mask_data = mask_8_8_center
mask_cor = signal.correlate2d(mask_data, mask_data)
mask = ~np.array(mask_cor, dtype=bool)


#### Static code ####
#Cov map radii
rad_map = np.array([[np.sqrt(i**2 + j**2) for j in np.arange(-7, 8)] for i in np.arange(-7, 8)])
rad_map_boolean = np.zeros((8, 15, 15))
# boolean mask matrix
for i in np.arange(0, 8):
    rad_max = ma.masked_where(rad_map <=  i + 0.5, rad_map)
    rad_min = ma.masked_where(rad_map >  i - 0.5, rad_map)
    rad = rad_max.mask & rad_min.mask
    rad_map_boolean[i-1] = rad
rad_map_inv = np.ones(rad_map_boolean.shape) - rad_map_boolean

#########################
#### Estimator Class ####
#########################

class Estimate_simple(object):
    # background subtraction variables
    bg_sub = False # by default, no bg subs
    sub_len = 200
    med_sub = False
    move_sub = False
    sig_clip = True
    # params
    clstr_method = "meanshift"
    detect_clp = 4
    sdv_comp = 3
    sdv_cnd = 0    # by default no condensing
    n_filter = 0   # by default no number filtering
    
    def __init__(self, file, **kwargs):
        self.out_fits = file
        self.data = Correlator("", "", "", f_file = file)
        self.name = self.data.name
        self.date = self.data.date
        self.hz = fits.open(file)[1].header["FSAMPLE"] # might error for injections?
        self.est_len = self.data.tmax # length of the estimation
        # these will change each run
        self.table = pd.DataFrame()
        # auto correlation 
        self.a_dtct = None
        self.a_spds = None
        self.a_dirs = None
        self.a_xs = None
        self.a_ys = None
        self.a_cstr = None
        # cross correlation
        self.x_dtct = None
        self.x_spds = None
        self.x_dirs = None
        self.x_xs = None
        self.x_ys = None
        self.x_cstr = None
        # estimation params
        self.update_params(**kwargs)
    
    def update_params(self, **kwargs):
        """
        Updates the tuneable parameters based on args given
        """
        if "clstr_method" in kwargs: 
            self.clstr_method = kwargs.get("clstr_method")
        if "detect_clp" in kwargs: 
            self.detect_clp = kwargs.get("detect_clp")
        if "sdv_comp" in kwargs: 
            self.sdv_comp = kwargs.get("sdv_comp")
        if "sdv_cnd" in kwargs: 
            self.sdv_cnd = kwargs.get("sdv_cnd")
        if "n_filter" in kwargs: 
            self.n_filter = kwargs.get("n_filter")
        if "hz" in kwargs: 
            self.hz = kwargs.get("hz")
        if "est_len" in kwargs: 
            self.est_len = kwargs.get("est_len")
        if "bg_sub" in kwargs: 
            self.bg_sub = kwargs.get("bg_sub")
        if "sub_len" in kwargs: 
            self.sub_len = kwargs.get("sub_len")
        if "med_sub" in kwargs: 
            self.med_sub = kwargs.get("med_sub")
        if "move_sub" in kwargs: 
            self.move_sub = kwargs.get("move_sub")
        if "sig_clip" in kwargs: 
            self.sig_clip = kwargs.get("sig_clip")
            
    def update_wfs(self, a_wfs):
        self.data.active_wfs = a_wfs
        return True
    
    def acor_map(self):
        data = self.data
        x_acor, y_acor = data.data_get_ac(sub_len=self.sub_len, bg_sub=self.bg_sub, med_sub=self.med_sub, mov_sub=self.move_sub, sig_clip=self.sig_clip)
        avg_wfs = np.average((x_acor + y_acor)/2, axis=0)
        ## add in 
        return avg_wfs[:self.est_len]
        
    def xcor_map(self):
        data = self.data
        x_cor, y_cor = data.data_get_cc_f_all(sub_len=self.sub_len, bg_sub=self.bg_sub, med_sub=self.med_sub, mov_sub=self.move_sub, sig_clip=self.sig_clip)
        avg_cor = (x_cor + y_cor)/2
        wfs_use = data.active_wfs
        avg_xcor = np.zeros_like(avg_cor[0][0])
        count = 0
        for i, row in enumerate(avg_cor):
            for j in range(len(row)):
                if i < j and wfs_use[i] and wfs_use[j]:
                    avg_xcor = avg_xcor + row[j]
                    count = count + 1
        avg_xcor = np.divide(avg_xcor, count)
        return avg_xcor

###########################
###### MAIN FUNCTION ######
###########################

    def run(self, xcor=False, **kwargs):
        """
        Does the 3 part run for radial estimation
        
        Optional Inputs
        --------
        xcor: Boolean
        if true runs xcor, if false runs acor (step1)
        
        detect_clp: sdv
        cut on detectins to determine significance (step2)
        
        clstr_method: string
        type of clustering method, options: meanshift, dbscan, affprop (step3)
        
        Outputs
        ----------
        spds: list
        list of detected speeds in m/s
        
        dirs: list
        list of detected directions in radians
        """
        self.update_params(**kwargs)  ## update_params
        #  saving parameters
        detect_clp = self.detect_clp # default 4
        clstr_method = self.clstr_method # default meanshift
        clstr_fn = dispatcher[clstr_method] # translating string to function
        #  choosing between auto or crodd correlations 
        if not xcor:
            avg_awfs = self.acor_map()
        if xcor:    
            avg_awfs = self.xcor_map()
        avg_awfs_sub = np.subtract(avg_awfs, np.average(avg_awfs, axis = 0)) #  subtracting averages pixel value from frame
        
        ### PART 1 : determining detection levels
        #detect_lvl = detect_map(avg_awfs_sub)
        detect_lvl = detect_map(avg_awfs)
        ### PART 2 : Speed Map
        dtct, spds, dirs, xs, ys = speed_map(detect_lvl, detect_clp, self.hz)
        ### PART 3: Clustering
        clusters, yhat = clstr_fn(self.data, spds, xs, ys)
        # saving these parts based on a or x cor     
        if not xcor:
            self.a_dtct = dtct
            self.a_spds, self.a_dirs = spds, dirs
            self.a_xs, self.a_ys= xs, ys
            self.a_cstr = clusters
            self.a_yhat = yhat
            self.a_d_lvl = detect_lvl
        if xcor: 
            self.x_dtct = dtct
            self.x_spds, self.x_dirs = spds, dirs
            self.x_xs, self.x_ys = xs, ys
            self.x_cstr = clusters
            self.x_yhat = yhat
            self.x_d_lvl = detect_lvl
        return spds, dirs
    
    def return_table(self, **kwargs):
        """
        Returns a df summary of clusters detected with parameters. 
        If no outputs, will call run function on acor and xcor
        
        Optional Inputs
        --------
        sdv_cnd: int
        if 0, no condensing. Otherwise, stdv for condensing peaks
        
        n_filter: int
        min number of counts for a cluster to be considered
        
        Outputs
        ----------
        df: DataFrame
        a pandas table summarizing clusters
        """
        self.update_params(**kwargs)  ## update_params
        # Checking to see if system has been run before
        if self.a_dtct is None:
            self.run()
        if self.x_dtct is None:
            self.run(xcor = True)
        # Returning summary of cluster
        asum = self.summary_cluster()
        xsum = self.summary_cluster(xcor = True)
        # Classifying peaks
        aclass, xclass = self.classify_peaks(asum, xsum)
        # df for acor
        dfa = pd.DataFrame(asum, columns =['dir', 'dir_std', 'spd', 'spd_std', 'count'], dtype = float)
        dfa['class'] = aclass
        dfa['xcor'] = 0 
        # df for xcor
        dfx = pd.DataFrame(xsum, columns =['dir', 'dir_std', 'spd', 'spd_std', 'count'], dtype = float)
        dfx['class'] = xclass
        dfx['xcor'] = 1 
        df = pd.concat([dfa, dfx], ignore_index=True)
        df['name'] = self.name
        # saving parameters:
        df['clstr_method'] = self.clstr_method
        df['detect_clp'] = self.detect_clp
        df['sdv_comp'] = self.sdv_comp
        df['sdv_cnd'] = self.sdv_cnd
        df['n_filter'] = self.n_filter
        # save table
        self.table = df
        return df
    
    def summary_cluster(self, xcor=False, **kwargs):
        """
        Returns a list summary of clusters detected. 
        If no outputs, will call run function on acor and xcor
        
        Optional Inputs
        --------
        sdv_cnd: int
        if 0, no condensing. Otherwise, stdv for condensing peaks
        
        n_filter: int
        min number of counts for a cluster to be considered
        
        Outputs
        ----------
        summary: list
        a lists of list containing values of ['dir', 'dir_std', 'spd', 'spd_std', 'count']
        """
        self.update_params(**kwargs)  ## update_params
        sdv_cnd = self.sdv_cnd
        n_filter = self.n_filter
        # based off of cluster output
        # include counts per cluster     
        # choosing saved variables base on xcor/acor
        if not xcor and self.a_dtct is None:
            self.run()
        if xcor and self.x_dtct is None:
            self.run(xcor = True)
        # pulling method vars
        dtct = self.x_dtct if xcor else self.a_dtct
        spds = self.x_spds if xcor else self.a_spds
        dirs = self.x_dirs if xcor else self.a_dirs
        xs = self.x_xs if xcor else self.a_xs
        ys = self.x_ys if xcor else self.a_ys
        clusters = self.x_cstr if xcor else self.a_cstr
        yhat = self.x_yhat if xcor else self.a_yhat 
        summary = [] # for each cluster return summary of cluster
        for cluster in clusters:
            row_ix = where(yhat == cluster) # get row indexes for cluster samples
            ## Standardized return 5 element structure
            mean_x = np.mean(xs[row_ix])
            mean_y = np.mean(ys[row_ix])
            mean_dir = to_degrees(dir_c_proj(mean_x, mean_y, 0))
            std_dir = to_degrees(dir_std(dirs[row_ix]))
            summary.append([mean_dir,
                            std_dir,
                            np.mean(spds[row_ix]),
                            np.std(spds[row_ix]),
                            row_ix[0].shape[0]])
        summary = self.condense_peaks(summary, sdv_cnd) # condensing peaks
        summary = self.filter_peaks(summary, n_filter) # filtering peaks
        return summary
    
    def classify_peaks(self, acor_sum, xcor_sum):
        # standard deviation and means
        # mean within a standard deviation
        # setting up flag arrays
        acor_clas = np.array(["FL" for peak in acor_sum])
        xcor_clas = np.array(["NA" for peak in xcor_sum])
        #iteratively comparing between acor and xcor
        for a in range(len(acor_sum)):
            for x in range(len(xcor_sum)):
                if self.compare_peaks(acor_sum[a], xcor_sum[x]):
                    xcor_clas[x] = "GL"
                    acor_clas[a] = "GL"
        return acor_clas, xcor_clas
    
    def condense_peaks(self, summary, sdv_cnd):
        # takes in a summary based on func summary_cluster
        # looks at peaks within n standard deviations (sdv_cnd)
        if sdv_cnd == 0:
            return summary
        # compare all the peaks, and condense those that are closer
        # sort by count
        summary = sorted(summary, key=lambda x: x[4], reverse=True)
        sum_new = []
        for a1 in range(len(summary)):
            # seed with first array element
            if a1 == 0:
                sum_new = [summary[a1]]
            else:
                avged = 0
                for a2 in range(len(sum_new)):
                    if avged == 0 and self.compare_peaks(summary[a1], sum_new[a2], sd = sdv_cnd):
                        # average peaks
                        pk = self.avg_peaks(summary[a1], sum_new[a2])
                        # change sum_new
                        sum_new[a2] = pk
                        avged = 1
                        break
                if avged == 0 :
                    sum_new.append(summary[a1])
        return sum_new

    def avg_peaks(self, pk1, pk2):
        # combines two peaks
        dir_avg = np.average([pk1[0], pk2[0]])
        spds_avg = np.average([pk1[2], pk2[2]])
        dir_sdv = np.sqrt(pk1[1]**2 + pk2[1]**2)
        spds_sdv = np.sqrt(pk1[3]**2 + pk2[3]**2)
        counts = pk1[4] + pk2[4]
        return [dir_avg, dir_sdv, spds_avg, spds_sdv, counts]
    
    def filter_peaks(self, summary, n_filter):
        # only returns peaks that are greater or eaqual to filter number
        return [pk for pk in summary if pk[4] >= n_filter]
    
    def compare_peaks(self, peak1, peak2, sd=3):
        #returns if two peaks are similar
        spd_diff = np.abs(peak1[0] - peak2[0])
        dir_diff = np.abs(peak1[2] - peak2[2])
        #HARDCODED looking in range of 3 times stdev
        spd_stdev = sd*np.max([peak1[1], peak2[1]])
        dir_stdev = sd*np.max([peak1[3], peak2[3]])
        #boolean values if these match
        spd_q = spd_diff <= spd_stdev
        dir_q = dir_diff <= dir_stdev
        return spd_q & dir_q
    
    ### Saving
    
    def est_file_gen(self, plot_dir, plot_type, file_type):
        """
        Generates a file name for a plot, renumbers if name already taken
            args: plot_type (varies by plot func), file_type ("png" or "gif"), Graphing specifications
            retuns: out_file (descriptive string)
        """
        out_dir = self.data.out_dir + "est/" + plot_dir + "/"
        out_file = self.data.name 
        # generating file structure more finely
        if self.data.s_sub and self.data.tt_sub:
            out_file = out_file + "_stt"
        elif self.data.tt_sub:
            out_file = out_file + "_tt"
        elif self.data.s_sub:
            out_file = out_file + "_s"
        else:
            out_file = out_file + "_raw"
        out_file = out_file + plot_type
        #check if path exists
        path = self.data.path_check(out_dir, loud=False)
        # Renumber if repeats 
        out_base = out_dir + out_file
        i = 0
        while os.path.exists(out_base + "_%s.%s" % (i, file_type)):
            i += 1
        out = out_base + "_%s.%s" % (i, file_type)
        return out
    
    def save_text(self, df, txt_type):
        """Takes in the type of txt and sends to the proper save location

        Inputs
        --------
        df: DataFrame
        the pandas dataframe to save
        
        txt_type: string
        the type of data being saved, put into file name
        
        Outputs
        ----------
        plot dir: string
        location plot saved to
        """
        # saving a text file
        plot_dir = "txt"
        file_type = "txt"
        plot_dir =  self.est_file_gen(plot_dir, txt_type, file_type)
        
        # saving a Pandas df to txt
        df.to_csv(plot_dir, header=True, index=False, sep='\t', mode='a')
        return plot_dir
    
    def save_plot(self, plot, **kwargs):
        """Takes in the type of plot and sends to the proper save location

        Inputs
        --------
        plot: string
        the type of plot desired
        
        Optional Inputs
        -------------------
        kwargs: dictionary
        arguments for plotting function

        Outputs
        ----------
        plot dir: string
        location plot saved to
        """
        #determining placement from 
        plot_dir =  self.est_file_gen(plot, "_"+plot, "png")
        if plot == "spds":
            ## speed plot
            self.plot_spds_detect(plot_dir=plot_dir, **kwargs)
        elif plot == "clstr":
            ## cluster plot
            self.plot_clstr(plot_dir=plot_dir, **kwargs)
        elif plot == "cor_prob":
            ## prob cor plot
            self.plot_prob_hist(plot_dir=plot_dir, **kwargs)
        elif plot == "lyr_prob":
            pass
        return plot_dir
    
    ### Plotting
        
    def plot_spds_detect(self, **kwargs):
        """Plots the speeds and directions above the detection limit

        Optional Inputs
        -------------------
        xcor : boolean
        plotting xcor values
       
        acor : boolean
        plotting acor values
        """
        # pipeline iterating          
        acor = kwargs.get("acor") if 'acor' in kwargs else True
        xcor = kwargs.get("xcor") if 'xcor' in kwargs else True
        # starting plot
        fig = plt.figure(figsize=(10,5))
        plt.xlim([0, 360])
        plt.ylim([0, 50])       
        if acor and self.a_dirs is not None:
            plt.scatter(np.degrees(self.a_dirs), self.a_spds, c="green",  alpha = .3, label="a")
        if xcor and self.x_dirs is not None:
            plt.scatter(np.degrees(self.x_dirs), self.x_spds, c="red",  alpha = .3, label="x")            
        param_title = "Detection Clip: "+ str(self.detect_clp)+ ", \n"
        if self.bg_sub:
            param_title = param_title + "avg sub len "+ str(self.sub_len)
            if self.move_sub:
                param_title = param_title + ", window"
            if self.sig_clip:
                param_title = param_title+ ", clipped"
        plt.title(self.name + ' Detections with ' + param_title)
        plt.ylabel('wind speed (m/s)')
        plt.xlabel('wind dir (degrees)')
        plt.legend()
        # starting saves / returns
        if 'plot_dir' in kwargs:
            fig.savefig(kwargs.get('plot_dir'))
            plt.close()
            return True
        else:    
            return fig
        
    def plot_clstr(self, **kwargs):
        """Plots the cluster centers and 3*std dev on error

        Optional Inputs
        -------------------
        table : df
        plot your own cluster table
        """
        table = kwargs.get("table") if 'table' in kwargs else self.table
        if table.empty:
            return False
        
        fig = plt.figure(figsize=(10,5))
        plt.xlim([0, 360])
        plt.ylim([0, 50])
        ax = plt.gca()

        for index, row in table.iterrows():
            clr  = det_color(row["class"])
            #alph = det_alph(row["count"])
            alph= 0.5
            plt.scatter(row['dir'], row['spd'], c=clr,  alpha = 1)
            spread = matplotlib.patches.Ellipse((row['dir'], row['spd']), 3*row['dir_std'], 3*row['spd_std'], angle=0, color=clr, alpha=alph)
            ax.add_patch(spread)

        param_title = "Method: " + self.clstr_method + ", N filter: " + str(self.n_filter) + ", sdv comp: " + str(self.sdv_cnd) 
        plt.title(self.name + ' Estimation Clusters \n ' + param_title)
        plt.ylabel('wind speed')
        plt.xlabel('wind dir')   
        
        if 'plot_dir' in kwargs:
            fig.savefig(kwargs.get('plot_dir'))
            plt.close()
            return True  
        return fig
    
    # Plotting histogram:
    def plot_prob_hist(self, **kwargs):
        # inputs: minimum velocity, maximum velocity, bin size, label increment
        
        #TODO: assumes you've run both acor and xcor already
        #settig args
        rmax = kwargs.get("rmax") if 'rmax' in kwargs else 25
        rmin = kwargs.get("rmin") if 'rmin' in kwargs else -rmax
        rinc = kwargs.get("rinc") if 'rinc' in kwargs else 2
        rlab = kwargs.get("rlab") if 'rlab' in kwargs else 5
        
        # generating normalization
        #TODO: assumes size of cor map
        full = speed_map(np.ones([201, 15, 15]), 0.5, self.hz)
        v_x, v_y = proj_xy(full[2], full[1])

        # we get speeds and directions out of  run
        av_x, av_y = proj_xy(self.a_dirs, self.a_spds)
        xv_x, xv_y = proj_xy(self.x_dirs, self.x_spds)

        param_title = "Detection Clip: "+ str(self.detect_clp) 
        title = 'UNNormed wind map, %  \n' +self.name +', ' + param_title
        #plt.xlabel('$v_x$ (m/s)')
        #plt.ylabel('$v_y$ (m/s)') 

        fig, (ax, ax2, ax3, cax) = plt.subplots(ncols=4,figsize=(9,3), 
                      gridspec_kw={"width_ratios":[1,1,1, 0.05]})
        fig.subplots_adjust(wspace=0.3)
        fig.suptitle(str(title))
        
        bin_range = np.arange(rmin, rmax, rinc)

        nH, xedges, yedges = np.histogram2d(v_x, v_y, bins=[bin_range, bin_range])
        aH, xedges, yedges = np.histogram2d(av_x, av_y, bins=[bin_range, bin_range])
        xH, xedges, yedges = np.histogram2d(xv_x, xv_y, bins=[bin_range, bin_range])

        am = ax.imshow(aH, interpolation='nearest', origin='low',
                    extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]], cmap=plt.cm.Reds, vmin=0, vmax=50)
        xm = ax2.imshow(xH, interpolation='nearest', origin='low',
                    extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]], cmap=plt.cm.Reds, vmin=0, vmax=50)
        m = ax3.imshow(aH- xH, interpolation='nearest', origin='low',
                    extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]], cmap=plt.cm.Reds, vmin=0, vmax=50)

        ax.set_xlabel("acor")
        ax2.set_xlabel("xcor")
        ax3.set_xlabel("sub")

        fig.colorbar(m, cax=cax)
         # starting saves / returns
        if 'plot_dir' in kwargs:
            fig.savefig(kwargs.get('plot_dir'))
            plt.close()
            return True
        else:    
            return fig
        
############################    
### Generic Plotting
############################

def plot_clstr_table(table, name="Data", min_count=0):
    fig = plt.figure(figsize=(10,5))
    plt.xlim([0, 360])
    plt.ylim([0, 50])
    ax = plt.gca()

    for index, row in table[table['count']>min_count].iterrows():
        clr  = det_color(row["class"])
        #alph = det_alph(row["count"])
        alph= 0.3
        plt.scatter(row['dir'], row['spd'], c=clr,  alpha = 1)
        spread = matplotlib.patches.Ellipse((row['dir'], row['spd']), 3*row['dir_std'], 3*row['spd_std'], angle=0, color=clr, alpha=alph)
        ax.add_patch(spread)

    param_title = "Method: " + row['clstr_method'] + ", N filter: " + str(row['n_filter']) + ", sdv cnd: " + str(row['sdv_cnd'])
    plt.title(name +' Estimation Clusters, all \n ' + param_title)
    plt.ylabel('wind speed (m/s)')
    plt.xlabel('wind dir (degrees)')   
    return fig

def det_color(c):
    if c == "NA":
        return 'black'
    elif c == "GL":
        return 'blue'
    elif c == "FL":
        return 'orange'
        
        
# static functions
 
##############
### Part 1 ###
##############
def detect_map(avg_wf):
    #detect_lvl = np.zeros_like(avg_wfs)
    t_mean = np.zeros_like(avg_wf)
    t_stdev = np.zeros_like(avg_wf)

    for t in range(avg_wf.shape[0]):
        # 1.1.  For each pixel in the cov map assign a radius
        t_slice = avg_wf[t]
        sigma = 3
        # 1.2. find the mean (m) and stddev (sd) of the pixels 
        # in a radial annuli one pixel wide (r from 1 to 7 +/- 0.5.  
        # I do a sigma clipping on this with a clip=3 si
        for r in np.arange(8):
            rad_vals = np.multiply(t_slice, rad_map_boolean[r])
            mean, median, std = sigma_clipped_stats(rad_vals, sigma=sigma, mask=rad_map_inv[r])
            # TRY REDOING MASKING => try to flatten?
            #print(t, r, mean, std)
            t_mean[t] = t_mean[t] + mean*rad_map_boolean[r]
            t_stdev[t] = t_stdev[t] + std*rad_map_boolean[r]
    # 1.3. For each pixel I assign a "detection level in sigma over the mean"
    # = (im[x,y] - m)/sd.   
    # This results in a detection map for each cov map time slice. 
    detect_map = np.divide(np.subtract(avg_wf, t_mean), t_stdev)
    return detect_map

def stdv_map(avg_wf):
    #detect_lvl = np.zeros_like(avg_wfs)
    t_mean = np.zeros_like(avg_wf)
    t_stdev = np.zeros_like(avg_wf)

    for t in range(avg_wf.shape[0]):
        # 1.1.  For each pixel in the cov map assign a radius
        t_slice = avg_wf[t]
        sigma = 3
        # 1.2. find the mean (m) and stddev (sd) of the pixels 
        # in a radial annuli one pixel wide (r from 1 to 7 +/- 0.5.  
        # I do a sigma clipping on this with a clip=3 si
        for r in np.arange(8):
            rad_vals = np.multiply(t_slice, rad_map_boolean[r])
            mean, median, std = sigma_clipped_stats(rad_vals, sigma=sigma, mask=rad_map_inv[r])
            # TRY REDOING MASKING => try to flatten?
            #print(t, r, mean, std)
            t_mean[t] = t_mean[t] + mean*rad_map_boolean[r]
            t_stdev[t] = t_stdev[t] + std*rad_map_boolean[r]
    # 1.3. For each pixel I assign a "detection level in sigma over the mean"
    # = (im[x,y] - m)/sd.   
    # This results in a detection map for each cov map time slice. 
    #detect_map = np.divide(np.subtract(avg_wf, t_mean), t_stdev)
    return t_stdev
    
    
##############
### Part 2 ###
##############

# For each detection greater than some value I calculate a speed (vspd), direction (vdir). 
# This results in a speed and a direction "detection" map for each time slice.

def dist_c(x, y, c):
    #BUG: hard-coded 2.2 meter telescope
    pix_to_m = 2.2/15.0
    return np.sqrt(np.square(x-c) + np.square(y-c))

# issues with wind projections

def dir_wind(x,y,c):
    # this projects from (0,7) as north, clockwise angle
    return 0

def dir_c(x,y,c):
    return np.arctan2(x-c, c-y) + np.pi

def dir_c_proj(xp,yp,c):
    return np.arctan2(xp, yp)

def proj_x(dir_rad, c):
    return c*np.sin(dir_rad)

def proj_y(dir_rad, c):
    return c*np.cos(dir_rad)

def proj_xy(dir_rad, c):
    return c*np.sin(dir_rad), c*np.cos(dir_rad)

def dir_std(dirs):
    dirs_shift = (dirs + np.pi) % (2*np.pi)
    return np.min([np.std(dirs), np.std(dirs_shift)]) 

def to_degrees(rads):
    return (np.degrees(rads) + 360)%360

    
##############
### Part 3 ###
##############

# I also calculate the x,y location at the edge of the "circle" where the vdir projects to.
# Namely, xedgei = r0 * cos(thetai) 
# This makes the clustering analysis easier as it does away with the 2pi wrapping
    
def speed_map(detection_map, detect_val, hz):
    radius = 7
    
    #maps of the distance and direction for each pixel
    dist_map = np.array([[dist_c(x,y, radius) for y in range(15)] for x in range(15)])
    dir_map = np.array([[dir_c(x,y, radius) for y in range(15)] for x in range(15)])
    
    #storing arrays
    detect_list = np.array([])
    speed_list = np.array([])
    dir_list = np.array([])
    x_proj = np.array([])
    y_proj = np.array([])
    
    #values above detection value => make a map for detections
    i = 1
    for t_slice in detection_map:
        t_slice = ma.masked_invalid(t_slice)
        detect_t = ma.masked_where(t_slice <= detect_val, t_slice)
        #BUG: positive detections only
        
        t = i/hz
        speed_t = np.divide(dist_map, t*(8/2.2))
        
        #detections
        detect_det = t_slice[detect_t.mask == False]
        detect_list = np.append(detect_list, detect_det)
        # speeds 
        speed_det = speed_t[detect_t.mask == False]
        speed_list = np.append(speed_list, speed_det)
        # directions 
        dir_det = dir_map[detect_t.mask == False]
        dir_list = np.append(dir_list, dir_det)
        # x-y projections
        x_t = proj_x(dir_det, radius)
        y_t = proj_y(dir_det, radius)
        x_proj = np.append(x_proj, x_t)
        y_proj = np.append(y_proj, y_t)
        
        i+=1
        #if i>100: break
        
    return [detect_list, speed_list, dir_list, x_proj, y_proj]

##############
### Part 4 ###
##############
# various clustering algorithms
# see cluster.py
    

################# Running estimator on a CSV
def est_csv(file, ttsub = True, **kwargs):
    tic = time.perf_counter()
    if not os.path.isfile(file):
        print("could not find file")
    df_main = pd.read_csv(file) # read in csv
    rslt_df = df_main[df_main['ttsub'] == ttsub]
    rslt_df = rslt_df.drop_duplicates(subset=['dataname']).sort_values(by=['dataname'])
    fits = rslt_df["outfits"]
    fits_in = fits.values
    for i, ft in enumerate(fits_in):
        est_file(ft, text = False, plot = False, prob_plot = True, **kwargs)
        toc = time.perf_counter()
        print(f"File {i} in {toc - tic:0.4f} seconds")
    return

################# running estimator for one file
def est_file(cor_f, method="meanshift",  text = True, plot = False, prob_plot = True, **kwargs):
    """
    Takes in arguments and access the Estimator. Passes kwargs to the function.  
    
    """
    ## TODO: set up kwargs
    ## TODO: set up multiplot access
    detect_cpl = kwargs.get("rlab") if 'rlab' in kwargs else 5
    
    # set up corr object
    er_pipe = Estimate_simple(cor_f)
    if 'detect_clp' in kwargs:
        er_pipe.detect_clp = kwargs.get("detect_clp")
    # run estimation and return table
    table = er_pipe.return_table(clstr_method = method, sdv_cnd = 2, n_filter = 20)
    
    # deciciding where to store locations
    if text:
        # saving a text file
        plot_type = "_est_simple"
        # for a new directory with estimations
        er_pipe.save_text(table, plot_type, **kwargs)
        #saving table to this directory:   
    if plot:
        # Typical plots
        # Detect plot
        #plot_dir = "spds"
        #plot_type = "_spds"
        #fig = er_pipe.plot_spds_detect() # generate figure
        #er_pipe.save_plot(fig) # saving a fig
        f = er_pipe.save_plot("spds", **kwargs)
        # Cluster plot
        #plot_dir = "clstr"
        #plot_type = "_clstr"
        #fig = plot_clstr_table(table, name=er_pipe.name, min_count=er_pipe.n_filter) # generate figure
        f = er_pipe.save_plot("clstr", **kwargs)
        plt.close() 
    if prob_plot:
        # generate figure
        #plot_dir = "prob"
        #plot_type = "prob"
        #fig = er_pipe.plot_prob_hist()
        #er_pipe.save_plot(fig, plot_dir, plot_type)
        f = er_pipe.save_plot("cor_prob", **kwargs)
        print(f)
        plt.close() 
        
    return table


#### 