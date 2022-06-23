# t_const_code.py
# Eden McEwen 
# August 11, 2021
# Code used to generate fits for correlation time constants

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
#import pipeline.code.data_table as d_t
import pipeline.code.graph_code as gc

import pipeline.code.Correlator as Cor

# old fitting code:
import scipy.optimize as opt

# Linear Fits
def Linear_tot(x, y):
    a, b = Linear_Fit(x, y)
    a_err, b_err = Linear_err(x, y, a, b)
    return np.array([a, a_err, b, b_err])

def Linear(x, a, b):
    return a * x + b

def Linear_Fit(x, y):
    # returns [slope, intercept]
    # the optimized curve returns the optimized abc values for Gaussian
    # [height, position of peak center (mean), standard deviation]
    popt, pcov = opt.curve_fit(Linear, x, y)
    return popt[0], popt[1]
    
def Linear_fit(x, y):
    # returns fit
    return opt.curve_fit(Linear, x, y)
    
def Linear_est(x, a, b):
    #returns y estimates
    y_est = [Linear(xi, a, b) for xi in x]
    return y_est

#### Ply fit:

def Poly(x, a, b):
    return b*x**a

def Poly_est(x, a, b):
    #returns y estimates
    y_est = [Poly(xi, a, b) for xi in x]
    return y_est

def Poly_fit(x, y):
    # returns [slope, intercept]
    # the optimized curve returns the optimized abc values for Gaussian
    # [height, position of peak center (mean), standard deviation]
    return opt.curve_fit(Poly, x, y)

##################### Applying Fits #####################
#########################################################

# Using per WFS fits
def WFS_noise_linear(p_file):
    tmax = 15
    mWFS = 3
    
    # pulling and averaging WFS data
    pull_data = Cor.Correlator("", "", "", f_file = p_file)
    t = np.arange(tmax )
    avg_acor = (pull_data.acor_x + pull_data.acor_y)/2
    center_peak_avg = avg_acor[:mWFS, 0:tmax, 7, 7]
    
    # averaging 
    center_avg = np.average(center_peak_avg, axis=0)
    
    #linear fit
    c_min = 1
    c_max = 4
    dt0_diff = []
    
    # fit for each wfs individually 
    # find difference between fit and actual
    for wfs in range(mWFS):
        w_center_peak_avg = center_peak_avg[wfs]
        f1_w, f2_w = Linear_fit(t[c_min:c_max], w_center_peak_avg[c_min:c_max])
        lin_fit_w = Linear_est(t[0:c_max], f1_w[0], f1_w[1])
        dt0_diff.append(w_center_peak_avg[0]-lin_fit_w[0])
        
    return dt0_diff

# Distinct for each WFS
def center_peak_pull(p_file, wfs=3):
    pull_data = Cor.Correlator("", "", "", f_file = p_file)
    
    avg_acor = (pull_data.acor_x + pull_data.acor_y)/2
    center_peak_avg = avg_acor[:3, 2:, 7, 7]
    
    return center_peak_avg

# Average over each WFS
def center_avg_pull(p_file):
    pull_data = Cor.Correlator("", "", "", f_file = p_file)
    
    avg_acor = (pull_data.acor_x + pull_data.acor_y)/2
    center_peak_avg = avg_acor[:3, :, 7, 7]
    
    center_avg = np.average(center_peak_avg[:2], axis=0)
    center_avg[center_avg < 0] = 0.000000000000000000000000001
    
    return center_avg

####################### WFS avg decays ##################

# re-writing
def decay_fn(p_file, tmax_p, c_min = 100, c_max = 500):
    # return the center peaks avg, and fit variables
    ## setting variables
    pull_data = Cor.Correlator("", "", "", f_file = p_file)
    t = np.arange(tmax_p)
    center_avg = center_avg_pull(p_file)
    
    # for fit
    f1, f2 = Linear_fit(t[c_min:c_max], np.log(center_avg)[c_min:c_max])
    
    return center_avg, t, f1


# re-writing to better account for outliers
def decay_fn_2(p_file, tmax_p, c_min = 100, c_max = 500):
    # return the center peaks avg, and fit variables
    ## setting variables
    pull_data = Cor.Correlator("", "", "", f_file = p_file)
    t = np.arange(tmax_p)
    # no longer "fixes" values less than 0
    center_avg = np.average(center_peak_pull(p_file)[:2], axis=0)
    
    #now masks those from t as well, just for fitting
    t_e = t[c_min:c_max]
    m_lavg = ma.log(center_avg[c_min:c_max])
    
    # for fit
    f1, f2 = Linear_fit(t_e[m_lavg.mask == False], m_lavg[m_lavg.mask == False])
    
    return center_avg, t, f1


def decay_fn_var(p_file, tmax_p, c_min = 100, c_max = 500):
    # return the center peaks avg, and fit variables
    ## setting variables
    center_avg, t, f1 = decay_fn_2(p_file, tmax_p, c_min = c_min, c_max = c_max)
    pull_data = Cor.Correlator("", "", "", f_file = p_file)
    
    # rounded fit variables
    a = np.around(1/f1[0] / pull_data.hz_pull(), 4) 
    b = np.around(f1[1] / pull_data.hz_pull(), 4) 
    
    return -a, b

#################### Fits per WFS #####################


def decay_fn_wfs(p_file, tmax_p, wfs=0, c_min = 100, c_max = 500):
    # return the center peaks avg, and fit variables
    ## setting variables
    pull_data = Cor.Correlator("", "", "", f_file = p_file)
    t = np.arange(tmax_p)
    center_avg = center_peak_pull(p_file)[wfs]
    
    # for fit
    f1, f2 = Linear_fit(t[c_min:c_max], np.log(center_avg)[c_min:c_max])
    
    return center_avg, t, f1


### Better handling of invalid values
def decay_fn_wfs_2(p_file, tmax_p, wfs=0, c_min = 100, c_max = 500):
    # return the center peaks avg, and fit variables
    ## setting variables
    pull_data = Cor.Correlator("", "", "", f_file = p_file)
    t = np.arange(tmax_p)
    center_avg = center_peak_pull(p_file)[wfs]
    
    #now masks those from t as well, just for fitting
    t_e = t[c_min:c_max]
    m_lavg = ma.log(center_avg[c_min:c_max])
    
    # for fit
    f1, f2 = Linear_fit(t_e[m_lavg.mask == False], m_lavg[m_lavg.mask == False])
    
    return center_avg, t, f1, f2


def decay_fn_wfs_var(p_file, tmax_p, wfs=0, c_min = 100, c_max = 500):
    # return the center peaks avg, and fit variables
    ## setting variables
    center_avg, t, f1, __ = decay_fn_wfs_2(p_file, tmax_p, wfs=wfs, c_min = c_min, c_max = c_max)
    pull_data = Cor.Correlator("", "", "", f_file = p_file)
    
    # rounded fit variables
    a = np.around(1/f1[0] / pull_data.hz_pull(), 4) 
    b = np.around(f1[1] / pull_data.hz_pull(), 4) 
    
    return -a, b


####################### DECAY PLOTTING #########################

# re-writing
def decay_plot(p_file, tmax_p, c_min = 100, c_max = 500):
    
    ## Center Peak Plots
    pull_data = Cor.Correlator("", "", "", f_file = p_file)
    center_avg, t, f1 = decay_fn(p_file, tmax_p, c_min, c_max)
    center_peak_avg = center_peak_pull(p_file)
    dt = np.arange(tmax_p) #integers
    
    # rounded fit variables
    a = np.around(1/f1[0] / pull_data.hz_pull(), 4) 
    b = np.around(f1[1] / pull_data.hz_pull(), 4) 

    # Plotting center peaks
    fig, ax =  plt.subplots(2,1,figsize=(10,8))
    g_cmap = cm.get_cmap('Greens', 8)
    for i in range(3): 
        c_float = (i+2)/6
        ax[0].plot(t, center_peak_avg[i][dt], label = "WFS" + str(i), color = g_cmap(c_float))
        ax[1].plot(t, np.log(center_peak_avg[i][dt]),  color = g_cmap(c_float))

    # Log fit
    log_fit = Linear_est(t, f1[0], f1[1])
    ax[1].plot(t, log_fit, label = f"fit = {b} - t / {-a}")
    ax[0].plot(t, np.exp(log_fit), label = f"fit = exp({b} - t / {-a})")

    ax[0].axhline(0, color="black", alpha=0.3)
    ax[1].axvline(c_min, color="black", alpha=0.3)
    ax[1].axvline(c_max, color="black", alpha=0.3)

    ax[0].legend(loc = 'upper right')
    ax[0].set_ylabel('Intensity')
    ax[1].legend(loc = 'upper right')
    ax[1].set_ylabel('Intensity, log')
    ax[1].set_xlabel('frames')

    ax[0].set_title(f'WFS average Autocorrelation peak \n {pull_data.name},  length = {pull_data.x_slopes.shape[1]}, tmax = {pull_data.tmax}')
    fname = p_file.replace("fits/", "decay/").replace(".fits", ".png")
    plt.savefig(fname, dpi=300)

    
# EXPONENTIAL FUNCTION

# new exponential curve fit
from scipy.optimize import curve_fit

def exp_func(x, a, k, b):
    return a * np.exp(-k*x) + b

# re-writing
def decay_exp_fn(p_file, tmax_p, c_min = 100, c_max = 500):
    # return the center peaks avg, and fit variables
    pull_data = Cor.Correlator("", "", "", f_file = p_file)
    
    ## setting variables
    t = np.arange(tmax_p)
    avg_acor = (pull_data.acor_x + pull_data.acor_y)/2
    center_peak_avg = avg_acor[:3, 2:, 7, 7]
    
    # for fit
    center_avg = np.average(center_peak_avg[:2], axis=0)
    center_avg[center_avg < 0] = 0.000000000000000000000000001
    
    x = t[c_min:c_max]
    y = center_avg[c_min:c_max]
    
    # Initial guess:
    k0 = -np.log(y[0] - y[1])
    b0 = y[0] - 1
    p0 = (1, 0.001, b0) # starting search coefs
    
    # curve fit
    opt, pcov = curve_fit(exp_func, x, y, p0, maxfev=10000)
    print(opt)
    print(pcov)
    a, k, b = opt
    
    return center_peak_avg, a, k, b

# re-writing
def decay_exp_plot(p_file, tmax_p, c_min = 100, c_max = 500):
    center_peak_avg, a, k, b = decay_exp_fn(p_file, tmax_p, c_min, c_max)
    t = np.arange(1,tmax_p)
    exp_fit = exp_func(t, a, k, b)
    
    ## Center Peak Plots
    g_cmap = cm.get_cmap('Greens', 8)
    fig,ax =  plt.subplots(2, 1, figsize=(10,8))
    pull_data = Cor.Correlator("", "", "", f_file = p_file)

    # Plotting center peaks
    for i in range(3): 
        c_float = (i+2)/6
        ax[0].plot(t, center_peak_avg[i][t], label = "WFS" + str(i), color = g_cmap(c_float))
        ax[1].plot(t, np.log(center_peak_avg[i][t]),  color = g_cmap(c_float))

    # Log fit
    ax[0].plot(t, exp_fit, label = f"fit = {round(a, 3)} * exp(-{round(k, 3)} t) - {round(b, 3)}")
    ax[1].plot(t, np.log(exp_fit), label = f"fit = - t / {round(1/k, 3)}")

    ax[0].axhline(0, color="black", alpha=0.3)
    ax[0].axvline(c_min, color="black", alpha=0.3)
    ax[0].axvline(c_max, color="black", alpha=0.3)
    ax[1].axvline(c_min, color="black", alpha=0.3)
    ax[1].axvline(c_max, color="black", alpha=0.3)

    ax[0].legend(loc = 'upper right')
    ax[0].set_ylabel('Intensity')
    ax[1].legend(loc = 'upper right')
    ax[1].set_ylabel('Intensity, log')
    ax[1].set_xlabel('frames')
    #ax[1].set_xscale('log')
    #ax[0].set_xscale('log')

    ax[0].set_title(f'WFS average Autocorrelation peak \n {pull_data.name}, tmax = {pull_data.tmax}')
    fname = p_file.replace("fits/", "decay/").replace(".fits", ".png")
    plt.savefig(fname, dpi=300)
