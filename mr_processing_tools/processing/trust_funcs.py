#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Sky Jones
"""

import sys
sys.path.append('..')

import numpy as np
from scipy import stats
import scipy
import statsmodels.stats.api as sms
import pandas as pd
from scipy.optimize import curve_fit
from scipy.stats.distributions import t

import helpers.trust_helpers as trhe


def get_trust_diff(trust_subset, nte):
    '''
    Takes a single acquistion of TRUST and calculates the difference matrix

    Parameters
    ----------
    trust_subset : 3d np ndarray
    data for the acquisition. First two dimensions should spatial, while the third is the dynamic number.
    The dynamics should be set up such that each ctrl-label pair is after one another, and repeated measurements
    for each echo time are contiguous.
        For example, if the third dimension has 24 dynamics, and there are 4 echo times,
        then the dynamics should proceed as
            C1_alpha L1_alpha C1_beta L1_beta C1_gamma L1_gamma C2_alpha L2_alpha C2_beta L2_beta C2_gamma L2_gamma C3_alpha L3_alpha C3_beta L3_beta C3_gamma L3_gamma
        where C/L is control/label, greek subscripts are repeated measurements at a given echo time, and numerals indicate the echo time
        Note that the number of repeated measurements per echo time is inferred from the length of the 3rd dimension and nte
    nte : int
        the number of discrete echoe times

    Returns
    -------
    trust_diff : 4d np array
        The difference matrix for the acquistion. Its first two dimensions
        match the spatial dimensions of trust_subset, while the size of the 3rd dimension
        is nte and the 4th is this measurment index within each echo time. Thus, trust_diff[:,:,i,j] is the jth intensity difference at the ith
        echo value.

    '''
    
    trust_nrows, trust_ncolumns, trust_ndynamics = trust_subset.shape
    
    # Gather Parameters
    trust_measures = int(trust_ndynamics/(nte*2)) # Number of ctrl-label pairs per echo time

    # Order of Control & Label Images
    trust_nlabel = np.arange(1, trust_ndynamics+1, 2) # these are the indices of the labeled images
    trust_ncontrol = np.arange(0, trust_ndynamics+1, 2) # these are the indices of the control images

    # Average Across Measurements
    
    #trust_avglabel = np.zeros(shape=[trust_nrows, trust_ncolumns, trust_measures+1])
    #trust_avgcontrol = np.zeros_like(trust_avglabel)
    
    trust_labels = np.zeros(shape=[trust_nrows, trust_ncolumns, nte, trust_measures])
    trust_controls = np.zeros_like(trust_labels)
    counter = 0;
    for i in range(nte):
        interm_label = np.zeros(shape=[trust_nrows, trust_ncolumns, trust_measures+1])
        interm_control = np.zeros_like(interm_label)
        for j in range(trust_measures):
            interm_label[:,:,j] = trust_subset[:,:,trust_nlabel[counter]]
            interm_control[:,:,j] = trust_subset[:,:,trust_ncontrol[counter]]
            
            trust_labels[:,:,i,j] = interm_label[:,:,j]
            trust_controls[:,:,i,j] =  interm_control[:,:,j]
            
            counter += 1
        
        #trust_avglabel[:,:,i] = np.sum(interm_label[:,:,:], axis=2) / trust_measures
        #trust_avgcontrol[:,:,i] = np.sum(interm_control[:,:,:], axis=2) / trust_measures
        
    # subtract Label from Control
    #trust_diff = trust_avgcontrol - trust_avglabel
    trust_diff = trust_controls - trust_labels
    
    return trust_diff


def extract_sss_signal_difference(trust_subset, nte, auto=True, override_inversion=False):
    '''
    Takes a single TRUST acquisition and extracts the mean signal difference in the superior saggital sinus for
    each echo time

    Parameters
    ----------
    trust_subset : 3d np ndarray
        data for the acquisition. First two dimensions should spatial, while the third is the dynamic number.
        the dynamics should be set up such that each ctrl-label pair is after one another, and repeated measurements
        for each echo time are contiguous.
            For example, if the third dimension has 24 dynamics, and there are 4 echo times,
            then the dynamics should proceed as
                C1_alpha L1_alpha C1_beta L1_beta C1_gamma L1_gamma C2_alpha L2_alpha C2_beta L2_beta C2_gamma L2_gamma C3_alpha L3_alpha C3_beta L3_beta C3_gamma L3_gamma
            where C/L is control/label, greek subscripts are repeated measurements at a given echo time, and numerals indicate the echo time
            Note that the number of repeated measurements per echo time is inferred from the length of the 3rd dimension and the length of ete
            
    blood_t1 : float
        the blood T1 in seconds (?).
    nte : int
        The number of echo times in ms
    auto : bool
        If true, automatically picks out the SSS. Else uses manual ROI
    override_inversions : bool
        Specifies whether to override the automated data inversion of the difference map

    Returns
    -------
    trust_avg_sag
        the average signal difference bewtween label and control images in the SSS for each echo time
    trust_avg_sag_ci
        the 95% confidence interval of the averaged signal difference

    '''
    
    trust_diff_full = get_trust_diff(trust_subset, nte=nte)
    
    trust_diff = np.mean(trust_diff_full, axis=3) # average the n repeated measurments for each echo time (n is usually 3)
    
    
    
    #print('Evaluating inversion need')
    trust_slice = trust_diff[:,:,0] # we evaluate at the first echo because the signal difference is most extreme here
    trust_skewness = stats.skew(trust_slice, axis=None) # don't correct for bias, and get skewness of entire slice (not column-wise)
    #print(f'Imaging skewness: {trust_skewness} \n');
    
    if trust_skewness < 0: # detect automatically if we need to invert the TRUST diff
    # -skewness means the most of the data is at the + end of the range,
    # meaning the SSS (which constrasts with the background and is
    # always comprises the "light" tail of the distribution) is very negative and so
    # we need to invert
        trust_diff = -trust_diff
        trust_diff_full = -trust_diff_full
    
    
    if override_inversion: # just incase the automatic inversion doesn't work right, still allow user-specified inversion
        trust_diff = -trust_diff
        trust_diff_full = -trust_diff_full
        print('Additional manual inversion specified')
        
    # query for superior saggital sinus
    imslice = 0
    the_image = trust_diff[:,:,imslice]
    if not auto:
        trust_bw, px, py = trhe.manual_sss(the_image)
    else:
        trust_bw, px, py = trhe.auto_sss(the_image)
    
    # Isolate Sagittal Sinus

    trust_full_sag = np.zeros(shape=(trust_diff_full.shape[2], trust_diff_full.shape[3]))
    num_voxels = 4
    for n in range(trust_full_sag.shape[0]): # iterate through echo times
        # need to find the indices of the brightest averaged pixels for the echo time
        trust_sagsinus = trust_diff[:,:,n] * trust_bw
        ss_unraveled = np.ravel(trust_sagsinus)
        
        sort_indices = ss_unraveled.argsort() # indices that sort the array in ASCENDING order
        pick_vec_indices = sort_indices[-num_voxels:] # the indices of the $num_voxels biggest values
        
        # now that we know that indices of the brightest voxels for each averaged echo time, let's
        # extract the values used to make those means
        for m in range(trust_full_sag.shape[1]): # iterate through measurements within the echo time
        
            trust_sagsinus_rep = trust_diff_full[:,:,n,m] * trust_bw
            ss_unraveled_rep = np.ravel(trust_sagsinus_rep)
            
            voxels = np.array([ss_unraveled_rep[r] for r in pick_vec_indices])
            
            mean_sag = np.mean(voxels)
            trust_full_sag[n,m] = mean_sag
            # trust_full_sag[n,m] is the mean of the $num_voxels voxels
            # for the mth repetition of the nth echo time that at the spatial index
            # of the highest mean signal that that echo time within the SSS mask
            
            # note that this is different than the brightest voxels for the mth repetition =
            # of the nth echo time within the SSS mask
        
            
    trust_avg_sag = np.mean(trust_full_sag, axis=1)
    trust_avg_sag_ci = sms.DescrStatsW(trust_full_sag.T).tconfint_mean(alpha=0.05, alternative='two-sided')
    
    
    return trust_avg_sag, trust_avg_sag_ci
    

def calculate_t2_from_decay(trust_avg_sag, ete, blood_t1):
    """
    Calculates the blood T1 from a list of TRUST difference signals

    Parameters
    ----------
    trust_avg_sag : numpy ndarray
        1d array of the signal difference for each echo time in the SSS.
    ete : numpy ndarray
        1d array of echo times (in ms) corresponding to each index in trust_avg_sag
    blood_t1 : numpy ndarray
        the blood T1, in seconds (?).

    Returns
    -------
    blood_t2
        the estimated blood T2 in s
    fit
        tuple of the exponential decay fit as (exponent, coefficient)
    fit_ci
        tuple of the lower and upper 95% confidence intervals of each fitted parameter

    """

    # Fit to T2-Decay Curve
    X = np.array(ete) / 1000
    Y = np.array(trust_avg_sag)
    
    def fit_simp(ex, exp, coef):
        fit = coef*((np.e)**(ex*exp))
        return fit
    
    initial_guess = [-10, trust_avg_sag[0]]
    pars, pcov = curve_fit(fit_simp, X, Y, p0=initial_guess)
    
    alpha = 0.05
    en = len(Y)
    p = len(pars)
    
    dof = max(0, en - p) # number of degrees of freedom
    
    # student-t value for the dof and confidence level
    tval = t.ppf(1.0-alpha/2., dof) 
    
    fit_lower = []
    fit_upper = []
    for i, param,var in zip(range(en), pars, np.diag(pcov)):
        sigma = var**0.5
        param_lower = param - sigma*tval
        param_upper = param + sigma*tval
        
        fit_lower.append(param_lower)
        fit_upper.append(param_upper)

    trust_c = pars[0] # exp
    trust_s0 = pars[1] # coef
    
    fit = trust_c, trust_s0
    fit_ci = fit_lower, fit_upper
    
    # Compute Blood T2
    blood_t2 = 1 / ((1 / blood_t1) - trust_c)
    return blood_t2, fit, fit_ci
    


def t2_from_trust(trust_data, blood_t1, ete=[0, 40, 80, 160], auto=True, override_inversions=None):
    '''
    
    Takes the intensity readings from a TRUST scan and converts it to a T2 relaxation value
    

    Parameters
    ----------
    trust_data : np ndarray
        4d scaled array of trust data.
            the first two indices should be spatial indices (x,y)
            the third index should be dynamics (time)
            the fourth index should be the acquisition number
            
        Within each acquisition (i.e., trust_data[:,:,:,i]) the dynamics should
        progress as a control image followed by a label image for each echo time, and there can be
        an arbitrary number of repeated measures for each echo time as long as the repeats
        are contiguous. The number of repeated measures will be inferred.
        
        
    blood_t1 : float
        the blood t1 in seconds (?)
    ete : list of floats or ints
        The effective echo times in ms. The default is [0, 40, 80, 160].
    auto : bool
        If true, automatically picks out the SSS. Else uses manual ROI
    override_inversions : list of bools, len=2
        Specifies whether to override the automated data inversion for each repeat
        e.g., if override_inversion[0] is True, then the data for the first acquisition
            will be inverted if the automated decision was to not invert, but if
            the automated decision was to invert then the data will not be inverted.
            This is useful if the automated inversion choice is incorrect for
            an acquisition, though generally it is not needed.
        

    Returns
    -------
    trust_t2
        np array of the estimated T2 relaxation time in seconds as a float for each acquisition
    trust_signals
        list of pandas dataframes of the mean extracted signal at each echo time for each acquisition with 95% confidence intervals
        note that though ete is given in ms, the returned fit is calibrated using seconds
    trust_fits
        list of np arrays giving the parameters of the decay curve for each acquisition as with 95% confidence intervals
        y = coefficient * x**exponent
    
    '''
    
    trust_nrows, trust_ncolumns, trust_ndynamics, trust_nrepeats = trust_data.shape
    # note that in the original code, dynamics here were called echoes, and echoes here were called dynamics
    
    trust_scan_indices = np.arange(0, trust_nrepeats)
    #trust_t2 = np.zeros_like(trust_scan_indices)
    trust_t2 = []
    trust_signals = []
    trust_fits = []
    
    
    if override_inversions is None:
        override_inversions = [False for i in trust_scan_indices]
        
    for numscan, override_inversion in zip(trust_scan_indices, override_inversions):
        
        fit_cols = ['exponent', 'coefficient']
        fit_params = pd.DataFrame(columns=fit_cols)
        
        sig_cols = [str(i) for i,val in enumerate(ete)]
        signals = pd.DataFrame(columns=sig_cols)
        
        # Compute Blood T2
        trust_subset = trust_data[:,:,:,numscan]
        
        sss_signal_diffs, sss_signal_diffs_ci = \
            extract_sss_signal_difference(trust_subset, len(ete), auto=auto, override_inversion=override_inversion)
            
        signals.loc[0] = ete
        signals.loc[1] = sss_signal_diffs
        signals.loc[2] = sss_signal_diffs_ci[0]
        signals.loc[3] = sss_signal_diffs_ci[1]
        
        signals['signal'] = ['echo_time', 'mean', 'lower_95', 'upper_95']
        signals = signals.set_index('signal')
        
        blood_t2, fit, fit_ci = calculate_t2_from_decay(sss_signal_diffs, ete, blood_t1)
        
        fit_params.loc[0] = fit
        fit_params.loc[1] = fit_ci[0]
        fit_params.loc[2] = fit_ci[1]
        fit_params['fit'] = ['value', 'lower_95', 'upper_95']
        fit_params = fit_params.set_index('fit')
        
        trust_t2.append(blood_t2)
        trust_signals.append(signals)
        trust_fits.append(fit_params)
        
        
    return np.asarray(trust_t2), trust_signals, trust_fits


