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
from scipy import signal
import statsmodels.stats.api as sms
import pandas as pd
from scipy.optimize import curve_fit
from scipy.stats.distributions import  t
from skimage.filters import threshold_otsu

import helpers.general as ge

def quantify_cbf_from_asl(asl_data, m0_data, pld, ld):
    
    # pCASL ONLY
    # asl_data should have dimensions of (x, y, z, dyn, asltype)
    # legacy has (x, y, asltype, z, dyn)
    
    # m0_data should have dimensions of dimensions of (x, y, z)
    # true of legacy too
    
    
    # sliding window average (e.g. mean(C1,C2) - L1)
    # dmb[:,:,:,i] is the sum of the ith and i+1th control image, divided by the ith label image
    # dmb necessarily has one less dynamic than asl_data
    dmb = np.zeros_like(asl_data[:,:,:,:-1,0])
    for j in range(dmb.shape[3]-1):
        dmb[:,:,:,j] = ((asl_data[:,:,:,j,0] + asl_data[:,:,:,j+1,0])) / 2 - asl_data[:,:,:,j,1]
        
    control_images = asl_data[:,:,:,:,0]
    control_mean = control_images.mean(axis=3) # this is essentially the 3d baseline intensity map
    
    msk = ge.mask_image_by_cornervals(control_mean, (15,15), 1.5)
    
    impulse_matrix = np.array([[0.1, 0.3, 0.1], [0.3, 1, 0.3], [0.1, 0.3, 0.1]]) / (0.1*4+0.3*4+1)
    m0_blur = ge.iterative_2d_convolutions(m0_data, impulse_matrix) * msk
    mdmm0b = (ge.iterative_2d_convolutions(np.mean(dmb,axis=3), impulse_matrix) / m0_blur) * msk

    # paramaterization
    fmatb = np.zeros_like(m0_data)

    #eflagb = fmatb
    slicetime = 30 # (ms)
    slicetimesec = slicetime/1000

    #dmm0 = 0.005 # delta M / M0
    cfact = 6000 # convert from ml/g/s to ml/100g/min
    gdelta_a = 0.5 # the arterial transit time (s)
    gdelta = 1.5 # the tissue transit time (s)
    # gw = 2.0 # the post-labeling delay (s) note: this should be updated by slice changed from 1.8 on 01/19/2018
    # gtau = 1.8 # the labeling duration (s) changed from 1.65 on 01/19/2018  changed again to 1.8 on 08/24/2021
    gtau = ld
    gt1ab = 1.65  # t1 of blood water (s)
    #gt1aa = 1.44  # t1 of blood water during carbogen (s)
    #lamb = 0.9 # ml/g. note the lambda (lowercase) is a python keyword for a lambda function
    #a = 0.85 # labeling efficiency
    wz = np.arange(pld, (pld+slicetimesec*m0_data.shape[2])+slicetimesec, slicetimesec) # the corrected PLD for slice timing
    
    # Initial parameters for fmincon  =  =  =  =  =  =  =  =  =  =  =  =  = 
    x0 = [50/cfact] # initial guess for cbf (ml/g/s)
    LB = 0
    UB = 200/cfact
    bounds = scipy.optimize.Bounds(LB, UB, keep_feasible=True)


    for nz in range(m0_data.shape[2]):
        for nx in range(m0_data.shape[0]):
            for ny in range(m0_data.shape[1]):
                print(f'{nx,ny,nz} / {m0_data.shape}')
                if msk[nx,ny,nz]==1:
                    
                    
                    res = scipy.optimize.minimize(fun=ge.wang_pcasl_func, x0=x0,
                                                        args=(mdmm0b[nx,ny,nz], gdelta_a, gdelta, wz[nz], gtau, gt1ab), # gw+slicetimesec*(nz) -> wz
                                                        bounds=bounds,
                                                        method='SLSQP'
                                                        )
                    '''
                    [xb,fvalb,exitflagb] = fmincon(@ge.wang_pcasl_func,
                                                   x0,
                                                   A,
                                                   b,
                                                   Aeq,
                                                   beq,
                                                   LB,
                                                   UB,
                                                   [],
                                                   ops,
                                                   mdmm0b[nx,ny,nz],
                                                   gdelta_a,
                                                   gdelta,
                                                   gw+0.030*(nz-1),
                                                   gtau,
                                                   gt1ab) # corrects for slice timing
                    '''
                    
                    voxel_cbf = res.x
                    fmatb[nx,ny,nz] = voxel_cbf[0]
                    #eflagb[nx,ny,nz] = exitflagb


                
    adjusted_fmatb = fmatb * cfact
    
    # not neccessary to manually threshold the results if the optimization method is compatible with bounds
    #adjusted_fmatb[adjusted_fmatb < LB*cfact] = LB*cfact
    #adjusted_fmatb[adjusted_fmatb > UB*cfact] = UB*cfact
        

    # fmatb is the cbf
    # need to reverse the y dim
    
    return adjusted_fmatb

