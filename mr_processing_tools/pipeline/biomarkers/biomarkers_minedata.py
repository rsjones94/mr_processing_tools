#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Sky Jones
"""

help_info = """
This script takes the data for a processed Biomarkers participant and
collates the data into a CSV suitable for pushing to REDCap

input:
    -f / --folder : the path to the master folder
    -h / --help : brings up this helpful information. does not take an argument
"""


import os
import sys
import getopt
import shutil
import pathlib
two_up = str((pathlib.Path(__file__) / ".." / "..").resolve())
three_up = str((pathlib.Path(__file__) / ".." / ".." / "..").resolve())
four_up = str((pathlib.Path(__file__) / ".." / ".." / ".." / "..").resolve())
sys.path.append(two_up)


import pandas as pd
import numpy as np
import nibabel as nib
from fsl.wrappers import fslreorient2std, flirt, bet, epi_reg, fast, LOAD, invxfm
import redcap
import matplotlib.pyplot as plt
from nibabel import processing


####

inp = sys.argv
bash_input = inp[1:]
options, remainder = getopt.getopt(bash_input, "f:h", ["folder=", 'help'])


for opt, arg in options:
    if opt in ('-f', '--folder'):
        master_folder = arg
    elif opt in ('-h', '--help'):
        print(help_info)
        sys.exit()

        
assert os.path.isdir(os.path.normpath(master_folder))      

# we need to extract a few things:
    # the GM and WM volume from FreeSurfer
    # the GM and WM CBF from ASL
    # mean T2 from trust, and the Yv+OEF from the HbAA model
    # GM and WM Yv from ASE
    
mining_folder = os.path.join(master_folder, 'biomarkers_mining')
if not os.path.exists(mining_folder):
    os.mkdir(mining_folder)
    
ase_folder = os.path.join(master_folder, 'ase')
asl_folder = os.path.join(master_folder, 'asl')
trust_folder = os.path.join(master_folder, 'trust')
vols_folder = os.path.join(master_folder, 'vols')
t1space_folder = os.path.join(master_folder, 't1space')


variable_list = ['wm_vol', 'gm_vol', 'wm_cbf', 'gm_cbf', 'trust_t2', 'hbaa_yv', 'hbaa_oef', 'wm_yv', 'gm_yv']
variable_vals = {i:np.nan for i in variable_list}

# get the tissue volumes
# we will also resample the masks (not register) from the FS T1 space to native T1 sampling
if os.path.isdir(vols_folder):
    original_t1_loc = os.path.join(master_folder, 't1.nii.gz')
    original_t1_in = nib.load(original_t1_loc)
    # we can just pull the values straight from the vol_vals spreadsheet
    volvals_loc = os.path.join(vols_folder, 'vol_vals.xlsx')
    volvals = pd.read_excel(volvals_loc).set_index('mask')
    
    for ttype in ['wm', 'gm']:
        variable_vals[f'{ttype}_vol'] = volvals.loc[ttype]['volume_mm3'] / 1000 # convert to mL
        
        mask_loc = os.path.join(vols_folder, f'fs_mask_{ttype}.nii.gz')
        mask_in = nib.load(mask_loc)
        mask_resampled = processing.resample_from_to(mask_in, original_t1_in, order=0)
        resampled_mask_loc = os.path.join(mining_folder, f'fs_mask_{ttype}_resampled.nii.gz')
        nib.save(mask_resampled, resampled_mask_loc)
        
    
# get the CBF in each tissue type. have to coregister CBF to T1, then resample
if os.path.isdir(asl_folder) and os.path.isdir(vols_folder):
    cbf_loc = os.path.join(asl_folder, 'cbf.nii.gz')
    cbf_in = nib.load(cbf_loc)
    
    asl_to_t1_omat_loc = os.path.join(master_folder, 'aslm0_to_t1space.omat')
    
    target_t1_loc = os.path.join(master_folder, 't1.nii.gz')
    target_t1_in = nib.load(target_t1_loc)
    
    cbf_in_t1 = flirt(cbf_in, ref=target_t1_in, out=LOAD,
      applyxfm=1, init=asl_to_t1_omat_loc, interp='trilinear')['out']
    
    cbf_in_t1 = nib.Nifti1Image(cbf_in_t1.get_fdata(), target_t1_in.affine, cbf_in_t1.header)
    
    cbf_in_t1_loc = os.path.join(mining_folder, 'cbf_t1space.nii.gz')
    nib.save(cbf_in_t1, cbf_in_t1_loc)
    
    cbf = cbf_in_t1.get_fdata()
    
    # now we can overlay the masks and cbf
    for ttype in ['wm', 'gm']:
        mask_loc = os.path.join(mining_folder, f'fs_mask_{ttype}_resampled.nii.gz')
        mask_in = nib.load(mask_loc)
        mask_data = mask_in.get_fdata().astype(int)
        
        cbf_masked = cbf.copy()
        cbf_masked[mask_data!=1] == np.nan
        meaned = np.nanmean(cbf_masked)
        
        variable_vals[f'{ttype}_cbf'] = meaned
    
    
    
    
    
    
    
    
    
    
    
    
    