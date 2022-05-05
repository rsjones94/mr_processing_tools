#!/usr/bin/env python3
# -*- coding: utf-8 -*-

help_info = """
This script uses takes a T1-weighted MR image of the brain and calls FreeSurfer
to segment its tissue types. It then uses the resulting data to make a simple
tabulation of some key values.

The script will create:
    1) a folder in the T1's directory called 'vols' which will contain the entirety of the FreeSurfer results as a subfolder called 'fs'
    2) an excel spreadsheet in this folder called 'vol_vals.xlsx' that contains the simple tabulation.
    3) The following NiFTis in Talairach space
        fs_t1.nii.gz: the T1
        fs_brain.nii.gz: the brain extracted T1
        fs_mask_gm.nii.gz: gray matter mask
        fs_mask_wm.nii.gz: white matter mask
    4) a .omat file specifying the linear transformation from the native T1 space to Talairach space

input:
    -t / --t1 : the path to the t1, which is a .nii.gz file.
        This should be a raw T1, not a brain-extracted T1
    -s / --subjid : subject id. optional
        default is 'default_subjid'
        does not change anything other than the name of the temporary FreeSurfer folder
        in $SUBJECTS_DIR
    -f / --freesurfer : whether to run FreeSurfer. 0 (do not run) or 1 (run)
        default is 1
        If FreeSurfer is not run, it is assumed that this script has already been used to run FS and
        thus the folder created by recon-all has been moved to the directory with the T1 already
        and is called 'fs'
    -h / --help : brings up this helpful information. does not take an argument
"""


import os
import sys
import getopt
import glob
import shutil
import pathlib
two_up = str((pathlib.Path(__file__) / ".." / "..").resolve())
three_up = str((pathlib.Path(__file__) / ".." / ".." / "..").resolve())
four_up = str((pathlib.Path(__file__) / ".." / ".." / ".." / "..").resolve())
sys.path.append(three_up)


import pandas as pd
import numpy as np
import nibabel as nib
from fsl.wrappers import fslreorient2std, flirt, bet, epi_reg, fast, LOAD

import helpers.registration as regi
from helpers.conversion import parrec_to_nifti, unpack_dims
from helpers.general import read_xfm
####

SUBJECTS_DIR = os.environ['SUBJECTS_DIR']
subjid = 'default_subjid'
run_fs = 1

inp = sys.argv
bash_input = inp[1:]
options, remainder = getopt.getopt(bash_input, "t:s:f:h", ["t1=", 'subjid=', 'freesurfer=', 'help'])


for opt, arg in options:
    if opt in ('-t', '--t1'):
        t1_file = arg
    if opt in ('-s', '--subjid'):
        subjid = arg
    if opt in ('-f', '--freesurfer'):
        run_fs = int(arg)
    elif opt in ('-h', '--help'):
        print(help_info)
        sys.exit()

assert os.path.isfile(t1_file)
    

t1_basename = os.path.basename(os.path.normpath(t1_file))
t1_parent = os.path.dirname(os.path.normpath(t1_file))

vol_loc = os.path.join(t1_parent, 'vols')
if not os.path.exists(vol_loc):
    os.mkdir(vol_loc)


fs_orig_loc = os.path.join(SUBJECTS_DIR, subjid)
fs_out_loc = os.path.join(vol_loc, 'fs')

seg_loc = os.path.join(fs_out_loc, 'mri', f'aparc+aseg.mgz')
mask_gm_loc = os.path.join(fs_out_loc, 'mri', f'mask_gm.mgz')
mask_wm_loc = os.path.join(fs_out_loc, 'mri', f'mask_wm.mgz')

call_initialize = f'recon-all -i {t1_file} -subjid {subjid}'
call_process = f'recon-all -all -subjid {subjid}'
call_gmmask = f'mri_binarize --i {seg_loc} --gm --o {mask_gm_loc}'
call_wmmask = f'mri_binarize --i {seg_loc} --all-wm --o {mask_wm_loc}'

if run_fs:
    print(f'Initializing recon-all')
    print(f'\t{call_initialize}')
    os.system(call_initialize)
    
    print(f'Executing recon-all')
    print(f'\t{call_process}')
    os.system(call_process)
    
    print(f'Moving FS data')
    shutil.move(fs_orig_loc, fs_out_loc)


print('Converting Talairach-space images from mgz to nii and binarizing masks')

# need to convert the intensity-normalized T1, the skullstripped T1
# also need to convert the transformation matrix to standard .omat


os.system(call_gmmask)
os.system(call_wmmask)

fnames = ['brain', 't1', 'mask_gm', 'mask_wm']
for fname in fnames:

    mgz_loc = os.path.join(fs_out_loc, 'mri', f'{fname}.mgz')
    nii_loc = os.path.join(vol_loc, f'fs_{fname}.nii.gz'.lower())
    
    loaded_mgz = nib.load(mgz_loc)
    loaded_nii = nib.Nifti1Image(loaded_mgz.get_fdata(), loaded_mgz.affine, loaded_mgz.header)
    loaded_nii = fslreorient2std(loaded_nii, output=LOAD)['output']
    
    nib.save(loaded_nii, nii_loc)

transformation_xfm_loc = os.path.join(vol_loc, 'fs', 'mri', 'transforms', 'talairach.xfm')
transformation_omat_loc = os.path.join(vol_loc, 't1_to_talairach.omat')

xfm = read_xfm(transformation_xfm_loc, out_loc=transformation_omat_loc)


print(f'Tabulating volumes')
mask_globber = os.path.join(vol_loc, 'fs_mask_*.nii.gz')
mask_locs = glob.glob(mask_globber, recursive=False)

tabulation_loc = os.path.join(vol_loc, 'vol_vals.xlsx')
cols = ['mask', 'volume_ml', 'n_voxels', 'voxelvol_ml']
tabulation_df = pd.DataFrame(columns=cols)

for mask_loc in mask_locs:
    
    mask_ser = pd.Series(dtype=np.float64)
    
    mask_base = os.path.basename(os.path.normpath(mask_loc))
    mask_wo_ext = mask_base.split('.')[0]
    mask_name = mask_wo_ext.split('_')[-1]
    
    loaded = nib.load(mask_loc)
    zooms = loaded.header.get_zooms()
    voxel_vol = np.product(zooms) # in mL
    
    data = loaded.get_fdata()
    n_voxels = data.sum()
    total_vol = n_voxels * voxel_vol # in mL
    
    mask_ser['mask'] = mask_name
    mask_ser['volume_ml'] = total_vol
    mask_ser['n_voxels'] = n_voxels
    mask_ser['voxelvol_ml'] = voxel_vol
    
    tabulation_df = tabulation_df.append(mask_ser, ignore_index=True)
    
    
tabulation_df.to_excel(tabulation_loc, index=False)

    