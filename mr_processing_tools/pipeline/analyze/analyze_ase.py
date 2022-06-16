#!/usr/bin/env python3
# -*- coding: utf-8 -*-

help_info = """
This script uses takes an ASE image and creates a voxelwise OEF map
Note that the ASE processing is natively done in MATLAB

The script will create:
    1) a folder in the ASE image directory called 'ase'
    2) four files in the above folder:
        ase_R2prime.nii.gz
        ase_rOEF.nii.gz
        ase_Rsquared.nii.gz
        ase_rvCBV.nii.gz

required input:
    -a / --ase : the path to the ASE image, which is a .nii.gz file
        there is no need for vlabels file, as the volume labels are assumed
    -c / --hct : the hematocrit as a float between 0 and 1
        default is 0.42
    -b / --artox : the arterial oxygen saturation as a float between 0 and 1
        default is 0.98
    -x / --tr : the repetition time in ms
        default is 4400
    -y / --te1 : the first echo time in ms
        default is 64
    -z / --te2 : the second echo time in ms
        default is 107
    -r / --rfon : whether the scan was acquired with RFon
        default is 1
    -v / --vasoase : whether the scan is VasoASE
        default is 0
    -h / --help : brings up this helpful information. does not take an argument
"""


import os
import sys
import getopt
import shutil
import pathlib
import subprocess

import numpy as np
import matplotlib.pyplot as plt
import nibabel as nib
from scipy.stats import mode

two_up = str((pathlib.Path(__file__) / ".." / "..").resolve())
three_up = str((pathlib.Path(__file__) / ".." / ".." / "..").resolve())
four_up = str((pathlib.Path(__file__) / ".." / ".." / ".." / "..").resolve())
sys.path.append(three_up)
####

ase_funcs_loc = os.path.join(three_up, 'processing', 'SLW_ASEproc_v4')
sys.path.append(ase_funcs_loc)



inp = sys.argv
bash_input = inp[1:]
options, remainder = getopt.getopt(bash_input, "a:c:b:x:y:z:r:v:m:h",
                                   ["ase=", "hct=", "artox=", "tr=", "te1=", "te2=", "rfon=", "vasoase=", "matlab=" 'help'])

TR = 4400
TE1 = 64
TE2 = 107

hctUncorr = 0.42
artox=0.98

RFon = 1
vasoase=0

matlab_spec_file = os.path.join(four_up, 'bin', 'matlab_location.txt')

#matlab_loc = '/Applications/MATLAB_R2020b.app/'
matlab_loc = open(matlab_spec_file).read()

for opt, arg in options:
    if opt in ('-a', '--ase'):
        ase_file = arg
    if opt in ('-c', '--hct'):
        hctUncorr = float(arg)
    if opt in ('-b', '--artox'):
        artox = float(arg)
    if opt in ('-x', '--tr'):
        TR = float(arg)
    if opt in ('-y', '--te1'):
        TE1 = float(arg)
    if opt in ('-z', '--te2'):
        TE2 = float(arg)
    if opt in ('-r', '--rfon'):
        RFon = int(arg)
    if opt in ('-v', '--vasoase'):
        vasoase = int(arg)
    if opt in ('-m', '--matlab'):
        matlab_loc = arg
    elif opt in ('-h', '--help'):
        print(help_info)
        sys.exit()

assert os.path.isfile(ase_file)
assert type(artox) == float

# create a mock pt folder for the MATLAB function to work with
temp_folder = os.path.join(ase_funcs_loc, 'DataDir', 'tempproc')
acq_folder = os.path.join(temp_folder, 'Acquired')

orig_dir = os.path.dirname(os.path.normpath(ase_file))

ase_base = 'ASE.nii.gz'
ase_base_noext = ase_base.split('.')[0]

temp_ase = os.path.join(acq_folder, ase_base)


if os.path.exists(temp_folder):
    shutil.rmtree(temp_folder)
    
os.mkdir(temp_folder)
os.mkdir(acq_folder)
shutil.copy(ase_file, temp_ase)

ase_loaded = nib.load(temp_ase)
ase_dims = str(list(ase_loaded.header['pixdim'][1:4])).replace(',', '')


pt_dir = f"'{temp_folder}'"
#hctUncorr = 0.42 # Uncorrected hematocrit. Estimated Value 0.42
slicegrab = 0
slice2grab = 0
echo2Exclude = 0
tauvec = '[19 20 6.5 8.5 1 16.5 1.5 14.5 3.5 13 9 19.5 17 2.5 4.5 12.5 2 9.5 6 14 7.5 5.5 12 8 11 0.5 16 10.5 13.5 5 10 0 18.5 17.5 3 11.5 7 18 15.5 15 4]'
#TR = 4400 # ms
#TE1 = 64 # ms
#TE2 = 107
excludedyns_1echo = []
excludedyns = []
B0 = 3 # Tesla
#ASEresolution = '[1.72, 1.72, 3]'
ASEresolution = f'{ase_dims}'
#RFon=1
#vasoase=0
#function [ ASE_doneflag ] = ASEcheck_v7_fn_fornifti(patdir,hctUncorr,sliceGrab,slice2grab,echo2Exclude,tauvec,TR,TE1,TE2,excludedyns_1echo,excludedyns,B0,ASEresolution,RFon,vasoase)

func_args = [pt_dir,
            hctUncorr,
            slicegrab,
            slice2grab,
            echo2Exclude,
            tauvec,
            TR,
            TE1,
            TE2,
            excludedyns_1echo,
            excludedyns,
            B0,
            ASEresolution,
            RFon,
            vasoase
            ]

func_args_asstr = ', '.join([str(i) for i in func_args])

processing_input = (
                    f'''{matlab_loc}/bin/matlab -nodesktop -nosplash -r'''
                    f''' "cd('{ase_funcs_loc}'); ASEcheck_v7_fn_fornifti({func_args_asstr})"'''
                    )



print(f'Calling MATLAB: {processing_input}')
os.system(processing_input)


output_suffixes = ['R2prime', 'rOEF', 'Rsquared', 'rvCBV']
cmaps = ['inferno', 'viridis', 'Reds', 'plasma']

target_dir = os.path.join(orig_dir, 'ase')
if not os.path.exists(target_dir):
    os.mkdir(target_dir)

# the MATLAB script saves files with no spatial information, so we need to reassociate it
# we can also make our qc plot at the same time
qc_plot_loc = os.path.join(target_dir, 'ase_qualitycontrol.png')
fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(12,12))
axs = np.ravel(axs)
for suf,ax,cmap in zip(output_suffixes,axs,cmaps):
    file_base = f'{ase_base_noext}_{suf}.nii.gz'
    orig = os.path.join(temp_folder, 'Results', file_base)
    target_file = os.path.join(target_dir, file_base)
    
    orig_loaded = nib.load(orig)
    the_data = orig_loaded.get_fdata()
    if suf=='rOEF':
        # the OEF was calculated assuming Ya=1. This is not true, so we need to adjust the data
        yv = 1 - the_data
        oef = (artox - yv) / artox
        the_data = oef
        
        # we can also use this to make a Yv map
        yv_im = nib.Nifti1Image(yv, ase_loaded.affine, header=ase_loaded.header)
        yv_base = f'{ase_base_noext}_yv.nii.gz'
        target_yv = os.path.join(target_dir, yv_base)
        nib.save(yv_im, target_yv)
        
    
    new_im = nib.Nifti1Image(the_data, ase_loaded.affine, header=ase_loaded.header)
    
    bgval = mode(np.ravel(the_data))[0][0]
    
    the_data_nanned = the_data.copy()
    the_data_nanned[np.isclose(the_data,bgval)] = np.nan
    
    minner = np.nanpercentile(the_data_nanned,5)
    maxxer = np.nanpercentile(the_data_nanned,95)
    
    data_shape = the_data.shape
    z_ind = int(data_shape[2]/2)
    the_data_slice = the_data[:,:,z_ind]
    
    pos = ax.imshow(the_data_slice, cmap=cmap, vmin=minner, vmax=maxxer)
    ax.set_title(suf)
    
    cbar = fig.colorbar(pos, ax=ax)
    #cbar.set_label(suf)
    
    nib.save(new_im, target_file)

fig.tight_layout()
fig.savefig(qc_plot_loc, dpi=400)


