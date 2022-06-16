#!/usr/bin/env python3
# -*- coding: utf-8 -*-

help_info = """
This script uses takes a source pCASL image and m0 of the brain and applies a single-PLD,
two-compartment kinetic model to the subtraction to quantify CBF in mL/100g/min

The script will create:
    1) a folder in the pCASL's directory called asl which will contain
        all of the following results
    2) a nifti called cbf.nii.gz containing CBF in ml/100g/min

input:
    -a / --asl : the path to the asl, which is a .nii.gz file
        There should be an associated vlabels xlsx
    -m / --aslm0 : the path to the m0, which is a .nii.gz file
    -p - / --pld : the post-labeling delay in ms
    -l - / --ld : the labeling duration in ms
    -c - / --hct : hematocrit
        optional. used to calculate blood t1
        if not specified, blood t1 is assumed to be 1650ms (hct~=0.4345)
    -e - / --labeff : labeling efficiency
        optional. if not specified, assumed to be 0.91
        For BOLD processing 0.85 is standard
        For SCD, 0.72
        Otherwise 0.91
    -t - / --ttt : the tissue transit time in seconds
        optional. if not specified assumed to be 1.4
        For BOLD processing 1.500 is standard
        For SCD, 1.290
        Otherwise 1.400
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
import matplotlib.pyplot as plt

import helpers.registration as regi
from helpers.conversion import parrec_to_nifti, unpack_dims
from helpers.general import read_xfm, calculate_blood_t1
from processing.asl_funcs import quantify_cbf_from_asl
####
inp = sys.argv
bash_input = inp[1:]
options, remainder = getopt.getopt(bash_input, "a:m:p:l:c:e:t:h", ["asl=", 'aslm0=', 'pld=', 'ld=', 'hct=', 'labeff=', 'ttt=', 'help'])

blood_t1_1 = 1.650
labeling_efficiency = 0.91
tissue_transit_time = 1.4
for opt, arg in options:
    if opt in ('-a', '--asl'):
        asl_file = arg
    if opt in ('-m', '--aslm0'):
        m0_file = arg
    if opt in ('-p', '--pld'):
        pld = float(arg)
        pld_s = pld/1000
    if opt in ('-l', '--ld'):
        ld = float(arg)
        ld_s = ld/1000
    if opt in ('-c', '--hct'):
        hct = float(arg)
        blood_t1_s = calculate_blood_t1(hct)
    if opt in ('-t', '--ttt'):
        tissue_transit_time = float(arg)
    if opt in ('-e', '--labeff'):
        labeling_efficiency = float(arg)
    elif opt in ('-h', '--help'):
        print(help_info)
        sys.exit()

assert os.path.isfile(asl_file)
assert os.path.isfile(m0_file)
    

asl_basename = os.path.basename(os.path.normpath(asl_file))
asl_parent = os.path.dirname(os.path.normpath(asl_file))

asl_loc = os.path.join(asl_parent, 'asl')
if not os.path.exists(asl_loc):
    os.mkdir(asl_loc)
    
asl_in = nib.load(asl_file)
m0_in = nib.load(m0_file)

asl_in_data = asl_in.get_fdata()
m0_in_data = m0_in.get_fdata()
 
asl_core = asl_basename.split('.')[0]
asl_vlabels_base = f'{asl_core}_vlabels.xlsx'
asl_vlabels_loc = os.path.join(asl_parent, asl_vlabels_base)

vlabels_loaded = pd.read_excel(asl_vlabels_loc, index_col=0)
ordering = ['dynamic scan number', 'label type']
vlabels_ordered = vlabels_loaded[ordering]

asl_unpacked = unpack_dims(asl_in_data, vlabels_ordered)

print('Quantifying CBF....')
cbf_map = quantify_cbf_from_asl(asl_unpacked, m0_in_data, pld_s, ld_s,
                                blood_t1=blood_t1_s, ttt=tissue_transit_time,
                                label_eff=labeling_efficiency)
print('Finished quantification')

cbf_out_loc = os.path.join(asl_loc, 'cbf.nii.gz')
cbf_out = nib.Nifti1Image(cbf_map, asl_in.affine, asl_in.header)

nib.save(cbf_out, cbf_out_loc)

### qc
qc_out_loc = os.path.join(asl_loc, 'asl_qualitycontrol.png')
fig, axs = plt.subplots(nrows=3, ncols=1, figsize=(6,12))

pdims = asl_in.header['pixdim'][1:4] # relevant for setting pixel aspect ratio for dims

cbf_shape = cbf_map.shape
x_half = int(cbf_shape[0]/2)
y_half = int(cbf_shape[1]/2)
z_half = int(cbf_shape[2]/2)

ax_slice = cbf_map[:,:,z_half]
cor_slice = cbf_map[:,y_half,:]
sag_slice = cbf_map[x_half,:,:]

ax_aspect = pdims[0] / pdims[1]
cor_aspect = pdims[0] / pdims[2]
sag_aspect = pdims[1] / pdims[2]

slices = [ax_slice, sag_slice, cor_slice]
slice_names = ['ax', 'sag', 'cor']
aspects = [ax_aspect, sag_aspect, cor_aspect]

minner = 0
maxer = np.percentile(cbf_map, 98)

for sli, sln, ax, asp in zip(slices, slice_names, axs, aspects):
    pos = ax.imshow(sli, cmap='magma', aspect=asp, vmin=minner, vmax=maxer)
    ax.set_title(sln)

cbar = fig.colorbar(pos, ax=axs[0])
cbar.set_label('CBF (ml/100g/min)')

fig.tight_layout()
fig.savefig(qc_out_loc, dpi=400)