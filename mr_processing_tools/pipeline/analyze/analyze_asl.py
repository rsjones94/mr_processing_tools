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
    -m / --m0 : the path to the m0, which is a .nii.gz file
    -p - / --pld : the post-labeling delay in ms
    -l - / --ld : the labeling duration in ms
    -c - / --hct : hematocrit
        optional. used to calculate blood t1
        if not specified, blood t1 is assumed to be 1650ms (hct~=0.4345)
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
from helpers.general import read_xfm, calculate_blood_t1
from processing.asl_funcs import quantify_cbf_from_asl
####


inp = sys.argv
bash_input = inp[1:]
options, remainder = getopt.getopt(bash_input, "a:m:p:l:t:h", ["asl=", 'm0=', 'pld=', 'ld=', 'tr=', 'help'])

blood_t1 = 1650
for opt, arg in options:
    if opt in ('-a', '--asl'):
        asl_file = arg
    if opt in ('-m', '--m0'):
        m0_file = arg
    if opt in ('-p', '--pld'):
        pld = float(arg)
        pld_s = pld/1000
    if opt in ('-l', '--ld'):
        ld = float(arg)
        ld_s = pld/1000
    if opt in ('-c', '--hct'):
        hct = float(arg)
        blood_t1 = calculate_blood_t1(hct)
    elif opt in ('-h', '--help'):
        print(help_info)
        sys.exit()

blood_t1_s = blood_t1/1000

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
cbf_map = quantify_cbf_from_asl(asl_unpacked, m0_in_data, pld_s, ld_s, blood_t1_s)
print('Finished quantification')

cbf_out_loc = os.path.join(asl_loc, 'cbf.nii.gz')
cbf_out = nib.Nifti1Image(cbf_map, asl_in.affine, asl_in.header)

nib.save(cbf_out, cbf_out_loc)
