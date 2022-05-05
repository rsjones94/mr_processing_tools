#!/usr/bin/env python3
# -*- coding: utf-8 -*-

help_info = """
This script uses takes a TRUST image and extracts the estimated blood T2 from
the superior sagittal sinus. It then applies a number of models to quantify
venous oxygenation.

The script will create:
    1) a folder in the TRUST images directory called 'trust'
    2) an excel spreadsheet in this folder called 't2_vals.xlsx' that contains the T2 for each acquisition along with uncertainty info
    3) an excel spreadsheet in this folder called 'yv_vals.xlsx' that contains the venous oxygenation using mean T2 applied to various quantification models

input:
    -t / --trust : the path to the trust, which is a .nii.gz file
        this file should have an associated vlabels file in the same folder
        (i.e., if your image is trust.nii.gz, there should be a trust_vlabels.xlsx in the same folder)
    -c / --hct : the hematocrit as a float between 0 and 1
        optional, but if not supplied then only T2 will be quantified, not Yv
    -s / --hbs : the hemoglobin s fraction as a float between 0 and 1
        optional, but if not supplied then some Yv models can not be used
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
options, remainder = getopt.getopt(bash_input, "t:c:s:h", ["trust=", 'hct=', 'hbs=', 'help'])


hct = None
hbs = None
for opt, arg in options:
    if opt in ('-t', '--trust'):
        trust_file = arg
    if opt in ('-c', '--hct'):
        hct = float(arg)
    if opt in ('-s', '--hbs'):
        hbs = float(arg)
    elif opt in ('-h', '--help'):
        print(help_info)
        sys.exit()

assert os.path.isfile(trust_file)

trust_base = os.path.basename(os.path.normpath(trust_file))   
trust_dir = os.path.dirname(os.path.normpath(trust_file))   
trust_core = trust_base.split('.')[0]
trust_vlabels_base = f'{trust_core}_vlabels.xlsx'
trust_vlabels_loc = os.path.join(trust_dir, trust_vlabels_base)

 
trust_loaded = nib.load(trust_file)
vlabels_loaded = pd.read_excel(trust_vlabels_loc)

trust_unpacked = XXX
