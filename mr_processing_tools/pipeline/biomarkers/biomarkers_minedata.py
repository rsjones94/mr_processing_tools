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
from fsl.wrappers import fslreorient2std, flirt, bet, epi_reg, fast, LOAD
import redcap
import matplotlib.pyplot as plt

import helpers.registration as regi
from helpers.conversion import parrec_to_nifti, unpack_dims

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

        
        
        
        
        
        
        
    