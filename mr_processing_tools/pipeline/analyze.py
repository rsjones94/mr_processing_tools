#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Sky Jones
"""

help_info = """
This script uses an execution specificiation file (especfile) to parameterize the processing of MR data.

Please see the espec_readme.md file for information on how to parameterize the pspecfile.
For the vast majority of processing needs, the especfile generated by
prepare.py is sufficient.


input:
    -e / --especfile : the path to the especfile, which is a .xlsx file.
        The especfile should be in the same folder as the acquired scans. 
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

import helpers.registration as regi
from helpers.conversion import parrec_to_nifti, unpack_dims
####

inp = sys.argv
bash_input = inp[1:]
options, remainder = getopt.getopt(bash_input, "e:h", ["especfile=", 'help'])

especfile = None

for opt, arg in options:
    if opt in ('-e', '--especfile'):
        especfile = arg
    elif opt in ('-h', '--help'):
        print(help_info)
        sys.exit()
    
assert os.path.isfile(especfile)
    
    
            
        
        
    