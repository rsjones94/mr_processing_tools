#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Sky Jones
"""

help_info = """
This script takes the spreadsheet output by biomarkers_minedata.py and pushes it to REDCap

input:
    -i / --infile : the path to the minedata.csv output by biomarkers_minedata.py
    -m / --mrid : the participant's MR ID in the Biomarker's database (e.g., Jordan_1050 or an accession number, depending on where the scan was)
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
import redcap

####

inp = sys.argv
bash_input = inp[1:]
options, remainder = getopt.getopt(bash_input, "i:m:h", ["infile=", 'mrid=', 'help'])


for opt, arg in options:
    if opt in ('-i', '--infile'):
        infile = arg
    if opt in ('-m', '--mrid'):
        mrid = arg
    elif opt in ('-h', '--help'):
        print(help_info)
        sys.exit()

in_data = pd.read_csv(infile)
cols_in_question = list(in_data['Unnamed: 0'])
new_data = list(in_data['values'])
        
keys_folder = os.path.join(four_up, '_secret_keys') # if you downloaded this repository from github you need to make this folder yourself and add the keys
# it is not under version control because it is meant to store REDCap API keys
biomarkers_token_loc = os.path.join(keys_folder, 'biomarkers.txt')

api_url = 'https://redcap.vanderbilt.edu/api/'
token = open(biomarkers_token_loc).read()

project = redcap.Project(api_url, token)
project_data_raw = project.export_records()
project_data = pd.DataFrame(project_data_raw)

id_cols = ['mri_scan_id', 'mri_scan_id_transp']

found_mrid = False
for idco_index,idco in enumerate(id_cols):
    n_matches = sum(project_data[idco]==mrid) 
    if n_matches == 1:
        found_mrid = True
        break
    elif n_matches > 1:
        raise Exception('There are more than one MRID matches for {mrid}. This should not happen. Check database.')
        
assert found_mrid


project_data_mrid = project_data.set_index(idco)
project_data_rowindex = list(project_data_mrid.index).index(mrid) # index in the original dataframe of our row of interest
project_row_mrid = project_data_mrid.loc[mrid]

sid = project_row_mrid['study_id']

old_data = list(project_row_mrid[cols_in_question])


print('\nHere is what will be updated:')
for c,old,new in zip(cols_in_question,old_data,new_data):
    if old=='':
        old='[NOTHING]'
    print(f'\t{c}: {old}  ----->  {new}')


    
do_proceed = False
valid_yes = ['y', 'yes']
valid_no = ['n',  'no']
while not do_proceed:
    feedback = input('\nIs this okay?\n\t')
    if feedback.lower() in valid_yes:
        do_proceed = True
    elif feedback.lower() in valid_no:
        raise Exception('User terminated processing')
    else:
        print('just yes or no, please')

# so this database is weird in that because the index is non-unique (repeated measures)
# we have to reupload the whole database as a list of dictionaries

project_data_new = project.export_records() # pull a fresh sheet. probably nothing has edited since last pulled, but we definitely don't want to overwrite anyone's work

# make sure that the row we're going to update corresponds to the correct SID/MRID. If not then something has gone wrong with our indexing
assert project_data_new[project_data_rowindex]['study_id'] == sid
assert project_data_new[project_data_rowindex][idco] == mrid

for c,val in zip(cols_in_question, new_data):
    project_data_new[project_data_rowindex][c] = val


np_out = project.import_records(project_data_new)


