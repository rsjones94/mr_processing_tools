#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Sky Jones
"""

help_info = """
This script finds the relevant lab values, etc. for a particular scan in the Biomarker's database
Note that it just pulls the data for one particular scan date for a patient, not the data for all scans that the patient has.

Writes the data to a file called redcap_pull.xlsx in the specified folder.

input:
    -m / mrid : the participant's MR ID in the Biomarker's database (e.g., Jordan_1050 or an accession number, depending on where the scan was)
    -f / --folder : the folder to write redcap_pull.xlsx to
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
options, remainder = getopt.getopt(bash_input, "f:s:h", ["folder=", 'studyid=', 'help'])


for opt, arg in options:
    if opt in ('-f', '--folder'):
        out_folder = arg
    if opt in ('-s', '--studyid'):
        mrid = arg
    elif opt in ('-h', '--help'):
        print(help_info)
        sys.exit()


common_columns = ['study_group', 'mh_heme_type']
common_columns_names = ['study_group', 'heme_type']
# the column columns will only have data for the first study_id

columns_of_interest_transp = ['studylab_hct_transp', 'studylab_hplcs_transp', 'spo2_trust_transp']
columns_of_interest_nontransp = ['hct_studylabs', 'hgb_s_studylabs', 'spo2_trust']
# the case and control columns are named slightly differently but should correspond to the same type of data. e.g., the hbs frac or transplant is under studylab_hplcs_transp, but for controls it's hgb_s_studylabs
# the columns_of_interest_names variablesi is used to unify the naming convention for ease of processing later
columns_of_interest_names = ['hct', 'hbs_frac', 'oxsat']
        
keys_folder = os.path.join(four_up, '_secret_keys') # if you downloaded this repository from github you need to make this folder yourself and add the keys
# it is not under version control because it is meant to store REDCap API keys
biomarkers_token_loc = os.path.join(keys_folder, 'biomarkers.txt')

api_url = 'https://redcap.vanderbilt.edu/api/'
token = open(biomarkers_token_loc).read()

project = redcap.Project(api_url, token)
project_data_raw = project.export_records()
project_data = pd.DataFrame(project_data_raw)

id_cols = ['mri_scan_id', 'mri_scan_id_transp']
interest_cols = [columns_of_interest_nontransp, columns_of_interest_transp]

found_mrid = False
for idco_index,(idco, coi) in enumerate(zip(id_cols, interest_cols)):
    n_matches = sum(project_data[idco]==mrid) 
    if n_matches == 1:
        found_mrid = True
        break
    elif n_matches > 1:
        raise Exception('There are more than one MRID matches for {mrid}. This should not happen. Check database.')
        
assert found_mrid


project_data_mrid = project_data.set_index(idco)
project_row_mrid = project_data_mrid.loc[mrid]

sid = project_row_mrid['study_id']

project_data_common = project_data.set_index('study_id')
project_row_common = project_data_common.loc[sid].iloc[0]
common_data = project_row_common[common_columns]
the_data_mrid = project_row_mrid[coi]

all_data = common_data.append(the_data_mrid)
all_renames = []
all_renames.extend(common_columns_names)
all_renames.extend(columns_of_interest_names)

the_data_renamed = all_data.rename({old:new for old,new in zip(list(all_data.index), all_renames)})


conv_to_frac = ['hct', 'hbs_frac', 'oxsat']
for cf in conv_to_frac:
    try:
        the_data_renamed[cf] = float(the_data_renamed[cf]) / 100
    except ValueError:
        pass


print('I found these values:')
for key,val in the_data_renamed.items():
    print(f'\t{key} : {val}')


write_loc = os.path.join(out_folder, 'redcap_pull.xlsx')
the_data_renamed.to_excel(write_loc)