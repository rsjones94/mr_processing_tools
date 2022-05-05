#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Sky Jones
"""

help_info = """
This script uses a pattern specificiation file (pspecfile) to search for specific images from
an acquisition, catalog them, convert them, standardize their spaces, and determine
what processing sequences can be executed.

Please see the pspec_readme.md file for information on how to parameterize the pspecfile.
For the vast majority of processing needs, the default pspecfile is sufficient.

The following steps are possible, and will be carried out in this order:
    
    1) find: search for each file type/pattern in the pspecfile and log it
    2) convert: convert each file to NiFTi if it is not a NiFTi already, and reorient to standard anatomical orientation
    3) warp: warp files to T1 space
    4) determine: determine what processing sequences can be executed and writes an execution specification file (especfile)

input:
    -p / --pspecfile : the path to the pspecfile, which is a .xlsx file.
        The pspecfile should be in the same folder as the acquired scans. 
        If the pspecfile does not exist, a template specfile will be written to the input path and used.
        If the argument passed is a folder, then the pspecfile name will be imputed as 'pattern_specs.xlsx'.
    -f / --find : whether to run the "find" step of processing. 0 (do not run) or 1 (run)
        Default is 1.
    -c / --convert : whether to run the "convert" step of processing. 0 (do not run) or 1 (run)
        Default is 1.
    -w / --warp : whether to run the "warp" step of processing. 0 (do not run) or 1 (run)
        Default is 1.
    -d / --determine : whether to run the "determine" step of processing. 0 (do not run) or 1 (run)
        Default is 1.
        the especfile is always saved as 'execution_specs.xlsx' in the same folder as the pspecfile
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
options, remainder = getopt.getopt(bash_input, "p:f:c:w:d:h", ["pspecfile=", "find=", 'convert=', 'warp=', 'determine=', 'help'])



seek = 1
convert = 1
reg_to_t1 = 1
determine = 1

for opt, arg in options:
    if opt in ('-s', '--specfile'):
        in_file = arg
    elif opt in ('-f', '--find'):
        seek = arg
    elif opt in ('-c', '--convert'):
        convert = arg
    elif opt in ('-w', '--warp'):
        reg_to_t1 = arg
    elif opt in ('-d', '--determine'):
        determine = arg
    elif opt in ('-h', '--help'):
        print(help_info)
        sys.exit()

        
if os.path.isdir(in_file):
    in_file = os.path.join(in_file, 'pattern_specs.xlsx')
    
try:
    assert os.path.isfile(in_file)
except AssertionError:
    print('WARNING: input pspecfile does not exist: a default pspecfile will be loaded instead')
    blank_specfile = os.path.join(three_up, 'bin', 'pattern_specs.xlsx')
    shutil.copyfile(blank_specfile, in_file)
    
    



####

pattern_file = in_file
data_folder = pathlib.Path(pattern_file).parent.absolute()
native_folder = data_folder.parent.absolute()


# finds files matching the patterns specified in your patterns.xlsx
if seek:
    print('\nRunning SEEK')
    loaded_patterns = pd.read_excel(pattern_file, dtype=str)
    match_cols = [j for j in loaded_patterns.columns if 'match' in j.lower()]
    exclude_cols = [j for j in loaded_patterns.columns if 'exclude' in j.lower()]
    
    add_cols = ['File found', 'N found']
    for ac in add_cols:
        loaded_patterns[ac] = None
    
    #sorted_patterns = loaded_patterns.sort_values(by='T1 Registration method')
    for i,row in loaded_patterns.iterrows():
        
        if row['Execute']!='1':
            continue
        
        match_expressions = [row[k] for k in match_cols if isinstance(row[k], str)]
        exclude_expressions = [row[k] for k in exclude_cols if isinstance(row[k], str)]
        
        matched, n, others = regi.find_image(base_folder=data_folder,
                                                     match=match_expressions,
                                                     exclude=exclude_expressions,
                                                     recursive=False)
        
        loaded_patterns.at[i,'File found'] = matched
        loaded_patterns.at[i,'n found'] = n
        
    loaded_patterns.to_excel(pattern_file, index=False)

# converts to NiFTi if not NiFTi already
# also assigns standard anatomical labels
if convert:
    print('\nRunning CONVERT')
    loaded_patterns = pd.read_excel(pattern_file, dtype=str)
    strips = ['.par', '.rec', '.nii', '.gz']
    for i,row in loaded_patterns.iterrows():
        if row['Execute']!='1':
            continue
        
        in_file = row['File found']
        scan_type = row['Scan name']
        if not isinstance(in_file, str):
            continue
        
        scale_type = row['Scaling']
        
        
        orig_bn = os.path.basename(os.path.normpath(in_file))
        bn = orig_bn
        for st in strips:
            bn = bn.replace(st.lower(), '')
            bn = bn.replace(st.upper(), '')
        
        new_bn = f'{scan_type}.nii.gz'
        vlabels = f'{scan_type}_vlabels.xlsx'
        
        new_file_loc = os.path.join(native_folder, new_bn)
        vlabels_loc = os.path.join(native_folder, vlabels)
        
        print(f'Converting and normalizing indices for {scan_type} ({orig_bn}) ({i+1} of {len(loaded_patterns)})')
        
        if '.nii' in orig_bn.lower():
            loaded = nib.load(in_file)
            reoriented = fslreorient2std(loaded, output=LOAD)['output']
            reoriented.to_filename(new_file_loc)
            # no vlabels for niftis
        elif '.par' in orig_bn.lower():
            loaded = nib.load(in_file, scaling=scale_type)
            if row['Has 4D'] != '1':
                vlabels_loc = None
            parrec_to_nifti(loaded, out_nifti=new_file_loc,
                             out_vlabels=vlabels_loc, reorient=True)
        else:
            raise NotImplementedError('Filetype not supported')
            
                    
# converts the files to T1 space (if specified)
if reg_to_t1:
    print('\nRunning WARP')
    loaded_patterns = pd.read_excel(pattern_file, dtype=str)
    sorted_patterns = loaded_patterns.sort_values(by='T1 Registration method')
    # sorting lets us register files that will use a prior registration's .omat
    # (rather than calculating their own) last
    
    t1_row = sorted_patterns[sorted_patterns['Is T1']=='1'].iloc[0]
    t1_name = t1_row['Scan name']
    f_val = float(t1_row['f'])
    t1_location = os.path.join(native_folder, f'{t1_name}.nii.gz')
    t1_location_betted = os.path.join(native_folder, f'{t1_name}_BET.nii.gz')
    fast_location_betted = os.path.join(native_folder, f'{t1_name}_BET_FAST')
    t1_wm_location = os.path.join(native_folder, f'{t1_name}_BET_FAST_pve_2.nii.gz')
    
    # first thing: BET the T1
    print('Primary: extracting brain from T1')
    bet(t1_location, t1_location_betted, f=f_val)
    print('Primary: FASTing T1')
    fast(t1_location_betted, out=fast_location_betted, n_classes=3)
    
    t1space_folder = os.path.join(native_folder, 't1space')
    if not os.path.exists(t1space_folder):
        os.mkdir(t1space_folder)
    
    
    for i,row in sorted_patterns.iterrows():
        
        reg_method = row['T1 Registration method']
        if row['Execute']!='1':
            continue
        if row['Is T1']=='1':
            # we already did what needs to be done to the T1
            continue
        in_file = row['File found']
        scan_type = row['Scan name']
        if not isinstance(in_file, str):
            continue
        
        print(f'\tSecondary: registering {scan_type}')
        
        scan_location = os.path.join(native_folder, f'{scan_type}.nii.gz')
        omat_location = os.path.join(native_folder, f'{scan_type}_to_t1space.omat')
        
        scan_location_in_t1space = os.path.join(t1space_folder, f'{scan_type}_t1space.nii.gz')
        
        if reg_method=='0': # file is already in T1 space
            continue
        elif reg_method=='1': # register using FLIRT
            scan_location_betted = os.path.join(native_folder, f'{scan_type}_BET.nii.gz')
            bet(scan_location, scan_location_betted)
            flirt(src=scan_location_betted, ref=t1_location_betted, omat=omat_location)
        elif reg_method=='2': # register using epi_reg
            base_out = omat_location.replace('.omat', '')
            mat_location = f'{base_out}.mat'
            epi_reg(epi=scan_location, t1=t1_location, t1brain=t1_location_betted, out=base_out, wmseg=t1_wm_location)
            os.rename(mat_location, omat_location) # for whatever reason epi_reg makes a .mat file, not .omat. their contents are equivalent
            
            to_remove = ['_fast_wmedge.nii.gz', '_fast_wmseg.nii.gz', '_init.mat', '.nii.gz']
            for tr in to_remove: # some intermediate files we want to clean up
                file_to_remove = os.path.join(native_folder, f'{base_out}{tr}')
                os.remove(file_to_remove)
        elif reg_method=='-1': # file should not be converted to T1 space
            continue
        else: 
            # if parameters are correctly specified, then this branch
            # handles cases where reg_method is the name of a scan
            # in this case, we use that scan's omat instead of generating a new one
            reference_omat = os.path.join(native_folder, f'{reg_method}_to_t1space.omat')
            shutil.copyfile(reference_omat, omat_location)
            
        
        flirt(scan_location, ref=t1_location, out=scan_location_in_t1space,
              applyxfm=1, init=omat_location, interp='trilinear')
    
    
if determine:
    print('\nRunning DETERMINE')
    espec_template = os.path.join(three_up, 'bin', 'execution_specs.xlsx')
    
    loaded_patterns = pd.read_excel(pattern_file, dtype=str)
    
    
            
        
        
    