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
    -s / --special : special processing params
        As of right now, the only special param is biomarkers,[mrid]
            This pulls data from the SCD biomarker study database using mrid as an indexer
            and also sets ASL labeling efficiency and tissue transit times appropriately
            depending on whether the participant is a control or has SCD
            
            Alternatively, biomarkers,control will use parameters for a generic control,
            and biomarkers,scd will use parameters for a generic participant with SCD.
            As of now, all that is changed is that biomarkers,scd will set the labeling efficiency to
            0.72 (from 0.91), the ttt to 1.29 (from 1.4) and the hbs to 0.5 (from 0)
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
options, remainder = getopt.getopt(bash_input, "p:f:c:w:d:s:h", ["pspecfile=", "find=", 'convert=', 'warp=', 'determine=', 'special=', 'help'])



seek = 1
convert = 1
reg_to_t1 = 1
determine = 1

special = ''
for opt, arg in options:
    if opt in ('-p', '--pspecfile'):
        in_file = arg
    elif opt in ('-f', '--find'):
        seek = int(arg)
    elif opt in ('-c', '--convert'):
        convert = int(arg)
    elif opt in ('-w', '--warp'):
        reg_to_t1 = int(arg)
    elif opt in ('-d', '--determine'):
        determine = int(arg)
    elif opt in ('-s', '--special'):
        special = arg
    elif opt in ('-h', '--help'):
        print(help_info)
        sys.exit()

is_biomarker_control = 0
is_biomarker_scd = 0
if 'biomarkers' in special:
    # pull some special SCD data
    splt = special.split(",")
    second = splt[1]
    
    if second == 'control':
        print('Using control parameters for Biomarkers study (note that these are identical to default parameters)')
        is_biomarker_control = 1
    elif second == 'scd':
        print('Using generic SCD parameters for Biomarkers study (shortened TTT, reduced labeling efficiency, guess of hbs=0.5)')
        is_biomarker_scd = 1
    else:
        print(f'Attempting to pull data from REDCap for {second}')
        raise NotImplementedError
        keys_folder = os.path.join(three_up, '_secret_keys') # if you downloaded this repository from github you need to make this folder yourself and add the keys
        # it is not under version control because it is meant to store REDCap API keys
        biomarkers_token_loc = os.path.join(keys_folder, 'biomarkers.txt')
        
        api_url = 'https://redcap.vanderbilt.edu/api/'
        token = open(biomarkers_token_loc).read()
        
        project = redcap.Project(api_url, token)
        project_data_raw = project.export_records()
        project_data = pd.DataFrame(project_data_raw)
        
        # as of right now I actually don't know the database structure, so we can't extract data just yet
        
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
native_folder = data_folder.absolute()

standard_folder = os.path.join(native_folder, 'standard')

if not os.path.exists(standard_folder):
    os.mkdir(standard_folder)


# finds files matching the patterns specified in your patterns.xlsx
if seek:
    print('\nRunning SEEK')
    loaded_patterns = pd.read_excel(pattern_file, dtype=str)
    match_cols = [j for j in loaded_patterns.columns if 'match' in j.lower()]
    exclude_cols = [j for j in loaded_patterns.columns if 'exclude' in j.lower()]
    
    add_cols = ['File found', 'n found']
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
        
        new_file_loc = os.path.join(standard_folder, new_bn)
        vlabels_loc = os.path.join(standard_folder, vlabels)
        
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
            # for anyone seeking to implement DICOM compatibility,
            # pydicom is a good utility to start with for conversion purposes 
            
                    
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
    t1_location = os.path.join(standard_folder, f'{t1_name}.nii.gz')
    t1_location_betted = os.path.join(standard_folder, f'{t1_name}_BET.nii.gz')
    fast_location_betted = os.path.join(standard_folder, f'{t1_name}_BET_FAST')
    t1_wm_location = os.path.join(standard_folder, f'{t1_name}_BET_FAST_pve_2.nii.gz')
    
    # first thing: BET the T1
    print('Primary: extracting brain from T1')
    bet(t1_location, t1_location_betted, f=f_val)
    print('Primary: FASTing T1')
    fast(t1_location_betted, out=fast_location_betted, n_classes=3)
    
    t1space_folder = os.path.join(standard_folder, 't1space')
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
        
        scan_location = os.path.join(standard_folder, f'{scan_type}.nii.gz')
        omat_location = os.path.join(standard_folder, f'{scan_type}_to_t1space.omat')
        
        scan_location_in_t1space = os.path.join(t1space_folder, f'{scan_type}_t1space.nii.gz')
        
        if reg_method=='0': # file is already in T1 space
            continue
        elif reg_method=='1': # register using FLIRT
            scan_location_betted = os.path.join(standard_folder, f'{scan_type}_BET.nii.gz')
            bet(scan_location, scan_location_betted)
            flirt(src=scan_location_betted, ref=t1_location_betted, omat=omat_location)
        elif reg_method=='2': # register using epi_reg
            base_out = omat_location.replace('.omat', '')
            mat_location = f'{base_out}.mat'
            epi_reg(epi=scan_location, t1=t1_location, t1brain=t1_location_betted, out=base_out, wmseg=t1_wm_location)
            os.rename(mat_location, omat_location) # for whatever reason epi_reg makes a .mat file, not .omat. their contents are equivalent
            
            to_remove = ['_fast_wmedge.nii.gz', '_fast_wmseg.nii.gz', '_init.mat', '.nii.gz']
            for tr in to_remove: # some intermediate files we want to clean up
                file_to_remove = os.path.join(standard_folder, f'{base_out}{tr}')
                os.remove(file_to_remove)
        elif reg_method=='-1': # file should not be converted to T1 space
            continue
        else: 
            # if parameters are correctly specified, then this branch
            # handles cases where reg_method is the name of a scan
            # in this case, we use that scan's omat instead of generating a new one
            reference_omat = os.path.join(standard_folder, f'{reg_method}_to_t1space.omat')
            shutil.copyfile(reference_omat, omat_location)
            
        
        flirt(scan_location, ref=t1_location, out=scan_location_in_t1space,
              applyxfm=1, init=omat_location, interp='trilinear')
        
        
        # qc
        plt.style.use('dark_background')
        qc_location = os.path.join(t1space_folder, f'{scan_type}_to_t1_qualitycontrol.png')
        
        t1_loaded = nib.load(t1_location_betted)
        t1_data = t1_loaded.get_fdata()
        brain_shape = t1_data.shape
        
        target_loaded = nib.load(scan_location_in_t1space)
        target_data = target_loaded.get_fdata()
        target_shape = target_data.shape
        
        n_extradims = len(target_shape) - len(brain_shape)
        
        allslice = [slice(None)] * 3
        slice_axes = [0,1,2]
        
        fig, axs = plt.subplots(nrows=3, ncols=3, figsize=(8,12))
        for i, (slax, axrow) in enumerate(zip(slice_axes, axs)):
            slc = allslice.copy()
            slc[slax] = int(brain_shape[slax]/2)+10
            
            brain_sliced = t1_data[slc]
            axrow[0].imshow(brain_sliced, cmap='gist_gray')
            axrow[1].imshow(brain_sliced, cmap='gist_gray')
            
            target_slc = slc.copy()
            for i in range(n_extradims):
                target_slc.append(1) # for each extra dimension we need to make sure to just take the 0th element
            target_sliced = target_data[target_slc]
            axrow[1].imshow(target_sliced, cmap='inferno', alpha=0.5)
            axrow[2].imshow(target_sliced, cmap='inferno')
            
            if i==0:
                axrow[0].set_title('BETted T1')
                axrow[1].set_title('Overlap')
                axrow[2].set_title(f'Registered {scan_type}')
                
                
        for ax in np.ravel(axs):
            ax.axis('off')
        
        plt.tight_layout()
        fig.savefig(qc_location, dpi=400)

    
if determine:
    print('\nRunning DETERMINE')
    espec_template = os.path.join(three_up, 'bin', 'execution_specs.xlsx')
    target_espec = os.path.join(standard_folder, 'execution_specs.xlsx')
    xls = pd.ExcelFile(espec_template)
    
    loaded_patterns = pd.read_excel(pattern_file)
    loaded_patterns = loaded_patterns.set_index('Scan name')

    # need custom code for each sheet to handle the processing input
    
    can_run = []
    cant_run = []
    guessed_vals = []
    updated_dfs = []
    # vols
    for sheet_name in xls.sheet_names:
        loaded_template = pd.read_excel(xls,
                                        sheet_name=sheet_name,
                                        dtype=str)
        loaded_template = loaded_template.set_index('arg')
        if sheet_name=='vols':
            if loaded_patterns.loc['t1']['n found']:
                loaded_template.at['t1', 'value'] = os.path.join(standard_folder, 't1.nii.gz')
                
        if sheet_name=='asl':
            if loaded_patterns.loc['asl']['n found']:
                loaded_template.at['asl', 'value'] = os.path.join(standard_folder, 'asl.nii.gz')
                
            if loaded_patterns.loc['aslm0']['n found']:
                loaded_template.at['aslm0', 'value'] = os.path.join(standard_folder, 'aslm0.nii.gz')
                
            hct = 0.42
            labeff = 0.91
            ttt = 1.4
            
            if is_biomarker_scd:
                labeff = 0.72
                ttt = 1.29
            
            loaded_template.at['hct', 'value'] = hct
            loaded_template.at['labeff', 'value'] = labeff
            loaded_template.at['ttt', 'value'] = ttt
            
            loaded_template.at['hct', 'guessed'] = 1
            loaded_template.at['labeff', 'guessed'] = 1
            loaded_template.at['ttt', 'guessed'] = 1
            
            filename_extractors = ['ld', 'pld']
            standard_vals = [1650, 1525]
            if loaded_patterns.loc['asl']['n found']:
                for fne, s_val in zip(filename_extractors,standard_vals):
                    fn = loaded_patterns.loc['asl']['File found'].lower()
                    split = fn.split('_')
                    has_fne = [i.replace(fne, '') for i in split if fne in i]
                    
                    loaded_template.at[fne, 'value'] = s_val
                    loaded_template.at[fne, 'guessed'] = 1
                    for hf in has_fne:
                        try:
                            val = float(hf)
                            loaded_template.at[fne, 'value'] = val
                            loaded_template.at[fne, 'guessed'] = np.nan
                            break
                        except ValueError: # if not a perfect match, then hf will have characters in it
                            # we will still keep looking to see if we can find the right one
                            pass
            
        if sheet_name=='trust':
            if loaded_patterns.loc['trust']['n found']:
                loaded_template.at['trust', 'value'] = os.path.join(standard_folder, 'trust.nii.gz')
                
            hct = 0.42
            ya = 0.98 # aka artox
            hbs = 0
            
            if is_biomarker_scd:
                hbs = 0.5
            
            loaded_template.at['hct', 'value'] = hct
            loaded_template.at['ya', 'value'] = ya
            loaded_template.at['hbs', 'value'] = hbs
            
            loaded_template.at['hct', 'guessed'] = 1
            loaded_template.at['ya', 'guessed'] = 1
            loaded_template.at['hbs', 'guessed'] = 1
        
        if sheet_name=='ase':
            if loaded_patterns.loc['ase']['n found']:
                loaded_template.at['ase', 'value'] = os.path.join(standard_folder, 'ase.nii.gz')
                
            hct = 0.42
            artox = 0.98 # aka ya
            
            loaded_template.at['hct', 'value'] = hct
            loaded_template.at['artox', 'value'] = artox
            
            loaded_template.at['hct', 'guessed'] = 1
            loaded_template.at['artox', 'guessed'] = 1
        
        if sheet_name=='bold':
            # not implemented
            updated_dfs.append(None)
            continue
        
        
        
        
        runner = 1
        for i,row in loaded_template.iterrows():
            if row['guessed']==1:
                guessed_vals.append(f'{sheet_name} : {i} = {row["value"]}')
                
            if float(row['required'])==1 and row['value'] in [np.nan, None, '']:
                # we don't have everything we need to run this
                #print(f'You do not have {row.index[0]} (it is {row["value"]})')
                runner = 0
                
        loaded_template.at['DORUN', 'value'] = runner
        updated_dfs.append(loaded_template)
        
        if runner:
            can_run.append(sheet_name)
        else:
            cant_run.append(sheet_name)
            
    
    print(f'\nYou can run the following processing sequences: {can_run}')
    print(f'You can NOT run the following processing sequences: {cant_run}')
    
    if guessed_vals:
        guessed = "\n\t".join(guessed_vals)
        print(f'\nAlso, I could not find the following values, and had to guess them:\n\t{guessed}\n')
        
        
    writer = pd.ExcelWriter(target_espec, engine='xlsxwriter')
    
    for sheet_name, df in zip(xls.sheet_names, updated_dfs):
        if df is None:
            continue
        df.to_excel(writer, sheet_name=sheet_name)
        
    writer.save()
        
        
        
        
        
        
        
        
        
        
    