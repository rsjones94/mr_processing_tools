#!/usr/bin/env python3
# -*- coding: utf-8 -*-

help_info = """
This script uses takes a TRUST image and extracts the estimated blood T2 from
the superior sagittal sinus. It then applies a number of models to quantify
venous oxygenation.

The script will create:
    1) a folder in the TRUST image directory called 'trust'
    2) an excel spreadsheet in this folder called 't2_vals.xlsx' that contains the T2 (s) for each acquisition
    3) an excel spreadhseet called 'fits.xlsx' that contains exponential fits for each acquisition with confidence intervals
    4) an excel spreadsheet called 'signals.xlsx' that gives the signal at each ete with confidence intervals
    5) an excel spreadsheet called 'models.xlsx' that contains the venous oxygenation using mean T2 applied to
        various quantification models as well as the OEF

required input:
    -t / --trust : the path to the trust, which is a .nii.gz file
        this file should have an associated vlabels file in the same folder
        (i.e., if your image is trust.nii.gz, there should be a trust_vlabels.xlsx in the same folder)
    -c / --hct : the hematocrit as a float between 0 and 1
essentially required input:
    -y / --ya : the arterial oxygen saturation as a float between 0 and 1
        assumed to be 0.98 if not specified
    -s / --hbs : the hemoglobin s fraction as a float between 0 and 1
        assumed to be 1 if not specified
optional input:
    -e / --ete : the echo times in ms as comma-separated values
        default is 0,40,80,160
    -a / --auto : whether to use the automatic superior sagittal sinus selection, or to select it yourself (0 or 1)
        default is 1
    -i / --invert : comma separated 0 or 1s to indicate whether to invert the intensities for each acquisition
        default is None, which will not invert anything
        note that the intensities are automatically inverted (or not) already to make the
        SSS bright against a dark background If this automatic inversion does not work, then
        manually specify which acquisitions needs to be inverted. For example, if you have
        two acquisitions, and want to invert the second, then pass --invert 0,1
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
from processing.trust_funcs import t2_from_trust
import helpers.oxsat_models as osm
####

inp = sys.argv
bash_input = inp[1:]
options, remainder = getopt.getopt(bash_input, "t:c:y:s:e:a:i:h",
                                   ["trust=", 'hct=', 'ya=', 'hbs=', 'etes=', 'auto=', 'invert=', 'help'])


hct = None
hbs = 1
ya = 0.98
auto = True
inverter = None
etes = [0,40,80,160]
'''

trust_file = '/Users/skyjones/Documents/legacy_samples/SCD_K045_02/TRUST.nii.gz'
hct = 0.5
hbs = 0.7'''


for opt, arg in options:
    if opt in ('-t', '--trust'):
        trust_file = arg
    if opt in ('-c', '--hct'):
        hct = float(arg)
    if opt in ('-y', '--ya'):
        ya = float(arg)
    if opt in ('-s', '--hbs'):
        hbs = float(arg)
    if opt in ('-a', '--etes'):
        ete_string = arg
        ete_split = ete_string.split(',')
        etes = [float(i) for i in ete_split]
    if opt in ('-a', '--auto'):
        auto = bool(int(arg))
    if opt in ('-i', '--invert'):
        invert_string = arg
        if invert_string.lower() == 'none' or invert_string=='0':
            inverter = None
        else:
            invert_split = invert_string.split(',')
            inverter = [int(i) for i in invert_split]
    elif opt in ('-h', '--help'):
        print(help_info)
        sys.exit()

assert os.path.isfile(trust_file)


blood_t1 = calculate_blood_t1(hct)

trust_base = os.path.basename(os.path.normpath(trust_file))   
trust_dir = os.path.dirname(os.path.normpath(trust_file))   
trust_core = trust_base.split('.')[0]
trust_vlabels_base = f'{trust_core}_vlabels.xlsx'
trust_vlabels_loc = os.path.join(trust_dir, trust_vlabels_base)


trust_loc = os.path.join(trust_dir, 'trust')
if not os.path.exists(trust_loc):
    os.mkdir(trust_loc)

 
trust_loaded = nib.load(trust_file)
trust_data = trust_loaded.get_fdata()

vlabels_loaded = pd.read_excel(trust_vlabels_loc, index_col=0)
ordering = ['dynamic scan number', 'label type']
vlabels_ordered = vlabels_loaded[ordering]

trust_unpacked = unpack_dims(trust_data, vlabels_ordered)
trust_clean = trust_unpacked[:,:,0,:,:] # dim 3 (z) should be a singleton


t2s, signals, fits = t2_from_trust(trust_clean,
                                   blood_t1,
                                   ete=etes, 
                                   auto=auto,
                                   override_inversions=inverter)

acqs = [i for i,val in enumerate(t2s)]

t2_file = os.path.join(trust_loc, 't2_vals.xlsx')
fit_file = os.path.join(trust_loc, 'fits.xlsx')
signals_file = os.path.join(trust_loc, 'signals.xlsx')
model_file = os.path.join(trust_loc, 'models.xlsx')

# write the t2 file
mean_t2 = np.mean(t2s)
t2_df = pd.DataFrame()


#t2_write = mean_t2.copy()
#t2_write.append(mean_t2)
t2_write = t2s.copy()
acq_write = acqs.copy()
#acq_write.append('mean')

t2_df['acq'] = acq_write
t2_df['t2_s'] = t2_write
t2_df.to_excel(t2_file)

# write the fit file
writer = pd.ExcelWriter(fit_file, engine='xlsxwriter')
for sdf, acq in zip(fits, acqs):
    sdf.to_excel(writer, sheet_name=f'acq_{acq}')
writer.save()

# write the uncertainty file
writer = pd.ExcelWriter(signals_file, engine='xlsxwriter')
for sdf, acq in zip(signals, acqs):
    sdf.to_excel(writer, sheet_name=f'acq_{acq}')
writer.save()

# write the model file
yvs = []
oefs = []

oxsat_funcs = [osm.bovine_oxsat, osm.aa_oxsat, osm.ss_oxsat, osm.f_oxsat, osm.mixture_oxsat]
oxsat_func_names = ['bovine', 'wood_aa', 'wood_ss', 'fetal', 'mixture']
for fun, name in zip(oxsat_funcs, oxsat_func_names):
    yv = fun(mean_t2, hct, hbs=hbs)
    oef = (ya - yv) / ya
    
    yvs.append(yv)
    oefs.append(oef)
    
model_df = pd.DataFrame()
model_df['model'] = oxsat_func_names
model_df['yv'] = yvs
model_df['ya'] = ya # note that ya is a scalar and so this column will be homogenous
model_df['oef'] = oefs

model_df.to_excel(model_file)


    
    
    
    
    
    
    
    




