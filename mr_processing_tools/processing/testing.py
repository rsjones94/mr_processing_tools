#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 15 13:49:07 2022

@author: skyjones
"""

import sys
sys.path.append('..')

import nibabel as nib
import scipy
import numpy as np
import matplotlib.pyplot as plt

import helpers.general as ge
from helpers.conversion import parrec_to_nifti, unpack_dims
import trust_funcs

in_parrec = '/Users/skyjones/Documents/repositories/legacy_processing/SCD/samples/SCD_K045_02/Acquired/SCD_K045_02_WIP_SOURCE_-_MJD_TRUST_VEIN_R57_FINAL_5_2.PAR'
out_nifti = '/Users/skyjones/Documents/repositories/legacy_processing/SCD/samples/SCD_K045_02/hold/sub/trust.nii.gz'
out_vlabels = '/Users/skyjones/Documents/repositories/legacy_processing/SCD/samples/SCD_K045_02/hold/sub/trust_vlabels.xlsx'

loaded_parrec = nib.load(in_parrec, scaling='fp')

renifti, relabels = parrec_to_nifti(loaded_parrec, out_nifti, out_vlabels)
redata = renifti.get_fdata()
unpacked = unpack_dims(redata, relabels)

reformat = unpacked[:,:,0,:,:]
reformat2 = np.swapaxes(reformat,2,3)

result = trust_funcs.t2_from_trust(reformat2, 10, auto=False)


#trust_source = 

"""
'''
pcasl_source = '/Users/skyjones/Documents/repositories/legacy_processing/SCD/samples/SCD_K045_02/Acquired/SCD_K045_02_WIP_SOURCE_-_pCASL_PLD1650_LD1525_REL5FIXED_WithBS_6_2.PAR'
m0_source = '/Users/skyjones/Documents/repositories/legacy_processing/SCD/samples/SCD_K045_02/Acquired/SCD_K045_02_WIP_MJD_SCDPEDS_pCASL_M0_R57_PATCH_7_1.PAR'

pcasl_in = nib.load(pcasl_source)
v_pcasl, vs_pcasl, *_ = ge.reorganize_parrec(pcasl_in)

m0_in = nib.load(m0_source)
v_m0, vs_m0, *_ = ge.reorganize_parrec(m0_in)


asl_data = vs_pcasl[:, :, :, :, 0, 0, 0, :]
m0_data = vs_m0[:, :, :, 0, 0, 0, 0, 0]
'''

t1_parrec_source = '/Users/skyjones/Documents/repositories/legacy_processing/SCD/samples/SCD_TRANSP_P025_01_hold/SCD_TRANSP_P025_01_04_01_09.08.03_(WIP_cs_T1W_3D_TFE).PAR'
t1_nifti_source = '/Users/skyjones/Documents/repositories/legacy_processing/SCD/samples/SCD_TRANSP_P025_01_hold/hold_WIPcs_T1W_3D_TFEtransp4.2_19770703150928_4.nii.gz'

t1_template_source = '/Users/skyjones/Documents/repositories/legacy_processing/SCD/samples/SCD_K045_02/hold/SCD_K045_02_WIP_MJD_SCDPEDS_3DT1_4_1.PAR'
t1_template_source_nii = '/Users/skyjones/Documents/repositories/legacy_processing/SCD/samples/SCD_K045_02/hold/hold_WIPMJD_SCDPEDS_3DT1WOCONTRAST_19770703150928_4.nii.gz'


par_in = nib.load(t1_parrec_source, scaling='fp')
nii_in = nib.load(t1_nifti_source)
template_in = nib.load(t1_template_source)
template_in_nii = nib.load(t1_template_source_nii)
trust_in = nib.load(trust_source)

raw, par_data, *_ = ge.reorganize_parrec(par_in)
nii_data = nii_in.get_fdata()
_, template_data, *_ = ge.reorganize_parrec(template_in)
trust_raw, trust_data, ax, bx = ge.reorganize_parrec(trust_in)

idex = 100

template_slice = template_in.get_fdata()[:,idex,:]
template_slice_nii = template_in_nii.get_fdata()[:,idex,:]

fig, axs = plt.subplots(1,2)

axs[0].imshow(template_slice)
axs[1].imshow(template_slice_nii)

axs[0].set_title('Template PAR')
axs[1].set_title('Template NiFTi')

print('PAR')
print(par_in.affine.round())
print(nib.aff2axcodes(par_in.affine))
print()

print('NiFTi')
print(nii_in.affine.round())
print(nib.aff2axcodes(nii_in.affine))
print()

print('Template')
print(template_in.affine.round())
print(nib.aff2axcodes(template_in.affine))

print('Template NiFTi')
print(template_in_nii.affine.round())
print(nib.aff2axcodes(template_in_nii.affine))
"""