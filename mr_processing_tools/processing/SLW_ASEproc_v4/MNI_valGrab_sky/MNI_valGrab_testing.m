% addpath(genpath('/Users/manusdonahue/Documents/Sky/MNI_valGrab/')); % this should be the folder that this file is in


%% testing
inFile = '/Users/spencerwaddle/Documents/SLWtools/ASEcheck/MNI_valGrab_sky/examples/rOEF_Rdelete.nii.gz';
roidir = '/Users/spencerwaddle/Documents/SLWtools/ASEcheck/MNI_valGrab_sky/MNI_valGrab_masks';
fxnHandles = {@mean};
csvname = '/Users/spencerwaddle/Documents/SLWtools/ASEcheck/MNI_valGrab_sky/rOEF_valgrab_results.csv';

[result, strArray] = MNI_valGrab_v2(inFile, roidir, fxnHandles, csvname );


% f1 = '/Users/manusdonahue/Documents/Sky/MNI_valGrab/examples/rOEF.nii.gz';
% f2 = '/Users/manusdonahue/Documents/Sky/MNI_valGrab/examples/Donahue_241665_08_01_17.25.27_\(WIP_ms_ASE_tauShuffle0step0.5to22.5_randomTau_TE64_107_RFoff\)_rvCBV_coreg.nii.gz';
% f3 = '/Users/manusdonahue/Documents/Sky/MNI_valGrab/MNI_valGrab_masks/temp_store/MNI_2mm_temporal_R.nii.gz';
% f4 = '/Users/manusdonahue/Desktop/Projects/SCD/Data/SCD_P011_01/Processed/SCD_P011_01_3D_brain_pve_1_T2.nii.gz';
% f5 = '/Users/manusdonahue/Desktop/Projects/SCD/Data/SCD_P011_01/Processed/SCD_P011_01_CBF.nii.gz';
% f6 = '/Users/manusdonahue/Documents/Sky/MNI_valGrab/MNI_valGrab_masks/MNI_2mm_temporal_R.nii.gz';