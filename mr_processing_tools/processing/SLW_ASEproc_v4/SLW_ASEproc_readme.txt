The main script the user will interact with is ASEmaster_Bioclinica.m. I did my best to make this the only script that you really need
to interact with to change ASE processing details. TR, TE1, TE2, tau shifting values, hematocrit, resolution, all of these are 
specified in this master file. pCASL and TRUST processing are also in this file, and you can comment them out if you do not want to process
those on a particular patient.

The majority of ASE processing is done in ASEcheck_v7_fn.m. Some other parameters are defined here, though these parameters shouldn't change.
For example, dChi0 is the measured susceptibility difference between oxygenated and deoxygenated blood. Smoothing and fitting all occur in this file

ASEmodelScripts are all of the files used in ASE fitting.

slw_resources contains all the third-party resources the code uses. I included this to try to make it easier for someone
to pick up the code without having to hunt these other resources down.

TRUST_standalone is the TRUST processing directory. Everything needed for processing TRUST should be in that one directory.

To process ASE data
1. Change standaloneDir and datadir1 to point to the SLW_ASEproc_v4 and data folders on your computer, respectively
2. Input the name of each patient directory as a structure in patstruct
3. Change any ASE parameters you have altered in acquisition. These should all be in the first 10-40 lines of code of ASEmaster_washu.m

Additional notes: 
- Coregistration uses fsl tools. The easiest way to access these if you don’t have them is to install them from
https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/ and open Matlab from terminal. If you encounter this problem, you will get an 
Error that says “$FSLOUTDIR not set”

- Please note that the ASE processing script searches for files based on their name. For example, Asymmetric Spin Echo
Has the keyword ASE, and RFoff or RFon (without a refocusing pulse, and with a refocusing pulse between the first
And second ASE echoes, respectively). TRUST is identified with TRUST_VEIN. T1 is identified with 3D_T1.
pCASL with pCASL and the associated M0 file is pCASL_M0.

- The approximate time shown during processing is for a single ASE file. Processing multiple files can take a day or more.

