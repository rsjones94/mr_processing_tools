# MR processing tools for the Donahue lab

This repository is a collection of Python and MATLAB scripts for processing MRI data, as well as complete pipelines for the automated analysis of standard acquisitions (e.g., adult sickle cell, moyamoya BOLD, etc.)

Do **not** put participant or patient data in the folder, even if they are deidentified. This repository is publicly visible and doing so may expose PHI.


After downloading this repository, you will need to:

	1) Make a text file called matlob_loc.txt in the /bin folder where the content is your MATLAB app location (e.g., /Applications/MATLAB_R2020b.app/) (otherwise ASE processing will not work). This file is not under version control

	2) Make a folder called _secret_keys and place a file called 'biomarkers.txt' in it that contains the API key for the SCD Biomarkers database. This folder is automatically excluded from version control by the .gitignore. (if you do not do this then you will not be able to call '--special biomarkers,[mrid]' in process.py) DO NOT PLACE THE API KEY IN ANY OTHER FOLDER, OR YOU WILL VERSION CONTROL IT AND EXPOSE IT TO THE PUBLIC

	3) For ease of use, add /path/to/this/mr_processing_tools/mr_processing_tools/pipeline/ to your PATH so that you can call the pipeline functions (e.g., prepare.py) directly without having to specify the full path