# MR processing tools for the Donahue lab

This repository is a collection of Python and MATLAB scripts for processing MRI data, as well as complete pipelines for the automated analysis of standard acquisitions (e.g., adult sickle cell, moyamoya BOLD, etc.)
For an overview on how to operate the pipeline, check the README in /mr_processing_tools/pipeline

Do **not** put participant or patient data in the folder, even if they are deidentified. This repository is publicly visible and doing so may expose PHI.


After downloading this repository, you will need to:

	- pip install pycap and fslpy (you also should have an Anaconda installation of Python)

	- install FreeSurfer and set your $SUBJECTS_DIR

	- Make a text file called matlab_location.txt in the /bin folder where the content is your MATLAB app location (e.g., /Applications/MATLAB_R2020b.app/) (otherwise ASE processing will not work). This file is not under version control

	- Make a folder called _secret_keys in the top level of this repository and place a file called 'biomarkers.txt' in it that contains the API key for the SCD Biomarkers database. This folder is automatically excluded from version control by the .gitignore. (if you do not do this then you will not be able to call '--special biomarkers,[mrid]' in process.py) DO NOT PLACE THE API KEY IN ANY OTHER FOLDER, OR YOU WILL VERSION CONTROL IT AND EXPOSE IT TO THE PUBLIC

	- For ease of use, add /path/to/this/mr_processing_tools/mr_processing_tools/pipeline/ to your PATH so that you can call the pipeline functions (e.g., prepare.py) directly without having to specify the full path. You may also want to add the other subfolders in this directory (e.g., /analyze, /biomarkers) to your PATH, depending on your intended usecase.