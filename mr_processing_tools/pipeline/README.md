
# Overview

There are two primary pipeline tools: prepare.py and analyze.py. A courtesy deidentification tool (deidentify.py) is also provided. It merely removes the participants name from filenames and file bodies, and so is not HIPAA compliant. You should still run it prior to preparation and analysis. Additional SCD Biomarkers study tools are provided in the /biomarkers subfolder; see the README for more info. All scripts in this folder and subfolders have help strings that can be accessed by calling the function and passing -h or --help.

## prepare.py
Takes a folder of scans and prepares them for analysis by converting them to NiFTi, normalizing their indices, coregistering them, and making an educated guess as to what processing sequences can be run. This script will output a .xlsx file specifying the processing sequences to run and any processing parameters necessary. Some processing parameters may guessed but when possible they will be pulled from available metadata. It is up to you to inspect this file (/standard/execution_specs.xlsx) and make any necessary changes before proceeding.
>Note that if you are processing a Biomarkers patient, then you can run the biomarkers_getredcap.py script, then pass -s usepull to automatically use this data in the execution specification file.

## analyze.py

Uses the execution_specs.xlsx file to call scripts from the /analyze folder using the parameters specified in each tab. Each tab *tabname* calls /analyze/analyze_*tabname*.py. Currently there are scripts in place to analyze ASE, ASL, TRUST and volumetrics via FreeSurfer. Each script makes a folder of the same name with the processed data.
>If you wish to add additional processing functionality, write a script (callable from the terminal) with the name analyze_*funcname*.py and place it in the analyze folder. Then add a tab *funcname* to the stock execution_specs.xlsx file in the top level /bin folder. In order to automatically fill out this tab, you will have to modify prepare.py to parse and add the required information as well.