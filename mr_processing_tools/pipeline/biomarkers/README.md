
# Biomarkers specific tools

These tools enable further automation of the Biomarkers processing pipeline.

## biomarkers_getredcap.py
Pulls relevant lab values, etc. from REDCap and writes it to a .xlsx file. Call this prior to prepare.py. This file should go into the same folder as your raw scans, and when you call prepare.py, pass *-s usepull* to incorporate this data into the execution specification file.

## biomarkers_minedata.py
Once you have run prepare.py and analyze.py, call this grab all of the information that goes into the REDCap database. This data is written to a file called mineddata.csv.

## biomarkers_pushdata.py
This takes mineddata.csv and uploads it to REDCap for you.