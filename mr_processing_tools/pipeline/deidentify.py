#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Sky Jones
"""

help_info = """
This script replaces strings in both the body and names of all files in a folder.
Be careful to not pass strings the represent non-patient data, as you may compromise the integrity of the file by replacing them
However, this script will create a /bak folder with bakups up the original files. Note that the backups will obviously
not be deidentified, so they should be deleted or moved before sharing.

input:
    -f / --folder : the folder to deidentify
    -t / --target : the target string to replace
        Default is AUTO, which will attempt to find the name as delimited in the filename by a double carat (^^)
        Otherwise, pass a string or you can pass multiple strings to replace by delimiting them with commas, e.g., str1,str2,str3
    -r / --replacement : the string to replace the target with
        default is nothing, i.e., ''
    -h / --help : brings up this helpful information. does not take an argument
"""


        
        
        
        
import os
import sys
import getopt
import shutil
import subprocess


####

inp = sys.argv
bash_input = inp[1:]
options, remainder = getopt.getopt(bash_input, "f:t:r:h", ["folder=", 'target=', 'replacement=', 'help'])

target = 'AUTO'
for opt, arg in options:
    if opt in ('-f', '--folder'):
        master_folder = arg
    if opt in ('-t', '--target'):
        target = arg
    if opt in ('-r', '--replacement'):
        replacement = arg
    elif opt in ('-h', '--help'):
        print(help_info)
        sys.exit()

onlyfiles = [f for f in os.listdir(master_folder) if os.path.isfile(os.path.join(master_folder, f))]
if target=='AUTO':
    delimit = '^^'
    has_dcarat = [i for i in onlyfiles if delimit in i]
    if len(has_dcarat) < 1:
        raise Exception(f'No files in this directory appear to have a {delimit} delimiter, so I cannot automatically find the name to replace.')
    sample = has_dcarat[0]
    name = sample.split(delimit)[0]
    target = f'{name}{delimit}'
    
    do_proceed = False
    valid_yes = ['y', 'yes']
    valid_no = ['n',  'no']
    while not do_proceed:
        feedback = input(f'\nI found {target} as the auto-selected name. Okay to replace this?\n\t')
        if feedback.lower() in valid_yes:
            do_proceed = True
        elif feedback.lower() in valid_no:
            raise Exception('User terminated processing')
        else:
            print('just yes or no, please')
            
    target = [target]
            
else:
    if ',' not in target:
        target = [target]
    else:
        target = target.split(',')
        
        
#make a bak folder and move the files there
bak_folder = os.path.join(master_folder, 'bak')
if not os.path.exists(bak_folder):
    os.mkdir(bak_folder)

for f in onlyfiles:
    if f == '.DS_Store':
        continue
    fromfile = os.path.join(master_folder, f)
    tofile = os.path.join(bak_folder, f)
    shutil.copy(fromfile, tofile)
    
#change dir command
os.chdir(master_folder)

os.environ['LC_CTYPE'] = "C"
os.environ['LANG'] = "C"



for t in target:
    
    
    #replace in file command (and .bak)
    
    rep_body = ['find',
                 '.',
                 '-type',
                 'f',
                 '-maxdepth',
                 '1',
                 '-name',
                 "'*.*'",
                 '-exec',
                 'sed',
                 '-i.bak',
                 '-e',
                 f'"s/{t}/{replacement}/g"',
                 '{}',
                 '+']
    rep_body = ' '.join(rep_body)
    
    bodystat = subprocess.call(rep_body,shell=True,executable='/bin/zsh')
    #print(rep_body)
    #print(' '.join(rep_body))
    #os.system(rep_body)
    
    #replace in filename command
    #rep_filename = f'for file in *.&; do mv -v "$file" "${{file/{t}/{replacement}}}"; done'
    #filenamestat = subprocess.call(rep_filename,shell=True,executable='/bin/zsh')
    
    # get rid of the .bak files
    onlyfilesnew = [f for f in os.listdir(master_folder) if os.path.isfile(os.path.join(master_folder, f))]
    bakfiles = [i for i in onlyfilesnew if '.bak' in i]
    
    
    for bf in bakfiles:
        thepath = os.path.join(master_folder, bf)
        os.remove(thepath)
    
    
    onlyfilescurrent = [f for f in os.listdir(master_folder) if os.path.isfile(os.path.join(master_folder, f))]
    replacefiles = [i for i in onlyfilescurrent if t in i]
    for rf in replacefiles:
        repname = rf.replace(t,replacement)
        
        
        fromfile = os.path.join(master_folder, rf)
        tofile = os.path.join(master_folder, repname)
        shutil.move(fromfile, tofile)
                    
    
