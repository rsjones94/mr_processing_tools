#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Sky Jones
"""

import os
import glob


def find_image(base_folder, match, exclude=None, recursive=False):
    '''
    GLOBs through a folder to find the path to an image matching the pattern
    
    The paths making up the final set of candidate images are the union of all
    the match patterns, followed by the subtraction of each exclude pattern
    
    

    Parameters
    ----------
    base_folder : str
        Folder to begin search in.
    match : list of st
        The patterns that the filename should satisfy. If there is only one match pattern
        it should still be passed as a list.
    exclude : list of str, optional
        A list of patterns that the filename should not satisfy. If there is only one
        exclude pattern it should still be passed as a list. If None, then
        no exclude patterns will be used. The default is None.
    recursive : bool, optional
        If True, searches all subfolders of base_folder. If false, only looks in
        base_folder. The default is False.

    Returns
    -------
    final : str
        The path to the matched image. If more than one match is found, the first is returned.
        If no match is found, None is returned
    n_potentials : int
        The number of potential matches found.
    potentials : list of str
        The list of all potential matches.

    '''
    
    
    in_between = ''
    if recursive:
        in_between = '**'
        
    all_matches = []
    for mat in match:
        incl_globber = os.path.join(base_folder, in_between, mat)
        matches = glob.glob(incl_globber, recursive=recursive)
        match_set = set(matches)
        potentials = match_set
        all_matches.append(potentials)
    
    potentials = set.intersection(*all_matches)
    
    for excl in exclude:
        excl_globber = os.path.join(base_folder, in_between, excl)
        excludes = glob.glob(excl_globber, recursive=recursive)
        exclude_set = set(excludes)
        potentials = potentials - exclude_set
        
    potentials = list(potentials)
    # PAR/RECs will count as double, so filter out RECs
    potentials = [i for i in potentials if '.rec' not in i.lower()]
    
    try:
        final = potentials[0]
    except IndexError:
        final = None
    
    n_potentials = len(potentials)
    
    return final, n_potentials, potentials