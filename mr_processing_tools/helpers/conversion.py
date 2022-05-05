#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Sky Jones
"""

import nibabel as nib
from fsl.wrappers import fslreorient2std, LOAD
import pandas as pd
import numpy as np

'''

Standard usage:
    
parrec_path = path/to/parrec
nifti_path = where/you/want/to/write/nifti.nii.gz
vlabel_path = where/you/want/to/write/vlabel.xlsx

loaded_parrec = nib.load(parrec_path, scaling='fp')
nifti, vlabels = parrec_to_nifti(loaded_parrec, nifti_path, vlabel_path, reorient=True)
nifti_data = renifti.get_fdata()
unpacked_data = unpack_dims(nifti_data, vlabels)


'''


def parrec_to_nifti(in_parrec, out_nifti=None, out_vlabels=None, reorient=True):
    '''
    Takes a PARREC object and converts it to a NiFTi image
    
    A note on scaling:
        Make sure your PARREC was loaded with the correct scaling option. You probably
        want to load with scaling='fp'
            'fp' is floating point value: FP = DV / (RS * SS)
            'dv' is value on console: DV = PV * RS + RI
        where:
            RS: rescale slope; RI: rescale intercept; SS: scale slope
    

    Parameters
    ----------
    in_parrec : Nibabel PARREC image
        The image to be converted
    out_nifti : str or None, optional
        The location to write the NiFTi to.
        If None, the file is not written (but is still retained in-memory).
        The default is None.
    out_vlabels : str or None, optional
        Location to write a companion Excel file specifying the volume labels
        of the 4th image dimension. If None, then the file will not be written
        but the labels are still retained as a pandas DataFrame. Note that purely
        spatial scans (e.g., T1, FLAIR) do not have any data in the 4th dimension
        and thus the vlabels are empty.
        The default is None.
    reorient : bool, optional
        Whether to reorient the NiFTI to standard anatomical alignment. The default is True.
        This may change anatomical labels and the order/direction of the spatial dimensions of
        the data matrix itself, but does not change the voxel locations in scanner space.

    Returns
    -------
    to_write
        the converted NiFTi image (may or may not actually be written to a file)
        as a Nibabel Nifti1image
    vlabels
        the volume labels as a pandas DataFrame (also may or may not be written to a file)

    '''
    
    
    vlabels_dict = in_parrec.header.get_volume_labels()
    vlabels = pd.DataFrame(vlabels_dict)
    
    orig_data = in_parrec.dataobj
    orig_affine = in_parrec.affine
    orig_header = in_parrec.header
    
    nifti = nib.Nifti1Image(orig_data, orig_affine, header=orig_header)
    nifti.set_data_dtype('<f4')
    
    if reorient:
        to_write = fslreorient2std(nifti, output=LOAD)['output']
    else:
        to_write = nifti
    
    
    if out_nifti:
        to_write.to_filename(out_nifti)
    if out_vlabels:
        vlabels.to_excel(out_vlabels)
    
    return to_write, vlabels
    

def unpack_dims(data, vlabels):
    """
    Unpacks an interleaved 4th dimension in an imaging data array

    Parameters
    ----------
    data : np array
        a numpy array of data. Should have 3 spatial dimensions followed by
        one nonspatial dimension of interleaved data
    vlabels : pandas DataFrame
        a dataframe indicating the label type for each slice in the 4th dimension.
        Each column specifies a label type, and each row gives the labeling combination for
        each index in the 4th dimension. It is assumed row 0 corresponds to index 0 in the 4th
        dimension, and so on. Additionally, each column should be gapless. In other words, there
        should be at least one entry for each integer between the min and max of a column. Otherwise
        you will get blank dimensions in the unpacked data.
        
        The vlabels are returned as part of parrec_to_nifti. They can also be found as
        an ordered dict with parrecobject.header.get_volume_labels() (but you must
        convert to an DataFrame). A NiFTi using parrec_to_nifti may also have the labels
        saved as an Excel file
        
        Note that you change the order of the column of vlabels,
        and this will change the order of the returned dimensions of the output

    Returns
    -------
    new_data : np array
        The data with the 4th dimension unpacked into n additional dimensions, where
        n is the number of columns in vlabels. For example, if relabels has three columns,
        then the 4th dimension of the original data will become the 4th, 5th and 6th dimensions
        in new_data. The 1st column in vlabels will be stored as the 4th dimension, the 2nd
        column in vlabels will be stored as the 5th dimension, and so on

    """
    
    adj_labels = vlabels.copy()
    for col in adj_labels.columns:
        # we need to make the labels zero indexed
        adj_labels[col] = adj_labels[col] - np.min(adj_labels[col])
        
    spatial_dim_maxes = data.shape[:3]
    extra_dim_maxes = [i for i in adj_labels.max()+1]
    
    dim_shape = []
    dim_shape.extend(spatial_dim_maxes)
    dim_shape.extend(extra_dim_maxes)
    
    new_data = np.zeros(dim_shape)
    for i, row in adj_labels.iterrows():
        sli = data[:,:,:,i]
        # construct an indexing tuple
        itup = [...]
        add_labels = list(row)
        itup.extend(add_labels)
        itup = tuple(itup)
        new_data[itup] = sli
        
    return new_data
    

        