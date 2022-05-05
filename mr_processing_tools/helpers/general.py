#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Sky Jones
"""

import os

import numpy as np
import nibabel as nib
from scipy import signal
from numba import jit, float64, int64


def calculate_blood_t1(hct):
    # hct as a float between 0 and 1
    #%s - Jordan et al.
    blood_t1 = 1./(0.52*hct+0.38)
    
    return blood_t1


def fread(fid, nelements, dtype=np.int16):
     if dtype is np.str:
         dt = np.uint8  # WARNING: assuming 8-bit ASCII for np.str!
     else:
         dt = dtype

     data_array = np.fromfile(fid, dt, nelements)
     data_array.shape = (nelements, 1)

     return data_array
 
         

def reorganize_parrec(loaded_parrec):
    '''
    Takes a parrec object and attempts to reformat and scale its data in a standardized way.
    
    Only tested with par version 4.2 and on, as this is when the ASL type header
    was added (called 'label type', last column in v4.2). However, it *should* work
    with at least 4.0 on

    Parameters
    ----------
    loaded_parrec : nibabel PARRECImage
        A PARREC loaded with the nibabel.load() function.

    Returns
    -------
    v : numpy ndarray
        The reformatted, unscaled image data.
    v_scaled : numpy ndarray
        The reformatted, scaled image data.
    dim_names : list of str
        The names corresponding to each dimension of the reformatted images
    dim_sizes : list of int
        The voxel resolution corresponding to each dimension of the reformatted images.
        It's up to you to know the units
        
    Both v and v_scaled have the following dimensions:
        [rows, n_cols, n_slices, n_echoes, n_dynamics, n_types, n_phases, n_asltypes]
        
    Note that pre-v4.2 PAR files do not record asltypes. For such PARRECs,
        this dimension will be imputed as a singleton.
        
    Singleton dimensions are retained for clarity.
    

    '''
    
    gen_info = loaded_parrec.header.general_info
    im_defs = loaded_parrec.header.image_defs
    
    phases_name = 'cardiac phase number'
    phases = np.unique(im_defs[phases_name])
    n_phases = len(phases)
    
    slices_name = 'slice number'
    slices_step_name = 'slice thickness'
    slices = np.unique(im_defs[slices_name])
    n_slices = len(slices)
    size_slices = im_defs[slices_step_name][0]
    
    dynamics_name = 'dynamic scan number'
    dynamics_step_name = 'dyn_scan_begin_time'
    dynamics = np.unique(im_defs[dynamics_name])
    n_dynamics = len(dynamics)
    dyn_times = np.unique(im_defs[dynamics_step_name]) # this should be ordered, and the spacing should be even
    dyn_diffs = np.diff(dyn_times)
    size_dynamics = np.mean(dyn_diffs)
    
    echoes_name = 'echo number'
    echoes_step_name = 'echo_time'
    echoes = np.unique(im_defs[echoes_name])
    n_echoes = len(echoes)
    size_echoes = im_defs[echoes_step_name][0]
    
    types_name = 'scanning sequence'
    types = np.unique(im_defs[types_name])
    n_types = len(types)
    
    # note that the asltypes is only recorded for parrecs v4.2 and later
    asltypes_name = 'label type'
    pre_asltypes_flag = False
    try:
        asltypes = np.unique(im_defs[asltypes_name]) # types of ASL. For TRUST, this is the number of acquisitions. The nomenclature is confusing
        n_asltypes = len(asltypes)
    except KeyError: # need to test this branch
        print(f'{asltypes_name} not found in parsed PAR. Assuming only 1 label type')
        asltypes = [1]
        n_asltypes = len(asltypes)
        pre_asltypes_flag = True
        
    
    resolution_name = 'recon resolution'
    resolution_step_name = 'pixel spacing'
    res = im_defs[resolution_name][0] # we're assuming here that every frame has the same number of rows and columns
    res_step = im_defs[resolution_step_name][0] # we're assuming here that every frame has the same pixel size
    n_rows, n_cols = res
    size_rows, size_cols = res_step
    
    data_binary = loaded_parrec.dataobj.file_like.fobj
    db_name = data_binary.name
    fid = open(db_name)
    
    rs_col = 'rescale slope' # scaling slope
    ri_col = 'rescale intercept' # scaling intercept
    ss_col = 'scale slope' # philips-specific scaling factor
    scale_cols = [rs_col, ri_col, ss_col]
    
    dim_names = ['rows', 'cols', 'slices', 'dynamics', 'echoes', 'types', 'phases', 'asltypes']
    dim_sizes = [size_rows, size_cols, size_slices, size_dynamics, size_echoes, np.nan, np.nan, np.nan]
    
    v = np.zeros([n_rows, n_cols, n_slices, n_dynamics, n_echoes, n_types, n_phases, n_asltypes])
    v_scaled = np.zeros_like(v)
    
    step = np.product(res)
    readin = np.fromfile(fid, dtype=np.int16)
    for i,row in enumerate(im_defs):
        i_from = step*i
        i_to = i_from + step
        dat = readin[i_from:i_to]
        dat_res = dat.reshape(n_rows, n_cols)
        
        scale_factors = row[scale_cols]
        dat_scaled = scale_parrec_data(dat_res, scale_factors)
        
        for j, (matrix, holder) in enumerate(zip([dat_res, dat_scaled], [v, v_scaled])):
            
            # note that python is 0-indexed, but some of these attributes are 1-indexed (e.g., slices start at 1)
            # because of this we need to subtract 1 from the appropriate attributes
            # argwhere entries dont need to be altered, as they are inherently pinned to python indexing even if the original attribute is not
            
            if pre_asltypes_flag:
                final_dim = 1
            else:
                final_dim = np.argwhere(asltypes==row[asltypes_name])[0][0]
            
            holder[:,
                   :,
                   row[slices_name] - 1,
                   np.argwhere(dynamics==row[dynamics_name])[0][0],
                   np.argwhere(echoes==row[echoes_name])[0][0], 
                   np.argwhere(types==row[types_name])[0][0],
                   row[phases_name]-1,
                   final_dim
                   ] = matrix.T
                
    # in the original matlab code, the return matrices are permuted as such in MATLAB:
        # v = permute(v, [2 1 3:length(size(v))]);
        # This rotates pi/2 clockwise and flips lr
    
    permuter = [1,0]
    permuter.extend(np.arange(2,len(v.shape)))
    v = np.transpose(v, permuter)
    v_scaled = np.transpose(v_scaled, permuter)
    
    # note that in the original MATLAB script, the dimension indices for dynamics and echoes
    # are swapped from what is shown here. The data is returned in this order becuase
    # dynamics here generally correspond to the time dimension, and having time be the 4th
    # dimension is standard for most MR image formats
                
    return v, v_scaled, dim_names, dim_sizes
    
    

def get_parrec_scaling(image_info):
    '''
    Gets the three philips scaling factors from parsed image info. Note that this
    function checks that the scaling factors are the same for every frame,
    and will throw an error if they are not.

    Parameters
    ----------
    image_info : numpy ndarray
        Parsed image info. If your loaded PARREC is X, then the image info
        can be retrieved with X.header.image_defs

    Returns
    -------
    rs : float
        rescaling slope.
    ri : float
        rescaling intercept.
    ss : float
        scaling slope.

    '''
    #image_cols = np.array(image_info.dtype.names)
    
    # get scaling factors for each frame
    rs_col = 'rescale slope' # scaling slope
    ri_col = 'rescale intercept' # scaling intercept
    ss_col = 'scale slope' # philips-specific scaling factor
    
    scale_cols = [rs_col, ri_col, ss_col]
    scale_vals = np.asarray([image_info[i] for i in scale_cols]) # returns the scaling values for each frame
    
    # we make the assumption that the scaling is the same for every frame
    # but let's enforce that assumption
    first_row = scale_vals[:,0]
    matches = np.asarray([first_row==scale_vals[:,i] for i in range(scale_vals.shape[1])])
    
    assert matches.all()
    
    # now we can actually rescale
    rs = first_row[0]
    ri = first_row[1]
    ss = first_row[2]
    
    return rs, ri, ss

    
def scale_parrec_data(unscaled_data, scaling_factors):
    '''
    Unfortunately nibabel does not scale parrec data correctly as it does not consider the philips-proprietary scale slope.
    This function scales raw parrec data per Philips specs.
    IMPORTANT: this function assumes the scaling is the same for every frame.

    Parameters
    ----------
    unscaled_data : numpy ndarray
        The raw data to be scaled. If your loaded PARREC is X, then the
        unscaled data can be retrieved with X.dataobj.get_unscaled()
    scaling_factors : list of length 3
        The rescaling slope, the rescaling intercept, and the scaling slope as floats.
        Can retrieve these by calling get_parrec_scaling() on X.header.image_defs

    Returns
    -------
    scaled_data : numpy ndarray
        The properly scaled data.

    '''
    
    rs, ri, ss = scaling_factors
    scaled_data = (unscaled_data * rs + ri) / (rs * ss) # this is just the way it is for philips
    
    return scaled_data



def fun_pcasl_wang_tissue(x, fmdata, pdelta_a, pdelta, pw, ptau, pt1a):
    # minimization function for pCASL
    # Wang J et al. MRM 48:242-254 (2002)
    # Note this assumes that blood water has reached the capillary exchange
    # site
    # note that fmdata is deltaM/M0tissue
    
    f = x[0] # cbf to solve for (ml/g/s)
    a = 0.85 # labeling efficiency
    r1t = 1/1.3 # r1 of unperfused brain tissue (s^-1)
    lam  = 0.9 # blood brain partition coefficient (ml/g). this is lambda
    
    # these variables can be passed in:
    fdelta = pdelta # the tissue transit time (s)
    fw = pw # the post-labeling delay (s)
    ftau = ptau # the labeling duration (s)
    ft1a = pt1a # t1 of arterial blood water (s)
    fdelta_a = pdelta_a # arterial arrial time (s)
    
    # calculate r1a of blood water
    r1a = 1/ft1a # r1 of arterial blood water (s^-1)
    
    # calculate r1app (tissue relaxation rate in presence of exchanged, labeled
    # blood water
    r1app = r1t+f/lam
    
    # the kinetic model
    term_coeff = (2*f*a/lam)
    term_tissue = (np.e**(-fdelta*r1a)/(r1app))*(np.e**(min(fdelta-fw,0)*r1app)-np.e**((fdelta-ftau-fw)*r1app))
    term_blood = (1/r1a)*(np.e**((min(fdelta_a-fw,0)-fdelta_a)*r1a)-np.e**((min(fdelta-fw,0)-fdelta)*r1a))
    kinmodel = term_coeff*(term_tissue+term_blood)
    
    dd = np.sqrt((fmdata-kinmodel)**2) # miniize the difference between the data and them odel
    return dd


def mask_image_by_cornervals(img, box_size, thresh_mod=1):
    '''
    Samples boxes in the corner of an MR image to determine a threshold that is then
    used to create a brain mask. Mask is generated slicewise but returned as a 3d image

    Parameters
    ----------
    img : np array with dimensions x,y,z
        The image to be masked. The z dimension should be the slices
    box_size : array or tuple of ints, len=2
        the x,y dimensions of the corner boxes to be sampled. the x,y dimensions should
        be big enough to sample a sizeable portion of the background noise, but small
        enough that the boxes do not extend into the brain
    thresh_mod : float
        the multiplier that will be applied to the mean of the sampled corners
        to set the masking threshold. Higher values create more conservative
        masks

    Returns
    -------
    a boolean np array with the same dimensions as img

    '''
    
    dx, dy = box_size
    im_dims = img.shape
    msk = np.zeros_like(img)
    
    for i in range(im_dims[2]):
        sli = img[:,:,i]
        
        bottom_left = sli[:dx,:dy]
        bottom_right = sli[-dx:,:dy]
        upper_left = sli[:dx,-dy:]
        upper_right = sli[-dx:,-dy:]
        
        corners_mean = np.mean(np.array([bottom_left, bottom_right, upper_left, upper_right]))
        thresh_val = corners_mean*thresh_mod
        
        masked_sli = sli > thresh_val
        msk[:,:,i] = masked_sli
        
        print(thresh_val)
        
    return msk



def iterative_2d_convolutions(img, impulse_matrix):
    '''
    Simple function that mocks behavior of Hanzhang Lu's 'smooth_filter_spatial'
    function in MATLAB. If impulse_matrix = 
        np.array([[0.1, 0.3, 0.1], [0.3, 1, 0.3], [0.1, 0.3, 0.1]]) / (0.1*4+0.3*4+1)
    then the equivalent mode is 1.
    
    Smooth a 3d image one slice at a time
    

    Parameters
    ----------
    img : np array
        should have three dimensions, and the last dimension is time or slice
    impulse_matrix : np array
        the matrix used to create the finiite impulse response filter.

    Returns
    -------
    the filtered image.

    '''
    
    convolved_img = np.zeros_like(img)
    im_dims = img.shape
    
    for i in range(im_dims[2]):
        sli = img[:,:,i]
        sli_blur = signal.convolve2d(sli, np.rot90(impulse_matrix), mode='same')
        convolved_img[:,:,i] = sli_blur
        
    return convolved_img


def wang_pcasl_func(x,fmdata,pdelta_a,pdelta,pw,ptau,pt1a):
    
    
    # Wang J et al. MRM 48:242-254 (2002)
    # Note this assumes that blood water has reached the capillary exchange
    # site
    # note that fmdata is deltaM/M0tissue
    
    f = x[0] # cbf to solve for (ml/g/s)
    a = 0.85 # labeling efficiency
    r1t = 1/1.3 # r1 of unperfused brain tissue (s^-1)
    lambd = 0.9 # blood brain partition coefficient (ml/g). lambda
    
    # these variables can be passed in:
    fdelta = pdelta # the tissue transit time (s)
    fw = pw # the post-labeling delay (s)
    ftau = ptau # the labeling duration (s)
    ft1a = pt1a # t1 of arterial blood water (s)
    fdelta_a = pdelta_a # arterial arrial time (s)
    
    # calculate r1a of blood water
    r1a = 1/ft1a # r1 of arterial blood water (s^-1)
    
    # calculate r1app (tissue relaxation rate in presence of exchanged, labeled
    # blood water
    r1app = r1t+f/lambd
    
    # the kinetic model
    term_coeff = (2*f*a/lambd)
    term_tissue = (np.exp(-fdelta*r1a)/(r1app))*(np.exp(min([fdelta-fw,0])*r1app)-np.exp((fdelta-ftau-fw)*r1app))
    term_blood = (1/r1a)*(np.exp((min([fdelta_a-fw,0])-fdelta_a)*r1a)-np.exp((min([fdelta-fw,0])-fdelta)*r1a))
    kinmodel = term_coeff*(term_tissue+term_blood)
    
    dd = np.sqrt((fmdata-kinmodel)**2) # minimize the difference between the data and the model
    
    return dd


def read_xfm(xfm, out_loc=None):
    # read in a FS xfm transformation and returns the 4x4 transformation matrix as a np array
    # last row is imputed as 0 0 0 1
    # if out_loc is specified, the matrix is written as a .omat file (FSL compatable) at that location
    
    with open(xfm) as readin:
        lines = readin.read().splitlines()
        spec = lines[4]
        if spec != 'Linear_Transform = ':
            raise Exception(f'read_xfm assumes the transformation is linear. This does not appear to be the case:\n\tspec={spec}')
        text_mat = lines[-3:] # the last three lines are the transformation matrix
        stripped = [j.strip(';').split(' ') for j in text_mat] #  the numbers are space delimited, and the last line has a semicolon that needs to be stripped
        num_mat = [[float(i) for i in j] for j in stripped] # convert each number to a float.
        last_line = [0,0,0,1]
        
        num_mat.append(last_line)
        
        final_mat = np.array(num_mat)
        
    
    if out_loc:
        save_mat = str(num_mat)
        save_mat = save_mat.replace('], ','\n')
        remove_chars = ['[', ']', ',', ';']
        for rc in remove_chars:
            save_mat = save_mat.replace(rc, '')
            
        with open(out_loc, 'w') as readout:
            readout.write(str(save_mat))
        
    return final_mat

        

    
    
    
    