# Structuring your processing code

With the addition of new sequences to standard acquisitions and updated processing/reporting requirements, it may be neccessary to modify processing code or write new processing modules. When doing so keep in the mind the **Ultimate and Most Supreme Rules of Writing Processing Code But Also Code in General**:

 1. Your functions and modules should cover the narrowest scope possible.
 2. The exact nature of the input and output of your functions and modules should be clearly defined.
 3. Your edits to functions and modules must not break other functions and modules, or pipelines. 

Consider a function *process_asl()*. 

	def process_asl(participant_id, is_pcasl=True):
		asl_location = find_participant_asl(patient_id) # get the text filepath of the asl
		t1_location = find_participant_t1(patient_id)
		
		asl_readin = read_mri(asl_location)
		t1_readin = read_mri(t1_location)

		asl_in_t1space = warp(asl_readin, t1_readin)
		
		if is_pcasl=True:
			some_function = func_a
		else:
			some_function = func_b
		
		asl_data = asl_in_t1space.data
		asl_data = gaussian_blur(asl_data)
		cbf_map = voxelwise_apply(asl_data, some_function)
		
		wm, gm = segment_from_t1(t1_readin.data)
		wm_cbf, gm_cbf = regional_analysis(cbf_map, [wm, gm])

		return wm_cbf, gm_cbf

This function takes in a participant's ID, uses that to search for a suitable ASL and T1 image, reads those in, co-registers them, smooths the data, calculates the CBF in T1 space, segments the WM and GM from the T1 and then uses that to report the CBF within those regions. This is a very common thing to do with our ASL data, but the scope of this function is more like that of a pipeline. It makes a lot of assumptions about how your data is structured, what format its in, and what you want out of it. These assumptions make it difficult to fit the function into a pipeline outside of what it was originally written for, and if you make modifications to the function to fit new requirements for one pipeline you may break something in another pipeline. Be specific and be narrow in scope when writing functions!
	
	def asl_to_cbf(asl_data, is_pcasl=True):
	    """
	    Applies a function to ASL source data to create a voxel-wise estimate of CBF.

	    Parameters
	    ----------
	    asl_matrix : 3d numpy array of floats
	        the m0 data
	    is_pcasl : bool
	        optional. If true, uses a function appropriate for pCASL for quantifying CBF. 
					        Otherwise assumes PASL and uses appropriate function
	    Returns
	    -------
	    3d numpy array of floats with same shape as asl_matrix,
	    giving the estimated CBF in each voxel
	    """
		if is_pcasl=True:
			some_function = func_a
		else:
			some_function = func_b
			
		cbf_map = voxelwise_apply(asl_data, some_function)

		return cbf_map

This function is significantly  reduced in scope. It just takes the ASL data (the voxel data itself—not a path to the file, or a data structure that includes both the voxel data and metadata!) and returns the CBF map, applying a different quantification function depending on if the supplied data is pCASL or not. By leaving it up to the user to read in the data from file and then register it to a new space or blur it (if desired),  this function is much more modular and can easily be a part of larger functions. Returning the simple CBF map has the same effect: we return the most general result possible. Users can then get regional statistics etc. from it if they want, but we don't make that decision for them.