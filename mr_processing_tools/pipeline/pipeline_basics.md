
# Pipeline basics

## Components

Most pipelines are run as a wrapper script from the command line, e.g.,

	process_scd_prebio.py -i some/path/to/a/folder -a args

This wrapper glues together of the processing pipeline, which generally consists of:

 1. Pre-processing (standardization) of images, such as co-registering images or quantifying functional imagery
 2. Post-processing (analysis) of images, such as finding the SSS in a TRUST image and extracting the T<sub>2</sub> relaxation or finding the mean CBF in a lobe
 3. Report generation
 4. Pushing results to REDCap (if appplicable)

Details on how to use each pipeline are available within the applicable pipeline package.