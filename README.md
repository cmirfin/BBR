# BBR
Multimodal non-rigid boundary based registration
## Synopsis
A novel two-stage non-rigid registration of MR images to histological sections is de- scribed. The method uses state-of-the art modality independent neighbourhood descriptors (MIND) to form a cost function which is placed within an efficient Gauss-Newton optimization scheme.
A second registration stage involves a modified boundary-based cost function which attempts to match tissue boundaries across modalities. Tissue boundaries are assumed to be independent structures which allows boundaries to deform and overlap. The final transformation is modelled as an elastic deformation and optimized through a fixed-point iteration scheme.

## Motivation

The registration of MR images to histological sections is crucial to the work of many current and future human brain studies.
Current registration methods do not deal with large scale tissue distortion such as tissue overlaps that unavoidably occur during histological imaging processes.

## Installation

Install in selected directory. 
Change 'filename1' and 'filename2' in full_script with appropriate paths to location of your images.
Note: nifti package is required to handle mri data.
Note: demos require user to addpath to above folder system.

## Contributors

Christopher J. Mirfin (christopher.mirfin@dtc.ox.ac.uk)
Mattias P. Heinrich for providing MIND code. 

## Acknowledgements

This work was funded by the ESPRC and MRC as part of the Doctoral Training Centre, Oxford, UK and undertaken at the Oxford University Centre for Functional MRI of the Brain, Oxford, UK.

## License

Please cite the appropriate author(s) if you use the code.
