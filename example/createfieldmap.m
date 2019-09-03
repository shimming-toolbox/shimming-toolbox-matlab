%PROCESSFIELDMAP
%
% e.g. script to create field map from GRE dicoms
%
% =========================================================================
% Author: ryan.topfer@polymtl.ca
% =========================================================================




% paths to DICOM folders:
pathToMag   = '/Users/ryan/Projects/Shimming/Acdc/20190511_Ismrm/data/acdc_70/22-realtime_fieldmap_deepBreath' ;
pathToPhase = '/Users/ryan/Projects/Shimming/Acdc/20190511_Ismrm/data/acdc_70/23-realtime_fieldmap_deepBreath' ;

%% -------
% Define a few processing parameters 
% for more info, see HELP FieldEval.mapfield

% relative to max magnitude: defines the binary mask for phase unwrapping
Params.threshold = 0.01 ;

%% -------
% Create FieldEval object:
Field = FieldEval( pathToMag, pathToPhase, Params ) ;

%% =========================================================================
% At this point, the user has access to a number of methods, e.g.: 
% see HELP for 
%
% Field.assessfielddistribution()
% Field.getacquisitiontime()
% Field.resliceimg()
%
% etc.

% e.g. interpolate field maps to match voxel positions of another image:
pathToT2w = '/Users/ryan/Projects/Shimming/Acdc/20190511_Ismrm/data/acdc_70/24-T2w_1mm' ;

T2w = MaRdI( pathToT2w ) ;


% e.g. create NifTI output
% 
% NOTE: 
% This calls MaRdI.write(), which creates temporary DICOMs, followed by a system call to dicm2niix 
% This has not been extensively tested and some of Hdr info may be incomplete:

pathToNii = './field/' ; % where to save the images
isSavingSingleNiis = true ;
Field.write( pathToNii, 'nii', isSavingSingleNiis ) ;
