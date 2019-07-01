%PROCESSFIELDMAP
%
% e.g. script to create field map from GRE dicoms
%
% =========================================================================
% Author: ryan.topfer@polymtl.ca
% =========================================================================

% paths to DICOM folders:
pathToMag   = '/Users/ryan/Projects/Shimming/Acdc/20190511_Ismrm/data/acdc_69/13-realtime_fieldmap' ;
pathToPhase = '/Users/ryan/Projects/Shimming/Acdc/20190511_Ismrm/data/acdc_69/14-realtime_fieldmap' ;

%% -------
% Define a few processing parameters 
% for more info, see HELP FieldEval.mapfield

% relative to max magnitude: defines the binary mask for phase unwrapping
Params.threshold = 0.1 ;

%% -------
% Create FieldEval object:
Field = FieldEval( pathToMag, pathToPhase, Params ) ;

%% =========================================================================
% At this point, the user has access to a number of methods, e.g.: 
% see HELP for 
%
% Field.assessfielddistribution()
% Field.getacquisitiontime()
%
% etc.

% e.g. create NifTI output
% 
% NOTE: 
% This calls MaRdI.write(), which creates temporary DICOMs, followed by a system call to dicm2niix 
% This has not been extensively tested and some of Hdr info may be incomplete:

pathToNii = './field/' ; % where to save the images
Field.write( pathToNii, 'nii' )
