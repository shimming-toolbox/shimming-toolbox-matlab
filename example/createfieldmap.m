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
Params.threshold = 0.01 ;

% assign the unwrapping method (must be 'sunwrap' for single 2d slice)
Params.unwrapper = 'sunwrap' ;

%% -------
% Create FieldEval object:

Field = FieldEval.mapfield( pathToMag, pathToPhase ) ;
 
