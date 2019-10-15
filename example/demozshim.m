%DEMOZSHIM
%
% e.g. script to create field map from GRE dicoms
%
% =========================================================================
% Author: ryan.topfer@polymtl.ca
% =========================================================================


%% -------
% load images and physio recording as FieldEval + ProbeTracking objects :

% paths to DICOM folders:
pathToMag   = '/Users/ryan/Projects/Shimming/Acdc/202004_Ismrm/data/acdc_72/07-realtime_fieldmap_deepBreath' ;
pathToPhase = '/Users/ryan/Projects/Shimming/Acdc/202004_Ismrm/data/acdc_72/08-realtime_fieldmap_deepBreath' ;

% field map time series
Fields = FieldEval( pathToMag, pathToPhase ) ; 

% Siemens PMU recording
Pmu   = ProbeTracking( '/Users/ryan/Projects/Shimming/Acdc/202004_Ismrm/data/acdc_72/PMU/acdc_72_02.resp' ) ;

%% -------
% link the two objects (interpolate PMU to the image times)
Fields.associateaux( Pmu ) ;

% modeled static + respiratory fields (in Field.img and Field.Model.Riro.img respectively)
Field = FieldEval.modelfield( Fields ) ;

%% -------
% define the shim system + target field
Shims = ShimOpt_IUGM_Prisma_fit( Field ) ;


%% -----

