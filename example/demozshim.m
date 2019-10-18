%DEMOZSHIM
%
% e.g. script to create field map from GRE dicoms
%
% =========================================================================
% Author: ryan.topfer@polymtl.ca
% =========================================================================

% load t2w magnitude for SCT segmentation
Mag = MaRdI( '/Users/ryan/Projects/Shimming/Acdc/202004_Ismrm/data/acdc_73/05-T2w_1mm' ) ;

%% -------
% load images and physio recording as FieldEval + ProbeTracking objects :

% paths to DICOM folders:
pathToMag   = '/Users/ryan/Projects/Shimming/Acdc/202004_Ismrm/data/acdc_73/08-realtime_fieldmap_normalBreath' ;
pathToPhase = '/Users/ryan/Projects/Shimming/Acdc/202004_Ismrm/data/acdc_73/09-realtime_fieldmap_normalBreath' ;

% /Users/ryan/Projects/Shimming/Acdc/202004_Ismrm/data/acdc_73/08-realtime_fieldmap_normalBreath
% /Users/ryan/Projects/Shimming/Acdc/202004_Ismrm/data/acdc_73/09-realtime_fieldmap_normalBreath

% field map time series
Fields = FieldEval( pathToMag, pathToPhase ) ; 

% Siemens PMU recording
Pmu   = ProbeTracking( '/Users/ryan/Projects/Shimming/Acdc/202004_Ismrm/data/acdc_73/Probes/acdc_73_01.resp' ) ;

%% -------
% link the two objects (interpolate PMU to the image times)
Fields.associateaux( Pmu ) ;

% modeled static + respiratory fields (in Field.img and Field.Model.Riro.img respectively)
Field = FieldEval.modelfield( Fields ) ;

%% -------
% define the shim system + target field
Shims = ShimOpt_IUGM_Prisma_fit( Field ) ;


%% -----

