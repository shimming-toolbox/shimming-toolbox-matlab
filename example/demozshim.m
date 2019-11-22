%DEMOZSHIM
%
% =========================================================================
% Author: ryan.topfer@polymtl.ca
% =========================================================================

%% -----
% Sort unorganized images from DICOM socket transfer 
% isCopying (boolean) is optional (move or copy the .IMA files)
MaRdI.sortimages( unsortedDicomDir, sortedDicomDir, isCopying )

% load t2w magnitude for SCT segmentation
Mag = MaRdI( '/Users/ryan/Projects/Shimming/Acdc/202004_Ismrm/data/acdc_73/05-T2w_1mm' ) ;

%% -------
% load images and physio recording as FieldEval + ProbeTracking objects :

% paths to DICOM folders:
pathToMag   = '/Users/ryan/Projects/Shimming/Acdc/202004_Ismrm/data/acdc_73/08-realtime_fieldmap_normalBreath' ;
pathToPhase = '/Users/ryan/Projects/Shimming/Acdc/202004_Ismrm/data/acdc_73/09-realtime_fieldmap_normalBreath' ;

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
% NAIVE Z-SHIM:

Gz = Shims.img(:,:,:,3) ;

zValues = -Field.img ./Gz ;
