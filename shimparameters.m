function [ Params] = shimparameters()
% =========================================================================
% Params for actual shim experiment/acquisition
% =========================================================================
%
% =========================================================================
% Updated::20171107::ryan.topfer@polymtl.ca
% =========================================================================
Params.projectDir = '/Users/ancha_admin/data/20180616_Ismrm' ;

% =========================================================================
%  
% =========================================================================
Params.dataLoadDir  = [Params.projectDir '/data/acdc_15/'] ;

% -------
% for phase unwrapping:
Params.threshold = 0.05 ; % as percent of max measured intensity.

% -------
% for shimming:
Params.pathToShimReferenceMaps = '/Users/ancha_admin/Documents/Acdc/Calibration/data/AcdcReferenceMaps20171107.mat' ;
Params.shimSystem='Acdc';

%Path to Matlab folder
Params.matlabPath='/Users/ancha_admin/Documents/Matlab';

%Command for SCT segmentation (Sct_propseg) with CSF segmentation :

%Params.command =sprintf('%s', 'sct_propseg -i ',Params.matlabPath,'/gre_field_mapping_shim0_ins.nii',' -c t1 ','-ofolder ',Params.matlabPath,' -CSF');
%Params.command2 =sprintf('%s', 'sct_propseg -i ',Params.matlabPath,'/gre_field_mapping_shim0_exp.nii',' -c t1 ','-ofolder ',Params.matlabPath,' -CSF');

%Command for SCT segmentation (Sct_deepseg_sc) :

Params.command =sprintf('%s', 'sct_deepseg_sc -i ',Params.matlabPath,'/gre_field_mapping_shim0_ins.nii',' -c t1 ','-ofolder ',Params.matlabPath);
Params.command2 =sprintf('%s', 'sct_deepseg_sc -i ',Params.matlabPath,'/gre_field_mapping_shim0_exp.nii',' -c t1 ','-ofolder ',Params.matlabPath);

%Command for SCT segmentation (Sct_get_centerline) :
Params.commandbis =sprintf('%s', 'sct_get_centerline -i ',Params.matlabPath,'/gre_field_mapping_shim0_ins.nii',' -c t1 ','-ofolder ',Params.matlabPath);
Params.commandbis2 =sprintf('%s', 'sct_get_centerline -i ',Params.matlabPath,'/gre_field_mapping_shim0_exp.nii',' -c t1 ','-ofolder ',Params.matlabPath);


%Command to call SortData.m in background and sort the folder -------------

Params.calltoSortdata=sprintf('%s','/Applications/MATLAB_R2016a.app/bin/matlab -nodesktop -nojvm -r ", SortData(''');

Params.ProbeSpecs    = [] ;
Params.ProbeSpecs.dt = 10 ; % sampling interval [units: ms]

Params.maxCurrentPerChannel = 2.2 ; % [units: A] 
Params.maxVoltagePerChannel = 2500 ; % [units: mV] Not use in Acdc project

Params.isSolvingAugmentedSystem    = false ;
Params.isPenalizingFieldDifference = false;
Params.regularizationParameter     = 0.01 ;

% -------
% for field mapping


Params.echoTimeDifference = 2.5;   % [units: ms]

Params.isFilteringField   = true ; 
Params.filterRadius       = 3 ; % [units: mm]

Params.maxAbsField        = 600 ; % [units: Hz]
Params.maxFieldDifference = 250 ; % [units: Hz]
% if the difference of measured field values for inspired vs expired cases
% is greater than this, the measurement is deemed unreliable 

% -------
% for real-time calibration:
Params.Inspired = [] ; 
Params.Expired  = [] ;
Params.nCalibrationScans    = 2 ;
Params.pressureLogFilenames = cell( Params.nCalibrationScans, 1) ;
Params.medianPressures      = zeros( Params.nCalibrationScans, 1 ) ;

% -------
Params.isFilteringMeasurements = true ;
Params.isClipping              = true ;
Params.nSamplesFilter          = 5;
Params.correctionOrder         = 1 ; % linear correction
Params.txDelay                 = 1000 ; % [units: ms]: approx. 50 ms acoustic delay from tube + UNKNOWN from serial communication.

%Params.feedbackcalibrationcoeff1=[0.65909, 0.65, 0.64547, 0.6499, 0.6591, 0.654, 0.65, 0.65];
%Params.feedbackcalibrationcoeff2 = [19.09, -10.908, -23.636, 20.91, 32.728, 11.308, -42.728, 34.546];   %Calibration coefficient for the Adc feedback


%fprintf('')
%Params;


end
