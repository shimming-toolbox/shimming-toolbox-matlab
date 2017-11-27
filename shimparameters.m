function [ Params] = shimparameters()
% =========================================================================
% Params for actual shim experiment/acquisition
% =========================================================================
%
% =========================================================================
% Updated::20171107::ryan.topfer@polymtl.ca
% =========================================================================
Params.projectDir = '/Users/ancha_admin/Documents/Acdc/Code' 

% =========================================================================
%  
% =========================================================================
Params.dataLoadDir  = [Params.projectDir '/20180616_Ismrm/data/acdc_15/'] ;

% -------
% for phase unwrapping:
Params.threshold = 0.05 ; % as percent of max measured intensity.

% -------
% for shimming:
Params.pathToShimReferenceMaps = [Params.projectDir '/Calibration/data/AcdcReferenceMaps20171107.mat'] ;

Params.ProbeSpecs    = [] ;
Params.ProbeSpecs.dt = 10 ; % sampling interval [units: ms]

Params.maxCurrentPerChannel = 0.4 ; % [units: A] 
% Params.maxVoltagePerChannel = 200 ; % [units: mV] 

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

fprintf('')
Params


fprintf('\n\n\n**** txDelay is UNKNOWN *****\n\n\n')

end
