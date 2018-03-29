function [ Params ] = declareshimparams2018ismrm()
% =========================================================================
% Params for actual shim experiment/acquisition
% =========================================================================
%
% =========================================================================
% Updated::20180328::ryan.topfer@polymtl.ca
% =========================================================================
Params.projectDir = '/home/vhosuser/Projects/Shimming/scripts/shimming_rri'

Params.isDebugging = true ;
% =========================================================================
%  
% =========================================================================
Params.dataLoadDir  = [Params.projectDir '/acdc_21/'] ;

% -------
% MRI 
Params.InstitutionName = 'IUGM' ;
Params.StationName     = 'MRC35049' ;

% -------
% for phase unwrapping:
Params.threshold = 0.05 ; % as percent of max measured intensity.

% -------
% for shimming:
Params.shimSystem = 'Acdc' ;
Params.pathToShimReferenceMaps = [Params.projectDir '/acdc_21/ShimReferenceMaps_Acdc_20180326.mat'] ;

Params.ProbeSpecs    = [] ;
Params.ProbeSpecs.dt = 10 ; % sampling interval [units: ms]

Params.maxCurrentPerChannel = 2 ; % [units: A] 
% Params.maxVoltagePerChannel = 200 ; % [units: mV] 

Params.isSolvingAugmentedSystem    = true ;
Params.isPenalizingFieldDifference = false;
Params.regularizationParameter     = 0.01 ;
Params.runMode='isGui';

Params.isAutoSegmenting = true ;
% -------
% for field mapping

Params.echoTimeDifference = (4.92 - 2.46) ; % [units: ms]

Params.isFilteringField   = false ; 
Params.filterRadius       = 3 ; % [units: mm]

Params.maxAbsField        = 600 ; % [units: Hz]
Params.maxFieldDifference = 250 ; % [units: Hz]
% if the difference of measured field values for inspired vs expired cases
% is greater than this, the measurement is deemed unreliable 

% -------
% for real-time calibration:
Params.Inspired = [] ; 
Params.Expired  = [] ;
Params.nTrainingFrames      = 2 ; % number of field maps to acquire for shim 'training'
Params.pressureLogFilenames = cell( Params.nTrainingFrames, 1) ;
Params.medianPressures      = zeros( Params.nTrainingFrames, 1 ) ;

% -------
Params.isFilteringMeasurements = true ;
Params.isClipping              = true ;
Params.nSamplesFilter          = 5;
Params.correctionOrder         = 1 ; % linear correction
Params.txDelay                 = 150 ; % [units: ms]: approx. 50 ms acoustic delay from tube + UNKNOWN from serial communication.

fprintf('')
Params

fprintf('\n\n\n**** txDelay is UNKNOWN *****\n\n\n')

end
