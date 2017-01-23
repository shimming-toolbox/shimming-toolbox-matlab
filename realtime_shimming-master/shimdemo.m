% SHIM DEMO 
%
% Script to demonstrate static and real-time shimming
% 
% NB: This serves as a guide for the relevant commands - the script can't
% simply be run automatically/on its own.
% 
% =========================================================================
% Updated::20161129::ryan.topfer@polymtl.ca
% =========================================================================

% =========================================================================
% DECLARE INITIAL VARIABLES 
% =========================================================================
 
Params.dataLoadDir = '/Users/ryan/Projects/Shimming/Dynamic/20170422_ISMRM/data/shim_029/';
mkdir(Params.dataLoadDir) ;

% -------
% for shimming:
Params.shimSystem               = 'Rri' ;
Params.pathToShimReferenceMaps  = '/Users/ryan/Projects/Shimming/Static/Calibration/Data/SpineShimReferenceMaps20161007.mat'
Params.ProbeSpecs               = [] ;
Params.ProbeSpecs.arduinoPeriod = 10 ; % [units: ms]
Params.maxCurrentPerChannel     = 2.5 ;

% -------
% for field mapping
Params.isFilteringField  = true ;

Params.maxAbsField        = 600 ; % [units: Hz]
Params.maxFieldDifference = 150 ; % [units: Hz]
% if the difference of measured field values for inspired vs expired cases
% is greater than this, the measurement is deemed unreliable 

% -------
% for phase unwrapping:
Params.threshold = 0.01 ; % as percent of max measured intensity.

% -------
% for real-time calibration:
Params.PressureLogs         = [] ;
Params.nCalibrationScans    = 2 ;
Params.pressureLogFilenames = cell( Params.nCalibrationScans, 1) ;
Params.medianPressures      = zeros( Params.nCalibrationScans, 1 ) ;

% -------
Params.isFilteringPressure  = false ;

Shims = ShimUse(  Params  ) ;

% =========================================================================
% RECORD PRESSURE LOGS 
% =========================================================================
RecParams.isSavingData        = true;
RecParams.isForcingOverwrite  = true ;
RecParams.runTime             = 5*60 ;
 
% -------
% Normal breathing spine shim OFF
Params.pressureLogFilename = [Params.dataLoadDir datestr(now,30) '-pressureLog-Breathing-ShimOff.bin'] ;
Params.sampleTimesFilename = [Params.dataLoadDir datestr(now,30) '-sampleTimes-Breathing-ShimOff.bin'] ;
[pressureLog, sampleTimes] =Shims.Opt.Probe.recordandplotpressurelog( Params ) ;

% -------
% Inspired condition
RecParams.pressureLogFilename = [Params.dataLoadDir datestr(now,30) '-pressureLog-INS.bin'] ;
RecParams.sampleTimesFilename = [Params.dataLoadDir datestr(now,30) '-sampleTimes-INS.bin'] ;

Params.pressureLogFilenames(1,1) = { RecParams.pressureLogFilename } ;
Params.Inspired.pressureLog = Shims.Opt.Probe.recordandplotpressurelog( RecParams ) ;

% -------
% Expired condition
RecParams.pressureLogFilename = [Params.dataLoadDir datestr(now,30) '-pressureLog-EXP.bin']
RecParams.sampleTimesFilename = [Params.dataLoadDir datestr(now,30) '-sampleTimes-EXP.bin']

Params.pressureLogFilenames(2,1) = { RecParams.pressureLogFilename } ;
Params.Expired.pressureLog = Shims.Opt.Probe.recordandplotpressurelog( RecParams ) ;

% =========================================================================
% PROCESS GRE FIELD MAPS
% =========================================================================
% Change directory names as necessary 

% -------
% Inspired 
Params.Path.Mag.echo1         = [ MaRdI.getfulldir( Params.dataLoadDir, 5 ) 'echo_4.92' ] ;
Params.Path.Mag.echo2         = [ MaRdI.getfulldir( Params.dataLoadDir, 5 ) 'echo_7.64' ] ;

Params.Path.Phase.echo1       = [ MaRdI.getfulldir( Params.dataLoadDir, 6 ) 'echo_4.92' ] ;
Params.Path.Phase.echo2       = [ MaRdI.getfulldir( Params.dataLoadDir, 6 ) 'echo_7.64' ] ;

[FieldInspired,Extras]        = ShimOpt.mapfield( Params ) ;

% -------
% Expired 
Params.Path.Mag.echo1         = [ MaRdI.getfulldir( Params.dataLoadDir, 7 ) 'echo_4.92' ] ;
Params.Path.Mag.echo2         = [ MaRdI.getfulldir( Params.dataLoadDir, 7 ) 'echo_7.64' ] ;

Params.Path.Phase.echo1       = [ MaRdI.getfulldir( Params.dataLoadDir, 8 ) 'echo_4.92' ] ;
Params.Path.Phase.echo2       = [ MaRdI.getfulldir( Params.dataLoadDir, 8 ) 'echo_7.64' ] ;

[FieldExpired,Extras] = ShimOpt.mapfield( Params ) ;

% =========================================================================
% INTERP everything to grid of field 
% =========================================================================
Shims.Opt.interpolatetoimggrid( FieldInspired ) ;

% =========================================================================
% DEFINE SHIM VOI 
% =========================================================================
mask = Shims.Opt.getvaliditymask( Params, FieldInspired, FieldExpired ) ;

% =========================================================================
% STATIC OPTIMIZATION
% =========================================================================

% -------
% Inspired 
Shims.Opt.setoriginalfield( FieldInspired ) ;
Shims.Opt.setshimvolumeofinterest( mask ) ;

Shims.Opt.optimizeshimcurrents( Params ) ;

Params.Inspired.currents = Shims.Opt.Model.currents ;

% -------
% Expired 
Shims.Opt.setoriginalfield( FieldExpired ) ;
Shims.Opt.setshimvolumeofinterest( mask ) ;

Shims.Opt.optimizeshimcurrents( Params ) ;

Params.Expired.currents = Shims.Opt.Model.currents ;


% =========================================================================
% STATIC SHIMMING
% =========================================================================

% -------
% Inspired 
Shims.Com.setandloadallshims( Params.Inspired.currents ) ;
% Params.pressureLogFilename = [Params.dataLoadDir datestr(now,30) '-pressureLog-INS-ShimOn.bin'] ; 
% Params.sampleTimesFilename = [Params.dataLoadDir datestr(now,30) '-sampleTimes-INS-ShimOn.bin'] ;
% [pressureLog, sampleTimes] =Shims.Opt.Probe.recordandplotpressurelog( Params ) ;

Shims.Com.resetallshims() ;

% -------
% Expired 
Shims.Com.setandloadallshims( Params.Expired.currents ) ;
% Params.pressureLogFilename = [Params.dataLoadDir datestr(now,30) '-pressureLog-EXP-ShimOn.bin'] ; 
% Params.sampleTimesFilename = [Params.dataLoadDir datestr(now,30) '-sampleTimes-EXP-ShimOn.bin'] ;
% [pressureLog, sampleTimes] =Shims.Opt.Probe.recordandplotpressurelog( Params ) ;

Shims.Com.resetallshims() ;

% =========================================================================
% CALIBRATE REAL-TIME UPDATES 
% =========================================================================
Shims.Opt.calibraterealtimeupdates( Params ) ;

% =========================================================================
% RUN REAL-TIME SHIMMING
% =========================================================================
% -------
% Normal breathing spine shim ON (real-time)
Params.maxCurrents = max( [Params.Inspired.currents Params.Expired.currents], [],2) ; 
Params.minCurrents = min( [Params.Inspired.currents Params.Expired.currents], [],2) ; 
Params.isFilteringPressure = false ;
Params.isClippingPressure  = true ;
Params.minClippingPressure = 210 ;
Params.maxClippingPressure = 260 ;
Params.pressureLogFilename = [Params.dataLoadDir datestr(now,30) '-pressureLog-Breathing-ShimOn-RT.bin'] ;
Params.sampleTimesFilename = [Params.dataLoadDir datestr(now,30) '-sampleTimes-Breathing-ShimOn-RT.bin'] ;
Shims.runrealtimeshim( Params ) ;

% =========================================================================
% TURN OFF SHIMS
% =========================================================================
