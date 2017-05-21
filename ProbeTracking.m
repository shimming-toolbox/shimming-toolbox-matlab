classdef ProbeTracking < Tracking
% PROBETRACKING - Respiratory probe for real-time shimming 
%
% Tracker = PROBETRACKING(  )
%
%   Tracker contains fields
%           
%       .ComPort    
%
%       .Data 
%
%       .Specs
%
% .......
%
%   Description
%
%
% =========================================================================
% Part of series of classes pertaining to shimming:
%
%    Tracking
%    ShimCal
%    ShimCom
%    ShimEval
%    ShimOpt
%    ShimSpecs
%    ShimTest 
%    ShimUse
%
% =========================================================================
% Updated::20170521::ryan.topfer@polymtl.ca
% =========================================================================

% *** TODO 
%
% ..... 
%
%   check/assert ComPort is in fact open before attempting read
%   (issue is that strcmp( ComPort.status, 'closed' ) is fairly slow (~ 3ms)
%   compared to the 1/10ms sample rate of the probe.      
% =========================================================================

properties   
    % Data ;
    % Specs ;
    ComPort ;
end

% =========================================================================
% =========================================================================
methods
% =========================================================================
function Tracker = ProbeTracking( Specs )
%PROBE - ProbeTracking  

if nargin < 1
    Specs = [] ;
end

[Tracker.ComPort, Tracker.Specs] = ProbeTracking.declareprobe( Specs ) ;

Tracker.Data.p = [];

end    
% =========================================================================
function [] = delete( Tracker )
%DELETE  
% DELETE( Tracker )
% 
% Destructor. Calls Tracker.deletecomport( ) 

if myisfield( Tracker, 'ComPort' ) && ~isempty( Tracker.ComPort ) 
    Tracker.deletecomport();
end

clear Tracker ;

end
% =========================================================================
function [Tracker] = deletecomport( Tracker )
%DELETECOMPORT  
% 
% DELETECOMPORT( Tracker )
% 
% Proper way to close + delete the serial port object

if myisfield( Tracker, 'ComPort' ) && ~isempty( Tracker.ComPort ) 
    fclose( Tracker.ComPort ) ;
    delete( Tracker.ComPort ) ;
    clear Tracker.ComPort  ;
else
    disp('Failed to delete .ComPort - it does not exist.');
end

end
% =========================================================================
function [isTracking] = begintracking( Tracker )
%BEGINTRACKING - Initialize (RS-232) Communication Port 
%
% [isTracking] = BEGINTRACKING( Tracker )
%
% Opens Tracker.ComPort, assigns sampling frequency (whereby probe begins sampling
% (filling internal buffer))
%
% Returns TRUE if successful

isTracking = false;

maxCommunicationAttempts = 3; 

fopen(Tracker.ComPort);

iAttempt = 1 ;
disp('Connecting to respiratory probe...')

while(~isTracking && iAttempt <= maxCommunicationAttempts )

    % Send refresh time 
    disp(['Attempt #' num2str(iAttempt)]);    
    
    firstWord = fscanf( Tracker.ComPort, '%f') 
    
    if( ~isempty(firstWord) && isnumeric(firstWord) )  
        isTracking = true;
    end

    iAttempt = iAttempt + 1;

end

if isTracking
    disp('Communication successful. Reading in from serial port...')
else
    Tracker.stoptracking() ;
    error('Communication to respiratory probe failed. Closing Tracker.ComPort.')
end

end
% =========================================================================
function [] = stoptracking( Tracker )
%STOPTRACKING - close (RS-232) Communication Port 
%
% [] = STOPTRACKING( Tracker )
%
% fclose( Tracker.ComPort )

fclose(Tracker.ComPort);

end
% =========================================================================
function [weightedAvg] = calculateweightedaverage( Tracker, normalizedWeights )
%CALCULATEWEIGHTEDAVERAGE
%
% weightedAvg = CALCULATEWEIGHTEDAVERAGE( Tracker, normalizedWeights )
% 
% weightedAvg = normalizedWeights .* ... 
%    Tracker.Data.p( (end-(numel(normalizedWeights)-1)) : end) ;

weightedAvg = normalizedWeights .* ...
    Tracker.Data.p( (end-(numel(normalizedWeights)-1)) : end) ;

end
% =========================================================================
function [p] = getupdate( Tracker )
%GETUPDATE 
%
% Reads a single (32bit) pressure measurement (p) in from open com port.
%
% p = GETUPDATE( Tracker )

assert( strcmp( Tracker.ComPort.Status, 'open' ), 'Error: Serial port is closed.' );

p = fscanf( Tracker.ComPort,'%f') ;

end  
% =========================================================================
function [pressureLog, sampleTimes] = recordandplotpressurelog( Tracker, Params ) 
%RECORDANDPLOTPRESSURELOG    
%
%   Description
%   
%   Reads pressure data from microcontroller 
%
%   Syntax
%
%   [pressureLog, sampleTimes] = RECORDANDPLOTPRESSURELOG( Tracker, Parameters )
%
%    .......................
%   
%   The following Parameters.fields are supported
%
%   .isSavingData
%       default = true
%
%   .isForcingOverwrite 
%       Will overwrite the log files if they already exist. 
%       default = false
%
%   .pressureLogFilename
%       default = ['./' datestr(now,30) '-pressureLog.bin' ] ; 
%
%   .sampleTimesFilename
%       default = ['./' datestr(now,30) '-sampleTimes.bin' ] ;
%
%   .runTime 
%       Total sampling time in seconds.
%       default = 30

assert( strcmp( Tracker.ComPort.Status, 'closed' ), 'Error: Serial port is open/in use.' );

DEFAULT_ISSAVINGDATA          = true ;
DEFAULT_ISFORCINGOVERWRITE    = false ;
DEFAULT_PRESSURELOGFILENAME   = ['./' datestr(now,30) '-pressureLog.bin' ] ;
DEFAULT_SAMPLETIMESFILENAME   = ['./' datestr(now,30) '-sampleTimes.bin' ] ;

if  nargin < 2 || isempty(Params)
    Params.dummy = [] ;
end
 
if  ~myisfield( Params, 'isSavingData' ) || isempty(Params.isSavingData)
    Params.isSavingData = DEFAULT_ISSAVINGDATA ;
end

if  ~myisfield( Params, 'isForcingOverwrite' ) || isempty(Params.isForcingOverwrite)
    Params.isForcingOverwrite = DEFAULT_ISFORCINGOVERWRITE ;
end

if  ~myisfield( Params, 'pressureLogFilename' ) || isempty(Params.pressureLogFilename)
    Params.pressureLogFilename = DEFAULT_PRESSURELOGFILENAME ;
end

if  ~myisfield( Params, 'sampleTimesFilename' ) || isempty(Params.sampleTimesFilename)
    Params.sampleTimesFilename  = DEFAULT_SAMPLETIMESFILENAME ; 
end

msg = input( ['\n Press [Enter] to begin recording pressure log.\n'], 's' );
assert( isempty(msg), 'Cancelled calibration.')

% ------- 
isUserSatisfied = false ;

while ~isUserSatisfied
    
    % ------- 
    Tracker.Data.p = 0 ;
    Tracker = Tracker.recordpressurelog( Params ) ;
    
    pressureLog = Tracker.Data.p ;    
    sampleTimes = (Tracker.Specs.dt/1000)*[0:(numel(Tracker.Data.p)-1)] ;
   
    % ------- 
    fprintf(['\n ----- \n Plotting pressure for user verification... \n']) ;
    figure ;
    plot( sampleTimes, Tracker.Data.p, '+' ) ;
    if Params.isSavingData
        title( Params.pressureLogFilename ) ;
    else
        title( 'Pressure log' )
    end
    xlabel('Time (s)');
    ylabel('Pressure (0.01 mbar)');
    
    response = input(['\n Is the current pressure log satisfactory? ' ...
        'Enter 0 to rerecord; 1 to continue. \n']) ;

    isUserSatisfied = logical(response) ;

end

% ------- 
if Params.isSavingData
    pressureLogFid = fopen( Params.pressureLogFilename, 'w+' ) ;
    fwrite( pressureLogFid, Tracker.Data.p, 'double' ) ;
    fclose( pressureLogFid );

    sampleTimesFid = fopen( Params.sampleTimesFilename, 'w+' ) ;
    fwrite( sampleTimesFid, sampleTimes, 'double' ) ;
    fclose( sampleTimesFid );
end

end
% =========================================================================
function [Tracker] = recordpressurelog( Tracker, Params )
%RECORDPRESSURELOG  
%
% Continuously tracks respiratory probe.
%
% Tracker = TRACKPROBE( Tracker, Params )
%
% Params
%   .runTime
%       [default : 30 s]

DEFAULT_RUNTIME = 30 ; % [units : s]

if  nargin < 2 || isempty(Params)
    Params.dummy = [] ;
end

if  ~myisfield( Params, 'runTime' ) || isempty(Params.runTime)
    Params.runTime = DEFAULT_RUNTIME ;
end
    
% ------- 
StopButton = stoploop({'Stop recording'}) ;
Tracker.Data.p = [] ;

nSamples = Params.runTime / (Tracker.Specs.dt/1000) ;
iSample  = 1 ; 

isTracking = Tracker.begintracking() ;

assert(isTracking, 'Could not begin tracking. Check device connection.')
while ( iSample < nSamples ) && ~StopButton.Stop()

    Tracker.Data.p(end+1) = Tracker.getupdate() ;
    iSample = iSample + 1 ;

end
Tracker.stoptracking();

StopButton.Clear() ;

end
% =========================================================================
end

% =========================================================================
% =========================================================================

% =========================================================================
% =========================================================================
methods(Static)
% =========================================================================
function [ComPort, TrackerSpecs] = declareprobe( TrackerSpecs )
%DECLAREPROBE Declares serial object for probe 
% 
% [ComPort,TrackerSpecs] = declareprobe( TrackerSpecs )
%
% TrackerSpecs can have the following fields
% 
% .dt
%   Interval between pressure samples [units: ms]
%   Should be a positive multiple of the minimum period of 10 ms.
%
% .portName 
%   Address of the probe-associated serial port within file system
%   default: 
%       if ismac 
%           portName = '/dev/tty.usbmodem*'
%       elseif isunix
%           portName = '/dev/ttyS100'

ShimUse.display( ['\n----- Pressure probe -----\n'] );

MIN_ARDUINOPERIOD     = 10 ; % [units: ms] 
DEFAULT_ARDUINOPERIOD = 10 ; % [units: ms] 

if nargin < 1 || isempty(TrackerSpecs)
    ShimUse.display('Default parameters will be used:')
    TrackerSpecs.dummy = [] ;
end

if  ~myisfield( TrackerSpecs, 'dt' ) || isempty(TrackerSpecs.dt) ...
    || (TrackerSpecs.dt < MIN_ARDUINOPERIOD) 

        TrackerSpecs.dt = DEFAULT_ARDUINOPERIOD ;
else
    
    skipSampleFactor = round( TrackerSpecs.dt / MIN_ARDUINOPERIOD ) ;
    TrackerSpecs.dt = skipSampleFactor * MIN_ARDUINOPERIOD ;

end


% -------------------------------------------------------------------------
% check and/or assign serial port
ComPort = [] ; 
isAssigningSerialPort = false ;

if  myisfield( TrackerSpecs, 'portName' ) && ~isempty(TrackerSpecs.portName)
    
    [fileDir,portname,fileExtension] = fileparts( TrackerSpecs.portName ) ;
    
    listOfDevices = dir( fileDir ) ;

    isPortDetected = false ;
    
    % Check for file existence. For some reason this doesn't work...
    for iFile = 1 : length(listOfDevices)-1
        if strcmp( [portname fileExtension], Tmp(iFile).name )
            isPortDetected = true ;
        end
    end

    if ~isPortDetected    
        disp(['Warning: Given port name (' TrackerSpecs.portName ' ) not found. ' ...
            'Checking default port names.']) ;
        TrackerSpecs.portName = [] ;
    end
end

if  ~myisfield( TrackerSpecs, 'portName' ) || isempty(TrackerSpecs.portName)

    if ismac
        % For some reason this doesn't work...
        % -------
        % if exist( '/dev/tty.usbmodem1411', 'file' ) ; % my macbook's left USB port
        %   TrackerSpecs.portName = '/dev/tty.usbmodem1411';
        % elseif exist( '/dev/tty.usbmodem1421', 'file' ) ;
        %     TrackerSpecs.portName = '/dev/tty.usbmodem1421' ; % my macbook's right USB port
        % else
        %     error(['Device file not found. ' ...
        %         'Is the microcontroller connected?']) ; 
        % ------
        listOfDevices = dir( '/dev/tty.usbmodem*' ) ;

        if isempty(listOfDevices) 
            ShimUse.display( 'Error: Device file not found. Check USB device is connected.' ) ;
        elseif length(listOfDevices) ~= 1, ...
            ShimUse.display( ['Error: Ambiguous device identifier.' ...
                'Enter device file name as function argument. See HELP TrackerTracking.declareprobe.'] ) ;
        else
            isAssigningSerialPort = true ;
            TrackerSpecs.portName = ['/dev/' listOfDevices.name] ;
        end

    elseif isunix
        if exist('/dev/ttyS100', 'file')  
            TrackerSpecs.portName = '/dev/ttyS100' ;
        else
            errorMsg = ['Device file ( ' portName ' ) not found.'  ...
            'Is the microcontroller connected?' ...
            'Ensure symbolic link exists between actual Arduino device' ...
            'file and the phantom-device file discoverable by MATLAB:' ...
            'In a terminal, type: ' ... 
            'ln /dev/ttyACM0 /dev/ttyS100' ] ; 
            error( errorMsg ) ;
        end

    else
        error( 'OS not supported' ) ;
    end
end


if isAssigningSerialPort
    ComPort = serial( TrackerSpecs.portName ) ;

    ShimUse.display( [ 'Arduino sampling frequency = ' ...
        num2str(1000/TrackerSpecs.dt) ' Hz'] )
    % ComPort = serial( TrackerSpecs.portName, ...
    %   'InputBufferSize', 8 ) ; % 2x 32-bit floats
else
    ComPort = [] ;
end

end
% =========================================================================
function [pressureLog, sampleTimes] = loadpressurelog( pressureLogFilename, sampleTimesFilename )
%LOADPRESSURELOG
% 
% Wraps to Tracking.loadmeasurementlog()
%
% pressureLog                = LOADPRESSURELOG( pressureLogFilename ) ;
% [pressureLog, sampleTimes] = LOADPRESSURELOG( pressureLogFilename, sampleTimesFilename )

if nargin < 1
    error( 'Insufficient arguments. Must provide full path to pressure log .bin file.' ) ;
else
    if nargin == 1
        [pressureLog] = Tracking.loadmeasurementlog( pressureLogFilename ) ;
    elseif nargin == 2
        [pressureLog, sampleTimes] = ...
            Tracking.loadmeasurementlog( pressureLogFilename, sampleTimesFilename )
    end
end

end
% =========================================================================
function [] = plotpressurelog( pressureLog, Params )
%PLOTPRESSURELOG
% 
% Wraps to Tracking.plotmeasurementlog()
%
% PLOTPRESSURELOG( pressureLog ) ;
% PLOTPRESSURELOG( pressureLog, Params )
%
% Supported fields to Params struct
%
%   .figureTitle
%       [default: 'Pressure log']
%
%   .sampleTimes
%       vector (length == length(pressureLog)) of sample times in seconds

DEFAULT_FIGURETITLE = 'Pressure log' ;
DEFAULT_YLABEL      = 'Pressure (kPa)' ;

if nargin < 1
    error( 'Insufficient arguments. Must provide measurement log vector.' ) ;
end

if nargin == 1 || isempty( Params ) 
    Params.dummy = [] ;
end

if ~myisfield( Params, 'figureTitle' ) || isempty( Params.figureTitle ) 
    Params.figureTitle = DEFAULT_FIGURETITLE ;
end

if ~myisfield( Params, 'yLabel' ) || isempty( Params.yLabel ) 
    Params.yLabel = DEFAULT_YLABEL ;
end

Tracking.plotmeasurementlog( pressureLog, Params ) ;

end
% =========================================================================
function [medianPressure] = userselectmedianpressure( pressureLog )
% USERSELECTMEDIANPRESSURE 
%
% Wraps to Tracking.userselectmedianmeasurement()
%
% medianPressure = USERSELECTMEDIANPRESSURE( pressureLog ) 
%
% Plots pressureLog and the user selects START and END (apnea) indices over
% which to calculate the median. The median pressure is superposed over the
% pressureLog graph and the user is asked if the result is satisfactory (or
% redo).

medianPressure = Tracking.userselectmedianmeasurement( pressureLog )

end
% =========================================================================

end
% =========================================================================
% =========================================================================

end
