classdef ProbeTracking < matlab.mixin.SetGet
% PROBETRACKING - Respiratory probe for real-time shimming 
%
% Probe = PROBETRACKING(  )
%
%   Probe contains fields
%           
%       .ComPort    
%
%       .Data 
%
%       .Specs
% .......
%
%   Description
%
%
% =========================================================================
% Part of series of classes pertaining to shimming:
%
%    ProbeTracking
%    ShimCal
%    ShimCom
%    ShimEval
%    ShimOpt
%    ShimSpecs
%    ShimTest 
%    ShimUse
%
% =========================================================================
% Updated::ryan.topfer@polymtl.ca::Mon 20 Jun 2016 18:15:05 EDT
% =========================================================================

% *** TODO 
%
% ..... 
%   check/assert ComPort is in fact open before attempting read
%   (issue is that strcmp( ComPort.status, 'closed' ) is fairly slow (~ 3ms)
%   compared to the 1/10ms sample rate of the probe.      
% =========================================================================

properties   
    ComPort ;
    Data ;
    Specs ;
end

% =========================================================================
% =========================================================================
methods
% =========================================================================
function Probe = ProbeTracking( Specs )
%PROBE - ProbeTracking  

if nargin < 1
    Specs = [] ;
end

[Probe.ComPort, Probe.Specs] = ProbeTracking.declareprobe( Specs ) ;

Probe.Data.pressure = [];

end    
% =========================================================================
function [] = delete( Probe )
%DELETE  
% DELETE( Probe )
% 
% Destructor. Calls Probe.deletecomport( ) 

if myisfield( Probe, 'ComPort' ) && ~isempty( Probe.ComPort ) 
    Probe.deletecomport();
end

clear Probe ;

end
% =========================================================================
function [Probe] = deletecomport( Probe )
%DELETECOMPORT  
% Probe = DELETECOMPORT( Probe )
% 
% Proper way to close + delete the serial port object

if myisfield( Probe, 'ComPort' ) && ~isempty( Probe.ComPort ) 
    fclose( Probe.ComPort ) ;
    delete( Probe.ComPort ) ;
    clear Probe.ComPort  ;
else
    disp('Failed to delete .ComPort - it does not exist.');
end

end
% =========================================================================
function Probe = initializecomport( Probe )
%INITIALIZECOMPORT - Initialize (RS-232) Communication Port 
%
% [Probe] = INITIALIZECOMPORT( Probe )
%
% Opens Probe.ComPort, assigns sampling frequency (whereby probe begins sampling
% (filling internal buffer))

maxCommunicationAttempts = 3; 

fopen(Probe.ComPort);

isSamplingFrequencyReceived = false;

iAttempt = 1 ;
disp('Attempting to assign probe sample frequency (may take about 15 s)...')

while(~isSamplingFrequencyReceived && iAttempt <= maxCommunicationAttempts )

    % Send refresh time 
    disp(['Attempt #' num2str(iAttempt)]);    
    fprintf( Probe.ComPort, '%d', Probe.Specs.arduinoPeriod);
    pause(0.1)
    firstWord = fscanf( Probe.ComPort, '%f') 
    
    if(~isempty(firstWord) && isnumeric(firstWord) && firstWord == Probe.Specs.arduinoPeriod/10)  
        isSamplingFrequencyReceived = true;
    end

    iAttempt = iAttempt + 1;

end

if isSamplingFrequencyReceived 
    disp('Communication successful. Reading in from serial port...')
else
    Probe.deletecomport() ;
    error('Communication to respiratory probe failed. Deleting serial object Probe.ComPort.')
end

end
% =========================================================================
function [weightedAvg] = calculateweightedaverage( Probe, normalizedWeights )
%CALCULATEWEIGHTEDAVERAGE
%
% weightedAvg = CALCULATEWEIGHTEDAVERAGE( Probe, normalizedWeights )
% 
% weightedAvg = normalizedWeights .* ... 
%    Probe.Data.pressure( (end-(numel(normalizedWeights)-1)) : end) ;

weightedAvg = normalizedWeights .* ...
    Probe.Data.pressure( (end-(numel(normalizedWeights)-1)) : end) ;

end
% =========================================================================
function [p] = readpressure( Probe )
%READPRESSURE 
%
% Reads one (4-byte) pressure measurement (p) in from open com port.
%
% p = LOGPRESSURE( Probe )

assert( strcmp( Probe.ComPort.Status, 'open' ), 'Error: Serial port is closed.' );

p = fscanf( Probe.ComPort,'%f') ;

end  
% =========================================================================
function [pressureLog, sampleTimes] = recordandplotpressurelog( Probe, Params ) 
%RECORDANDPLOTPRESSURELOG    
%
%   Description
%   
%   Reads pressure data from microcontroller 
%
%   Syntax
%
%   [pressureLog, sampleTimes] = RECORDANDPLOTPRESSURELOG( Probe, Parameters )
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

assert( strcmp( Probe.ComPort.Status, 'closed' ), 'Error: Serial port is open/in use.' );

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
    Probe.Data.pressure = 0 ;
    Probe = Probe.recordpressurelog( Params ) ;
    
    pressureLog = Probe.Data.pressure ;    
    sampleTimes = (Probe.Specs.arduinoPeriod/1000)*[0:(numel(Probe.Data.pressure)-1)] ;
   
    % ------- 
    fprintf(['\n ----- \n Plotting pressure for user verification... \n']) ;
    figure ;
    plot( sampleTimes, Probe.Data.pressure, '+' ) ;
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
    fwrite( pressureLogFid, Probe.Data.pressure, 'double' ) ;
    fclose( pressureLogFid );

    sampleTimesFid = fopen( Params.sampleTimesFilename, 'w+' ) ;
    fwrite( sampleTimesFid, sampleTimes, 'double' ) ;
    fclose( sampleTimesFid );
end

end
% =========================================================================
function [Probe] = recordpressurelog( Probe, Params )
%RECORDPRESSURELOG  
%
% Continuously tracks respiratory probe.
%
% Probe = TRACKPROBE( Probe, Params )
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
Probe.Data.pressure = [] ;

nSamples = Params.runTime / (Probe.Specs.arduinoPeriod/1000) ;
iSample  = 1 ; 

Probe = Probe.initializecomport() ;

while ( iSample < nSamples ) && ~StopButton.Stop()

    Probe.Data.pressure(end+1) = Probe.readpressure() ;
    iSample = iSample + 1 ;

end

fclose(Probe.ComPort)

StopButton.Clear() ;

end
% =========================================================================
end

% =========================================================================
% =========================================================================
methods(Static)
% =========================================================================
function [ComPort, ProbeSpecs] = declareprobe( ProbeSpecs )
%DECLAREPROBE Declares serial object for probe 
% 
% [ComPort,ProbeSpecs] = declareprobe( ProbeSpecs )
%
% ProbeSpecs can have the following fields
% 
% .arduinoPeriod
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

if nargin < 1 || isempty(ProbeSpecs)
    ShimUse.display('Default parameters will be used:')
    ProbeSpecs.dummy = [] ;
end

if  ~myisfield( ProbeSpecs, 'arduinoPeriod' ) || isempty(ProbeSpecs.arduinoPeriod) ...
    || (ProbeSpecs.arduinoPeriod < MIN_ARDUINOPERIOD) 

        ProbeSpecs.arduinoPeriod = DEFAULT_ARDUINOPERIOD ;
else
    
    skipSampleFactor = round( ProbeSpecs.arduinoPeriod / MIN_ARDUINOPERIOD ) ;
    ProbeSpecs.arduinoPeriod = skipSampleFactor * MIN_ARDUINOPERIOD ;

end


% -------------------------------------------------------------------------
% check and/or assign serial port
ComPort = [] ; 
isAssigningSerialPort = false ;

if  myisfield( ProbeSpecs, 'portName' ) && ~isempty(ProbeSpecs.portName)
    
    [fileDir,portname,fileExtension] = fileparts( ProbeSpecs.portName ) ;
    
    listOfDevices = dir( fileDir ) ;

    isPortDetected = false ;
    
    % Check for file existence. For some reason this doesn't work...
    for iFile = 1 : length(listOfDevices)-1
        if strcmp( [portname fileExtension], Tmp(iFile).name )
            isPortDetected = true ;
        end
    end

    if ~isPortDetected    
        disp(['Warning: Given port name (' ProbeSpecs.portName ' ) not found. ' ...
            'Checking default port names.']) ;
        ProbeSpecs.portName = [] ;
    end
end

if  ~myisfield( ProbeSpecs, 'portName' ) || isempty(ProbeSpecs.portName)

    if ismac
        % For some reason this doesn't work...
        % -------
        % if exist( '/dev/tty.usbmodem1411', 'file' ) ; % my macbook's left USB port
        %   ProbeSpecs.portName = '/dev/tty.usbmodem1411';
        % elseif exist( '/dev/tty.usbmodem1421', 'file' ) ;
        %     ProbeSpecs.portName = '/dev/tty.usbmodem1421' ; % my macbook's right USB port
        % else
        %     error(['Device file not found. ' ...
        %         'Is the microcontroller connected?']) ; 
        % ------
        listOfDevices = dir( '/dev/tty.usbmodem*' ) ;

        if isempty(listOfDevices) 
            ShimUse.display( 'Error: Device file not found. Check USB device is connected.' ) ;
        elseif length(listOfDevices) ~= 1, ...
            ShimUse.display( ['Error: Ambiguous device identifier.' ...
                'Enter device file name as function argument. See HELP ProbeTracking.declareprobe.'] ) ;
        else
            isAssigningSerialPort = true ;
            ProbeSpecs.portName = ['/dev/' listOfDevices.name] ;
        end

    elseif isunix
        if exist('/dev/ttyS100', 'file')  
            ProbeSpecs.portName = '/dev/ttyS100' ;
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
    ComPort = serial( ProbeSpecs.portName ) ;

    ShimUse.display( [ 'Arduino sampling frequency = ' ...
        num2str(1000/ProbeSpecs.arduinoPeriod) ' Hz'] )
    % ComPort = serial( ProbeSpecs.portName, ...
    %   'InputBufferSize', 8 ) ; % 2x 32-bit floats
else
    ComPort = [] ;
end

end
% =========================================================================
function [pressureLog, sampleTimes] = loadpressurelog( pressureLogFilename, sampleTimesFilename )
%LOADPRESSURELOG
%
% pressureLog                = LOADPRESSURELOG( pressureLogFilename ) ;
% [pressureLog, sampleTimes] = LOADPRESSURELOG( pressureLogFilename, sampleTimesFilename )

if nargin < 1
    error( 'Insufficient arguments. Must provide full path to pressure log .bin file.' ) ;

else
    if nargin >= 1
        pressureLogFid = fopen( pressureLogFilename, 'r' ) ;
        pressureLog    = fread( pressureLogFid, inf, 'double' ) ;
        fclose( pressureLogFid );
    end

    if nargin == 2 
        sampleTimesFid = fopen( sampleTimesFilename, 'r' ) ;
        sampleTimes    = fread( sampleTimesFid, inf, 'double' ) ;
        fclose( sampleTimesFid );
    end
end

end
% =========================================================================
function [] = plotpressurelog( pressureLog, Params )
%PLOTPRESSURELOG
%
% PLOTPRESSURELOG( pressureLog ) ;
% PLOTPRESSURELOG( pressureLog, Params )

DEFAULT_FIGURETITLE = 'Pressure log' ;

if nargin < 1
    error( 'Insufficient arguments. Must provide pressure log vector.' ) ;
end

if nargin == 1 || isempty( Params ) 
    Params.dummy = [] ;
end

if ~myisfield( Params, 'figureTitle' ) || isempty( Params.figureTitle ) 
    Params.figureTitle = DEFAULT_FIGURETITLE ;
end


% ------- 
figure 

if myisfield( Params, 'sampleTimes' ) && ~isempty( Params.sampleTimes ) 
    plot( Params.sampleTimes, pressureLog, '+' ) ;
    xlabel('Time (s)');
else
    plot( pressureLog, '+' ) ;
    xlabel('Sample index');
end
    
title( Params.figureTitle ) ;
ylabel('Pressure (Pa)');

end
% =========================================================================
function [medianPressure] = userselectmedianpressure( pressureLog )
% USERSELECTMEDIANPRESSURE - Set respiratory probe specifications
%
%   medianPressure = USERSELECTMEDIANPRESSURE( pressureLog ) 
%
%   Plots pressureLog and the user selects START and END (apnea) indices
%   over which to calculate the median. The median pressure is superposed
%   over the pressureLog graph and the user is asked if the result is 
%   satisfactory (or redo).

isUserSatisfied = false ;

while ~isUserSatisfied

    gcf ;
    plot( pressureLog, '+' ) ;
    title( 'Pressure Log' ) ;
    
    xlabel('Sample index');
    ylabel('Pressure (0.01 mbar)');
    
    apneaStartIndex = ...
        input( ['Identify sample index corresponding to beginning of apnea ' ...
            '([Enter] selects sample 1): '] ) ;
    
    if isempty(apneaStartIndex)
        apneaStartIndex = 1;
    end

    apneaEndIndex = ...
        input( ['Identify sample index corresponding to end of apnea ' ...
            '([Enter] selects the last recorded sample): '] ) ;

    if isempty(apneaEndIndex)
       medianPressure = ...
           median( pressureLog( apneaStartIndex : end ) ) ;
    else
       medianPressure = ...
           median( pressureLog( apneaStartIndex : apneaEndIndex ) ) ;
    end

    gcf; 
    plot( pressureLog, '+' );
    hold on;
    plot( medianPressure*ones( size( pressureLog ) ) ) ;
    title( 'Pressure Log' ) ;
    xlabel('Sample index');
    ylabel('Pressure (0.01 mbar)');
    legend('Pressure log','Median pressure over given interval');    
    hold off;

    response = input(['Is the current median estimate satisfactory?' ...
        '0 to re-enter the data range; 1 (or enter) to continue: ']) ;

     if ~isempty(response)
        isUserSatisfied = logical(response) ;
     else
         isUserSatisfied = true;
     end

end

end
% =========================================================================
end

% =========================================================================
% =========================================================================

end
