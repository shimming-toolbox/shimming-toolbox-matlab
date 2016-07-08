classdef ProbeTracking 
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
%     ProbeTracking
%     ShimCom
%     ShimOpt
%     ShimSpecs
%     ShimUse
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

Probe.Specs = ProbeTracking.setprobespecs( Specs ) ;

Probe.ComPort = ProbeTracking.initialisecomport( Probe.Specs ) ;

Probe.Data.pressure = [];

end    
% =========================================================================
function [] = delete( Probe )
%DELETE  
% DELETE( Probe )
% 
% Destructor. Calls Probe.deletecomport( ) 

Probe.deletecomport();

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
    disp('ComPort is not open.');
end

end
% =========================================================================
function [p] = readpressure( Probe )
%READPRESSURE 
%
% Reads one (4-byte) pressure measurement (p) in from com port.
%
% p = LOGPRESSURE( Probe )
% ------- 

p = fscanf( Probe.ComPort,'%f') ;
fprintf(['p = ' num2str(p) '\n'])
end  
% =========================================================================
function [Probe] = recordandplotpressurelog( Probe, Params ) 
%RECORDANDPLOTPRESSURELOG    
%
%
%

if  nargin < 2 || isempty(Params)
    Params.dummy = [] ;
end

if  ~myisfield( Params, 'filename' ) || isempty(Params.filename)
    isSavingPressureLog = false ;
    Params.filename = 'Pressure log (Unsaved)' ;
else
    isSavingPressureLog = true ;
end

isUserSatisfied = false ;

msg = input( ['\n Press [Enter] to begin recording pressure log.\n'], 's' );
assert(isempty(msg),'Cancelled calibration.')

while ~isUserSatisfied
    
    Probe.Data.pressure = 0 ;
    Probe = Probe.recordpressurelog( Params ) ;
    
    sampleTimes = (Probe.Specs.arduinoPeriod/1000)*[0:(numel(Probe.Data.pressure)-1)] ;
   
    % ------- 
    fprintf(['\n ----- \n Plotting pressure for user verification... \n']) ;
    figure ;
    plot( sampleTimes, Probe.Data.pressure, '+' ) ;
    title( Params.filename ) ;
    xlabel('Time (s)');
    ylabel('Pressure (0.01 mbar)');
    
    response = input(['\n Is the current pressure log satisfactory? ' ...
        'Enter 0 to rerecord; 1 to continue. \n']) ;

    isUserSatisfied = logical(response) ;

end

if isSavingPressureLog

    pressureLogFid = fopen( Params.filename, 'w+' ) ;
    fwrite( pressureLogFid, Probe.Data.pressure, 'double' ) ;
    fclose( pressureLogFid );

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
% ------- 

DEFAULT_RUNTIME = 30 ; % [units : s]

if  nargin < 2 || isempty(Params)
    Params.dummy = [] ;
end

if  ~myisfield( Params, 'runTime' ) || isempty(Params.runTime)
    Params.runTime = DEFAULT_RUNTIME ;
end
    
% ------- 
StopButton = stoploop({'Stop recording'}) ;

nSamples = Params.runTime / (Probe.Specs.arduinoPeriod/1000) ;
iSample  = 1 ; 

while ( iSample < nSamples ) && ~StopButton.Stop()
    
    Probe.Data.pressure(end+1) = Probe.readpressure() ;
    iSample = iSample + 1 
end

end
% =========================================================================
end

% =========================================================================
% =========================================================================
methods(Static)
% =========================================================================
function [ComPort] = initialisecomport( Specs )
%INITIALISECOMPORT - Initialise (RS-232) Communication Port 
% 

maxCommunicationAttempts = 3; 

ComPort = serial( Specs.portName ) ;
fopen(ComPort);

isSamplingFrequencyReceived = false;

iAttempt = 1 ;
disp('Attempting to assign probe sample frequency...')

while(~isSamplingFrequencyReceived && iAttempt <= maxCommunicationAttempts )

    % Send refresh time 
    disp(['Attempt #' num2str(iAttempt)]);    
    fprintf(ComPort, '%d', Specs.arduinoPeriod);
    pause(0.1)
    firstWord = fscanf(ComPort,'%f') ;

    if(~isempty(firstWord) && isnumeric(firstWord) && firstWord == Specs.arduinoPeriod/10)  
        isSamplingFrequencyReceived = true;
    end

    iAttempt = iAttempt + 1;

end

if isSamplingFrequencyReceived
    disp('Communication successful. Ready to read IO port.')

else
    error('Communication to respiratory probe failed.')
end

end
% =========================================================================
% =========================================================================
function [Specs] = setprobespecs( Specs )
% SETPROBESPECS - Set respiratory probe specifications

DEFAULT_ARDUINOPERIOD    = 10 ; % [units: ms] ;
DEFAULT_SKIPSAMPLEFACTOR = 1;

if nargin < 1 || isempty(Specs)
    disp('Default parameters will be used')
    Specs.dummy = [] ;
end

if  ~myisfield( Specs, 'skipSampleFactor' ) || isempty(Specs.skipSampleFactor)
    Specs.skipSampleFactor = DEFAULT_SKIPSAMPLEFACTOR ;
end

Specs.arduinoPeriod = Specs.skipSampleFactor * DEFAULT_ARDUINOPERIOD ;

% -------------------------------------------------------------------------
% check and/or assign serial port

if  myisfield( Specs, 'portName' ) && ~isempty(Specs.portName)

    if ~exist( Specs.portName, 'file' )    
        disp(['Warning: Given port name (' Specs.portName ' ) not found. ' ...
            '...Attempting alternatives']) ;
        Specs.portName = [] ;
    end
end

if  ~myisfield( Specs, 'portName' ) || isempty(Specs.portName)

    if ismac
        % if exist( '/dev/tty.usbmodem1411', 'file' ) ; % my macbook's left USB port
        %     Specs.portName = '/dev/tty.usbmodem1411'
        % elseif exist( '/dev/tty.usbmodem1421', 'file' ) ;
            Specs.portName = '/dev/tty.usbmodem1421' ; % my macbook's right USB port
        % else
        %     error(['Device file not found. ' ...
        %         'Is the microcontroller connected?']) ; 

    elseif isunix
        if exist('/dev/ttyS100', 'file')  
            portName = '/dev/ttyS100' ;
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
            '([Enter] selects sample 1)'] ) ;
    if isempty(apneaStartIndex)
        apneaStartIndex = 1;
    end

    apneaEndIndex = ...
        input( ['Identify sample index corresponding to end of apnea ' ...
            '([Enter] selects the last recorded sample)'] ) ;

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

    response = input(['Is the current pressure log satisfactory?' ...
        'Enter 0 to rerecord; 1 to continue']) ;

    isUserSatisfied = logical(response) ;

end

end
% =========================================================================
end

% =========================================================================
% =========================================================================

end
