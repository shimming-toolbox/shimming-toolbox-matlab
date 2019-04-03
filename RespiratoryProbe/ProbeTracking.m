classdef ProbeTracking < matlab.mixin.SetGet 
% PROBETRACKING - Respiratory probe for real-time shimming 
%
% Aux = PROBETRACKING(  )
%
%   Aux contains fields
%           
%       .Source    
%
%       .Data 
%
%       .Specs
%
% .......
%
%   Description
%
% =========================================================================
%
% =========================================================================
% Updated::20190402::ryan.topfer@polymtl.ca
% =========================================================================

properties   
    Source ;
    Data ;
    Specs ; % state = {active, inactive, inert, void}
end

% =========================================================================
% =========================================================================
methods
% =========================================================================
function Aux = ProbeTracking( varargin )
%PROBE - ProbeTracking  

Aux.Data.p    = []; % may be filtered & limited
Aux.Data.pRaw = []; % raw measurement

if nargin < 1 || isempty( varargin{1} )
    Specs = [] ;

elseif isstruct( varargin{1} )
    Specs = varargin{1} ;

elseif ischar( varargin{1} )
% Non-urgent TODO:
%   This form of ProbeTracking() initialization/construction is not to
%   be called by a user, but rather, from the ProbeTracking() constructor
%   itself, making it better suited as a 'private' constructor. Apparently
%   Matlab does not permit this. There is a work-around described here:
%   https://stackoverflow.com/questions/29671482/private-constructor-in-matlab-oop
    filename = varargin{1} ;
    load( filename ) ;
    Aux.beginrecordingdaemon() ; % runs continuously in background
    return;
end

if myisfield( Specs, 'state' ) && strcmp( Specs.state, 'inert' )
    Aux.Source = [] ;
    Aux.Specs   = Specs ;
else
    [ Aux.Source, Aux.Specs ] = ProbeTracking.declareprobe( Specs ) ;
    
    if ~isempty( Aux.Source )
        Aux.createbufferfile() ;
        Aux.createrecordingdaemon() ;
    end
end

end    
% =========================================================================
function [AuxCopy] = copy( Aux )
%COPY  
% 
% Aux = COPY( Aux )

Specs       = Aux.Specs ;
Specs.state = 'inert' ;

AuxCopy = ProbeTracking( Specs ) ; % or just allow that constructor to accept ProbeTracking as an input and copy the fields...

AuxCopy.Data = Aux.Data ;

end
% =========================================================================
function [] = delete( Aux )
%DELETE  
% DELETE( Aux )
% 
% Destructor. Calls Aux.deletecomport( ) 

if myisfield( Aux, 'Source' ) && ~isempty( Aux.Source ) 
    Aux.deletecomport();
end

clear Aux ;

end
% =========================================================================
function [Aux] = deletecomport( Aux )
%DELETECOMPORT  
% 
% DELETECOMPORT( Aux )
% 
% Proper way to close + delete the serial port object

if myisfield( Aux, 'Source' ) && ~isempty( Aux.Source ) 
    fclose( Aux.Source ) ;
    delete( Aux.Source ) ;
    clear Aux.Source  ;
else
    disp('Failed to delete .Source - it does not exist.');
end

end
% =========================================================================
function [isTracking] = begintracking( Aux )
%BEGINTRACKING - Initialize & open (RS-232) communication port 
%
% [isTracking] = BEGINTRACKING( Aux )
%
% Opens Aux.Source 
%
% Returns TRUE if successful

fopen(Aux.Source);
disp('Connecting to respiratory probe...')

iAttempt    = 1 ;
maxAttempts = 3 ; 

isTracking  = false ;

while( ~isTracking && iAttempt <= maxAttempts )

    disp(['Attempt #' num2str(iAttempt)]);    
   
    firstWord = fscanf( Aux.Source, '%u') ;
    
    if( ~isempty(firstWord) && isnumeric(firstWord) )  
        isTracking = true ;
    end
    
    iAttempt = iAttempt + 1;

end

if isTracking
    disp('Communication successful. Reading in from serial port...')
else
    Aux.stoptracking() ;
    error('Communication to respiratory probe failed. Closing Aux.Source.')
end

end
% =========================================================================
function [] = calibratelimiting( Aux, Params )
%CALIBRATELIMITING 
%
% [] = CALIBRATELIMITING( Aux )
%
% Record 1 min of signal to determine the theshold levels beyond which limiting
% will be applied. 
%
% 5 standard deviations either above or below the mean determines levels.
%
% Limits are saved in Aux.Specs.clipLimits

% reset limits:
Aux.Specs.clipLimits = [-Inf Inf] ;

signal = Aux.recordandplotphysiosignal( Params ) ;

Aux.Specs.clipLimits = [ ( mean(signal) - 5*std(signal) ), ...
                         ( mean(signal) + 5*std(signal) ) ] ;

end
% =========================================================================
function [] = stoptracking( Aux )
%STOPTRACKING - close (RS-232) Communication Port 
%
% [] = STOPTRACKING( Aux )
%
% fclose( Aux.Source )

fclose( Aux.Source ) ;

end
% =========================================================================
function [] = clearrecording( Aux )
%CLEARRECORDING  
%
% [] = CLEARRECORDING( Aux )
%
% Empties Aux.Data.t, Aux.Data.p, and  Aux.Data.pRaw  

Aux.Data.t    = [] ;
Aux.Data.p    = [] ;
Aux.Data.pRaw = [] ;

end
% =========================================================================
function [p,t] = getupdate( Aux )
%GETUPDATE 
%
% [p,t] = GETUPDATE( Aux )
%
% Reads a single (16-bit) measurement (p) and returns p as the typecasted double.
% t is the sample time in units of milliseconds.
%
% p is either read from the open com port, or from the temp file buffer.

t = 0;

switch Aux.Specs.readSource
    case 'Source' 
        
        assert( strcmp( Aux.Source.Status, 'open' ), 'Error: Serial port is closed.' );

        tmp = fscanf( Aux.Source, '%u', [1 1] ) ;
        Aux.Data.pRaw(end+1) = tmp(end) ;
        % %% No of measurements before the actual polyfit will begin
        % smooth_window = 5; % place in declare probe Params
        % timepoint     = length( Aux.Data.pRaw );
        %
        % time_data(timepoint) = timepoint * 0.1;
        % p = polyfit(time_data(smooth_window:timepoint),data_cprobe(smooth_window:timepoint),4);
        % y = polyval(p,time_data(smooth_window:timepoint));

        if ( Aux.Data.pRaw(end) > Aux.Specs.clipLimits(1) ) ...
                && ( Aux.Data.pRaw(end) < Aux.Specs.clipLimits(2) )

            Aux.Data.p(end+1) = Aux.Data.pRaw(end) ;
        else
            % replace with most recent unclipped/undistorted sample:
            Aux.Data.p(end+1) = Aux.Data.p(end) ;
        end

        p = Aux.Data.p(end) ;
    
    case 'fileBuffer'
        [p, t] = Aux.getupdatefromfilebuffer( ) ;
end

end  
% =========================================================================
function [p,t] = getupdatefromfilebuffer( Aux )
%GETUPDATEFROMFILEBUFFER

[p, t] = Aux.readfromfilebuffer( ) ;

Aux.Data.p(end+1) = p ;
Aux.Data.t(end+1) = t ;

end
% =========================================================================
function [] = killrecordingdaemon( Aux )
%KILLRECORDINGDAEMON
%
% Sends STOP byte to file buffer

Aux.Specs.Buffer.Data(1) = 0 ; 

end
% =========================================================================
function [physioSignal, sampleTimes] = recordandplotphysiosignal( Aux1, Params, Aux2 ) 
%RECORDANDPLOTPHYSIO    
%
%   Description
%   
%   Reads physiological data from microcontroller 
%
%   Syntax
%
%   [measurements, sampleTimes] = RECORDANDPLOTPHYSIO( Aux1, Parameters )
%   [measurements, sampleTimes] = RECORDANDPLOTPHYSIO( Aux1, Parameters, Aux2 )
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
%   .physioSignalFilename
%       default = ['./' datestr(now,30) '-physioSignal.bin' ] ; 
%
%   .runTime 
%       Total sampling time in seconds.
%       default = 60
%
%   .isPlottingInRealTime
%       [default : true ]
%
%   .refreshRate  
%       Rate at which the real-time display refreshes. 
%       Problems may arise if this is too fast!
%       [default : 4 Hz ]

assert( strcmp( Aux1.Source.Status, 'closed' ), 'Error: Serial port is open/in use.' );

if nargin == 3
    assert( strcmp( Aux2.Source.Status, 'closed' ), 'Error: Serial port is open/in use.' );
    isDualTracking = true ;
else
    isDualTracking = false ;
end

DEFAULT_ISSAVINGDATA          = true ;
DEFAULT_ISFORCINGOVERWRITE    = false ;
DEFAULT_PHYSIOSIGNALFILENAME  = ['./' datestr(now,30) '-physioSignal.bin' ] ;

if  nargin < 2 || isempty(Params)
    Params.dummy = [] ;
end
 
if  ~myisfield( Params, 'isSavingData' ) || isempty(Params.isSavingData)
    Params.isSavingData = DEFAULT_ISSAVINGDATA ;
end

if  ~myisfield( Params, 'isForcingOverwrite' ) || isempty(Params.isForcingOverwrite)
    Params.isForcingOverwrite = DEFAULT_ISFORCINGOVERWRITE ;
end

if  ~myisfield( Params, 'physioSignalFilename' ) || isempty(Params.physioSignalFilename)
    Params.physioSignalFilename = DEFAULT_PHYSIOSIGNALFILENAME ;
end

% ------- 
isUserSatisfied = false ;

while ~isUserSatisfied
   
   if ~isDualTracking 
        Aux1.clearrecording() ;
        Aux1.recordphysiosignal( Params ) ;
    else
        Aux1.clearrecording() ;
        Aux2.clearrecording() ;
        Aux1.recordphysiosignal( Params, Aux2 ) ;
    end
    
    physioSignal = Aux1.Data.p ;    
    sampleTimes = (Aux1.Specs.dt/1000)*[0:(numel(Aux1.Data.p)-1)] ;
   
    % ------- 
    fprintf(['\n ----- \n Plotting physio recording for user verification... \n']) ;
    figure ;
    if ~isDualTracking
        plot( sampleTimes, Aux1.Data.p, '+' ) ;

        if Params.isSavingData
            title( Params.physioSignalFilename ) ;
        else
            title( 'Physio signal' )
        end

        xlabel('Time (s)');
        ylabel('Amplitude (AU)');
    else
        subplot(211)
        plot( sampleTimes, Aux1.Data.p, '+' ) ;

        if Params.isSavingData
            title( Params.physioSignalFilename ) ;
        else
            title( 'Physio signal 1' )
        end

        xlabel('Time (s)');
        ylabel('Amplitude (AU)');
        
        subplot(212)
        plot( sampleTimes, Aux2.Data.p, '+' ) ;

        if Params.isSavingData
            title( [Params.physioSignalFilename '-2'] ) ;
        else
            title( 'Physio signal 2' )
        end

        xlabel('Time (s)');
        ylabel('Amplitude (AU)');
    end
    
    response = input(['\n Is the current physio recording satisfactory? ' ...
        'Enter 0 to rerecord; 1 to continue. \n']) ;

    isUserSatisfied = logical(response) ;

end

% ------- 
if Params.isSavingData

    Aux1.saverecording( Params.physioSignalFilename )
    
    if isDualTracking
        Aux2.saverecording( [ Params.physioSignalFilename '-2' ] )
    end

end

end
% =========================================================================
function [] = recordphysiosignal( Aux1, Params, Aux2 )
%RECORDPHYSIOSIGNAL  
%
% Continuously tracks respiratory probe.
%
% [] = TRACKPROBE( Aux, Params )

DEFAULT_RUNTIME              = 60 ; % [units : s]
DEFAULT_ISPLOTTINGINREALTIME = true ;
DEFAULT_REFRESHRATE          = 4 ; % [units : Hz]

if  nargin < 2 || isempty(Params)
    Params.dummy = [] ;
end

if  ~myisfield( Params, 'runTime' ) || isempty(Params.runTime)
    Params.runTime = DEFAULT_RUNTIME ;
end

if  ~myisfield( Params, 'isPlottingInRealTime' ) || isempty(Params.isPlottingInRealTime)
    Params.isPlottingInRealTime = DEFAULT_ISPLOTTINGINREALTIME ;
end

if  ~myisfield( Params, 'refreshRate' ) || isempty(Params.refreshRate)
    Params.refreshRate = DEFAULT_REFRESHRATE ;
end

% ------- 
StopButton     = stoploop({'Stop recording'}) ;
Aux1.clearrecording() ;

sampleIndices  = [] ;

nSamples = Params.runTime / (Aux1.Specs.dt/1000) ;
iSample  = 0 ; 
sampleTimes = [] ;
% t = 0 ; % sample time [units: ms]

if Params.isPlottingInRealTime

    % ------- 
    % figure
    figureHandle = figure('NumberTitle','off',...
        'Name','Physio signal',...
        'Color',[0 0 0],'Visible','off');
        
    % Set axes
    axesHandle = axes('Parent',figureHandle,...
        'YGrid','on',...
        'YColor',[0.9725 0.9725 0.9725],...
        'XGrid','on',...
        'XColor',[0.9725 0.9725 0.9725],...
        'Color',[0 0 0]);

    hold on;
    
        plotHandle = plot(axesHandle,0,0,'Marker','.','LineWidth',1,'Color',[1 0 0]);
            
        xlim(axesHandle,[0 nSamples]);
            
        title('Respiration Aux','FontSize',15,'Color',[1 1 0]);
        ylabel('Amplitude','FontWeight','bold','FontSize',14,'Color',[1 1 0]);
        xlabel('Time [ms]','FontWeight','bold','FontSize',14,'Color',[1 1 0]);
        % xlabel('Sample index','FontWeight','bold','FontSize',14,'Color',[1 1 0]);
    
        drawnow limitrate; 
    set(figureHandle,'Visible','on');

end

display('Reading probe measurements...');

% for real-time plotting, updating @the same rate as the samples
% (e.g. 100 Hz) poses a problem (computer can't seem to keep up with the incoming samples).
% solution is to update the display ~every so often~ (e.g. 4x per second seems OK)
nSamplesBetweenRefresh = (1/Params.refreshRate)/(Aux1.Specs.dt/1000) ;

while ( iSample < nSamples ) && ~StopButton.Stop()

    iSamplesBetweenRefresh = 0;

    for iSamplesBetweenRefresh = 1 : nSamplesBetweenRefresh 
        
        iSample = iSample + 1 ;
        [p,t] = Aux1.getupdate() ;
        sampleIndices(end+1) = iSample ;
        sampleTimes(end+1) = t ;

    end

    if Params.isPlottingInRealTime
            % set(plotHandle,'YData',Aux1.Data.p,'XData',sampleTimes);
            set(plotHandle,'YData',Aux1.Data.p,'XData',sampleIndices);
    end

end

Aux1.stoptracking();

StopButton.Clear() ;

end
% =========================================================================
function [] = saverecording( Aux, logFilename )
%SAVERECORDING
%
%   SAVERECORDING( Aux )
%   SAVERECORDING( Aux, logFilename )

if  nargin < 2 || ~ischar( logFilename ) 
    logFilename = ['./' datestr( now, 30 ) '-physioSignal.txt' ] ;
end

fid = fopen( logFilename, 'w') ;
fprintf( fid, '%d  %d\n', [uint32(Aux.Data.t); int16(Aux.Data.p)] ) ;
fclose(fid) ;

end
% =========================================================================

% =========================================================================
% =========================================================================
end

% =========================================================================
% =========================================================================

% =========================================================================
% =========================================================================
methods( Access =  private)
% =========================================================================
function [] = createbufferfile( Aux )
%CREATEBUFFERFILE
% 
% [] = CREATEBUFFERFILE( Aux )
% 
% Creates buffer file to which a daemon process writes, while a second,
% user-controlled process, reads. 
%
% The file is saved in the Matlab-defined tempdir() under
% 'capacitive_probe.buffer.dat' or 'pressure_probe_buffer.dat' depending on the
% case.
%
% The file is memory-mapped using memmapfile() with the return argument
% assigned to Aux.Specs.Buffer
%
% BYTE SIGNIFICANCE:
% The file is 24 bytes (char) in length and structured as following:
%
% [1] : STOP RECORDING FLAG 
% [2] : PAUSE READ/WRITE FLAG 
% [3] : nBytes sample time
% [4] : nBytes sample value
% [5:(4+[3])] : sample time [units: ms]
% [(5+[3]):4+[3]+[4])] : sample value [A.U.]

% Create or overwrite the file: 
filename = fullfile( tempdir, [ Aux.Specs.probeType '_probe_buffer.dat'] ) ;
[f, msg] = fopen(filename, 'w') ;

if f ~= -1
    fwrite(f, [1 1 zeros(1,22)], 'uint8') ; 
    fclose(f) ;
else
    error('MATLAB:demo:send:cannotOpenFile', ...
          'Cannot open file "%s": %s.', filename, msg) ;
end
 
% Memory map the file.
Aux.Specs.Buffer = memmapfile(filename, 'Writable', true, 'Format', 'uint8') ;

end
% =========================================================================
function [p,t] = readfromfilebuffer( Aux )
%READFROMFILEBUFFER
%
% BYTE SIGNIFICANCE:
%
% [1] : STOP RECORDING FLAG 
% [2] : PAUSE READ/WRITE FLAG 
% [3] : nBytes sample time
% [4] : nBytes sample value
% [5:(4+[3])] : sample time [units: ms]
% [(5+[3]):4+[3]+[4])] : sample value [A.U.]

% time-out condition occurs if no new measurement is registered after MAX_DWELL
DWELL_TIME  = 0.01 ; % [units: s]
MAX_DWELL   = 5 ;  
totalDwell  = 0 ;

% Wait until bytes [3] and [4] are non-zero (indicating a new measurement)
while (~Aux.Specs.Buffer.Data(3)) & (~Aux.Specs.Buffer.Data(4)) & (totalDwell < MAX_DWELL)
    pause( DWELL_TIME ) ;
    totalDwell = totalDwell + DWELL_TIME ; 
end

% [2] -> 0 : freeze values from overwrite
Aux.Specs.Buffer.Data(2) = 0 ;

if ( totalDwell >= MAX_DWELL )
    warning('Time-out occured before a new value was read. Returning previous measurement.')
    p = Aux.Data.p(end) ;
    t = Aux.Data.t(end) ;
    return;
end

% read sample time [units: ms]
t = str2double( char(Aux.Specs.Buffer.Data( 5:(4+double(Aux.Specs.Buffer.Data(3))) ))' ) ;

% byte index of the sample value most-significant digit 
iP1 = 5 + double(Aux.Specs.Buffer.Data(3)) ;

% read sample value [units: A.U.]
p = str2double( char( Aux.Specs.Buffer.Data( iP1:(iP1-1+double(Aux.Specs.Buffer.Data(4)) )))') ;

% indicates to next call to read() that the latest measurement has already been read
Aux.Specs.Buffer.Data(3) = 0 ;
Aux.Specs.Buffer.Data(4) = 0 ;

% [2] -> 1 : unfreeze and continue recording 
Aux.Specs.Buffer.Data(2) = 1 ;
 
end
% =========================================================================
function [] = writetofilebuffer( Aux, p, t )
%WRITETOFILEBUFFER
%
% BYTE SIGNIFICANCE:
%
% [1] : STOP RECORDING FLAG 
% [2] : PAUSE READ/WRITE FLAG 
% [3] : nBytes sample time
% [4] : nBytes sample value
% [5:(4+[3])] : sample time [units: ms]
% [(5+[3]):4+[3]+[4])] : sample value [A.U.]

DWELL_TIME  = 0.01 ; % [units: s]

% Wait until BYTE[2] == 1 (presence of 0 indicates file is being read)
while Aux.Specs.Buffer.Data(2) ~= 1
    pause( DWELL_TIME ) ;
end

t_str   = num2str( t ) ;
t_nChar = length( t_str ) ;

p_str   = num2str( p ) ;
p_nChar = length( p_str ) ;

% ------
% write sample time [units: ms]
Aux.Specs.Buffer.Data(5:(4+t_nChar)) = t_str ;

% ------
% write sample value [units: AU] 

% byte index of the sample value most-significant digit 
iP1 = 5 + t_nChar ;
Aux.Specs.Buffer.Data(iP1:(iP1-1+p_nChar)) = p_str ;

% 
Aux.Specs.Buffer.Data(3) = t_nChar ; 
Aux.Specs.Buffer.Data(4) = p_nChar ; 
   
end
% =========================================================================
function [] = createrecordingdaemon( Aux )
%CREATERECORDINGDAEMON
% 
% CREATERECORDINDAEMON( Aux )
%
% Save the instantiated Aux object, and launch a background (daemon) Matlab session
% to load the object and begin background recording.

% The daemon session reads directly from the USB (Com) port while the user
% session reads from a file buffer.
Aux.Specs.readSource = 'Source' ;

pathToAuxObject = [ tempdir 'Aux' ] ;
save( pathToAuxObject, 'Aux' ) ;

tmpDir = pwd ;
cd( '~/' ) ; % change to home folder so Matlab runs startup.m (defines path)
cmd = sprintf( '%s', 'matlab -r " ProbeTracking( ''', pathToAuxObject, ''' );" &') ;
unix( cmd ) ;
pause(5) ;
cd( tmpDir ) ;

Aux.Specs.readSource = 'fileBuffer' ;

end
% =========================================================================
function [] = beginrecordingdaemon( Aux )
%BEGINRECORDINGDAEMON
% 
% BEGINRECORDINDAEMON( Aux )

% if  nargin < 2 || ~ischar( logFilename ) 
%     logFilename = ['./' datestr( now, 30 ) '-physioSignal.txt' ] ;
% end

Aux.clearrecording() ;

t = 0 ; % sample time [units: ms]

% sample index
iSample  = 0 ; 

isTracking = Aux.begintracking() ;

assert( isTracking, 'Could not begin tracking. Check device connection.')

fprintf('\n\n Recording from probe...\n\n');

tic;
while Aux.Specs.Buffer.Data(1) 
    
    p = Aux.getupdate() ;
    Aux.Data.t(end+1) = 1000*toc;
    iSample = iSample + 1 ;
    % t = iSample*Aux.Specs.dt ; % sample time [units: ms]
    Aux.writetofilebuffer( p, Aux.Data.t(end) ) ;
    
end

Aux.stoptracking() ;

Aux.saverecording( ) ;

% ------
% close serial port, delete object, and quit matlab session 
Aux.delete() ;
quit ;

end
% =========================================================================

end

% =========================================================================
% =========================================================================
methods(Static)
% =========================================================================
function [Source, AuxSpecs] = declareprobe( AuxSpecs )
%DECLAREPROBE Declares serial object for probe 
% 
% [Source,AuxSpecs] = declareprobe( AuxSpecs )
%

% AuxSpecs can have the following fields 
%  
% NOTE 
%   these are not really optional/configurable... TODO: remove deprecated/false 'options'
%
% .dt
%   Interval between samples [units: ms]
%   Should be a positive multiple of the minimum period of 10 ms.
% 
% .baudRate
%   default : 115200
%
% .portName 
%   Address of the probe-associated serial port within file system
%   default: 
%       if ismac 
%           portName = '/dev/tty.usbmodem*'
%       elseif isunix
%           portName = '/dev/ttyS100'

DEFAULT_ARDUINOPERIOD = 50 ; % [units: ms] 
DEFAULT_TEENSYPERIOD  = 100 ; % [units: ms] 

DEFAULT_BAUDRATE      = 115200 ;
DEFAULT_CLIPLIMITS    = [-Inf Inf] ; % [units: probe signal] 
DEFAULT_PROBETYPE     = 'pressure' ; 

if nargin < 1 || isempty(AuxSpecs)
    ShimUse.customdisplay('Default parameters will be used:')
    AuxSpecs.dummy = [] ;
end

if ~myisfield( AuxSpecs, 'baudRate' ) || isempty(AuxSpecs.baudRate) ...
    AuxSpecs.baudRate = DEFAULT_BAUDRATE ;
end

if ~myisfield( AuxSpecs, 'clipLimits' ) || isempty(AuxSpecs.clipLimits) ...
    AuxSpecs.clipLimits = DEFAULT_CLIPLIMITS ;
end

if ~myisfield( AuxSpecs, 'probeType' ) || isempty(AuxSpecs.probeType) ...
    AuxSpecs.probeType = DEFAULT_PROBETYPE ;
end

% -------------------------------------------------------------------------
% check and/or assign serial port
Source = [] ; 
isAssigningSerialPort = false ;


if  myisfield( AuxSpecs, 'portName' ) && ~isempty(AuxSpecs.portName)
   
    [fileDir,portname,fileExtension] = fileparts( AuxSpecs.portName ) ;
    
    listOfDevices = dir( fileDir ) ;

    isPortDetected = false ;
    
    for iFile = 1 : length(listOfDevices)-1
        if strcmp( [fileDir '/' listOfDevices(iFile).name fileExtension], AuxSpecs.portName )
            isPortDetected = true ;
            isAssigningSerialPort = true ;
            if strcmp( AuxSpecs.portName, '/dev/tty.usbmodem4471890' )
                AuxSpecs.probeType = 'capacitive' ;
            end
       end
    end

    if ~isPortDetected    
        disp(['Warning: Given port name (' AuxSpecs.portName ' ) not found. ' ...
            'Checking default port names.']) ;
        AuxSpecs.portName = [] ;
    end
end

if  ~myisfield( AuxSpecs, 'portName' ) || isempty(AuxSpecs.portName)

    if ismac
        % For some reason this doesn't work...
        % -------
        % if exist( '/dev/tty.usbmodem1411', 'file' ) ; % my macbook's left USB port
        %   AuxSpecs.portName = '/dev/tty.usbmodem1411';
        % elseif exist( '/dev/tty.usbmodem1421', 'file' ) ;
        %     AuxSpecs.portName = '/dev/tty.usbmodem1421' ; % my macbook's right USB port
        % else
        %     error(['Device file not found. ' ...
        %         'Is the microcontroller connected?']) ; 
        % ------
        listOfDevices = dir( '/dev/tty.usbmodem*' ) ;

        if isempty(listOfDevices) 
            warning( 'Device file not found. Check USB device is connected.' ) ;
        elseif length(listOfDevices) == 1
            
            isAssigningSerialPort = true ;
            AuxSpecs.portName     = ['/dev/' listOfDevices(1).name ] ;
            
            switch listOfDevices(1).name
            
                case 'tty.usbmodem4471890'
                    AuxSpecs.probeType = 'capacitive' ;
                otherwise
                    AuxSpecs.probeType = 'pressure' ;
            end
        else
            warning('Ambiguous device identifier. Consider entering portName as argument. See HELP.');
            
            for iD = 1 : length( listOfDevices ) 
                if strcmp( listOfDevices(iD).name, 'tty.usbmodem4471890' ) ; % teensy 3.5 on RT macbook pro
                    isAssigningSerialPort = true ;
                    AuxSpecs.probeType    = 'capacitive' ;
                    AuxSpecs.portName     = ['/dev/' listOfDevices(iD).name ] ;
                end
            end
        end

    elseif isunix
        % if exist('/dev/ttyS100', 'file')  
            warning('ProbeTracking.declareprobe() assumes a serial port address/name that may be invalid!')
            AuxSpecs.portName = '/dev/ttyS100' ;
        % % else
        %     errorMsg = ['Device file ( ' portName ' ) not found.'  ...
        %     'Is the microcontroller connected?' ...
        %     'Ensure symbolic link exists between actual Arduino device' ...
        %     'file and the phantom-device file discoverable by MATLAB:' ...
        %     'In a terminal, type: ' ... 
        %     'sudo ln -f /dev/ttyACM0 /dev/ttyS100' ] ; 
        %     error( errorMsg ) ;
        % end

    else
        error( 'OS not supported' ) ;
    end
end

if isAssigningSerialPort
    Source = serial( AuxSpecs.portName, 'BaudRate',  AuxSpecs.baudRate ) ;
else
    Source = [] ;
    AuxSpecs.state = 'inert' ;
end

switch AuxSpecs.probeType
    case 'pressure'
        ShimUse.customdisplay( ['\n----- Pressure probe -----\n'] );
        AuxSpecs.dt = DEFAULT_ARDUINOPERIOD ; % [units: ms]
    case 'capacitive'
        ShimUse.customdisplay( ['\n----- Capacitive probe -----\n'] );
        AuxSpecs.dt = DEFAULT_TEENSYPERIOD ; 
end

ShimUse.customdisplay( [ 'Sampling frequency = ' num2str(1000/AuxSpecs.dt) ' Hz'] )

end
% =========================================================================
function [measurementLog, sampleTimes] = loadmeasurementlog( measurementLogFilename, sampleTimesFilename )
%LOADMEASUREMENTLOG
% 
% Reads binary file of data measurements (e.g. pressure recording) to return
% vector(s) of doubles.
%
% measurementLog                = LOADMEASUREMENTLOG( measurementLogFilename ) ;
% [measurementLog, sampleTimes] = LOADMEASUREMENTLOG( measurementLogFilename, sampleTimesFilename )

if nargin < 1
    error( 'Insufficient arguments. Must provide full path to measurement log .bin file.' ) ;

else
    if nargin >= 1
        [~,~,ext] = fileparts( measurementLogFilename ) ;

        if strcmp( ext, '.txt' )
            X = load( measurementLogFilename ) ;
            sampleTimes    = X(:,1) ;
            measurementLog = X(:,2) ;
            return;
        else % saved as binary
            measurementLogFid = fopen( measurementLogFilename, 'r' ) ;
            measurementLog    = fread( measurementLogFid, inf, 'double' ) ;
            fclose( measurementLogFid );
        end
    end
    
    if nargin == 2 
        sampleTimesFid = fopen( sampleTimesFilename, 'r' ) ;
        sampleTimes    = fread( sampleTimesFid, inf, 'double' ) ;
        fclose( sampleTimesFid );
    end
end

end
% =========================================================================
function [] = plotmeasurementlog( measurementLog, Params )
%PLOTMEASUREMENTLOG
%
% PLOTMEASUREMENTLOG( measurementLog ) ;
% PLOTMEASUREMENTLOG( measurementLog, Params )
%
% Supported fields to Params struct
%
%   .figureTitle
%       [default: 'Pressure log']
%
%   .sampleTimes
%       vector (length == length(measurementLog)) of sample times in seconds
%
%   .yLabel
%       [default: 'Pressure (kPa)']

DEFAULT_FIGURETITLE = 'Respiration' ;
DEFAULT_YLABEL      = 'Amplitude (AU)' ;

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

% ------- 
figure 

if myisfield( Params, 'sampleTimes' ) && ~isempty( Params.sampleTimes ) 
    plot( Params.sampleTimes, measurementLog, '+' ) ;
    xlabel('Time (s)');
else
    plot( measurementLog, '+' ) ;
    xlabel('Sample index');
end
    
title( Params.figureTitle ) ;
ylabel( Params.yLabel ) ;

end
% =========================================================================
function [iFlattest] = findflattest( measurementLog, nSamples )
%FINDFLATTEST 
%
% iFlattest = FINDFLATTEST( measurementLog, nSamples ) 
% 
% Calculates measurementLog variance over sliding window (nSamples long) and
% returns index (iFlattest) corresponding to start of the most constant segment
% (e.g. a breath-hold).

assert( nSamples > 0 )
assert( nSamples <= length(measurementLog) )

nVariances = length(measurementLog) - nSamples ;
variances  = zeros( nVariances, 1 );

for iFlattest =  1 : nVariances
   variances(iFlattest) = var( measurementLog( iFlattest:(iFlattest+nSamples) ) ) ;
end

[~, iFlattest] = min( variances ) ;

end
% =========================================================================
function [medianMeasure] = selectmedianmeasurement( measurementLog, nSamplesApnea, isUserSelectionEnabled )
% SELECTMEDIANMEASUREMENT
%
%   medianMeasure = SELECTMEDIANMEASUREMENT( measurementLog ) 
%   medianMeasure = SELECTMEDIANMEASUREMENT( measurementLog, nSamplesApnea ) 
%   medianMeasure = SELECTMEDIANMEASUREMENT( measurementLog, nSamplesApnea, isUserSelectionEnabled ) 
%
%   Plots measurementLog and the user selects START and END (apnea) indices
%   over which to calculate the median. The median measurement is superposed
%   over the measurementLog graph and the user is asked if the result is 
%   satisfactory (or redo).

DEFAULT_ISUSERSELECTIONENABLED = true ;

if ( nargin == 1 ) 
    nSamplesApnea = [] ;

elseif ~isempty( nSamplesApnea ) 
    assert( nSamplesApnea > 0 ) ;
    
end

if nargin < 3
    isUserSelectionEnabled = DEFAULT_ISUSERSELECTIONENABLED ;
end


isUserSatisfied = false ;

while ~isUserSatisfied

    gcf ; 
    clf ;
    plot( measurementLog, '+' ) ;
    title( 'Measure Log' ) ;
    
    xlabel('Sample index');
    ylabel('Amplitude');

    if ~isempty( nSamplesApnea )

        if nSamplesApnea < length( measurementLog )
            % Auto-selection of start/end breath-hold indices
            trainingFrameStartIndex = ProbeTracking.findflattest( measurementLog, nSamplesApnea ) ;
            trainingFrameEndIndex   = trainingFrameStartIndex + nSamplesApnea ;
        else
            trainingFrameStartIndex = 1 ;
            trainingFrameEndIndex   = length( measurementLog ) ;
        end
        
        medianMeasure = ...
           median( measurementLog( trainingFrameStartIndex : trainingFrameEndIndex ) ) ;

        gcf; 
        plot( measurementLog, '+' );
        hold on;
        plot( [trainingFrameStartIndex : trainingFrameEndIndex-1], medianMeasure*ones( 1, nSamplesApnea ), 'LineWidth',3 ) ;
        title( 'Measure Log' ) ;
        xlabel('Sample index');
        ylabel('Amplitude');
        legend('Measure log',['Median over interval of min. variance']);    
        hold off;
    
    end
    
    if ~isUserSelectionEnabled 
         isUserSatisfied = true ;
         return;
    
     else isUserSelectionEnabled 
        
        response = input(['Is the current median estimate satisfactory? ' ...
            '0 to manually specify the data range; 1 (or enter) to accept & continue: ']) ;

         if ~isempty(response)
            isUserSatisfied = logical(response) ;
         else
             isUserSatisfied = true;
         end
    


        if ~isUserSatisfied
            % Manual selection    
            trainingFrameStartIndex = ...
                input( ['Identify sample index corresponding to beginning of training frame ' ...
                    '([Enter] selects sample 1): '] ) ;
            
            if isempty(trainingFrameStartIndex)
                trainingFrameStartIndex = 1;
            end

            trainingFrameEndIndex = ...
                input( ['Identify sample index corresponding to end of training frame ' ...
                    '([Enter] selects the last recorded sample): '] ) ;

            if isempty(trainingFrameEndIndex)
               medianMeasure = ...
                   median( measurementLog( trainingFrameStartIndex : end ) ) ;
            else
               medianMeasure = ...
                   median( measurementLog( trainingFrameStartIndex : trainingFrameEndIndex ) ) ;
            end

            gcf; 
            plot( measurementLog, '+' );
            hold on;
            plot( medianMeasure*ones( size( measurementLog ) ) ) ;
            title( 'Measure Log' ) ;
            xlabel('Sample index');
            ylabel('Amplitude');
            legend('Measure log','Median measurement over given interval');    
            hold off;

            response = input(['Is the current median estimate satisfactory? ' ...
                '0 to re-enter the data range; 1 (or enter) to continue: ']) ;

             if ~isempty(response)
                isUserSatisfied = logical(response) ;
             else
                 isUserSatisfied = true;
             end
        end
    end
end

end
% =========================================================================1
% =========================================================================

end
% =========================================================================
% =========================================================================

end
