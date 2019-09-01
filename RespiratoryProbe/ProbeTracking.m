classdef ProbeTracking < matlab.mixin.SetGet
% PROBETRACKING - Respiratory probe
% 
% This class deals with respiratory sensor recording.
% 
% .......
%
% Usage
%
%   P = ProbeTracking() ;
%   P = ProbeTracking( Specs ) ;
%
%   When called without arguments, ProbeTracking will attempt to find a
%   connected probe device (in /dev/) based on a set of known device names
%   (operating system and station dependent). 
%
%       [TODO: Add KnownDevices.m ?]
%
%   If a device is discovered, by default and if properly configured (see
%   NOTE [1]), a second ("daemon") MATLAB session is launched in the background
%   to record continuously from the device. 
%   
%   To record and display the signal at any given moment, the user calls 
%   P.recordphysiosignal(), which reads back the live recording from the
%   daemon session through a shared memory-mapping. A graphical STOP button
%   appears which, when pressed, ends this windowed recording and saves it to file
%   without ending the background recording (so P.recordphysiosignal() can be
%   called again). 
%
%   For more options, type: help ProbeTracking.recordphysiosignal
%
%   P = ProbeTracking( Specs ) ;
%
%       Specs.isRecordingDaemonEnabled (bool)
%           true: device I/O (probe recording) and signal processing are 
%             performed in a secondary background/daemon MATLAB session.
%             User session can still monitor the live recording with Probe.recordphysiosignal()
%           false: device I/O and signal processing is done in the same 
%             session. This is useful for debugging. 
%
%   For dual recording (e.g., pressure + capacitive probe), do the
%   following:
%   Plug both probes, then launch the deamon for one probe:
%     P = ProbeTracking();
%   Then, create a new structure with portName of the 2nd probe:
%     AuxSpecs.portName = 'tty.usbmodem4873121' (find the proper port address under /dev after plugging the usb)
%   Create a new object (daemon):
%     C = ProbeTracking(AuxSpecs)
%   To record with both probes:
%     P.recordphysiosignal(C)
% 
% .......
%
%   NOTE [1]: The daemon configuration requires 3 things: 
%
%       1. The file '~/startup.m' must exist and you must add the
%       realtime_shimming_repository directory to the MATLAB path, e.g.
%
%           addpath( genpath( '~/Code/realtime_shimming_repository/' ) ) ;
%       
%       2. For #3 to work, matlab needs to exist in the system path, 
%       e.g. add the following lines (adapted to refer to your version of MATLAB)
%       to ~/.bash_profile (or ~/.bashrc)
%
%           # add MATLAB path
%           export PATH=$PATH:/Applications/MATLAB_R2015a.app/bin/ 
%
%       3. The user needs to start their MATLAB session from the command line,
%       and from their home folder, e.g. 
%           
%           user@polymtl ~ $ matlab &
%
% .......
%
%   P has properties:
%       .Data 
%       .Log
%       .Source    
%       .Specs
%
% =========================================================================
% Author: ryan.topfer@polymtl.ca
% =========================================================================

properties   
    Data ;
    Log ;
    Source ;
    Specs ; % state = {active, inactive, inert, void}
end

methods
    
%% ========================================================================
function Aux = ProbeTracking( varargin )
%PROBETRACKING  

Aux.Log            = []; % memmapfile object pertaining to recording log

Aux.Data.p         = []; % may be filtered & limited
Aux.Data.pRaw      = []; % raw measurement
Aux.Data.t         = []; % measurement time [units: ms]
Aux.Data.trigger   = []; % measurement time [units: ms]

Aux.Data.startTime = [] ;
Aux.Data.endTime   = [] ;

Aux.Data.iSample   = [] ;

if nargin < 1 || isempty( varargin{1} )
    Specs = [] ;

elseif isstruct( varargin{1} )
    Specs = varargin{1} ;

elseif ischar( varargin{1} )
    
    filename = varargin{1} ;
    load( filename ) ;

    if exist('Data')
    % loaded file is Data struct corresponding to an old recording
        Aux.Data    = Data ;
        Specs.state = 'inert' ;
    else
    % loaded file is itself a ProbeTracking() object to be managed by the daemon session 
    %
    % Nonurgent TODO:
    %   This form of ProbeTracking() initialization/construction is not to
    %   be called by a user, but rather, from the ProbeTracking() constructor
    %   itself, making it better suited as a 'private' constructor. Apparently
    %   Matlab does not permit this. There is a work-around described here:
    %   https://stackoverflow.com/questions/29671482/private-constructor-in-matlab-oop
        Aux.launchrecordingdaemon() ; % runs continuously in background
        return;
    end
end

if myisfield( Specs, 'state' ) && strcmp( Specs.state, 'inert' )
    Aux.Source = [] ;
    Aux.Specs  = Specs ;
else
    [ Aux.Source, Aux.Specs ] = ProbeTracking.declareprobe( Specs ) ;
    
    if isa( Aux.Source, 'serial' )    
        Aux.createlogfile() ;
        if Aux.Specs.isRecordingDaemonEnabled
            Aux.createrecordingdaemon() ;
        end
    end
end

end

%% ========================================================================
function [AuxCopy] = copy( Aux )
%COPY  
% 
% Aux = COPY( Aux )

Specs       = Aux.Specs ;
Specs.state = 'inert' ;

AuxCopy = ProbeTracking( Specs ) ; 

AuxCopy.Data = Aux.Data ;

end

%% ========================================================================
function [] = delete( Aux )
%DELETE  
%
% DELETE( Aux )

if isa( Aux.Source, 'serial' ) 
    fclose( Aux.Source ) ;
    delete( Aux.Source ) ;
    clear Aux.Source  ;
else
    % TODO 
    %   prompt user to save latest recording?
end

clear Aux ;

end

%% ========================================================================
function [isRecording] = beginrecording( Aux )
%BEGINRECORDING - Initialize & open (RS-232) communication port 
%
% [isRecording] = BEGINRECORDING( Aux )
%
% Opens Aux.Source 
%
% Returns TRUE if successful

isRecording  = false ;

if strcmp( Aux.Specs.state, 'inert' )
    warning('Invalid command: Aux object is inert and cannot record.')
    return ;
end

if isa( Aux.Source, 'serial' )
    
    fopen( Aux.Source ) ;
    disp('Connecting to respiratory probe...')
    
    % wait for Arduino reset then assign sampling period 
    pause(2) 
    fprintf( Aux.Source, '%d\n', Aux.Specs.dt ) ;
    
    iAttempt    = 1 ;
    nAttempsMax = 3 ; 

    while( ~isRecording && iAttempt <= nAttempsMax )

        disp(['Attempt #' num2str(iAttempt)]);    
        
        tmp = fgetl( Aux.Source ) ;
        
        if ~isempty(tmp) 
            str2num( tmp ) ;
            isRecording = true ;
        end
        
        iAttempt = iAttempt + 1;

    end

    if isRecording
        disp('Communication successful. Reading in from serial port...')
        Aux.Log.Data.startTime = str2num( datestr( now, 'yyyymmddHHMMSS.FFF') ) ;  % to read it while debugging, use vpa(Aux.Log.Data.startTime)
        Aux.Data.startTime     = Aux.Log.Data.startTime ;
        Aux.Log.Data.endTime   = Inf ; 
        Aux.Data.endTime       = Aux.Log.Data.endTime ; 

        Aux.Specs.state = 'active' ;
    else
        Aux.stoprecording() ;
        error('Communication to respiratory probe failed. Closing Aux.Source.')
    end

else
    % Check daemon session is still recording
    if Aux.Log.Data.isLogging 
        isRecording = true ;
        Aux.Data.startTime = str2num( datestr( now, 'yyyymmddHHMMSS.FFF') ) ; 
    end
end

end

%% ========================================================================
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

%% ========================================================================
function [] = stoprecording( Aux )
%STOPRECORDING 
% 
% Closes communication port/source + marks recording end time
%
% [] = STOPRECORDING( Aux )
    
Aux.Data.endTime = str2num( datestr( now, 'yyyymmddHHMMSS.FFF') ) ; 

if isa( Aux.Source, 'serial' ) 
    fclose( Aux.Source ) ;
    Aux.Specs.state = 'inactive' ;
end

end

%% ========================================================================
function [] = clearrecording( Aux )
%CLEARRECORDING  
%
% Empties subfields of Aux.Data
%
% [] = CLEARRECORDING( Aux )

Aux.Data.pRaw      = [] ;
Aux.Data.p         = [] ;
Aux.Data.trigger   = [] ;

Aux.Data.t         = [] ;
Aux.Data.startTime = [] ;
Aux.Data.endTime   = [] ;

Aux.Data.iSample   = [] ; 

end

%% ========================================================================
function [pRaw, p, t] = getupdate( Aux )
%GETUPDATE 
%
% [pRaw, p] = GETUPDATE( Aux )

% Reads a single (16-bit) measurement (p) and returns p as the typecasted double.
% t is the sample time in units of milliseconds.
%
% p is either read from the open com port, or from the temp file buffer.

% Reading from open com port
if isa( Aux.Source, 'serial' ) 
    assert( strcmp( Aux.Source.Status, 'open' ), 'Error: Serial port is closed.' );

    tmp = str2num( fgetl( Aux.Source ) ) ;

    Aux.Data.trigger(end+1) = logical( tmp(1) ) ;
    Aux.Data.pRaw(end+1)    = tmp(2:end) ;
    
    p = filter_signal( Aux.Specs.probeType, Aux.Data.pRaw ) ;
    p = p(end) ; 
    Aux.Data.p(end+1) = p ;
    
    Aux.writeupdatetologfile( ) ;
    
    t = Aux.Log.Data.nSamples * Aux.Specs.dt ;
    Aux.Data.t(end+1) = t ;

% Reading from temp file buffer
elseif ischar( Aux.Source ) 

    [iSample, pRaw, p, t, trigger] = Aux.readupdatefromlogfile( ) ;

    Aux.Data.iSample        = iSample ;
    Aux.Data.pRaw(end+1)    = pRaw ;
    Aux.Data.p(end+1)       = p ;
    Aux.Data.t(end+1)       = t ;
    Aux.Data.trigger(end+1) = trigger ;

end

end

% % =========================================================================
% function [p] = processupdate( Aux )
% %PROCESSUPDATE
% % 
% % Filters/detrends raw recording in Aux.Data.pRaw
% %
% % [p] = PROCESSRECORDING( Aux )
% 
% % % limiting
% % if length( Aux.Data.pRaw ) > 100
% %
% %     pMean = mean( Aux.Data.pRaw(end-100:end) ) ;
% %     pStd  = std( Aux.Data.pRaw(end-100:end) ) ;
% %
% %     pRaw( abs( Aux.Data.pRaw(end) - pMean ) > 5*pStd ) = [] ;
% %
% % end
% 
% switch Aux.Specs.probeType
%     case 'capacitive'
%         
%         pRaw      = Aux.Data.pRaw ;
%       nSamples    = length(pRaw) ;
%       t           = Aux.Specs.dt*[0:nSamples-1] ;
%       nSamplesMin = 5 ;
% 
%       if nSamples < nSamplesMin
%           p = detrend( pRaw, 'constant' ) ;
%       else
%           y = polyfit( t, pRaw, 4 ) ;
%           p = pRaw - polyval( y, t ) ;
%       end
% 
%     case 'pressure'
%         % local windowing
%         if length( Aux.Data.pRaw ) > 3000
%             pRaw = Aux.Data.pRaw((end-3000):end) ;
%         else
%             pRaw = Aux.Data.pRaw ;
%         end
%       
%         p = detrend( pRaw, 'constant' ) ;
% end
%         
% end

%% ========================================================================
function [] = killrecordingdaemon( Aux )
%KILLRECORDINGDAEMON
%
% Sends STOP byte to daemon session

Aux.Log.Data.isLogging = uint64(0) ; 

end
%% ========================================================================
function [] = recordandplotphysiosignal( Aux, Params ) 
%RECORDANDPLOTPHYSIOSIGNAL
%   
% Calls ProbeTracking.recordphysiosignal() to record probe data, after which
% user is prompted whether to proceed (e.g. save + return) or re-record.
%
% [] = RECORDANDPLOTPHYSIOSIGNAL( Aux, Parameters )
%
% See HELP ProbeTracking.recordphysiosignal() for accepted Parameters

DEFAULT_ISSAVINGDATA          = true ;
DEFAULT_FILENAME  = [] ; % ProbeTracking.saverecording() will name it

if  nargin < 2 || isempty(Params)
    Params.dummy = [] ;
end

if  ~myisfield( Params, 'isSavingData' ) || isempty(Params.isSavingData)
    Params.isSavingData = DEFAULT_ISSAVINGDATA ;
else
    isSavingData = Params.isSavingData ;
    % Save only the final (user-approved) recording:
    Params.isSavingData = false ; 
end

if  ~myisfield( Params, 'filename' ) || isempty(Params.filename)
    Params.filename = DEFAULT_FILENAME ;
end

% ------- 
isUserSatisfied = false ;

while ~isUserSatisfied
    % boo hoo </3
    Aux.recordphysiosignal( Params ) ;

    % ------- 
    fprintf(['\n ----- \n Plotting physio recording for user verification... \n']) ;
    figure ;
    plot( Aux.Data.t/1000, Aux.Data.p, '+' ) ;
    xlabel('Time (s)');
    ylabel('Amplitude (AU)');
    title( 'Physio recording' ) ;
    
    response = input(['\n Is the current physio recording satisfactory? ' ...
        'Enter 0 to rerecord; 1 to continue. \n']) ;

    isUserSatisfied = logical(response) ;

end

% ------- 
if isSavingData
    Aux.saverecording( Params.filename )
end

end
%% ========================================================================
function [] = recordphysiosignal( Aux1, varargin )
%RECORDPHYSIOSIGNAL  
%
% Continuously tracks respiratory probe.
%
% Syntax
%
% [] = RECORDPHYSIOSIGNAL( Aux, Params )
% [] = RECORDPHYSIOSIGNAL( Aux, Aux2, Params )
%
%
% RECORDPHYSIOSIGNAL( Aux, Params ) will begin a new recording (and, optionally, real-time plotting)
% from the probe 'Aux'. 
% 
% With an additional probe input 'Aux2', dual recording is performed 
% NOTE: Aux and Aux2 must possess the same sampling period (i.e. Aux.Specs.dt == Aux2.Specs.dt)
%
%  .......................
%   
%   The following Params.fields are supported
%
%   .isSavingData
%       default = true
%
%   .filename
%       default: Named by ProbeTracking.saverecording() 
%
%   .runTime 
%       Total sampling time in seconds.
%       default = 15*60 [i.e. 15 minutes] 
%
%   .isPlottingInRealTime
%       [default : true ]
%
%   .refreshRate  
%       Rate at which the real-time display refreshes. 
%       Problems may arise if this is too fast!
%       [default : 4 Hz ]

DEFAULT_ISSAVINGDATA          = true ;
DEFAULT_FILENAME  = [] ; % ProbeTracking.saverecording() will name it
DEFAULT_RUNTIME               = 15*60 ; % [units : s]
DEFAULT_ISPLOTTINGINREALTIME  = true ;
DEFAULT_REFRESHRATE           = 4 ; % [units : Hz]

for iIn = 1 : length( varargin )
    
    if isstruct( varargin{iIn} )
        Params = varargin{iIn} ;

    elseif isa( varargin{iIn}, 'ProbeTracking' )
        Aux2 = varargin{iIn} ;
    end
end

if exist( 'Aux2' )
    assert( Aux1.Specs.dt == Aux2.Specs.dt, 'For dual tracking, the 2 probes must have the same sampling period.' )
    isDualTracking = true ;
    nProbes = 2 ;
else
    isDualTracking = false ;
    nProbes = 1 ;
end

if  nargin < 2 || ~exist( 'Params' )
    Params.dummy = [] ;
end

if  ~myisfield( Params, 'isSavingData' ) || isempty(Params.isSavingData)
    Params.isSavingData = DEFAULT_ISSAVINGDATA ;
end

if  ~myisfield( Params, 'filename' ) || isempty(Params.filename)
    Params.filename = DEFAULT_FILENAME ;
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

if isDualTracking
    Aux2.clearrecording() ;
end

sampleIndices  = [] ;
nSamples = Params.runTime / (Aux1.Specs.dt/1000) ;
iSample  = 0 ; 

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
   
    plotHandle1 = plot(axesHandle,0,0,'Marker','.','LineWidth',1,'Color',[1 0 0]);
    
    title(['RECORDING: ' Aux1.Specs.probeType ' probe'], 'FontSize',15,'Color',[1 1 0]);
    ylabel('Amplitude','FontWeight','bold','FontSize',14,'Color',[1 1 0]);
    xlabel('Time [s]','FontWeight','bold','FontSize',14,'Color',[1 1 0]);
    % xlabel('Sample index','FontWeight','bold','FontSize',14,'Color',[1 1 0]);

    drawnow limitrate; 
    set(figureHandle,'Visible','on');
   
    % -----
    % TTL/Trigger figure 
    figureHandleTrigger = figure('NumberTitle','off',...
        'Name','Trigger',...
        'Color',[0 0 0],'Visible','off');
        
    % Set axes
    axesHandleTrigger = axes('Parent',figureHandleTrigger,...
        'YGrid','on',...
        'YColor',[0.9725 0.9725 0.9725],...
        'XGrid','on',...
        'XColor',[0.9725 0.9725 0.9725],...
        'Color',[0 0 0]);

    hold on;
   
    plotHandleTrigger = plot(axesHandleTrigger,0,0,'Marker','.','LineWidth',1,'Color',[1 1 0]);
    
    title(['Recording: Trigger'],'FontSize',15,'Color',[1 1 0]);
    ylabel('Trigger','FontWeight','bold','FontSize',14,'Color',[1 1 0]);
    xlabel('Time [s]','FontWeight','bold','FontSize',14,'Color',[1 1 0]);

    drawnow limitrate; 
    set( figureHandleTrigger,'Visible','on');
   
    if isDualTracking 
        figureHandle2 = figure('NumberTitle','off',...
            'Name','Physio signal 2',...
            'Color',[0 0 0],'Visible','off');
            
        % Set axes
        axesHandle2 = axes('Parent',figureHandle2,...
            'YGrid','on',...
            'YColor',[0.9725 0.9725 0.9725],...
            'XGrid','on',...
            'XColor',[0.9725 0.9725 0.9725],...
            'Color',[0 0 0]);

        hold on;
       
        plotHandle2 = plot(axesHandle2,0,0,'Marker','.','LineWidth',1,'Color',[1 0 0]);
        
        title(['RECORDING: ' Aux2.Specs.probeType ' probe'],'FontSize',15,'Color',[1 1 0]);
        ylabel('Amplitude','FontWeight','bold','FontSize',14,'Color',[1 1 0]);
        xlabel('Time [s]','FontWeight','bold','FontSize',14,'Color',[1 1 0]);
        % xlabel('Sample index','FontWeight','bold','FontSize',14,'Color',[1 1 0]);

        drawnow limitrate; 
        set(figureHandle2,'Visible','on');
    end

end

isRecording = Aux1.beginrecording() ;

if isRecording
    display( ['Recording from probe 1/' num2str(nProbes) '...'] ) ;
else
    error(['Failed to record from probe 1']);
end

if isDualTracking
    isRecording = Aux2.beginrecording() ;
    
    if isRecording
        display( ['Recording from probe 2/' num2str(nProbes) '...'] ) ;
    else
        error(['Failed to record from probe 2']);
    end
end

% For real-time plotting, updating @the same rate as the samples
% (e.g. 100 Hz) poses a problem (computer can't seem to keep up with the incoming samples).
% solution is to update the display ~every so often~ (e.g. 4x per second seems OK)
nSamplesBetweenRefresh = (1/Params.refreshRate)/(Aux1.Specs.dt/1000) ;

while ( iSample < nSamples ) && ~StopButton.Stop()

    iSamplesBetweenRefresh = 0;

    for iSamplesBetweenRefresh = 1 : nSamplesBetweenRefresh 
        
        iSample = iSample + 1 ;
        Aux1.getupdate() ;
        
        if isDualTracking
            Aux2.getupdate() ;
        end
        
        sampleIndices(end+1) = iSample ;

    end

    if Params.isPlottingInRealTime
        
        set( plotHandle1,'XData', Aux1.Data.t/1000 ,'YData', Aux1.Data.p ); 
        set( plotHandleTrigger,'XData', Aux1.Data.t/1000 ,'YData', Aux1.Data.trigger ); 

        if isDualTracking
            set( plotHandle2,'XData', Aux2.Data.t/1000 ,'YData', Aux2.Data.p ); 
        end
   end

end

StopButton.Clear() ;

Aux1.stoprecording();

if isDualTracking
    Aux2.stoprecording() ;
end

% ------- 
if Params.isSavingData
    Aux1.saverecording( Params.filename )

    if isDualTracking
        Aux2.saverecording( )
    end
end

end

%% ========================================================================
function [] = saverecording( Aux, logFilename )
%SAVERECORDING
%
%   SAVERECORDING( Aux )
%   SAVERECORDING( Aux, logFilename )

if  nargin < 2 || ~ischar( logFilename ) 
    logFilename = ['./' datestr( now, 30 ) '_' Aux.Specs.probeType 'ProbeRec' ] ;
end
            
ShimUse.customdisplay( ['\n----- Saving physio recording -----\n'] );
ShimUse.customdisplay( ['Filename:  ' logFilename '\n'] );

Data = Aux.Data ;

Data.probeType = Aux.Specs.probeType ;

save( logFilename, 'Data' ) ;

end
%% ========================================================================
function [] = resetdaemon( Aux )
%RESETDAEMON
% 
% Resets the daemon recording: i.e. Preceding recording is retained in Aux.Log
% but the signal processing will begin afresh. (Use when patient position
% changes, for example.)

Aux.clearrecording() ;

% Changing .endTime from Inf signals daemon session to call this function &
% clear its current recording.
Aux.Log.Data.endTime = 0 ; 

end

end

% =========================================================================
methods( Access =  private)
%% ========================================================================
function [] = createlogfile( Aux )
%CREATELOGFILE
% 
% [] = CREATELOGFILE( Aux )

% Create or overwrite the file: 
filename = fullfile( tempdir, [ Aux.Specs.probeType '_probe_log.dat'] ) ;
[f, msg] = fopen( filename, 'w' ) ;

nSamplesMax = 21600000 ;

if f ~= -1
    fwrite( f, [ zeros(nSamplesMax+4, 3);], 'double' ) ; 
    fclose( f ) ;
else
    error('MATLAB:demo:send:cannotOpenFile', ...
          'Cannot open file "%s": %s.', filename, msg) ;
end

% Memory map the file.
Aux.Log = memmapfile( filename, 'Writable', true, 'Format', ...
    { 'uint64', [1 1], 'nSamples' ;
      'uint64', [1 1], 'isLogging' ;
      'double', [1 1], 'startTime' ;
      'double', [1 1], 'endTime' ;
      'uint8', [nSamplesMax 1], 'trigger' ;
      'double', [nSamplesMax 1], 'pRaw' ;
      'double', [nSamplesMax 1], 'p' } ) ;

Aux.Log.Data.isLogging = uint64(true) ;
Aux.Log.Data.endTime   = Inf ;

end

%% ========================================================================
function [nSamples, pRaw, p, t, trigger] = readupdatefromlogfile( Aux )
%READUPDATEFROMLOGFILE

DWELL_TIME     = 0.01 ; % [units: s]
MAX_DWELL_TIME = 1 ;

totalDwellTime = 0 ;

while ( Aux.Data.iSample == Aux.Log.Data.nSamples ) & ( totalDwellTime < MAX_DWELL_TIME )
    % Value has already been read. Wait for a new one.
    pause( DWELL_TIME ) ;
    totalDwellTime = totalDwellTime + DWELL_TIME ; 
end

% copy to new variable to ensure the value remains constant during reading
% NOTE: For some reason (?) typecasting as a double copies by value. Otherwise
% it is copied by reference,
nSamples = double( Aux.Log.Data.nSamples ) ;

pRaw    = Aux.Log.Data.pRaw( nSamples ) ;
p       = Aux.Log.Data.p( nSamples ) ;
t       = length( Aux.Data.p ) * Aux.Specs.dt ;
trigger = Aux.Log.Data.trigger( nSamples ) ;
 
end
%% ========================================================================
function [] = writeupdatetologfile( Aux )
%WRITEUPDATETOLOGFILE

% update N samples logged
Aux.Log.Data.nSamples = Aux.Log.Data.nSamples + 1 ;

Aux.Log.Data.pRaw( Aux.Log.Data.nSamples+1 )    = Aux.Data.pRaw(end) ;
Aux.Log.Data.p( Aux.Log.Data.nSamples+1 )       = Aux.Data.p(end) ;
Aux.Log.Data.trigger( Aux.Log.Data.nSamples+1 ) = Aux.Data.trigger(end) ;
   
end
%% ========================================================================
function [] = createrecordingdaemon( Aux )
%CREATERECORDINGDAEMON
% 
% CREATERECORDINDAEMON( Aux )
%
% Saves the instantiated Aux object, and launches a background (daemon) Matlab
% session to load the object and begin background recording.
%
% The daemon session reads directly from the USB (Com) port while the user
% session reads from a file buffer.

pathToAuxObject = [ tempdir 'Aux' ] ;
save( pathToAuxObject, 'Aux' ) ;

tmpDir = pwd ;
cd( '~/' ) ; % change to home folder so Matlab runs startup.m (defines path)
cmd = sprintf( '%s', 'matlab -r " ProbeTracking( ''', pathToAuxObject, ''' );" &') ;

ShimUse.customdisplay( ['\n----- Launching daemon session -----\n'] );
unix( cmd ) ;
pause(5) ;
cd( tmpDir ) ;

Aux.Source = 'logFile' ;

end
%% ========================================================================
function [] = launchrecordingdaemon( Aux )
%LAUNCHRECORDINGDAEMON
% 
% [] = LAUNCHRECORDINDAEMON( Aux )

Aux.clearrecording() ;
Aux.beginrecording() ;

while Aux.Log.Data.isLogging

    if Aux.Log.Data.endTime ~= Inf 
        Aux.resetdaemon() ;
        Aux.Log.Data.endTime = Inf ;
    end

    Aux.getupdate() ;
end

Aux.stoprecording() ;
Aux.saverecording( ) ;

% ------
% close serial port, delete object, and quit matlab session 
Aux.delete() ;
quit ;

end
% =========================================================================

end

% =========================================================================
methods(Static)

%% ========================================================================
function [Source, AuxSpecs] = declareprobe( AuxSpecs )
%DECLAREPROBE Declares serial object for probe 
% 
% [Source, AuxSpecs] = declareprobe( AuxSpecs )
%
% AuxSpecs can have the following fields 
% 
% .isRecordingDaemonEnabled 
%   Perform device I/O (probe recording) and signal processing in a 
%   secondary background/daemon MATLAB session.
%   User session can monitor the live recording with Probe.recordphysiosignal()
%   [default = true]
%
% .portName 
%   Address of the probe-associated serial port within file system
%   e.g. AuxSpecs.portName = '/dev/tty.usbmodem487312' ;

% NOTE 
%   System specs are hardcoded into the probe microcontroller. 
%
%   TODO: rather than hardcode the values here as well, there should be 
%   a parameter file specifying values for baudrate and sampling period.
DEFAULT_ARDUINOPERIOD = 100 ; % [units: ms] 
DEFAULT_TEENSYPERIOD  = 100 ; % [units: ms] 
DEFAULT_BAUDRATE      = 115200 ;

DEFAULT_ISRECORDINGDAEMONENABLED = true ;

if nargin < 1 || isempty( AuxSpecs )
    AuxSpecs.dummy = [] ;
end

if ~myisfield( AuxSpecs, 'isRecordingDaemonEnabled' ) || isempty(AuxSpecs.isRecordingDaemonEnabled)
    AuxSpecs.isRecordingDaemonEnabled = DEFAULT_ISRECORDINGDAEMONENABLED ;
end

AuxSpecs.clipLimits = [-Inf Inf] ; % [units: probe signal] 

% Check the presence of the device 
isDeviceFound = false ;
% If user has specified the portName:
if myisfield( AuxSpecs, 'portName' ) && ~isempty(AuxSpecs.portName)
   
    [fileDir,portname, fileExtension] = fileparts( AuxSpecs.portName ) ;
    
    listOfDevices = dir( fileDir ) ;
    
    for iFile = 1 : length(listOfDevices)-1
        if strcmp( [fileDir '/' listOfDevices(iFile).name fileExtension], AuxSpecs.portName )
            isDeviceFound = true ;
       end
    end

    if ~isDeviceFound    
        disp(['Warning: Given port name [ ' AuxSpecs.portName ' ] not found. Checking default port names.']) ;
        AuxSpecs.portName = [] ;
    end
end

% If portName unspecified by user, or if user-specified portName was not discovered:
if ~myisfield( AuxSpecs, 'portName' ) || isempty(AuxSpecs.portName)
    if ismac  % macOS
        
        listOfDevices = dir( '/dev/tty.usbmodem*' ) ;

        if length(listOfDevices) == 0
            warning( 'Respiratory Probe: Device file not found. Check USB device is connected.' ) ;
        elseif length(listOfDevices) == 1
            isDeviceFound     = true ;
            AuxSpecs.portName = ['/dev/' listOfDevices(1).name ] ;
        else
            warning('Ambiguous device identifier. Consider entering portName as argument. See HELP.') ;
        end

    elseif isunix  % Linux
        % if exist('/dev/ttyS100', 'file')  
            AuxSpecs.portName = '/dev/ttyS100' ;
            warning( [ 'ProbeTracking.declareprobe() currently assumes a device address of ' AuxSpecs.portName ' which may be invalid!' ] ) ;
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

% Configure the device
if isDeviceFound
    
    AuxSpecs.state = 'inactive' ;
    
    % Check device type 
    % NOTE : names probably need to change computer-to-computer!
    if ~isempty(strfind(AuxSpecs.portName, '/dev/tty.usbmodem487312')) || ~isempty(strfind(AuxSpecs.portName, '/dev/tty.usbmodem447189'))
        AuxSpecs.probeType = 'capacitive' ;
    elseif ~isempty(strfind(AuxSpecs.portName, '/dev/tty.usbmodem141')) || ~isempty(strfind(AuxSpecs.portName, '/dev/tty.usbmodem142'))
        AuxSpecs.probeType = 'pressure' ;
    else
        error('Probe was not found. Please check your /dev folder and see if you identify the probe (should start with tty.usbmodem and add this case in the code');
    end

    switch AuxSpecs.probeType
        case 'pressure'
            ShimUse.customdisplay( ['\n----- Pressure probe -----\n'] );
            AuxSpecs.dt       = DEFAULT_ARDUINOPERIOD ; % [units: ms]
            AuxSpecs.baudRate = DEFAULT_BAUDRATE ; % [units: ms]
        case 'capacitive'
            ShimUse.customdisplay( ['\n----- Capacitive probe -----\n'] );
            AuxSpecs.dt       = DEFAULT_TEENSYPERIOD ; 
            AuxSpecs.baudRate = DEFAULT_BAUDRATE ; % [units: ms]
    end
    
    Source = serial( AuxSpecs.portName, 'BaudRate',  AuxSpecs.baudRate ) ;

    samplingFrequency = 1000/AuxSpecs.dt ; % [units: Hz]
    
    ShimUse.customdisplay( [ 'Sampling frequency = ' num2str(samplingFrequency) ' Hz'] )
    
    % [b,a] = butter( 4, 0.064/(samplingFrequency/2) ) ;
    %
    % AuxSpecs.Filter = [] ;
    % AuxSpecs.Filter.Lowpass.order  = 4 ;
    % AuxSpecs.Filter.Lowpass.cutoff = 0.0001 ;
    % AuxSpecs.Filter.Lowpass.Coefficients = [] ;
    % AuxSpecs.Filter.Lowpass.Coefficients.numerator   = b ;
    % AuxSpecs.Filter.Lowpass.Coefficients.denominator = a ;

else
    Source = [] ;
    AuxSpecs.state = 'inert' ;
end

end

%% ========================================================================
function [p, sampleTimes] = loadmeasurementlog( measurementLogFilename, sampleTimesFilename )
%LOADMEASUREMENTLOG
% 
% Reads binary file of data measurements (e.g. pressure recording) to return
% vector(s) of doubles.
%
% measurementLog                = LOADMEASUREMENTLOG( measurementLogFilename ) ;
% [measurementLog, sampleTimes] = LOADMEASUREMENTLOG( measurementLogFilename, sampleTimesFilename )

if nargin < 1
    error( 'Insufficient arguments. Must provide full path to measurement log file.' ) ;

else
    if nargin >= 1
        [~,~,ext] = fileparts( measurementLogFilename ) ;

        if strcmp( ext, '.txt' )
            X = load( measurementLogFilename ) ;
            sampleTimes    = X(:,1) ;
            p = X(:,2) ;
            return;
        elseif strcmp( ext, '.bin' )
            measurementLogFid = fopen( measurementLogFilename, 'r' ) ;
            p    = fread( measurementLogFid, inf, 'double' ) ;
            fclose( measurementLogFid );
        end
    end
    
    if nargin == 2 % saved as binary
        sampleTimesFid = fopen( sampleTimesFilename, 'r' ) ;
        sampleTimes    = fread( sampleTimesFid, inf, 'double' ) ;
        fclose( sampleTimesFid );
    end
end

end

%% ========================================================================
function [] = plotmeasurementlog( measurementLog, Params )
%PLOTMEASUREMENTLOG
%
% PLOTMEASUREMENTLOG( measurementLog ) ;
% PLOTMEASUREMENTLOG( measurementLog, Params )
%
% Supported fields to Params struct
%
%   .figureTitle
%       [default: 'Respiration']
%
%   .sampleTimes
%       vector (length == length(measurementLog)) of sample times in seconds
%
%   .yLabel
%       [default: 'Amplitude (AU)']

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

%% ========================================================================
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
% =========================================================================
% =========================================================================

end
% =========================================================================
% =========================================================================

end
