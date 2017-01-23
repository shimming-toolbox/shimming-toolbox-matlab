classdef ShimUse 
%SHIMUSE - Shim Use
%
% .......
% 
% Description
% 
% ShimUse is a high-level user interface to operate the shim system. 
%
% .......
%   
% Usage
%
% Shim = ShimUse( Params )
% 
% Params
%   
%   .shimSystem 
%       'Acdc' or 'Rri' [default]
%
%   .pathToShimReferenceMaps
%   
%   .ProbeSpecs
%       .arduinoPeriod
%
%
%   Shim contains fields
%
%       .Com
%           Object of type ShimCom
%
%       .Opt
%           Object of type ShimOpt
%
%
% =========================================================================
% Notes 
% 
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
% Updated::20161129::ryan.topfer@polymtl.ca
% =========================================================================

% =========================================================================
% *** TODO 
%
% ..... 
% ()
%
%
% ..... 
% =========================================================================

properties   
    Com;
    Opt;
end


% =========================================================================
% =========================================================================
methods
% =========================================================================
function Shim = ShimUse( Params )
%SHIMUSE   

DEFAULT_SHIMSYSTEM = 'Rri' ; 
DEFAULT_RUNMODE    = 'isCmdLine' ;% vs. 'isGui'

if nargin < 1
    Params.dummy = [];
end

if ~myisfield( Params, 'runMode' ) || isempty( Params.runMode ) 
   Params.runMode = DEFAULT_RUNMODE ;
end

if ~myisfield( Params, 'shimSystem' ) || isempty( Params.shimSystem ) 
   Params.shimSystem = DEFAULT_SHIMSYSTEM ;
end

switch Params.shimSystem
    
    case 'Rri'
        Shim.Opt = ShimOptRri( Params ) ;
        Shim.Com = ShimComRri( ) ;

    case 'Acdc'
        Shim.Opt = ShimOptAcdc( Params ) ;
        Shim.Com = ShimComAcdc( ) ;

    otherwise
        error([ Params.shimSystem 'is an invalid or unimplemented shim system. See HELP ShimUse().' ])

end

if ~strcmp( Params.runMode,  'isDebugging' ) ;
    Shim.testshimconnection() ;
end

end
% =========================================================================
function [] = delete( Shim )
%DELETE  (custom helper function)
% 
% DELETE( Shim )
% 
% Destructor method: calls Shim.Com.delete( ) and Shim.Opt.delete( ) 

Shim.Opt.delete();
Shim.Com.delete();
clear Shim ;

end
% =========================================================================
function [] = runrealtimeshim( Shim, Params )
%RUNREALTIMESHIMMING
% 
% Compute and set optimal shim current update 
%
% [] = RUNREALTIMESHIM( Shim, Params  )

assert( strcmp( Shim.Opt.Probe.ComPort.Status, 'closed' ), ...
    'Error: Serial port is open/in use.' );

DEFAULT_ISSAVINGDATA          = true ;
DEFAULT_ISFORCINGOVERWRITE    = false ;
DEFAULT_PRESSURELOGFILENAME   = ['./' datestr(now,30) '-pressureLog.bin' ] ;
DEFAULT_SAMPLETIMESFILENAME   = ['./' datestr(now,30) '-sampleTimes.bin' ] ;
DEFAULT_RUNTIME               = 5*60 ; % [units: s]
DEFAULT_EXTRAPOLATIONORDER    = 0 ;
DEFAULT_EXTRAPOLATIONDELAY    = 0 ;

DEFAULT_ISFILTERINGPRESSURE   = false ;
DEFAULT_ISPLOTTINGINREALTIME  = true ;

DEFAULT_ISCLIPPINGPRESSURE    = true ;
DEFAULT_MINCLIPPINGPRESSURE   = 1 ;
DEFAULT_MAXCLIPPINGPRESSURE   = 1000 ;

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

[pathStr,name,ext] = fileparts( Params.pressureLogFilename ) ;
Params.rawPressureLogFilename = [pathStr name '_raw' ext] ;

if  ~myisfield( Params, 'sampleTimesFilename' ) || isempty(Params.sampleTimesFilename)
    Params.sampleTimesFilename  = DEFAULT_SAMPLETIMESFILENAME ; 
end

if  ~myisfield( Params, 'runTime' ) || isempty(Params.runTime)
    Params.runTime = DEFAULT_RUNTIME ;
end

if  ~myisfield( Params, 'extrapolationOrder' ) || isempty(Params.extrapolationOrder)
    Params.extrapolationOrder  = DEFAULT_EXTRAPOLATIONORDER ; 
end

if  ~myisfield( Params, 'extrapolationDelay' ) || isempty(Params.extrapolationDelay)
    Params.extrapolationDelay  = DEFAULT_EXTRAPOLATIONDELAY ; 
end

if  ~myisfield( Params, 'isPlottingInRealTime' ) || isempty(Params.isPlottingInRealTime)
    Params.isPlottingInRealTime  = DEFAULT_ISPLOTTINGINREALTIME ; 
end

if  ~myisfield( Params, 'maxCurrents' ) || isempty(Params.maxCurrents)
    Params.maxCurrents  = ones( Shim.Com.Specs.Amp.nActiveChannels, 1) ; 
end

if  ~myisfield( Params, 'minCurrents' ) || isempty(Params.minCurrents)
    Params.minCurrents  = -ones( Shim.Com.Specs.Amp.nActiveChannels, 1) ; 
end

if  ~myisfield( Params, 'isFilteringPressure' ) || isempty(Params.isFilteringPressure)
    Params.isFilteringPressure  = DEFAULT_ISFILTERINGPRESSURE ; 
end

if  ~myisfield( Params, 'isClippingPressure' ) || isempty(Params.isClippingPressure)
    Params.isClippingPressure  = DEFAULT_ISCLIPPINGPRESSURE ; 
end

if  ~myisfield( Params, 'minClippingPressure' ) || isempty(Params.minClippingPressure)
    Params.minClippingPressure  = DEFAULT_MINCLIPPINGPRESSURE ; 
end

if  ~myisfield( Params, 'maxClippingPressure' ) || isempty(Params.maxClippingPressure)
    Params.maxClippingPressure  = DEFAULT_MAXCLIPPINGPRESSURE ; 
end

Params.nSamplesFilter = 10 ;
Params.filterScaling  = 4 ; 

filterWeights = (-Params.filterScaling/Params.nSamplesFilter)*[0:Params.nSamplesFilter-1]' ;
filterWeights = flipud( filterWeights ) ;
filterWeights = exp( filterWeights ) ;

summedWeights = sum( filterWeights ) ;
filterWeights = filterWeights' / summedWeights ;

    
normCurrents = [];

% ------- 
Shim.Opt.Probe.Data.pressure = [] ;

nSamples = Params.runTime / (Shim.Opt.Probe.Specs.arduinoPeriod/1000) ;
iSample  = 0 ; 

nSamplesBetweenUpdates = ... % CHANGE THIS LINE FOR ACDC
    Shim.Com.Specs.Com.mxdUpdatePeriod/(Shim.Opt.Probe.Specs.arduinoPeriod/1000)

% ------- 
% figure
figureHandle = figure('NumberTitle','off',...
    'Name','Sonde de respiration',...
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
    
title('Pressure log','FontSize',15,'Color',[1 1 0]);
xlabel('Sample index','FontWeight','bold','FontSize',14,'Color',[1 1 0]);
ylabel('Pressure (Pa)','FontWeight','bold','FontSize',14,'Color',[1 1 0]);

drawnow limitrate; 
 
if Params.isPlottingInRealTime
    set(figureHandle,'Visible','on');
end

sampleIndices = [] ;

% ------- 
StopButton = stoploop({'Stop recording'}) ;
Shim.Opt.Probe = Shim.Opt.Probe.initializecomport() ;
ShimUse.display( 'Issuing real-time shim updates. Begin scanning!' )
ShimUse.display( 'Warning: 1st shim update based on raw pressure recording.' )

rawPressureLog = [] ;

while ( iSample < nSamples ) && ~StopButton.Stop()
    
    iSamplesBetweenUpdates = 0;

    for iSamplesBetweenUpdates = 1 : nSamplesBetweenUpdates 

        iSample = iSample + 1 ;
        
        sampleIndices(iSample) = iSample ;
       
        rawPressureLog(iSample) = Shim.Opt.Probe.readpressure() ;
        
        if Params.isFilteringPressure  && ( iSample > Params.nSamplesFilter )
            Shim.Opt.Probe.Data.pressure(iSample) = rawPressureLog(iSample) ;
            Shim.Opt.Probe.Data.pressure(end) = ...
                sum( filterWeights .* Shim.Opt.Probe.Data.pressure( (end - (Params.nSamplesFilter-1)) : end ) ) ;

        else
            
            Shim.Opt.Probe.Data.pressure(iSample) = rawPressureLog(iSample) ;
        
        end
        
        if Params.isClippingPressure 
            
            Shim.Opt.Probe.Data.pressure(iSample) = clippressure( Shim.Opt.Probe.Data.pressure(iSample) ) ;
        
        end

    end

    currents = Shim.Opt.computerealtimeupdate( Params ) 
    % currents = limitcurrents( currents ) ;
    currentsNorm = norm(currents) 
    
    Shim.Com.setandloadallshims( currents ) ;
    
    if Params.isPlottingInRealTime
        set(plotHandle,'YData',Shim.Opt.Probe.Data.pressure,'XData',sampleIndices);
    end

end

sampleTimes = (Shim.Opt.Probe.Specs.arduinoPeriod/1000)*[0:(numel(Shim.Opt.Probe.Data.pressure)-1)] ;

fclose( Shim.Opt.Probe.ComPort )

StopButton.Clear() ;


Shim.resetallshims() ;

% ------- 
if Params.isSavingData
    pressureLogFid = fopen( Params.pressureLogFilename, 'w+' ) ;
    fwrite( pressureLogFid, Shim.Opt.Probe.Data.pressure, 'double' ) ;
    fclose( pressureLogFid );
    
    pressureLogFid = fopen( Params.rawPressureLogFilename, 'w+' ) ;
    fwrite( pressureLogFid, rawPressureLog, 'double' ) ;
    fclose( pressureLogFid );

    sampleTimesFid = fopen( Params.sampleTimesFilename, 'w+' ) ;
    fwrite( sampleTimesFid, sampleTimes, 'double' ) ;
    fclose( sampleTimesFid );
end

function pressureSample = clippressure( pressureSample )
%CLIPPRESSURE
% 
% pressureSample = clippressure( pressureSample )

if pressureSample < Params.minClippingPressure
    pressureSample = Params.minClippingPressure ;
elseif pressureSample > Params.maxClippingPressure ;
    pressureSample = Params.maxClippingPressure ;
end

end

% function currents = limitcurrents( currents )
% %limitcurrents
%
%     isClipping = false ;
%
%     iChannelsOver = currents > Params.maxCurrents ; 
%     currents( iChannelsOver ) = Params.maxCurrents( iChannelsOver )  ;
%
%     iChannelsUnder = currents < Params.minCurrents ; 
%     currents( iChannelsUnder ) = Params.minCurrents( iChannelsUnder )  ;
%
%     if (nnz(iChannelsOver) + nnz(iChannelsUnder)) > 0
%         ShimUse.display('CLIPPING');
%     end
% end

end
% =========================================================================
function isAckReceived = testshimconnection( Shim )
%TESTSHIMCONNECTION
% 
% isAckReceived = TESTSHIMCONNECTION( Shim )
% 
% Queries shim amp for response (1 or 0)

isAckReceived = Shim.Com.getsystemheartbeat ; % should respond with 'ACK'

if isAckReceived
    msg = [ '-----\n'...
        'Communication to shim amplifier successful. \n -----'] ;
    ShimUse.display( msg ) ;
else
    msg = [ '-----\n'...
            'Communication to shim amplifier failed. Check device connections.' ... 
            '\n -----'] ;
    ShimUse.display( msg ) ;
end

end
% =========================================================================

end
% =========================================================================
% =========================================================================
methods(Static) 
% =========================================================================
function [systemResponses, meanings] = definesystemresponses( )
%DEFINESYSTEMRESPONSES
%
% [systemResponses, meanings] = DEFINESYSTEMRESPONSES( )
%
%   systemResponses is a string array of HEX messages
% 
%   meanings is a cell array where each entry gives a brief definition of the
%   corresponding entry in systemResponses.
% 
% -------------------------------------------------------------------------
% According to RRI HEX specification protocol, system responses from the 
% MXD/DSU can be one of the following:  
%
% i.   0x00: ACK , command has been received and processed
% ii.  0x01: Unimplemented command
% iii. 0x02: A command parameter is out of range
% iv.  0x03: An error is present that prevents execution of the command
% v.   0x08: Invalid command
% vi.  0xFF: NACK, unrecognized command or bad CRC

systemResponse = ['0x00'; '0x01'; '0x02'; '0x03'; '0x08'; '0xFF'] ;

meanings    = cell{length(systemResponse), 1} ;
meanings{1} = 'ACK , command has been received and processed' ;
meanings{2} = 'Unimplemented command' ;
meanings{3} = 'A command parameters is out of range' ;
meanings{4} = 'An error is present that prevents command execution' ;
meanings{5} = 'Invalid command' ;
meanings{6} = 'NACK, unrecognized command or bad CRC' ;

end
% =========================================================================
function [] = display( msg )
%DISPLAY
%

if nargin < 1 || isempty(msg)
    fprintf('\n')    ;
    % SystemInfo = Shim.getsysteminformation( )
    % msg = SystemInfo ;

else 
    assert( isstr(msg), 'Given message is not a string.' ) ;

    % switch Shim.Params.runMode 
        % case 'isCmdLine'
            fprintf(['\n' msg '\n']) ;
        % case 'isGui'
        %     fprintf(['\n' 'Error: GUI not yet supported!' '\n\n']) ;
    % end
end

end
% % =========================================================================
% function [] = listcommands( option )
% %LISTCOMMANDS
% %
% %   Displays all the implemented HEX commands for MXD/DSU.
% %   
% %   [] = LISTCOMMANDS( )
%
% Cmd = ShimCom.definecommands( ) ;
%
% fprintf(['\n ---------------------------------------------------------\n'])
% fprintf(['\n MXD Commands:\n\n'])
% commands = fieldnames( Cmd.Mxd ) ;
%
% for iCmd = 1 : length(commands)
%     fprintf([ commands{iCmd} '\n'])
% end
%
% fprintf(['\n ---------------------------------------------------------\n'])
% fprintf(['\n DSU Commands:\n\n'])
% commands = fieldnames( Cmd.Dsu ) ;
%
% for iCmd = 1 : length(commands)
%     fprintf([commands{iCmd} '\n'])
% end
%
% fprintf(['\n\n'])
%
% end
% =========================================================================

end
% =========================================================================
% =========================================================================

end
