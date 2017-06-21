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
%   .TrackerSpecs
%       .dt
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

Shim.Opt = [];
Shim.Com = [];

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

% assert( strcmp( Shim.Opt.Tracker.ComPort.Status, 'closed' ), ...
%     'Error: Serial port is open/in use.' );

DEFAULT_ISSAVINGDATA            = true ;
DEFAULT_ISFORCINGOVERWRITE      = false ;
DEFAULT_MEASUREMENTLOGFILENAME  = ['./' datestr(now,30) '-pressureLog.bin' ] ;
DEFAULT_SAMPLETIMESFILENAME     = ['./' datestr(now,30) '-sampleTimes.bin' ] ;
DEFAULT_UPDATETIMESFILENAME     = ['./' datestr(now,30) '-updateTimes.bin' ] ;

DEFAULT_RUNTIME                 = 5*60 ; % [units: s]
DEFAULT_EXTRAPOLATIONORDER      = 0 ;
DEFAULT_EXTRAPOLATIONDELAY      = 0 ;

DEFAULT_ISFILTERINGMEASUREMENTS = false ; % Tracker measurements
DEFAULT_ISPLOTTINGINREALTIME    = true ;

DEFAULT_ISCLIPPING            = true ;
DEFAULT_MINCLIPTHRESHOLD  = 1 ;
DEFAULT_MAXCLIPTHRESHOLD  = 1000 ;


if  nargin < 2 || isempty(Params)
    Params.dummy = [] ;
end

if  ~myisfield( Params, 'isSavingData' ) || isempty(Params.isSavingData)
    Params.isSavingData = DEFAULT_ISSAVINGDATA ;
end

if  ~myisfield( Params, 'isForcingOverwrite' ) || isempty(Params.isForcingOverwrite)
    Params.isForcingOverwrite = DEFAULT_ISFORCINGOVERWRITE ;
end

if  ~myisfield( Params, 'measurementLogFilename' ) || isempty(Params.measurementLogFilename)
    Params.measurementLogFilename = DEFAULT_MEASUREMENTLOGFILENAME ;
end

[pathStr,name,ext] = fileparts( Params.measurementLogFilename ) ;
Params.rawMeasurementLogFilename = [pathStr '/' name '_raw' ext] ;

if  ~myisfield( Params, 'sampleTimesFilename' ) || isempty(Params.sampleTimesFilename)
    Params.sampleTimesFilename  = DEFAULT_SAMPLETIMESFILENAME ; 
end

if  ~myisfield( Params, 'updateTimesFilename' ) || isempty(Params.updateTimesFilename)
    Params.updateTimesFilename  = DEFAULT_UPDATETIMESFILENAME ; 
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

if  ~myisfield( Params, 'isFilteringMeasurements' ) || isempty(Params.isFilteringMeasurements)
    Params.isFilteringMeasurements  = DEFAULT_ISFILTERINGMEASUREMENTS ; 
end

if  ~myisfield( Params, 'isClipping' ) || isempty(Params.isClipping)
    Params.isClipping  = DEFAULT_ISCLIPPING ; 
end

if  ~myisfield( Params, 'minClipThreshold' ) || isempty(Params.minClipThreshold)
    Params.minClipThreshold  = DEFAULT_MINCLIPTHRESHOLD ; 
end

if  ~myisfield( Params, 'maxClipThreshold' ) || isempty(Params.maxClipThreshold)
    Params.maxClipThreshold  = DEFAULT_MAXCLIPTHRESHOLD ; 
end





Params.nSamplesFilter  = 5 ; % Window length -> 5 samples * 10 ms/sample = 50 ms
Params.nSamplesHalfWindow = (Params.nSamplesFilter + 1)/2 - 1 ; 

Params.correctionOrder = 1;  % linear correction

[~,Params.filterWeights] = sgolay(Params.correctionOrder, Params.nSamplesFilter);   % Calculate S-G coefficients

% delay = tx/transmission delay plus the time shift involved from the filter
delay = ( Params.txDelay + Shim.Opt.Tracker.Specs.dt*(Params.nSamplesFilter-1)/2)/1000 ; % [units: s]



% ------- 
Shim.Opt.Tracker.Data.p = [] ;

nSamples = Params.runTime / (Shim.Opt.Tracker.Specs.dt/1000) ;
iSample  = 0 ; 

iSamplesBetweenUpdates = 0;
nSamplesBetweenUpdates = ... % CHANGE THIS LINE FOR ACDC --REFERS TO RRI/MXD
    Shim.Com.Specs.Com.mxdUpdatePeriod/(Shim.Opt.Tracker.Specs.dt/1000)
updatePeriod = Shim.Com.Specs.Com.mxdUpdatePeriod ;


    rawMeasurementLog    = [] ;
    Shims.Data.Tracker.p = [] ;

    sampleTimes          = [] ; %
    updateTimes          = [] ; % when shim updates occur

    sampleIndices        = [] ;
    updateIndices        = [] ; % corresponding to when shim updates occur
    iUpdate              = 0 ;

    currentsNorm         = 0 ;



if Params.isPlottingInRealTime
    
    close all

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
        
    xlim(axesHandle,[0 Params.runTime]);
        
    title('Respiration Tracker','FontSize',15,'Color',[1 1 0]);
    xlabel('Time (s)','FontWeight','bold','FontSize',14,'Color',[1 1 0]);
    ylabel('Amplitude','FontWeight','bold','FontSize',14,'Color',[1 1 0]);

    drawnow limitrate; 
    set(figureHandle,'Visible','on');

end


% ------- 
StopButton = stoploop({'Stop recording'}) ;

isTracking = Shim.Opt.Tracker.begintracking(); 

if isTracking

    ShimUse.display( 'Issuing real-time shim updates. Begin scanning!' )

    while ( iSample < nSamples ) && ~StopButton.Stop()
        
        iSamplesBetweenUpdates = 0;

        % acquire batch of respiration samples
        for iSamplesBetweenUpdates = 1 : nSamplesBetweenUpdates 

            iSample = iSample + 1 ;
            
            sampleIndices(iSample)     = iSample ;
            
            rawMeasurementLog(iSample) = Shim.Opt.Tracker.getupdate() ;

        end

        iUpdate = iUpdate + 1;        

        updateIndices(iUpdate) = iUpdate ;
        
        updateTimes(iUpdate) = updatePeriod * iUpdate ; 

        % lowpass + delay correction of respiratory signal  
        if Params.isFilteringMeasurements  && ( iSample > Params.nSamplesFilter )
                
             % 0th order corr (weighted avg)
             p = dot( Params.filterWeights(:,1), rawMeasurementLog( (iSample - Params.nSamplesFilter + 1) : iSample ) ) ;

             % if Params.isPredictingMeasurement
             % 1st order corr
             p = p + delay*dot( Params.filterWeights(:,2), rawMeasurementLog( (iSample - Params.nSamplesFilter + 1) : iSample ) ) ;
             %end
             p = clipvalue( p ) ;

            Shim.Opt.Tracker.Data.p(iUpdate) = p ;

            currents = Shim.Opt.computerealtimeupdate( ) 
        else

            Shim.Opt.Tracker.Data.p(iUpdate) = rawMeasurementLog(iSample) ;

            currents = Shim.Opt.computerealtimeupdate(  ) ;
        end
        
        % currents = limitcurrents( currents ) ;
        currentsNorm = norm(currents) ; 
        
        Shim.Com.setandloadallshims( currents ) ;
        
        if Params.isPlottingInRealTime
            set(plotHandle,'YData',Shim.Opt.Tracker.Data.p,'XData',updateTimes);
        end

    end

    Shim.Opt.Tracker.stoptracking() ;
    Shim.Com.resetallshims() ;

    sampleTimes = (Shim.Opt.Tracker.Specs.dt/1000)*sampleIndices ;

    % ------- 
    if Params.isSavingData
        measurementLogFid = fopen( Params.measurementLogFilename, 'w+' ) ;
        fwrite( measurementLogFid, Shim.Opt.Tracker.Data.p, 'double' ) ;
        fclose( measurementLogFid );
        
        measurementLogFid = fopen( Params.rawMeasurementLogFilename, 'w+' ) ;
        fwrite( measurementLogFid, rawMeasurementLog, 'double' ) ;
        fclose( measurementLogFid );

        sampleTimesFid = fopen( Params.sampleTimesFilename, 'w+' ) ;
        fwrite( sampleTimesFid, sampleTimes, 'double' ) ;
        fclose( sampleTimesFid );
        
        updateTimesFid = fopen( Params.updateTimesFilename, 'w+' ) ;
        fwrite( updateTimesFid, updateTimes, 'double' ) ;
        fclose( updateTimesFid );
    end

end

StopButton.Clear() ;



function p = clipvalue( p )
%CLIPVALUE
% 
% pClipped = CLIPVALUE( p )

if p < Params.minClipThreshold
    p = Params.minClipThreshold ;
elseif p > Params.maxClipThreshold ;
    p = Params.maxClipThreshold ;
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
