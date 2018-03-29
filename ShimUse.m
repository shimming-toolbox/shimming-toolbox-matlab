classdef ShimUse < matlab.mixin.SetGet
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
%       .Data
%           .Img
%
%               A cell array of MaRdI-type image data objects used for shim
%               training (i.e. these are the GRE images): 
%               Rows correspond to echo number
%               1st column = Magnitude  
%               2nd column = Phase (e.g. inter-echo difference) 
%               
%               3rd dimension corresponds to "frame" (e.g.  .Data{:,:,1} may be
%               the "inspired" data, and .Data{:,:,2} the "expired"
%               acquisition.) 
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
% Updated::20180328::ryan.topfer@polymtl.ca
% =========================================================================

% =========================================================================
% *** TODO 
%
% ..... 
% UPDATEPARAMS()/setshimparam()...
%   Shim should have its own .Params field,
%   but at various stages, the user may wish to change the parameters...
% ..... 
% SAVESHIMEXPERIMENT()
%   saves data + params as they were during the experiment
%
%
% ..... 
% =========================================================================

properties   
    Com;
    Opt;
    Data;
    Params;
end

% =========================================================================
% =========================================================================
methods
% =========================================================================
function Shim = ShimUse( Params )
%SHIMUSE   

Shim.Opt      = [];
Shim.Com      = [];
Shim.Data.Img = [];
Shim.Data.Aux = [];
Shim.Params   = [];

if nargin < 1
    Params.dummy = [];
end

Shim.Params = ShimUse.assigndefaultparameters( Params ) ;

Shim.uiconfirmdataloaddir( ) ;

switch Shim.Params.shimSystem
    
    case 'Rri'
        Shim.Opt = ShimOptRri( Shim.Params ) ;
        Shim.Com = ShimComRri( ) ;

    case 'Acdc'
        Shim.Opt = ShimOpt_Acdc( Shim.Params ) ;
        Shim.Com = ShimComAcdc( ) ;
    
    case 'UnfPrisma'
        Shim.Opt = ShimOpt_IUGM_Prisma_fit( Shim.Params ) ;
        Shim.Com = [] ;

    otherwise
        error([ Shim.Params.shimSystem 'is an invalid or unimplemented shim system. See HELP ShimUse().' ])

end

if ~Params.isDebugging ;
    Shim.testshimconnection() ;
end

if Shim.Params.isLoggingCommands
    diary( [Shim.Params.dataLoadDir Shim.Params.commandLogFilename] ) ;
end

Shim.Data.Aux.Tracker = cell( 1, 1, Shim.Params.nTrainingFrames + 1 ) ; % to record (in this order): (1) normal breathing; (2) breath-hold insp; (3) breath-hold exp

if strcmp( Shim.Params.runMode, 'isGui' ) 
    ShimGUI( Shim ) ;
end

end
% =========================================================================
function [] = delete( Shim )
%DELETE  (custom helper function)
% 
% DELETE( Shim )
% 
% Destructor method: 
% calls Shim.Com.delete( ) and Shim.Opt.delete( ) 
% (also suspends command log (i.e. 'diary') ).

% SAVE PARAMS  
Params = Shim.Params ;
save( [Shim.Params.dataLoadDir datestr(now,30) '-Params' ], 'Params' ) ;

if Shim.Params.isLoggingCommands
    diary off
end

Shim.Opt.delete();
Shim.Com.delete();
clear Shim ;

end
% =========================================================================
function [] = acquiretrainingdata( Shim, iTrainingFrame )
%ACQUIRETRAININGDATA  
%
% Wrapper to ProbeTracking.recordandplotpressurelog()
% 
% ACQUIRETRAININGDATA( Shim )


Params.isSavingData        = true;
Params.isForcingOverwrite  = true ;
Params.runTime             = 3*60 ; % [units: s], max runTime = 3 min.

if nargin < 2

    iTrainingFrame = ...
        input( ['Enter number corresponding to training acquisition/recording type: (1) Inspired breath-hold; (2) Expired breath-hold; (3) Breathing. \n'] ) ;

    if isempty(iTrainingFrame)
        iTrainingFrame = 3;
    end

end

Params.pressureLogFilename = [Shim.Params.dataLoadDir datestr(now,30) '-pressureLog-Training' num2str(iTrainingFrame) '.bin'] ;
Params.sampleTimesFilename = [Shim.Params.dataLoadDir datestr(now,30) '-sampleTimes-Training' num2str(iTrainingFrame) '.bin'] ;

% -------
% begin (pressure) tracking
Shim.Opt.Tracker.recordandplotpressurelog( Params ) ;

Shim.Params.pressureLogFilenames{iTrainingFrame} = Params.pressureLogFilename ;

Shim.Data.Aux.Tracker{iTrainingFrame} = Shim.Opt.Tracker.copyinert() ;


% TODO 
%
% sort field map image files and associate pressure reading with corresponding FieldEval object
%
% each MaRdI obj. has Aux. field (i.e. Mag.Aux.Tracker )
% ---> Tracker recordings should be saved there 
%(but with a reference retained in Shim.Data.Aux, which may have more tracking recordings than images?)

end
% =========================================================================
function [Fields] = getprocessedfieldmaps( Shim )
%GETPROCESSEDFIELDMAPS

assert( size( Shim.Data.Img, 2 ) == 3 ) ;

Fields = cell(Shim.Params.nTrainingFrames, 1) ; 

for iFrame = 1 : Shim.Params.nTrainingFrames
 
    Fields{iFrame} = Shim.Data.Img{1, 3, iFrame} ;

end

end
% =========================================================================
function [] = loadandprocesstrainingdata( Shim, imgDirectories )
%LOADANDPROCESSTRAININGDATA  (Wrapper function)
% 
%  LOADANDPROCESSTRAININGDATA( Shim  ) 
% 
%  Wrapper function that calls: 
%   Shim.loadtrainingdata( ) 
%   Shim.processtrainingdata( )

Shim.loadtrainingdata( ) ;
Shim.processtrainingdata( ) ;

end
% =========================================================================
function [imgDirectories] = loadtrainingdata( Shim, imgDirectories )
%LOADTRAININGDATA  
% 
% imgDirectories = LOADTRAININGDATA( Shim ) 
% imgDirectories = LOADTRAININGDATA( Shim, imgDirectories ) 
%
% imgDirectories is a 2-column cell array containing the paths to the DICOM
% containing folders of the GRE training data 
%
% e.g. 
%   For dual-echo GRE_FIELD_MAPPING, 2 rows, 2 columns : 
%
%   { MAG_DIRECTORY(1st echo) } { PHASE_DIRECTORY(difference) } 
%   { MAG_DIRECTORY(2nd echo) } { [] } 
% 
% The 3rd dimension corresponds to the index of the training data acquisition
%
% e.g.  
%   imgDirectories{:,:,1} --> gre_field_mapping-INSPIRED  
%   imgDirectories{:,:,2} --> gre_field_mapping-EXPIRED 
%
% If imgDirectories is not provided, LOADTRAININGDATA() calls uigetdir() and
% the user selects the directories manually.

if nargin < 2 || isempty( imgDirectories )
    
    isManualSelection = true ;

    if Shim.Params.isGrePhaseDifference && ( Shim.Params.nEchoes == 2 )
        imgDirectories = cell( 1, 2, Shim.Params.nTrainingFrames ) ;
    else
        error('Unimplemented feature')
        imgDirectories = cell( Shim.Params.nEchoes, 2, Shim.Params.nTrainingFrames ) ;
    end
else
    assert( ( size( imgDirectories, 3 ) == Shim.Params.nTrainingFrames ) ...
            & size( imgDirectories, 2 ) == 2 ) ;
end

% check for expected number of echo subdirectories
errMsg = [ 'Did not find the expected number of echo_* subfolders. Check .Params.nEchoes is correct'] ;

for iFrame = 1 : Shim.Params.nTrainingFrames

    for iImg = 1 : 2 

        if iImg == 1
            imgType = 'MAGNITUDE' ;
        elseif iImg == 2
            imgType = 'PHASE' ;
        end

        if isempty( imgDirectories{ 1, iImg, iFrame } )

            uiBoxTitle = ['Select the ' imgType ' training data parent folder containing the echo_* DICOM subfolders.' ...
                ' (gre_field_mapping training data set ' num2str(iFrame) ' of ' num2str(Shim.Params.nTrainingFrames) ')'] ;
            
            ShimUse.customdisplay(uiBoxTitle) ;
            
            parentDir     = [ uigetdir( Shim.Params.dataLoadDir, uiBoxTitle ) '/']; 
            if ~parentDir % user cancelled
                return;
            end
            dicomSubDirs  = dir( [parentDir 'echo*/'] ) ;
            nDicomSubDirs = length( dicomSubDirs ) ;
             
            switch imgType % check for expected number of echo subdirectories
                case 'MAGNITUDE'
                    assert( nDicomSubDirs == Shim.Params.nEchoes, errMsg ) ;
                case 'PHASE'
                    if Shim.Params.isGrePhaseDifference
                        assert( ( nDicomSubDirs + 1 ) == Shim.Params.nEchoes, errMsg ) ;
                    else
                        assert( nDicomSubDirs == Shim.Params.nEchoes, errMsg ) ;
                    end
            end
            
            % load the images 
            for iDicomSubDir = 1 : nDicomSubDirs
                imgDirectories{ iDicomSubDir, iImg, iFrame } = [ parentDir dicomSubDirs(iDicomSubDir).name '/' ] ;
                Shim.Data.Img{ iDicomSubDir, iImg, iFrame }  = MaRdI( imgDirectories{ iDicomSubDir, iImg, iFrame } ) ;
            end 
        end

    end
end 

end
% =========================================================================
function [] = processtrainingdata( Shim )
%PROCESSTRAININGDATA  
%
% PROCESSTRAININGDATA( Shim )
% 
% Shim.Data.Img gains a 3rd column with cells containing the processed b0 field
% map
% i.e. 
%   Columns [1 | 2 | 3] = [MAG | PHASE | FIELD]
%
% If the corresponding [respiratory] tracker measurements for a given training
% frame (iFrame) are available in:
% 
% Shim.Data.Aux.Tracker{ iFrame } 
% 
% then the user is prompted to extract the median (most representative) value
% corresponding to the actual breath-hold (generally the full recording extends
% beyond the breath-hold itself). 
%
% This scalar value is then entered in: 
%
% Shim.Data.Img{1,3,iFrame}.Aux.Tracker.Data.p
Shim.Data.Img(1,3,:) = cell(1,1) ; % 3rd column for field maps

for iFrame = 1 : Shim.Params.nTrainingFrames

    ImgArray         = cell( 1, 2 ) ;
    ImgArray{ 1, 1 } = Shim.Data.Img{ 1, 1, iFrame } ; % mag
    ImgArray{ 1, 2 } = Shim.Data.Img{ 1, 2, iFrame } ; % phase
 
    % -------
    Shim.Data.Img{ 1, 3, iFrame } = FieldEval.mapfield( ImgArray, Shim.Params ) ;

    % link field measurements to appropriate tracker values 
    if ~isempty( Shim.Data.Aux.Tracker{ 1, 1, iFrame } )
        Shim.Data.Img{ 1, 3, iFrame }.Aux.Tracker = [] ;
        Shim.Data.Img{ 1, 3, iFrame }.Aux.Tracker = Shim.Data.Aux.Tracker{ 1, 1, iFrame }.copy() ;

        % ------
        % extract a single scalar
        Shim.Data.Img{ 1, 3, iFrame }.Aux.Tracker.Data.p = ...
            ProbeTracked.userselectmedianmeasurement( Shim.Data.Img{ 1, 3, iFrame }.Aux.Tracker.Data.p ) ;
    end
    
end

if ( Shim.Params.nTrainingFrames > 1 ) && ~isempty( Shim.Data.Aux.Tracker{ 1, 1, Shim.Params.nTrainingFrames + 1 } )
   
    % this is assumed to be a free breathing recording
    Shim.Params.pDc  = median( Shim.Data.Aux.Tracker{ 1,1, Shim.Params.nTrainingFrames + 1 }.Data.p ) ;
    % TODO extract min + max pressures across ALL recordings instead?
    Shim.Params.pMin = min( Shim.Data.Aux.Tracker{ 1,1, Shim.Params.nTrainingFrames + 1 }.Data.p ) ;
    Shim.Params.pMax = max( Shim.Data.Aux.Tracker{ 1,1, Shim.Params.nTrainingFrames + 1 }.Data.p ) ;
    
else
    error('Free breathing measurements not available')
end

if Shim.Params.isAutoSegmenting == true
    
    Shim.Params.cordVoi     = true( Shim.Data.Img{ 1, 1 , 1}.getgridsize( ) ) ;
    Shim.Params.dataWeights = zeros( Shim.Data.Img{ 1, 1 , 1}.getgridsize( ) ) ;

    % call spinal cord toolbox 
    for iFrame = 1 : 1 % Shim.Params.nTrainingFrames TODO : make it faster ? takes too long with multiple training frames
        TmpParams.dataSaveDir = [ Shim.Params.dataLoadDir ] ;
        TmpParams.isUsingPropsegCsf = true ;
        
        [cordVoi, sctWeights] = Shim.Data.Img{ 1, 1, iFrame }.segmentspinalcanal( TmpParams ) ;

        % retain the intersection with previous VOI:
        Shim.Params.cordVoi = Shim.Params.cordVoi & cordVoi ;
        
        t2sWeights = Shim.Opt.derivedataweights( ...
            { Shim.Data.Img{ 1, 1, iFrame }  ; Shim.Data.Img{ 2, 1, iFrame } }, 16, Shim.Params.cordVoi ) ;

        Shim.Params.dataWeights = Shim.Params.dataWeights + sctWeights + t2sWeights ;
    end
   
    % averaging
    Shim.Params.dataWeights = Shim.Params.dataWeights/Shim.Params.nTrainingFrames ;

end

% Field to be shimmed gets defined in Shim.Opt.Field
Shim.Opt.setoriginalfield( FieldEval.modelfield( Shim.Data.Img(1,3,:), Shim.Params ) ) ;

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

DEFAULT_ISFILTERINGMEASUREMENTS = true ; % Tracker measurements
DEFAULT_ISPLOTTINGINREALTIME    = true ;

DEFAULT_ISCLIPPING              = true ;
DEFAULT_MINCLIPTHRESHOLD        = 1 ;
DEFAULT_MAXCLIPTHRESHOLD        = 1000 ;


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
    Params.minClipThreshold  = Shim.Params.pMin ; % DEFAULT_MINCLIPTHRESHOLD ; 
end

if  ~myisfield( Params, 'maxClipThreshold' ) || isempty(Params.maxClipThreshold)
    Params.maxClipThreshold  = Shim.Params.pMax ; % DEFAULT_MAXCLIPTHRESHOLD ; 
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

if strcmp( class( Shim.Com ), 'ShimComRri' )
    nSamplesBetweenUpdates = ... % CHANGE THIS LINE FOR ACDC --REFERS TO RRI/MXD
        Shim.Com.Specs.Com.mxdUpdatePeriod/(Shim.Opt.Tracker.Specs.dt/1000)

    updatePeriod = Shim.Com.Specs.Com.mxdUpdatePeriod ;

else
   
   updatePeriod  = Shim.Com.Specs.updatePeriod ;
    
   nSamplesBetweenUpdates = updatePeriod/( Shim.Opt.Tracker.Specs.dt/1000 ) ;

end

rawMeasurementLog    = [] ;
Shim.Data.Tracker.p = [] ;

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

    ShimUse.customdisplay( 'Issuing real-time shim updates. Begin scanning!' )

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

            Shim.Opt.Tracker.Data.p(iUpdate) = p ;

        else

            Shim.Opt.Tracker.Data.p(iUpdate) = rawMeasurementLog(iSample) ;

        end
        
         if Params.isClipping
             Shim.Opt.Tracker.Data.p(iUpdate) = clipvalue( Shim.Opt.Tracker.Data.p(iUpdate) ) ;
         end

        currents = Shim.Opt.computerealtimeupdate(  ) ;
        
        % currents = limitcurrents( currents ) ;
        currentsNorm = norm(currents) 
        
        Shim.Com.setandloadallshims( currents ) 
        
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
%         ShimUse.customdisplay('CLIPPING');
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
    ShimUse.customdisplay( msg ) ;
else
    msg = [ '-----\n'...
            'Communication to shim amplifier failed. Check device connections.' ... 
            '\n -----'] ;
    ShimUse.customdisplay( msg ) ;
end

end
% =========================================================================
function [] = uiconfirmdataloaddir( Shim )
%UICONFIRMDATALOADDIR
% 
% isAckReceived = TESTSHIMCONNECTION( Shim )
% 
% Queries shim amp for response (1 or 0)

uiBoxTitle = ['Select primary DATA folder.'] ;

Shim.Params.dataLoadDir  = [ uigetdir( Shim.Params.dataLoadDir, uiBoxTitle ) '/']; 

end
% =========================================================================

end
% =========================================================================
% =========================================================================
methods(Static) 
% =========================================================================
function [ Params ] = assigndefaultparameters( Params )
%ASSIGNDEFAULTPARAMETERS  
% 
% Params = ASSIGNDEFAULTPARAMETERS( Params )
% 
% Add default parameters fields to Params without replacing values (unless empty)
%
% DEFAULT_SHIMSYSTEM = 'Rri' ; 
% DEFAULT_RUNMODE    = 'isCmdLine' ;% vs. 'isGui'
%
% Re: 'Training' protocol :
% 
% DEFAULT_NECHOES = 2 ; 
%   Current protocol consists of 2 gre_field_mapping acquisitions, hence:
% DEFAULT_NTRAININGFRAMES = 2 ; % training time-points (e.g. inspired + expired field maps = 2 frames)
%   and the Siemens default is to only save the inter-echo phase difference:
% DEFAULT_ISGREPHASEDIFFERENCE = true ;
%
% DEFAULT_ISDEBUGGING = false ; 

DEFAULT_ISDEBUGGING = false ; 

DEFAULT_RUNMODE     = 'isCmdLine' ;% vs. 'isGui'

DEFAULT_SHIMSYSTEM  = 'Rri' ; 

DEFAULT_ISLOGGINGCOMMANDS  = true ; 
DEFAULT_COMMANDLOGFILENAME = ['commandLog_' datestr(now,30)] ;

% Re: 'Training' protocol
DEFAULT_NTRAININGFRAMES      = 2 ; % training time-points (e.g. inspired + expired field maps = 2 frames)
DEFAULT_ISGREPHASEDIFFERENCE = true ; % phase images are inter-echo phase difference images
DEFAULT_NECHOES              = 2 ;

if ~myisfield( Params, 'isDebugging' ) || isempty( Params.isDebugging ) 
   Params.isDebugging = DEFAULT_ISDEBUGGING ;
end

if ~myisfield( Params, 'runMode' ) || isempty( Params.runMode ) 
   Params.runMode = DEFAULT_RUNMODE ;
end

if ~myisfield( Params, 'isLoggingCommands' ) || isempty( Params.isLoggingCommands ) 
   Params.isLoggingCommands = DEFAULT_ISLOGGINGCOMMANDS ;
end

if ~myisfield( Params, 'commandLogFilename' ) || isempty( Params.commandLogFilename ) 
   Params.commandLogFilename = DEFAULT_COMMANDLOGFILENAME ;
end

if ~myisfield( Params, 'shimSystem' ) || isempty( Params.shimSystem ) 
   Params.shimSystem = DEFAULT_SHIMSYSTEM ;
end

if ~myisfield( Params, 'nTrainingFrames' ) || isempty( Params.nTrainingFrames ) 
   Params.nTrainingFrames = DEFAULT_NTRAININGFRAMES ;
end

if ~myisfield( Params, 'isGrePhaseDifference' ) || isempty( Params.isGrePhaseDifference ) 
   Params.isGrePhaseDifference = DEFAULT_ISGREPHASEDIFFERENCE ;
end

if ~myisfield( Params, 'nEchoes' ) || isempty( Params.nEchoes ) 
   Params.nEchoes = DEFAULT_NECHOES ;
end



end
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
function [] = customdisplay( msg )
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
