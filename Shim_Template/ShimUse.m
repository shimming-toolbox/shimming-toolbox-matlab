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
%       'Greg' or 'Rri' [default]
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
% Updated::20181030::ryan.topfer@polymtl.ca
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

if Shim.Params.isConfirmingDataLoadDir
    Shim.uiconfirmdataloaddir( ) ;
end

if ~Params.isConnectingToShim
    switch Shim.Params.shimSystem
        
        case 'Rriyan'
            Shim.Opt = ShimOptRri( Shim.Params ) ;

        case 'Greg'
            Shim.Opt = ShimOpt_Greg( Shim.Params ) ;
        
        case 'UnfPrisma'
            Shim.Opt = ShimOpt_IUGM_Prisma_fit( Shim.Params ) ;
        
        case 'Des' % i.e. shim design, a virtual shim system
            Shim.Opt = ShimOpt_Des( Shim.Params ) ;

        otherwise
            error([ Shim.Params.shimSystem 'is an invalid or unimplemented shim system. See HELP ShimUse().' ])

    end
else
    switch Shim.Params.shimSystem
        
        case 'Rriyan'
            Shim.Opt = ShimOptRri( Shim.Params ) ;
            Shim.Com = ShimComRri( ) ;

        case 'Greg'
            Shim.Opt = ShimOpt_Greg( Shim.Params ) ;
            Shim.Com = ShimCom_Greg( ) ;
        
        case 'UnfPrisma'
            Shim.Opt = ShimOpt_IUGM_Prisma_fit( Shim.Params ) ;
            Shim.Com = [] ;
        
        case 'Des' % i.e. shim design, a virtual shim system
            Shim.Opt = ShimOpt_Des( Shim.Params ) ;
            Shim.Com = [] ;

        otherwise
            error([ Shim.Params.shimSystem 'is an invalid or unimplemented shim system. See HELP ShimUse().' ])

    end
end

% if ~Params.isConnectingToShim ;
%     Shim.testshimconnection() ;
% end

if Shim.Params.isLoggingCommands
    diary( [Shim.Params.dataLoadDir Shim.Params.commandLogFilename] ) ;
end

Shim.Data.Aux.Tracker = cell( 1, 1, Shim.Params.nTrainingFrames + 1 ) ; % to record (in this order): (1) normal breathing; (2) breath-hold insp; (3) breath-hold exp

if strcmp( Shim.Params.uiMode, 'isGui' ) 
    ShimGui( Shim ) ;
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
% Wrapper to ProbeTracking.recordandplotphysiosignal( )
% 
% ACQUIRETRAININGDATA( Shim )


Params.isSavingData        = true;
Params.isForcingOverwrite  = true ;
Params.runTime             = 5*60 ; % [units: s], max runTime = 5 min.

if nargin < 2

    iTrainingFrame = ...
        input( ['Enter number corresponding to training acquisition/recording type: (1) Inspired breath-hold; (2) Expired breath-hold; (3) Breathing. \n'] ) ;

    if isempty(iTrainingFrame)
        iTrainingFrame = 3;
    end

end

Params.physioSignalFilename = [Shim.Params.dataLoadDir datestr(now,30) '-physioSignal-Training' num2str(iTrainingFrame) '.bin'] ;
Params.sampleTimesFilename  = [Shim.Params.dataLoadDir datestr(now,30) '-sampleTimes-Training' num2str(iTrainingFrame) '.bin'] ;

% -------
% begin physio tracking
if iTrainingFrame == 3
    Shim.Opt.Tracker.calibratelimiting( Params ) ;
else
    Shim.Opt.Tracker.recordandplotphysiosignal( Params ) ;
end

Shim.Params.physioSignalFilenames{iTrainingFrame} = Params.physioSignalFilename ;

Shim.Data.Aux.Tracker{iTrainingFrame} = Shim.Opt.Tracker.copy() ;


% TODO 
%
% sort field map image files and associate physio signal reading with corresponding FieldEval object
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

if nargin == 1 || isempty( imgDirectories )
    imgDirectories = [];
end

Shim.loadtrainingdata( imgDirectories ) ;
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

        if ~isempty( imgDirectories{ 1, iImg, iFrame } )
            
            nDicomSubDirs = size( imgDirectories( :, iImg, iFrame ), 1 ) ;
            
            for iDicomSubDir = 1 : nDicomSubDirs
                if ~isempty( imgDirectories{ iDicomSubDir, iImg, iFrame } )
                    Shim.Data.Img{ iDicomSubDir, iImg, iFrame }  = MaRdI( imgDirectories{ iDicomSubDir, iImg, iFrame } ) ;
                end
            
            end

        else % User prompt via GUI 
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
% Shim.Data.Img gains a 3rd column with cells containing the b0 field map
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

if ~myisfield( Shim.Params, 'isUserSelectionEnabled' ) || isempty( Shim.Params.isUserSelectionEnabled )
    isUserSelectionEnabled = true ; % default
else
    assert( islogical( Shim.Params.isUserSelectionEnabled ) ) ;
    isUserSelectionEnabled = Shim.Params.isUserSelectionEnabled ;
end

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

        % duration of breath-hold in number of samples 
        % (enables auto-estimation of the period corresponding to the breath-hold from
        % the complete physio log)
        nSamplesApnea = (1000*Shim.Data.Img{1,1,iFrame}.Hdr.MrProt.lTotalScanTimeSec)/Shim.Opt.Tracker.Specs.dt ;
        
        % ------
        % extract a single scalar
        Shim.Data.Img{ 1, 3, iFrame }.Aux.Tracker.Data.p = ...
            ProbeTracking.selectmedianmeasurement( Shim.Data.Img{ 1, 3, iFrame }.Aux.Tracker.Data.p, nSamplesApnea, isUserSelectionEnabled ) ;
    end
    
end

if ( Shim.Params.nTrainingFrames > 1 ) && ~isempty( Shim.Data.Aux.Tracker{ 1, 1, Shim.Params.nTrainingFrames + 1 } )
   
    % this is assumed to be a free breathing recording
    % i.e. trainingFrame 1 = inspired; 2 = expired.
    Shim.Params.pBreathing = Shim.Data.Aux.Tracker{ 1,1, Shim.Params.nTrainingFrames + 1 }.Data.p ;

    Shim.Params.pDc  = median( Shim.Params.pBreathing ) ;
    % TODO extract min + max amplitudes across ALL recordings instead?
    Shim.Params.pMin = min( Shim.Params.pBreathing ) ;
    Shim.Params.pMax = max( Shim.Params.pBreathing ) ;
else
    error('Free breathing measurements not available')
end
    
% -------
% derive t2star reliability weights:
t2sWeights = zeros( [ Shim.Data.Img{ 1, 1 , 1}.getgridsize( ) Shim.Params.nTrainingFrames ] ) ;
% t2sMasks: false wherever the t2sWeights are clearly unreliable (i.e. Mag(TE2) >= Mag(TE1) )
t2sMasks   = false( [ Shim.Data.Img{ 1, 1 , 1}.getgridsize( ) Shim.Params.nTrainingFrames ] ) ;

for iFrame = 1 : Shim.Params.nTrainingFrames 
    
    t2sWeights(:,:,:, iFrame) = Shim.Opt.derivedataweights( ...
        { Shim.Data.Img{ 1, 1, iFrame }  ; Shim.Data.Img{ 2, 1, iFrame } }, 5 ) ;

    t2sMasks(:,:,:, iFrame) = t2sWeights(:,:,:, iFrame) ~= 0 ;

end

% averaging
t2sWeights = sum( t2sMasks .* t2sWeights, 4 ) ./ sum( t2sMasks, 4 ) ;

% filtering
t2sWeights( ~( sum( t2sMasks, 4 ) ) ) = NaN ; % exclude from medfilt3( )
t2sWeights = medfilt3( t2sWeights, [3 3 1] ) ; 
assert( ~any( isnan( t2sWeights(:) ) ), 'Mapped t2sWeights include NaN values...' ) ;

% -------
% call SCT to segment spinal canal

if Shim.Params.isAutoSegmenting == true
    
    Shim.Params.cordVoi     = true( Shim.Data.Img{ 1, 1 , 1}.getgridsize( ) ) ;
    Shim.Params.dataWeights = ones( Shim.Data.Img{ 1, 1 , 1}.getgridsize( ) ) ;

    for iFrame = 1 : Shim.Params.nTrainingFrames 

        TmpParams.dataSaveDir       = [ Shim.Params.dataLoadDir ] ;
        TmpParams.isUsingPropsegCsf = true ;
        
        % call spinal cord toolbox 
        [cordVoi, sctWeights] = Shim.Data.Img{ 1, 1, iFrame }.segmentspinalcanal( TmpParams ) ;
        
        % retain the intersection with previous iteration:
        Shim.Params.cordVoi = Shim.Params.cordVoi & cordVoi ;
        
        Shim.Params.dataWeights = Shim.Params.dataWeights .* sctWeights  ;
    end
    
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

if ~myisfield( Shim.Params, 'dataLoadDir' )
    Params.dataSaveDir = './' ;
else
    Params.dataSaveDir = Shim.Params.dataLoadDir ;
end

DEFAULT_ISSAVINGDATA            = true ;
DEFAULT_ISFORCINGOVERWRITE      = false ;

DEFAULT_PHYSIOSIGNALFILENAME = [ Params.dataSaveDir datestr(now,30) '-physioSignal.bin' ] ;
DEFAULT_SAMPLETIMESFILENAME  = [ Params.dataSaveDir datestr(now,30) '-sampleTimes.bin' ] ;
DEFAULT_UPDATETIMESFILENAME  = [ Params.dataSaveDir datestr(now,30) '-updateTimes.bin' ] ;

DEFAULT_RUNTIME                 = 10*60 ; % [units: s]
DEFAULT_EXTRAPOLATIONORDER      = 0 ;
DEFAULT_EXTRAPOLATIONDELAY      = 0 ;

DEFAULT_ISFILTERINGMEASUREMENTS = true ; % Tracker measurements
DEFAULT_ISPLOTTINGINREALTIME    = true ;

DEFAULT_ISCLIPPING              = true ;
DEFAULT_MINCLIPTHRESHOLD        = 1 ;
DEFAULT_MAXCLIPTHRESHOLD        = 1000 ;

DEFAULT_TXDELAY                 = 100 ; % [units: ms]

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

[pathStr,name,ext] = fileparts( Params.physioSignalFilename ) ;
Params.rawPhysioSignalFilename = [pathStr '/' name '_raw' ext] ;

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

if  ~myisfield( Params, 'txDelay' ) || isempty(Params.txDelay)
    Params.txDelay = DEFAULT_TXDELAY ;
end

if ~myisfield( Shim.Params, 'pDc' ) || isempty( Shim.Params.pDc )
    error('Requires the DC offset of the respiratory signal: Shim.Params.pDc')
end


Params.nSamplesFilter  = 5 ; % Window length -> 5 samples * 10 ms/sample = 50 ms
Params.nSamplesHalfWindow = (Params.nSamplesFilter + 1)/2 - 1 ; 

Params.correctionOrder = 1;  % linear correction

[~,Params.filterWeights] = sgolay(Params.correctionOrder, Params.nSamplesFilter);   % Calculate S-G coefficients

% delay = tx/transmission delay plus the time shift involved from the filter
delay = ( Params.txDelay + Shim.Opt.Tracker.Specs.dt*(Params.nSamplesFilter-1)/2)/1000 ; % [units: s]


% ------- 
Shim.Opt.Tracker.clearrecording() ;
Shim.Opt.Tracker.Data.pRaw = 0;
Shim.Opt.Tracker.Data.p    = 0;


nSamples = Params.runTime / (Shim.Opt.Tracker.Specs.dt/1000) ;
iSample  = 0 ; 

iSamplesBetweenUpdates = 0;

if strcmp( class( Shim.Com ), 'ShimComRri' )
    nSamplesBetweenUpdates = ... % CHANGE THIS LINE FOR ACDC --REFERS TO RRI/MXD
        Shim.Com.Specs.Com.mxdUpdatePeriod/(Shim.Opt.Tracker.Specs.dt/1000)

    updatePeriod = Shim.Com.Specs.Com.mxdUpdatePeriod ;

else
   
   updatePeriod  = Shim.Com.Specs.Com.updatePeriod ;
    
   nSamplesBetweenUpdates = updatePeriod/( Shim.Opt.Tracker.Specs.dt/1000 ) ;

end

phys =[] ;

sampleTimes          = [] ; %
updateTimes          = [] ; % when shim updates occur

updateIndices        = [] ; % corresponding to when shim updates occur
iUpdate              = 0 ;

currentsNorm         = 0 ;

currentsIssued       = zeros( Shim.Com.Specs.Amp.nChannels, nSamples ) ;

if Params.isPlottingInRealTime
    
    close all

    % ------- 
    % figure
    figureHandle = figure('NumberTitle','off',...
        'Name','Respiration',...
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

Params.isPlottingInRealTime = false;
% -----
display('RAMPING SHIM CURRENTS TO DC') 
Shim.Com.setandrampallshims( Shim.Opt.Model.currents ) ;
pause(1)
Shim.Com.getallchanneloutputs

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
            
            Shim.Opt.Tracker.getupdate() ;
            % rawPhysioSignal(iSample) = Shim.Opt.Tracker.getupdate() ;

        end

        iUpdate = iUpdate + 1;        

        updateIndices(iUpdate) = iUpdate ;
        
        updateTimes(iUpdate) = updatePeriod * iUpdate ; 

        % % lowpass + delay correction of respiratory signal  
        % if Params.isFilteringMeasurements  && ( iSample > Params.nSamplesFilter )
        %         
        %      % 0th order corr (weighted avg)
        %      p = dot( Params.filterWeights(:,1), rawPhysioSignal( (iSample - Params.nSamplesFilter + 1) : iSample ) ) ;
        %
        %      % if Params.isPredictingMeasurement
        %      % 1st order corr
        %      p = p + delay*dot( Params.filterWeights(:,2), rawPhysioSignal( (iSample - Params.nSamplesFilter + 1) : iSample ) ) ;
        %      %end
        %
        %     Shim.Opt.Tracker.Data.p(iUpdate) = p ;
        %
        % else

            % Shim.Opt.Tracker.Data.p(iUpdate) = rawPhysioSignal(iSample) ;

        % end
        
         % if Params.isClipping
         %     Shim.Opt.Tracker.Data.p(iUpdate) = clipvalue( Shim.Opt.Tracker.Data.p(iUpdate) ) ;
         % end
         phys(iUpdate) = mean( Shim.Opt.Tracker.Data.p(end-1:end) ) ; % avg 2 samples
         pS            = phys(iUpdate) - Shim.Params.pDc ; % debiased measurement
         % pS = Shim.Opt.Tracker.Data.p(end) - Shim.Params.pDc ; % debiased measurement

        currents = Shim.Opt.computerealtimeupdate( pS ) ;

        
        % currents = limitcurrents( currents ) ;
        
        Shim.Com.setandloadallshims( currents ) ;
        
        currentsNorm = norm(currents) 
       
        % currentsIssued(:, iUpdate ) = currents ;

        if Params.isPlottingInRealTime
            set(plotHandle,'YData',phys,'XData',updateTimes);
            % set(plotHandle,'YData',Shim.Opt.Tracker.Data.p);
            % set(plotHandle,'YData',Shim.Opt.Tracker.Data.p,'XData',updateTimes);
        end

    end

    Shim.Opt.Tracker.stoptracking() ;
    Shim.Com.resetallshims() ;

    sampleTimes = (Shim.Opt.Tracker.Specs.dt/1000)*ones(size(Shim.Opt.Tracker.Data.p)) ;

    % ------- 
    if Params.isSavingData
        physioSignalFid = fopen( Params.physioSignalFilename, 'w+' ) ;
        fwrite( physioSignalFid, Shim.Opt.Tracker.Data.p, 'double' ) ;
        fclose( physioSignalFid );
        
        physioSignalFid = fopen( Params.rawPhysioSignalFilename, 'w+' ) ;
        fwrite( physioSignalFid, Shim.Opt.Tracker.Data.pRaw, 'double' ) ;
        fclose( physioSignalFid );

        sampleTimesFid = fopen( Params.sampleTimesFilename, 'w+' ) ;
        fwrite( sampleTimesFid, sampleTimes, 'double' ) ;
        fclose( sampleTimesFid );
        
        updateTimesFid = fopen( Params.updateTimesFilename, 'w+' ) ;
        fwrite( updateTimesFid, updateTimes, 'double' ) ;
        fclose( updateTimesFid );
        
        save( [Params.dataSaveDir datestr(now,30) '-currentsIssued'], currentsIssued ) ;
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

display('Select primary DATA folder.')
uiBoxTitle = ['Select primary DATA folder.'] ;

Shim.Params.dataLoadDir = [ uigetdir( Shim.Params.dataLoadDir, uiBoxTitle ) '/']; 

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
% DEFAULT_ISCONNECTINGTOSHIM = true ; 

DEFAULT_ISCONNECTINGTOSHIM = true ; 

DEFAULT_ISCONFIRMINGDATALOADDIR = true ; 

DEFAULT_UIMODE      = 'isCmdLine' ;% vs. 'isGui'

DEFAULT_SHIMSYSTEM  = 'Rriyan' ; 

DEFAULT_ISLOGGINGCOMMANDS  = true ; 
DEFAULT_COMMANDLOGFILENAME = ['commandLog_' datestr(now,30)] ;

% Re: 'Training' protocol
DEFAULT_NTRAININGFRAMES      = 2 ; % training time-points (e.g. inspired + expired field maps = 2 frames)
DEFAULT_ISGREPHASEDIFFERENCE = true ; % phase images are inter-echo phase difference images
DEFAULT_NECHOES              = 2 ;

if ~myisfield( Params, 'isConnectingToShim' ) || isempty( Params.isConnectingToShim ) 
   Params.isConnectingToShim = DEFAULT_ISCONNECTINGTOSHIM ;
end

if ~myisfield( Params, 'isConfirmingDataLoadDir' ) || isempty( Params.isConfirmingDataLoadDir ) 
   Params.isConfirmingDataLoadDir = DEFAULT_ISCONFIRMINGDATALOADDIR ;
end

if ~myisfield( Params, 'uiMode' ) || isempty( Params.uiMode ) 
   Params.uiMode = DEFAULT_UIMODE ;
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
% =========================================================================
function [] = writeadjvalidatebat( outputDir, shimValues )
%WRITEADJVALIDATEBAT
%
% [] = WRITEADJVALIDATEBAT( outputDir, shimValues )
% 
% Accepts a 8-element vector of shimValues in multipole units and writes
% a batch file to the outputDir to run the Siemens AdjValidate command for
% updating the gradient + 2nd order shim settings.
%    
% shimValues(1:3) are the x,y, and z gradient terms [units : micro-T/m]
% shimValues(4:8) are the 2nd order shim terms [units : micro-T/m^2]
%

%TODO
%   Move the method to ShimCom_IUGM_Prisma_fit() 

assert( length(shimValues) == 8 )

if exist( outputDir ) ~= 7
    warning(['Output directory ' outputDir ' does not exist. Creating it.']);
    mkdir( outputDir );
end

cmd = 'AdjValidate -shim -set -mp %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f \n';
fid = fopen([ outputDir '/adjshim_' datestr(now,30) '.bat' ],'w') ;
fprintf(fid, cmd, shimValues);
fclose(fid);

end
% =========================================================================
% =========================================================================
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
