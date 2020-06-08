classdef ShimUse < matlab.mixin.SetGet
%SHIMUSE - Shim Use
%
% A high-level user interface to operate the shim system. 
% 
% .......
%   
% Usage
%
% Shims = ShimUse(  )
% Shims = ShimUse( ShimParams )
% Shims = ShimUse( ShimParams, ProbeParams )
%
% Shims contains public properties:
%
%    .Aux
%       Object of type ProbeTracking (for tracking respiration (proxy for dynamic delta B0))
%
%    .Com
%       Object of type ShimCom (for communication with shim amplifiers)
%
%    .Data
%       .Aux
%           Cell array of 'inert' ProbeTracking objects (i.e. copies of what Shims.Aux has already recorded).
%
%       .Img  
%           Cell array of MaRdI-type objects (i.e. 'training data' for shim optimization)
%
%    .Opt
%       Object of type ShimOpt (for optimizing shim currents)
%
%    .Params
%       Struct of misellaneous parameters 
%
% =========================================================================
% Author::ryan.topfer@polymtl.ca
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
    Aux;
    Com;
    Data;
    Opt;
    Params;
end

% =========================================================================
% =========================================================================
methods
% =========================================================================
function Shim = ShimUse( ShimParams, ProbeParams )
%SHIMUSE   

Shim.Opt      = [];
Shim.Com      = [];
Shim.Data.Img = [];
Shim.Data.Aux = [];
Shim.Params   = [];

if nargin < 2
    ProbeParams.dummy = [];
    
    if nargin < 1
        ShimParams.dummy  = [];
    end
end

Shim.Params = ShimUse.assigndefaultparameters( ShimParams ) ;

Shim.Aux    = ProbeTracking( ProbeParams )

switch Shim.Params.trainingMode
    case 'breath-hold'
        % to record (in this order): (1) normal breathing; (2) breath-hold insp; (3) breath-hold exp
        Shim.Data.Aux = cell( Shim.Params.nTrainingSeries + 1, 1 ) ; 
        Shim.Params.physioSignalFilenames = cell( Shim.Params.nTrainingSeries + 1, 1 ) ;
    case 'free breathing'
        Shim.Data.Aux = cell( Shim.Params.nTrainingSeries, 1 ) ; 
        Shim.Params.physioSignalFilenames = cell( Shim.Params.nTrainingSeries, 1 ) ;
end

if Shim.Params.isConfirmingDataLoadDir
    Shim.uiconfirmdataloaddir( ) ;
end

if ~Shim.Params.isConnectingToShim
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
%
% if Shim.Params.isLoggingCommands
%     diary( [Shim.Params.dataLoadDir Shim.Params.commandLogFilename] ) ;
% end
%
%
% if strcmp( Shim.Params.uiMode, 'isGui' ) 
%     ShimGui( Shim ) ;
% end

end
% =========================================================================
function [] = delete( Shim )
%DELETE  (custom helper function)
% 
% DELETE( Shim )
% 
% Destructor method: 
% calls Shim.Aux.delete(), Shim.Com.delete( ), and  Shim.Opt.delete( ) 
% (also suspends command log (i.e. 'diary') ).

% SAVE PARAMS  
Params = Shim.Params ;
save( [Shim.Params.dataLoadDir datestr(now,30) '-Params' ], 'Params' ) ;

if Shim.Params.isLoggingCommands
    diary off
end

Shim.Aux.delete();
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

switch Shim.Params.trainingMode
    case 'breath-hold'
        if nargin < 2
            
            disp('Enter number corresponding to training acquisition/recording type: ')
            iTrainingFrame = input( ['(1) Inspired breath-hold; (2) Expired breath-hold; (3) Breathing. \n'] ) ;

            if isempty(iTrainingFrame)
                iTrainingFrame = 3 ;
            end
        end
        
    case 'free breathing'
        iTrainingFrame = 1 ;

end

Params.physioSignalFilename = [ Shim.Params.dataLoadDir datestr(now,30) '-physioSignal-Training' num2str(iTrainingFrame) ] ;

% -------
% begin physio tracking
% if iTrainingFrame == 3
%     Shim.Aux.calibratelimiting( Params ) ;
% else
    Shim.Aux.recordandplotphysiosignal( Params ) ;
% end

Shim.Params.physioSignalFilenames{iTrainingFrame} = Params.physioSignalFilename ;

Shim.Data.Aux{iTrainingFrame} = Shim.Aux.copy() ;


% TODO 
%
% sort field map image files and associate physio signal reading with corresponding FieldEval object
%
% each MaRdI obj. has Aux. field (i.e. Mag.Aux )
% ---> Aux recordings should be saved there 
%(but with a reference retained in Shim.Data.Aux, which may have more tracking recordings than images?)

end
% =========================================================================
function [Fields] = getprocessedfieldmaps( Shim )
%GETPROCESSEDFIELDMAPS

assert( size( Shim.Data.Img, 2 ) == 3, 'Could not find processed field maps in Shim.Data.Img' ) ;

Fields = cell(Shim.Params.nTrainingSeries, 1) ; 

for iSeries = 1 : Shim.Params.nTrainingSeries
    Fields{iSeries} = Shim.Data.Img{1, 3, iSeries}.copy() ;
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

if nargin == 1 
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
% containing folders of the GRE training data. The 3rd dimension corresponds to
% the index of the training data acquisition
%
% If imgDirectories is not provided, LOADTRAININGDATA() calls uigetdir() and
% the user selects the directories manually.
%
% e.g. 
%   For dual-echo GRE_FIELD_MAPPING, 2 rows, 2 columns : 
%
%   { MAG_DIRECTORY(1st echo) } { PHASE_DIRECTORY(difference) } 
%   { MAG_DIRECTORY(2nd echo) } { [] } 
%
%   imgDirectories{:,:,1} --> gre_field_mapping-INSPIRED  
%   imgDirectories{:,:,2} --> gre_field_mapping-EXPIRED 

if nargin == 2 && ~isempty( imgDirectories )
    
    assert( ( size( imgDirectories, 3 ) == Shim.Params.nTrainingSeries ) ...
            & size( imgDirectories, 2 ) == 2 ) ;

    for iSeries = 1 : Shim.Params.nTrainingSeries
        for iImg = 1 : 2 
            if ~isempty( imgDirectories{ 1, iImg, iSeries } )
                
                nDicomSubDirs = size( imgDirectories( :, iImg, iSeries ), 1 ) ;
                
                for iDicomSubDir = 1 : nDicomSubDirs
                    if ~isempty( imgDirectories{ iDicomSubDir, iImg, iSeries } )
                        Shim.Data.Img{ iDicomSubDir, iImg, iSeries }  = MaRdI( imgDirectories{ iDicomSubDir, iImg, iSeries } ) ;
                    end
                end
            end
        end
    end
    return;
end

% else, manually select the directories using GUI
imgDirectories = cell( 1, 2, Shim.Params.nTrainingSeries ) ;
imgType        = { 'MAGNITUDE' ; 'PHASE' } ;

errMsg         = [ 'Did not find the expected number of echo_* subfolders. Ensure the chosen directory and Shim.Params.trainingMode are correct.'] ;

for iSeries = 1 : Shim.Params.nTrainingSeries

    for iImg = 1 : 2 

        % if strcmp( Shim.Params.trainingMode, 'free breathing' ) && ( iSeries == 1 )
        %     uiBoxTitle = ['Select the ' imgType{iImg} ' training data folder containing the DICOM images.' ...
        %             ' (gre_field_mapping training data set ' num2str(iSeries) ' of ' num2str(Shim.Params.nTrainingSeries) ')'] ;
        %
        % else
            % uiBoxTitle = ['Select the ' imgType{iImg} ' training data parent folder containing the echo_* DICOM subfolders.' ...
            %         ' (gre_field_mapping training data set ' num2str(iSeries) ' of ' num2str(Shim.Params.nTrainingSeries) ')'] ;
        % end
        uiBoxTitle = ['Select the ' imgType{iImg} ' training data parent folder containing the echo_* DICOM subfolders.' ...
                    ' (gre_field_mapping training data set ' num2str(iSeries) ' of ' num2str(Shim.Params.nTrainingSeries) ')'] ;
        
        ShimUse.customdisplay(uiBoxTitle) ;
        
        parentDir  = [ uigetdir( Shim.Params.dataLoadDir, uiBoxTitle ) '/']; 

        if ~parentDir % user cancelled
            return;
        end
        
        dicomSubDirs  = dir( [parentDir 'echo*/'] ) ;
        nDicomSubDirs = length( dicomSubDirs ) ;
        
        if nDicomSubDirs == 0 
            error( errMsg ) ;
            % assert( strcmp( Shim.Params.trainingMode, 'free breathing' ) && (iSeries == 2), errMsg ) ;
            %
            % imgDirectories{ 1, iImg, iSeries } = parentDir ;
            % Shim.Data.Img{ 1, iImg, iSeries }  = MaRdI( parentDir ) ;

        else
            for iDicomSubDir = 1 : nDicomSubDirs
       
                imgDirectories{ iDicomSubDir, iImg, iSeries } = [ parentDir dicomSubDirs(iDicomSubDir).name '/' ] ;
                Shim.Data.Img{ iDicomSubDir, iImg, iSeries }  = MaRdI( imgDirectories{ iDicomSubDir, iImg, iSeries } ) ;
            
                switch imgType{iImg} % check for expected number of echo subdirectories
                    
                    case 'MAGNITUDE'
                        assert( nDicomSubDirs == Shim.Data.Img{ 1, iImg, iSeries }.Hdr.MrProt.lContrasts, errMsg ) ;
                    case 'PHASE'
                        % if phase subdirectories are phase difference images,
                        % then there will be one fewer subdir than the number
                        % of echoes (contrasts).  if the subdirectories are the
                        % phase images themselves, then the number of
                        % subdirectories == nEchoes.  therefore, this assertion
                        % is generic:
                        assert( ( nDicomSubDirs + 1) >= Shim.Data.Img{ 1, iImg, iSeries }.Hdr.MrProt.lContrasts, errMsg ) ;
                end
            end 
        end

    end
end 
    

end
% =========================================================================
function [] = associatetrackerdata( Shim, iSeries )
%ASSOCIATETRACKERDATA
%
error('DEPRECATED')
Shim.Data.Img{ 1, 3, iSeries }.Aux = Shim.Data.Aux{ iSeries }.copy() ;

% number of acquisitions
nAcq = size( Shim.Data.Img{ 1, 3, iSeries }.img, 4 ) ;

if  nAcq == 1
    % This is a single field map -- assumed to correspond to a single tracker measurement:

    % number of samples corresponding to breath-hold
    % (enables auto-estimation of the time window corresponding to the breath-hold):
    nSamplesApnea = (1000*Shim.Data.Img{1,1,iSeries}.Hdr.MrProt.lTotalScanTimeSec)/Shim.Aux.Specs.dt ;

    % extract a single scalar
    Shim.Data.Img{ 1, 3, iSeries }.Aux.Data.p = ...
        ProbeTracking.selectmedianmeasurement( Shim.Data.Img{ 1, 3, iSeries }.Aux.Data.p, nSamplesApnea ) ;

else
    % ------
    % Correct for system delays before associating tracker measurements
    Field = Shim.Data.Img{ 1, 3, iSeries } ;

    tAcq = Field.getacquisitiontime() ;
    tAcq = tAcq - tAcq(1) ;

    % time between images 
    dtAcq = median( diff( tAcq ) ) ;

    pt      = Shim.Data.Img{ 1, 3, iSeries }.Aux.Data.p ;
    tGlobal = Shim.Data.Img{ 1, 3, iSeries }.Aux.Data.t ;

    % time interval between respiratory samples:
    dt = 0.001*round( ( tGlobal(end) - tGlobal(1) )/length(tGlobal) ) ; % [units: s]

    tGlobal = 0.001*tGlobal ; % convert to [units: s]
    
    gridSize = Field.getgridsize() ;
    
    % % ------
    % % bandpass filter images?
    % 
    % frequencyCutoff = [ 0 0.6 ] ;
    %
    % for iRow = 1 : gridSize(1) 
    %     for iColumn = 1 : gridSize(2)  
    %         for iSlice = 1 : gridSize(3)
    %         Field.img( iRow, iColumn, 1, : ) = ...
    %             bpfilt( squeeze( Field.img(iRow, iColumn, iSlice,:) ), frequencyCutoff(1), frequencyCutoff(2), 1/dtAcq, false ) ;
    %         end
    %     end
    % end
        
    % ------
    % Interpolate field across time to match the interval of the physio recording 
    tInterp = [ 0: dt : tAcq(end) ]' ;

    imgInterp = zeros( [ gridSize length(tInterp) ] ) ;

    for iRow = 1 : gridSize(1) 
        for iColumn = 1 : gridSize(2)  
            for iSlice = 1 : gridSize(3)
                imgInterp( iRow, iColumn, iSlice, : ) = interp1( tAcq, squeeze( Field.img( iRow, iColumn, iSlice, :) ), tInterp, 'linear' ) ; 
            end
        end
    end
    
    % Assume the entire dB0(tInterp) series occurs somewhere within the tGlobal window.
    % given that we assign tInterp(1) = 0, 
    % the delay time (tInterp' = tGlobal' - t_delay) is necessarily positive, 
    % and less than max(tGlobal) - max(tInterp)

    maxLag = round(( max( tGlobal ) - max( tInterp ) )/dt) ;

    T_delay = nan( gridSize ) ;
    C       = zeros( [gridSize 1] ) ;
    
    % reliable region:
    mask  = ( sum( Field.Hdr.MaskingImage, 4 ) == nAcq ) ;

    ShimUse.customdisplay( 'Estimating probe-to-image delay...')
    % ------
    % X-corr
    for iRow = 1 : Field.Hdr.Rows 
        for iColumn = 1 : Field.Hdr.Columns  
             if mask( iRow, iColumn )

                [c, lags] = xcorr( pt, imgInterp(iRow,iColumn,:), maxLag ) ;

                % truncate to positive delay times
                c    = c( lags >=0 ) ;
                lags = lags( lags >=0 ) ;
                
                [maxC, iMax] = max(abs(c)) ;
                t_xcorr = dt*lags ;
                t_delay = t_xcorr( iMax ) ;

                C(iRow, iColumn, 1, 1:length(c)) = c ;
                T_delay(iRow, iColumn) = t_delay ;
            end
        end
    end

    % the median should be more robust than the average for determining the delay
    % once noise + uncorrelated voxels are taken into consideration
    t_medianDelay = median( T_delay( mask ) ) ;

    ShimUse.customdisplay( ['Estimated delay: ' num2str( t_medianDelay )] ) ;
        
    response = input(['Is the current estimate satisfactory? ' ...
            '0 to manually specify delay; 1 (or enter) to accept & continue: ']) ;

    if ~response
        
        figure
        imagesc( mask .* median( imgInterp, 4 ) ) ;
        title('Median field') ;
        ShimUse.customdisplay( 'Choose a voxel to plot its field time-course: ' ) ;

        isVoxelInReliableRegion = false ;

        while ~isVoxelInReliableRegion
            iRow = input( ['Enter the voxel image row: ' ] ) 
            iCol = input( ['Enter the voxel image column: ' ] )
            
            isVoxelInReliableRegion = mask( iRow, iCol ) ;
            
            if ~isVoxelInReliableRegion
                ShimUse.customdisplay( 'Selected voxel lies outside the reliable region. Choose again' ) ;
            end
        end
        
        figure; 
        subplot(121); plot( tGlobal, pt ); 
        xlabel('Time (s)') ; ylabel('Amplitude (a.u.)') ; 
        title('Physio signal') ;
        
        subplot(122); plot( tInterp, squeeze(imgInterp(iRow,iCol,:)) ) ;
        xlabel('Time (s)') ; ylabel('Field (Hz)') ; 
        title(['Field at voxel [' num2str(iRow) ',' num2str(iCol) ']' ]) ;

        isDelayValid = false ;
        while ~isDelayValid
            t_delay = input(['Enter the delay (between 0-' num2str(maxLag) ' s) in units of seconds: '] ) 
            
            isDelayValid = ( t_delay >= 0 ) & ( t_delay <= maxLag ) ;
            
            if ~isDelayValid
                ShimUse.customdisplay( 'Given delay exceeds the expected range. Choose again' ) ;
            end
        end 
    else
        t_delay = t_medianDelay ;
    end

    %  % % figure; subplot(121); plot( squeeze( Field.img(50,40,1,:) ) ); subplot(122); plot( Field.Aux.Data.p )
    %  % figure; subplot(121); plot( tGlobal, pt ); subplot(122); plot( tAcq, squeeze( Field.img(50,30,:)) )
    %  figure; subplot(121); plot( tGlobal, pt ); subplot(122); plot( tAcq, squeeze( Field.img(50,40,:)) )
    %  hold on; plot( tInterp, squeeze(imgInterp(50,40,:)) )
    %  % figure; subplot(121); plot( tGlobal, pt ); subplot(122); plot( tAcq, squeeze( Field.img(40,40,:)) )
    %  % figure; subplot(121); plot( tGlobal, pt ); subplot(122); plot( tAcq, squeeze( Field.img(30,40,:)) )
    %
    %  meanMag = repmat(mean(Shim.Data.Img{1,1,1}.img,4) ,[1 1 1 200]) ;
    %  stdMag  = repmat(std(Shim.Data.Img{1,1,1}.img,[],4) ,[1 1 1 200]) ;
    % zScore = ( Shim.Data.Img{1,1,1}.img - meanMag ) ./ stdMag ;
    % zMask = ( sum( abs(zScore) < 1, 4 ) >=199 ) ;    
    %
    % t_medianDelay = median( T_delay( mask & zMask ) ) 

    [~,iT] = min( abs( tGlobal - t_delay ) ) ;

    % Finally, interpolate pt to the acquired image times:
    tp = tGlobal( iT:iT+length(tInterp)-1 ) -tGlobal(iT) ;

    % NOTE: tp(end) can be slightly < tAcq(end)  
    % Using 'extrap' to avoid an interpolated p terminating with NaN's
    Shim.Data.Img{ 1, 3, iSeries }.Aux.Data.p = interp1( tp, pt( iT:iT+length(tInterp)-1 ), tAcq, 'linear', 'extrap' ) ;
    Shim.Data.Img{ 1, 3, iSeries }.Aux.Data.t = tAcq ;

    if any( isnan(Shim.Data.Img{ 1, 3, iSeries }.Aux.Data.p) )
        warning('Interpolated probe recording contains nonnumeric values. Check inputs. (E.g. Shim.Data.Img{ }: Entries should be ordered by acquisition time)' ) ;
    end

    % used in nonlinear shim optimization:
    Shim.Params.pMax     = max( Shim.Data.Img{ 1, 3, iSeries }.Aux.Data.p ) ;
    Shim.Params.pMin     = min( Shim.Data.Img{ 1, 3, iSeries }.Aux.Data.p ) ;
    
    % for posterity (might be useful?):
    Shim.Params.imgDelay = t_delay ;    

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
% frame (iSeries) are available in:
% 
% Shim.Data.Aux{ iSeries } 
% 
% then the user is prompted to extract the median (most representative) value
% corresponding to the actual breath-hold (generally the full recording extends
% beyond the breath-hold itself). 
%
% This scalar value is then entered in: 
%
% Shim.Data.Img{1,3,iSeries}.Aux.Data.p

Shim.Data.Img(1,3,:) = cell(1,1) ; % 3rd column for field maps

for iSeries = 1 : Shim.Params.nTrainingSeries

    Field = FieldEval.mapfield( Shim.Data.Img{ 1, 1, iSeries }, Shim.Data.Img{ 1, 2, iSeries }, Shim.Params ) ;

    % associate field measurements to tracker values 
    if ~isempty( Shim.Data.Aux{ iSeries } )
        Field.associateaux( Shim.Data.Aux{ iSeries }, Shim.Params ) ;
    end
    
    Shim.Data.Img{ 1, 3, iSeries } = Field.copy() ;
end

switch Shim.Params.trainingMode
    case 'free breathing'
        % Shim.Opt.setoriginalfield( FieldEval.modelfield( Shim.Data.Img{1,3,1}, Shim.Params ) ) ;
        return;
    
    case 'breath-hold'
        if ~isempty( Shim.Data.Aux{ Shim.Params.nTrainingSeries + 1 } )
           
            % this is assumed to be a free breathing recording
            % i.e. trainingFrame 1 = inspired; 2 = expired.
            Shim.Params.pBreathing = Shim.Data.Aux{ Shim.Params.nTrainingSeries + 1 }.Data.p ;

            Shim.Params.pDc  = median( Shim.Params.pBreathing ) ;
            % TODO extract min + max amplitudes across ALL recordings instead?
            Shim.Params.pMin = min( Shim.Params.pBreathing ) ;
            Shim.Params.pMax = max( Shim.Params.pBreathing ) ;
        else
            error('Free breathing measurements not available')
        end
end

% -------
% derive t2star reliability weights:
Shim.Params.t2sWeights = zeros( [ Shim.Data.Img{ 1, 1 , 1}.getgridsize( ) Shim.Params.nTrainingSeries ] ) ;
% t2sMasks: false wherever the Shim.Params.t2sWeights are clearly unreliable (i.e. Mag(TE2) >= Mag(TE1) )
t2sMasks   = false( [ Shim.Data.Img{ 1, 1 , 1}.getgridsize( ) Shim.Params.nTrainingSeries ] ) ;

for iSeries = 1 : Shim.Params.nTrainingSeries 
    
    Shim.Params.t2sWeights(:,:,:, iSeries) = Shim.Opt.derivedataweights( ...
        { Shim.Data.Img{ 1, 1, iSeries }  ; Shim.Data.Img{ 2, 1, iSeries } }, 5 ) ;

    t2sMasks(:,:,:, iSeries) = Shim.Params.t2sWeights(:,:,:, iSeries) ~= 0 ;

end

% averaging
Shim.Params.t2sWeights = sum( t2sMasks .* Shim.Params.t2sWeights, 4 ) ./ sum( t2sMasks, 4 ) ;

% filtering
Shim.Params.t2sWeights( ~( sum( t2sMasks, 4 ) ) ) = NaN ; % exclude from medfilt3( )
Shim.Params.t2sWeights = medfilt3( Shim.Params.t2sWeights, [3 3 1] ) ; 
Shim.Params.t2sWeights( isnan( Shim.Params.t2sWeights ) ) = 0 ;

% -------
% call SCT to segment spinal canal

if Shim.Params.isAutoSegmenting == true
    
    Shim.Params.cordVoi     = true( Shim.Data.Img{ 1, 1 , 1}.getgridsize( ) ) ;
    Shim.Params.dataWeights = ones( Shim.Data.Img{ 1, 1 , 1}.getgridsize( ) ) ;

    for iSeries = 1 : Shim.Params.nTrainingSeries 

        TmpParams.dataSaveDir       = [ Shim.Params.dataLoadDir ] ;
        TmpParams.isUsingPropsegCsf = true ;
        
        % call spinal cord toolbox 
        [cordVoi, sctWeights] = Shim.Data.Img{ 1, 1, iSeries }.segmentspinalcanal( TmpParams ) ;
        
        % retain the intersection with previous iteration:
        Shim.Params.cordVoi = Shim.Params.cordVoi & cordVoi ;
        
        Shim.Params.dataWeights = Shim.Params.dataWeights .* sctWeights  ;
    end
    
    Shim.Params.dataWeights = Shim.Params.dataWeights/Shim.Params.nTrainingSeries ;

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

% assert( strcmp( Shim.Aux.ComPort.Status, 'closed' ), ...
%     'Error: Serial port is open/in use.' );

if ~myisfield( Shim.Params, 'dataLoadDir' )
    Params.dataSaveDir = './' ;
else
    Params.dataSaveDir = Shim.Params.dataLoadDir ;
end

DEFAULT_ISSAVINGDATA            = true ;
DEFAULT_ISFORCINGOVERWRITE      = false ;

DEFAULT_PHYSIOSIGNALFILENAME = [ Params.dataSaveDir datestr(now,30) '-physioSignal-RTShimming' ] ;
DEFAULT_SAMPLETIMESFILENAME  = [ Params.dataSaveDir datestr(now,30) '-sampleTimes' ] ;
DEFAULT_UPDATETIMESFILENAME  = [ Params.dataSaveDir datestr(now,30) '-updateTimes' ] ;

DEFAULT_RUNTIME                 = 10*60 ; % [units: s]

DEFAULT_ISPLOTTINGINREALTIME    = true ;

DEFAULT_ISCLIPPING              = true ;

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

if  ~myisfield( Params, 'isPlottingInRealTime' ) || isempty(Params.isPlottingInRealTime)
    Params.isPlottingInRealTime  = DEFAULT_ISPLOTTINGINREALTIME ; 
end

if  ~myisfield( Params, 'maxCurrents' ) || isempty(Params.maxCurrents)
    Params.maxCurrents  = ones( Shim.Com.Specs.Amp.nActiveChannels, 1) ; 
end

if  ~myisfield( Params, 'minCurrents' ) || isempty(Params.minCurrents)
    Params.minCurrents  = -ones( Shim.Com.Specs.Amp.nActiveChannels, 1) ; 
end

if  ~myisfield( Params, 'isClipping' ) || isempty(Params.isClipping)
    Params.isClipping  = DEFAULT_ISCLIPPING ; 
end

if  ~myisfield( Params, 'txDelay' ) || isempty(Params.txDelay)
    Params.txDelay = DEFAULT_TXDELAY ;
end

% ------- 
Shim.Aux.clearrecording() ;

nSamples = Params.runTime / (Shim.Aux.Specs.dt/1000) ;
iSample  = 0 ; 

iSamplesBetweenUpdates = 0;

if strcmp( class( Shim.Com ), 'ShimComRri' )
    nSamplesBetweenUpdates = ... 
        Shim.Com.Specs.Com.mxdUpdatePeriod/(Shim.Aux.Specs.dt/1000)

    updatePeriod = Shim.Com.Specs.Com.mxdUpdatePeriod ;

else
   
   updatePeriod  = Shim.Com.Specs.Com.updatePeriod ;
    
   nSamplesBetweenUpdates = updatePeriod/( Shim.Aux.Specs.dt/1000 ) ;

end

sampleTimes          = [] ; %
updateTimes          = [] ; % when shim updates occur

updateIndices        = [] ; % corresponding to when shim updates occur
iUpdate              = 0 ;

currentsNorm         = 0 ;
currentsIssued       = zeros( Shim.Com.Specs.Amp.nChannels, 1 ) ;

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
        
    title('Physio signal','FontSize',15,'Color',[1 1 0]);
    xlabel('Time (s)','FontWeight','bold','FontSize',14,'Color',[1 1 0]);
    ylabel('Amplitude','FontWeight','bold','FontSize',14,'Color',[1 1 0]);

    drawnow limitrate; 
    set(figureHandle,'Visible','on');

end

Params.isPlottingInRealTime = true;
% -----
display('RAMPING SHIM CURRENTS TO DC') 
Shim.Com.setandrampallshims( Shim.Opt.Model.currents ) ;
pause(1)
Shim.Com.getallchanneloutputs

% ------- 
StopButton = stoploop({'Stop recording'}) ;

isTracking = Shim.Aux.beginrecording(); 

if isTracking

    ShimUse.customdisplay( 'Issuing real-time shim updates. Begin scanning!' )

    while ( iSample < nSamples ) && ~StopButton.Stop()
        
        iSamplesBetweenUpdates = 0;

        % acquire batch of respiration samples
        for iSamplesBetweenUpdates = 1 : nSamplesBetweenUpdates 

            iSample = iSample + 1 ;
            Shim.Aux.getupdate() ;

        end

        iUpdate                = iUpdate + 1;
        updateIndices(iUpdate) = iUpdate ;
        updateTimes(iUpdate)   = updatePeriod * iUpdate ;

        currents = Shim.Opt.computerealtimeupdate( Shim.Aux.Data.p(end) ) ;
        
        % currents = limitcurrents( currents ) ;
        
        Shim.Com.setandloadallshims( currents ) ;
        
        currentsIssued(:, iUpdate ) = currents ;
       
        currentsNorm(iUpdate) = norm(currents) ;
        
        if Params.isPlottingInRealTime
            % set(plotHandle,'YData',Shim.Aux.Data.p,'XData',Shim.Aux.Data.t);
            set(plotHandle,'XData',updateTimes, 'YData', currentsNorm);
        end

    end

    Shim.Aux.stoprecording() ;
    Shim.Com.resetallshims() ;

    sampleTimes = (Shim.Aux.Specs.dt/1000)*ones(size(Shim.Aux.Data.p)) ;

    % ------- 
    if Params.isSavingData
        baseFilename = [ Params.dataSaveDir datestr(now,30) '-physioSignal-RTShimming' ] ;
        p = Shim.Aux.Data.p ;
        t = Shim.Aux.Data.t ;
        save( baseFilename, 'currentsIssued', 'updateTimes', 'p', 't' )
    end

end

StopButton.Clear() ;

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
% DEFAULT_TRAININGMODE = ['breath-hold'] ; vs. 'free breathing'
%
% DEFAULT_ISCONNECTINGTOSHIM = true ; 

DEFAULT_ISCONNECTINGTOSHIM      = true ; 

DEFAULT_ISCONFIRMINGDATALOADDIR = true ; 

DEFAULT_UIMODE      = 'isCmdLine' ;% vs. 'isGui'

DEFAULT_SHIMSYSTEM  = 'Rriyan' ; 

DEFAULT_ISLOGGINGCOMMANDS  = false ; 
DEFAULT_COMMANDLOGFILENAME = ['commandLog_' datestr(now,30)] ;

% Re: 'Training' protocol
DEFAULT_TRAININGMODE         = 'breath-hold' ;

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

if ~myisfield( Params, 'trainingMode' ) || isempty( Params.trainingMode ) 
   Params.trainingMode = DEFAULT_TRAININGMODE ;
end

switch Params.trainingMode
    case 'breath-hold' 
    % training time-points (2 series: inspired + expired field maps)
        Params.nTrainingSeries = 2 ;
        % not sure if these are needed anymore:
        Params.Inspired        = [] ; 
        Params.Expired         = [] ;

    case 'free breathing'
        Params.nTrainingSeries = 1 ;
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

meanings    = cell(length(systemResponse), 1) ;
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
%   Move method to ShimCom_IUGM_Prisma_fit() 

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
