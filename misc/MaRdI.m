classdef MaRdI < MrdiUtil & MrdiProt & MrdiProc & MrdiMgmt & MrdiInfo & matlab.mixin.SetGet 
%MaRdI Ma(t)-R-dI(com)
%
%   Dicom into Matlab for Siemens MRI data
%
% .......
%   
% Usage
%
%   Img = MaRdI( imgPath )
% 
%   where imgPath is the path to a single .dcm image OR a directory containing
%   the .dcm or .IMA images.
%
% Img contains public properties:
%
%   .img
%       Array of images: 
%       If multi-echo, the echo index is along the 4th dimension;
%       If multiple volume repetitions, the measurement index is along the 5th dimension.
%
%   .Aux
%       Aux-objects: auxiliary measurements (e.g. respiratory ProbeTracking)
%       By default .Aux is empty. To fill it, call Img.associateaux( Aux ) with a valid
%       Aux object. For more info, type: help MaRdI.associateaux( )
% 
% In addition to the read-only properties
%
%    .Hdr 
%       The full Siemens DICOM header corresponding to Img.img(:,:,1,1,1) 
%
%    .Hdrs 
%       Cell array of (truncated) DICOM headers courtesy of dicominfo().
%       (One entry for every image)
%
% =========================================================================
% Author::ryan.topfer@polymtl.ca
% =========================================================================

%% ========================================================================
%
% ..... 
% WRITE()
%   Saving as DICOM (and, by extension, NifTI) not rigorously/fully tested
% 
%%

% =========================================================================
% =========================================================================
properties
    img ;
    Aux ;
end

properties(SetAccess={?MaRdI, ?MrdiUtil, ?MrdiProc})
    Hdr ; % full Siemens DICOM header of 1st img (i.e. Img.img(:,:,1) )
    Hdrs ; % cell array of (truncated) DICOM headers courtesy of dicominfo()
    Ref ; % Reference properties - prior to manipulation
end

% properties(SetAccess=protected, Hidden = true)
%     Ref ; % Reference properties - prior to manipulation
% end

% =========================================================================
% =========================================================================    
methods
% =========================================================================    
function Img = MaRdI( imgPath )

Img.img  = [] ;
Img.Hdr  = [] ;
Img.Hdrs = [] ;
Img.Aux  = [] ;

if nargin == 1 
    
    if ~exist( imgPath )
        error( 'DICOM images not found. Check input path is valid.' ) ;
        return;
    
    elseif ( exist( imgPath ) == 2 ) % input is a single file (e.g. DICOM) 
        
        [~, ~, ext] = fileparts( imgPath ) ;
        
        assert( strcmp( ext, '.dcm' ) || strcmp( ext, '.IMA' ), ...
            'Input must be a path string leading to a single dicom image OR to a dicom-containing folder' )

        Img.img = double( dicomread( imgPath ) ) ;
        Img.Hdr = MaRdI.dicominfosiemens( imgPath ) ;

        %  Add/replace a few Hdr entries:
        
        if ~myisfield( Img.Hdr, 'SpacingBetweenSlices' ) 
            Img.Hdr.SpacingBetweenSlices = Img.Hdr.SliceThickness ;
        end
       
        Img.setslicenormalvector() ;
    
    elseif ( exist( imgPath  ) == 7 ) % input is a directory (e.g. containing DICOMs)
        
        imgDir = [ imgPath '/' ] ; 

        imgList = MaRdI.findimages( imgDir ) ;
        nImages = length( imgList ) ;

        if ~nImages
            return;
        end
        
        % Read protocol info from 1st loaded image (arbitrary):
        SpecialHdr = MaRdI.dicominfosiemens( imgList{1} ) ;

        nRows      = SpecialHdr.Height ;
        nColumns   = SpecialHdr.Width ;
        
        if ~myisfield( SpecialHdr.MrProt, 'lRepetitions' ) 
            nVolumes = 1 ;
        else
            nVolumes = SpecialHdr.MrProt.lRepetitions + 1 ;
        end
        
        % containers for a few terms varying between dicoms:
        sliceLocations = zeros( nImages, 1 ) ; % [units: mm]
        echoTimes      = zeros( nImages, 1 ) ; % [units: ms]

        % Image headers, not yet organized:
        RawHdrs        = cell( nImages, 1 ) ; 
        
        % Read all the image headers and structure the MaRdI object accordingly:
        %
        % NOTE: It may well be possible to determine all the necessary info from
        % the complete Siemens header (e.g. SpecialHdr) of a single image;
        % practically, however, the following is easier since the header
        % information is abstrusely defined for some sequences.
        for iImg = 1 : nImages
            % using dicominfo() here as parse-siemens-shadow() takes much longer
            RawHdrs{ iImg }      = dicominfo( imgList{iImg} ) ;

            sliceLocations(iImg) = RawHdrs{ iImg }.SliceLocation ;
            echoTimes(iImg)      = RawHdrs{ iImg }.EchoTime ;
        end

        sliceLocations = sort( unique( sliceLocations ), 1, 'ascend' ) ; 
        nSlices        = length( sliceLocations ) ;
        
        echoTimes      = sort( unique( echoTimes ), 1, 'ascend' ) ; 
        nEchoes        = length( echoTimes ) ; 
        
        Img.img = zeros( nRows, nColumns, nSlices, nEchoes, nVolumes ) ;
        
        % copy of RawHdrs to be reorganized:
        Img.Hdrs = cell( nSlices, nEchoes, nVolumes ) ;
        
        for iImg = 1 : nImages

            iHdr    = RawHdrs{ iImg } ;

            iVolume = iHdr.AcquisitionNumber ;
            iSlice  = find( iHdr.SliceLocation == sliceLocations ) ;
            iEcho   = find( iHdr.EchoTime == echoTimes ) ;

            Img.Hdrs{ iSlice, iEcho, iVolume } = iHdr ;

            Img.img( :, :, iSlice, iEcho, iVolume ) = dicomread( imgList{iImg} ) ;

        end
        
        % Save the complete 1st Hdr
        Img.Hdr = MaRdI.dicominfosiemens( Img.Hdrs{1}.Filename ) ;
            
        Img.rescaleimg() ;

        if ~myisfield( Img.Hdr, 'SpacingBetweenSlices' ) 
            Img.Hdr.SpacingBetweenSlices = Img.Hdr.SliceThickness ;
        end
        
        if ~isempty( strfind( Img.Hdr.ImageType, 'MOSAIC' ) ) 
            Img = reshapemosaic( Img ) ;
        end

        Img.setslicenormalvector() ;

    end
end

end
% =========================================================================
function [] = associateaux( Img, Aux, Params )
%ASSOCIATEAUX - link image to corresponding auxiliary recording object
%
%  Usage 
%
%   [] = ASSOCIATEAUX( Img, Aux )
%   [] = ASSOCIATEAUX( Img, Aux, Params )
% 
% .......
%
% Description
%   
% Function compares an acquired image (typically, a time-series of multiple
% images) and an auxiliary recording (e.g. respiratory trace, ideally 
% featuring synchronization triggers) and returns (copies to Img.Aux) an
% estimate of the Aux recording corresponding to each image measurement.
%
% Cases:
%   
%   1. Single measurement Aux:
%       The value is copied across all image time points
%
%   2. Single measurement Img:
%       If Aux possesses multiple measurements, the image is presumed to have
%       been taken during a breath-hold. In this case, the Aux signal variance
%       is calculated over a shifting time-window (width = total image
%       acquisition time) and the median Aux value over the window pertaining
%       to the least variance is returned. (Length of Aux recording must be >=
%       total scan time).
%
%   3. Img and Aux are both time-series:
%       Length of Aux recording must be >= length of the image time-series.
%
%       3.0: Img and Aux time-points already correspond: 
%            Aux is simply copied. 
%
%       3.1: Aux recording possesses a single synchronization trigger:
%           The trigger is assumed to correspond to the first image in the time series.
%           Aux is cropped and linearly interpolated.
%
%       3.2: Aux recording possesses multiple synchronization triggers:
%           NOT IMPLEMENTED
%
%       3.3: Aux recording does not possess triggers:
%           DEPRECATED 
% .......
%
% Params - optional struct for which the following Params.fields are supported
% 
% Only used in Case 3 (image time series):
%
%   .interpolationMethod    [default = 'linear']
%       argument to INTERP1(), used to interpolate between Aux samples
%       (see INTERP1 for other options)
%
%   .auxDelay   [default = 0]
%       estimation of transmission delay inherent in the Aux recording process [units: ms]

if isempty( Aux ) || ~myisfield(Aux, 'Data') || ~myisfieldfilled(Aux.Data, 'p')  
    error('Aux recording is empty.')
end
 
DEFAULTS.interpolationMethod = 'linear' ;
DEFAULTS.auxDelay            = 0 ;

if nargin < 3 || isempty(Params)
    Params.dummy = [] ;
end

Params = assignifempty( Params, DEFAULTS ) ;

nSamples = length( Aux.Data.p ) ; % number of aux samples
tk0      = Img.estimatekorigintime() ;
tAcq     = Img.getacquisitiontime() ;
kDelay   = tk0(1) - tAcq(1) ; % delay between acquisition start and image ~content time~ around k-space origin

if size( tk0, 1 ) > 1
    warning( 'Img contains multiple slices. Returned Aux recording will be associated with the 1st slice.') 
elseif size( tk0, 2 ) > 1
    warning( 'Img contains multiple echoes. Returned Aux recording will be associated with the 1st echo.') 
end

tk0      = squeeze( tk0(1,1,:) )' ;
nAcq     = size( Img.img, 5 ) ; % number of acquisitions
Img.Aux  = Aux.copy() ;

%% -----
% check if image + aux recording times already coincide:
if ( nAcq == nSamples ) && ( all( Img.Aux.Data.t == tk0 ) )
    warning('Image and Aux recording times already coincide. Not performing interpolation.') ;
    return ;
end

%% -----
% trivial case
if nSamples == 1
    Img.Aux.Data.p       = Aux.Data.p*ones(1 , nAcq) ;
    Img.Aux.Data.pRaw    = Aux.Data.p*ones(1, nAcq) ;
    Img.Aux.Data.t       = tk0 ;
    Img.Aux.Data.trigger = zeros(1, nAcq) ;
    return ;
else
    assert( length(Aux.Data.t) == length(Aux.Data.p), ...
        'Function expects one sample time-point for each input sample.' )
end

% shifted image origin times (t = 0 will correspond to the first image):
tImg = tk0 - tk0(1) ;

% time between images 
dtImg = median( diff( tImg ) ) ;
    
% time between respiratory samples 
dtAux = round( ( Img.Aux.Data.t(end) - Img.Aux.Data.t(1) )/nSamples ) ; % [units: ms]
    
assert( max(Img.Aux.Data.t)>=max(tImg), 'Length of Aux recording should be greater or equal to the total image acquisition time.' )

%% ----- 
% assumed breath-hold case
if nAcq == 1
    % single image time-point assumed to correspond to a breath-hold and, thus, a single tracker value:

    % number of samples corresponding to breath-hold
    % (enables auto-estimation of the time window corresponding to the breath-hold):
    nSamplesApnea = (1000*Img.Hdr.MrProt.lTotalScanTimeSec)/dtAux ;

    % extract a single scalar
    p = ProbeTracking.selectmedianmeasurement( Img.Aux.Data.p, nSamplesApnea ) ;
    
    Img.Aux.Data.p       = p*ones(1 , nAcq) ;
    Img.Aux.Data.pRaw    = p*ones(1, nAcq) ;
    Img.Aux.Data.t       = tk0 ;
    Img.Aux.Data.trigger = zeros(1, nAcq) ;
    return ;
end

%% -----
% Determine if aux recording features synchronization (e.g. trigger pulses)
if myisfield( Img.Aux.Data, 'logStartMdhTime' ) && myisfield( Img.Aux.Data, 'logStopMdhTime' )
% Recording is from the Siemens PMU : find sample corresponding to scan start + insert trigger there
    tScan = Img.getacquisitiontime() ;
    
    % MDH time format i.e. "msecs since midnight", see https://cfn.upenn.edu/aguirre/wiki/public:pulse-oximetry_during_fmri_scanning
    if ( min(tScan(:)) < Img.Aux.Data.logStartMdhTime ) || ( max(tScan(:)) > Img.Aux.Data.logStopMdhTime  )
        error( ['Image acquisition times out of the range of the PMU sample times.'] ) ;
    end
    
    [~, iMin] = min( abs( Img.Aux.Data.t - tScan(1) ) ) ;
    Img.Aux.Data.trigger( iMin ) = 1 ; 

end

if ~isempty( Img.Aux.Data.trigger )
    nTriggers = nnz( Img.Aux.Data.trigger ) ;
else
    nTriggers = 0;
end

%% -----
if nTriggers > 0 
    % NOTE: 1st trigger assumed to correspond to 1st recorded value of the 1st
    % acquisition. However, the image 'content' (around k=0) will generally
    % occur sometime later (exceptions include spiral imaging)
     
    iTriggers = find( Img.Aux.Data.trigger==1 ) ;
    
    % % crop superfluous recording preceding 1st trigger
    % Img.Aux.Data.p       = Img.Aux.Data.p( iTriggers(1) : end ) ;
    % Img.Aux.Data.pRaw    = Img.Aux.Data.pRaw( iTriggers(1) : end ) ;
    % Img.Aux.Data.trigger = Img.Aux.Data.t( iTriggers(1) : end ) ;
    % Img.Aux.Data.t       = Img.Aux.Data.t( iTriggers(1) : end ) ;

    % shift in time s.t. t=0 corresponds to 1st image content time:
    Img.Aux.Data.t = Img.Aux.Data.t - Img.Aux.Data.t( iTriggers(1) ) - kDelay - Params.auxDelay ; 
    
    if nTriggers == 1
        %% -----
        % Interpolate aux recordings across time  
        tInterp              = [ 0: tImg(end)/(nAcq-1) : tImg(end) ]' ; % interp1 requires regular grid spacing
        Img.Aux.Data.p       = interp1( Img.Aux.Data.t, Img.Aux.Data.p, tInterp, Params.interpolationMethod )' ; 
        Img.Aux.Data.pRaw    = interp1( Img.Aux.Data.t, Img.Aux.Data.pRaw, tInterp, Params.interpolationMethod )' ; 
        Img.Aux.Data.t       = tk0 ;
        Img.Aux.Data.trigger = zeros(1, nAcq) ;
        return ;
    elseif nTriggers > 1
        error('Unimplemented feature: Interpolation for multi-trigger Aux recording. TODO')
    end

else
    error('Deprecated feature: Associating image and Aux recordings without synchronization trigger pulses via x-correlation. TODO')
    % % estimate via correlation
    % % determine intersystem (image acq. to respiratory tracker) delays
    %
    % gridSize = Img.getgridsize() ;
    %
    % % % ------
    % % % bandpass filter images?
    % % 
    % % frequencyCutoff = [ 0 0.6 ] ;
    % %
    % % for iRow = 1 : gridSize(1) 
    % %     for iColumn = 1 : gridSize(2)  
    % %         for iSlice = 1 : gridSize(3)
    % %         Field.img( iRow, iColumn, 1, : ) = ...
    % %             bpfilt( squeeze( Img.img(iRow, iColumn, iSlice,:) ), frequencyCutoff(1), frequencyCutoff(2), 1/dtAcq, false ) ;
    % %         end
    % %     end
    % % end
    %     
    % % ------
    % % Interpolate img across time to match the interval of the physio recording 
    % tInterp = [ 0: dtAux : tAcq(end) ]' ;
    %
    % imgInterp = zeros( [ gridSize length(tInterp) ] ) ;
    %
    % for iRow = 1 : gridSize(1) 
    %     for iColumn = 1 : gridSize(2)  
    %         for iSlice = 1 : gridSize(3)
    %             imgInterp( iRow, iColumn, iSlice, : ) = interp1( tAcq, squeeze( Img.img( iRow, iColumn, iSlice, :) ), tInterp, 'linear' ) ; 
    %         end
    %     end
    % end
    %
    % % Assume the entire img series occurs somewhere within the Aux.Data.t window.
    % % given that we assign tInterp(1) = 0, 
    % % the delay time (tInterp' = Aux.Data.t' - t_delay) is necessarily positive, 
    % % and less than max(Aux.Data.t) - max(tInterp)
    %
    % maxLag = round(( max( Aux.Data.t ) - max( tInterp ) )/dt) ;
    %
    % T_delay = nan( gridSize ) ;
    % C       = zeros( [gridSize 1] ) ;
    %
    % % reliable region:
    % mask  = ( sum( Field.Hdr.MaskingImage, 4 ) == nAcq ) ;
    %
    % ShimUse.customdisplay( 'Estimating probe-to-image delay...')
    % % ------
    % % X-corr
    % for iRow = 1 : Field.Hdr.Rows 
    %     for iColumn = 1 : Field.Hdr.Columns  
    %          if mask( iRow, iColumn )
    %
    %             [c, lags] = xcorr( Aux.Data.p, imgInterp(iRow,iColumn,:), maxLag ) ;
    %
    %             % truncate to positive delay times
    %             c    = c( lags >=0 ) ;
    %             lags = lags( lags >=0 ) ;
    %             
    %             [maxC, iMax] = max(abs(c)) ;
    %             t_xcorr = dt*lags ;
    %             t_delay = t_xcorr( iMax ) ;
    %
    %             C(iRow, iColumn, 1, 1:length(c)) = c ;
    %             T_delay(iRow, iColumn) = t_delay ;
    %         end
    %     end
    % end
    %
    % % the median should be more robust than the average for determining the delay
    % % once noise + uncorrelated voxels are taken into consideration
    % t_medianDelay = median( T_delay( mask ) ) ;
    %
    % ShimUse.customdisplay( ['Estimated delay: ' num2str( t_medianDelay )] ) ;
    %     
    % response = input(['Is the current estimate satisfactory? ' ...
    %         '0 to manually specify delay; 1 (or enter) to accept & continue: ']) ;
    %
    % if ~response
    %     
    %     figure
    %     imagesc( mask .* median( imgInterp, 4 ) ) ;
    %     title('Median field') ;
    %     ShimUse.customdisplay( 'Choose a voxel to plot its field time-course: ' ) ;
    %
    %     isVoxelInReliableRegion = false ;
    %
    %     while ~isVoxelInReliableRegion
    %         iRow = input( ['Enter the voxel image row: ' ] ) 
    %         iCol = input( ['Enter the voxel image column: ' ] )
    %         
    %         isVoxelInReliableRegion = mask( iRow, iCol ) ;
    %         
    %         if ~isVoxelInReliableRegion
    %             ShimUse.customdisplay( 'Selected voxel lies outside the reliable region. Choose again' ) ;
    %         end
    %     end
    %     
    %     figure; 
    %     subplot(121); plot( Aux.Data.t, Aux.Data.p ); 
    %     xlabel('Time (s)') ; ylabel('Amplitude (a.u.)') ; 
    %     title('Physio signal') ;
    %     
    %     subplot(122); plot( tInterp, squeeze(imgInterp(iRow,iCol,:)) ) ;
    %     xlabel('Time (s)') ; ylabel('Field (Hz)') ; 
    %     title(['Field at voxel [' num2str(iRow) ',' num2str(iCol) ']' ]) ;
    %
    %     isDelayValid = false ;
    %     while ~isDelayValid
    %         t_delay = input(['Enter the delay (between 0-' num2str(maxLag) ' s) in units of seconds: '] ) 
    %         
    %         isDelayValid = ( t_delay >= 0 ) & ( t_delay <= maxLag ) ;
    %         
    %         if ~isDelayValid
    %             ShimUse.customdisplay( 'Given delay exceeds the expected range. Choose again' ) ;
    %         end
    %     end 
    % else
    %     t_delay = t_medianDelay ;
    % end
    %
    % [~,iT] = min( abs( Aux.Data.t - t_delay ) ) ;
    %
    % % Finally, interpolate Aux.Data.p to the acquired image times:
    % tp = Aux.Data.t( iT:iT+length(tInterp)-1 ) -Aux.Data.t(iT) ;
    %
    % % NOTE: tp(end) can be slightly < tAcq(end)  
    % % Using 'extrap' to avoid an interpolated p terminating with NaN's
    % Shim.Data.Img{ 1, 3, iSeries }.Aux.Tracker.Data.p = interp1( tp, Aux.Data.p( iT:iT+length(tInterp)-1 ), tAcq, 'linear', 'extrap' ) ;
    % Shim.Data.Img{ 1, 3, iSeries }.Aux.Tracker.Data.t = tAcq ;
    %
    % if any( isnan(Shim.Data.Img{ 1, 3, iSeries }.Aux.Tracker.Data.p) )
    %     warning('Interpolated probe recording contains nonnumeric values. Check inputs. (E.g. Shim.Data.Img{ }: Entries should be ordered by acquisition time)' ) ;
    % end
    %
    % % used in nonlinear shim optimization:
    % Shim.Params.pMax     = max( Shim.Data.Img{ 1, 3, iSeries }.Aux.Tracker.Data.p ) ;
    % Shim.Params.pMin     = min( Shim.Data.Img{ 1, 3, iSeries }.Aux.Tracker.Data.p ) ;
    %
    % % for posterity (might be useful?):
    % Shim.Params.imgDelay = t_delay ;    

end

end
% =========================================================================
function [imgType] = getimagetype( Img ) 
%GETIMAGETYPE   Returns image type as string 
% 
% imgType = GETIMAGETYPE( Img )
%
% Returns either 'PHASE', 'MAGNITUDE', or 'UNKNOWN'

if strfind( Img.Hdr.ImageType, '\P\' )
    imgType = 'PHASE' ;
elseif strfind( Img.Hdr.ImageType, '\M\' )
    imgType = 'MAGNITUDE' ;
else
    imgType = 'UNKNOWN' ;
end

end
% =========================================================================
function mask = getreliabilitymask( Mag, threshold )
%GETRELIABILITYMASK
% 
%  mask = getreliabilitymask( Mag )
%  mask = getreliabilitymask( Mag, threshold )
% 
%  For each echo and each measurement, GETRELIABILITYMASK normalizes
%  Mag.img(:,:,:,iEcho,iMeasurement) and returns a logical mask wherein
%  elements are assigned TRUE whenever the corresponding normalized magnitude
%  voxel is > threshold
%
%  By default, threshold = 0.01

assert( Mag.ismagnitude(), 'Function expected magnitude image input' ) ;

DEFAULTS.threshold = 0.01 ;

if ( nargin == 1 ) || isempty( threshold )
    threshold = DEFAULTS.threshold ;
end

mask = false( size( Mag.img ) ) ;

for iVolume = 1 : size( mask, 5 ) 
    for iEcho = 1 : size( mask, 4 ) 
        mag = Mag.img(:,:,:,iEcho, iVolume) ;
        mask(:,:,:,iEcho,iVolume) = ( mag./max(mag(:)) ) > threshold ;
    end
end

end
% =========================================================================
function [isMag] = ismagnitude( Img )
%ISMAGNITUDE    Returns TRUE if Img is a magnitude image, FALSE otherwise.
% 
% isMag = ISMAGNITUDE( Img )

switch Img.getimagetype()
    case 'MAGNITUDE'
        isMag = true ;
    otherwise
        isMag = false ;
end

end
% =========================================================================
function [isPhase] = isphase( Img )
%ISPHASE    Returns TRUE if Img is a phase image, FALSE otherwise.
%
% isPhase = ISPHASE( Img )

switch Img.getimagetype()
    case 'PHASE'
        isPhase = true ;
    otherwise
        isPhase = false ;
end

end
% =========================================================================
function xyzIso = isocenter( Img )
%ISOCENTER
% 
% xyzIso = ISOCENTER( Img ) 
%
% Returns the 3-element vector of the x, y and z coordinates of the magnet
% isocenter in the patient coordinate system

xyzIso = Img.Hdr.Img.ImaRelTablePosition()' ;

assert( xyzIso(1) == 0, 'Table shifted in L/R direction?' ) ;
assert( xyzIso(2) == 0, 'Table shifted in A/P direction?' ) ;

end
% =========================================================================
function [mask, weights] = segmentspinalcanal( Img, Params )
%SEGMENTSPINALCANAL
% 
% segment T2* multiecho data using the Spinal Cord Toolbox (must be installed + in path)
%
% [ mask, weights ] = SEGMENTSPINALCANAL( Img, Params )
%
% Params
%
%   .dataLoadDir 
%       DICOM folder
%   
%   .dataSaveDir 
%
%   .isUsingPropsegCsf
%       [default = false]
%
% NOTE
%   The protocol is basically that of Topfer R, et al. Magn Reson Med, 2018. 
%   It hasn't been tested extensively for different acquisition prtocols/systems

mask = false ;

if nargin < 2 || isempty(Params)
    disp('Default parameters will be used')
    Params.dummy = [] ;
end

if  ~myisfield( Params, 'dataLoadDir' ) || isempty(Params.dataLoadDir)
    [Params.dataLoadDir,~,~] = fileparts( Img.Hdr.Filename ) ;
end

[mask, weights] = MaRdI.segmentspinalcanal_s( Params ) ;

end
% =========================================================================
function timeAverage = timeaverage( Img )
%TIMEAVERAGE
% 
% Img = TIMEAVERAGE( Img) 
% 
% Assumes 5th dimension of Img.img corresponds to time:
%   
%   timeAverage = mean( Img.img, 5 ) ;

timeAverage = mean( Img.img, 5 ) ;

end
% =========================================================================
function timeStd = timestd( Img )
%TIMESTD
% 
% standardDeviation = TIMESTD( Img ) 
% 
% Assumes 5th dimension of Img.img corresponds to time:
%   
%   standardDeviation = std( Img.img, 0, 5 ) ;

timeStd = std( Img.img, 0, 5 ) ;

end
% =========================================================================
function [] = write( Img, saveDirectory, imgFormat, isSavingSingleNiis )
%WRITE Ma(t)-R-dI(com)
% 
% Write MaRdI image object to DICOM and/or NifTI
%
%.....
%
%   Usage 
%
%   WRITE( Img )
%   WRITE( Img, saveDirectory )
%   WRITE( Img, saveDirectory, imgFormat ) 
%   WRITE( Img, saveDirectory, imgFormat, isSavingSingleNiis ) 
%
%   Inputs
%
%   default saveDirectory = './tmp'
%   
%   imgFormat can be:
%       'dcm' [default]
%       'nii' (creating temporary DICOMs which are deleted after the system call to dcm2niix)
%       'both' (does not delete the DICOMs)
%
%   isSavingSingleNiis (boolean):
%       false [default] : DICOMs are combined into single NifTI file
%       true : Separate .nii output for each image (passes '-s y' argument to dcm2niix)
%
%.....
%
% Adapted from dicom_write_volume.m (D.Kroon, University of Twente, 2009)
% https://www.mathworks.com/matlabcentral/fileexchange/27941-dicom-toolbox?focused=5189263&tab=function

DEFAULT_SAVEDIRECTORY      = './tmp' ;
DEFAULT_IMGFORMAT          = 'dcm' ;
DEFAULT_ISSAVINGSINGLENIIS = false ;

if nargin < 2 || isempty(saveDirectory)
    saveDirectory = DEFAULT_SAVEDIRECTORY ;
end

if nargin < 3 || isempty(imgFormat)
    imgFormat = DEFAULT_IMGFORMAT ;
end

if nargin < 4 || isempty(isSavingSingleNiis)
    isSavingSingleNiis = DEFAULT_ISSAVINGSINGLENIIS ;
end

fprintf(['\n Writing images to: ' saveDirectory ' ... \n'])

[~,~,~] = mkdir( saveDirectory ) ;

rescaleimg( Img, true ) ;

[X,Y,Z] = Img.getvoxelpositions() ;

%-------
% write Hdr
Hdr.NumberOfSlices          = size( Img.img, 3 ) ;

% Make random series number
SN                          = round(rand(1)*1000);
% Get date of today
today                       = [datestr(now,'yyyy') datestr(now,'mm') datestr(now,'dd')];
Hdr.SeriesNumber            = SN;
Hdr.AcquisitionNumber       = SN;
Hdr.StudyDate               = today;
Hdr.StudyID                 = num2str(SN);
Hdr.PatientID               = num2str(SN);
Hdr.AccessionNumber         = num2str(SN);

% copy from original
Hdr.ImageType               = Img.Hdr.ImageType ; 

Hdr.StudyDescription        = Img.Hdr.StudyDescription ;
Hdr.SeriesDescription       = Img.Hdr.SeriesDescription ;
Hdr.Manufacturer            = Img.Hdr.Manufacturer ;
Hdr.ScanningSequence        = Img.Hdr.ScanningSequence ;
Hdr.SequenceVariant         = Img.Hdr.SequenceVariant ;
Hdr.ScanOptions             = Img.Hdr.ScanOptions ;
Hdr.MRAcquisitionType       = Img.Hdr.MRAcquisitionType ;
Hdr.SliceThickness          = Img.Hdr.SliceThickness ;
Hdr.SpacingBetweenSlices    = Img.Hdr.SpacingBetweenSlices ;
Hdr.PatientPosition         = Img.Hdr.PatientPosition ;
Hdr.PixelSpacing            = Img.Hdr.PixelSpacing ;

Hdr.ImageOrientationPatient = Img.Hdr.ImageOrientationPatient ;
Hdr.SliceLocation           = Img.Hdr.SliceLocation ; 

Hdr.AcquisitionMatrix       = Img.Hdr.AcquisitionMatrix ; 
Hdr.RepetitionTime          = Img.Hdr.RepetitionTime ; 
Hdr.NumberOfAverages        = Img.Hdr.NumberOfAverages ; 
Hdr.PercentSampling         = Img.Hdr.PercentSampling ; 
Hdr.PercentPhaseFieldOfView = Img.Hdr.PercentPhaseFieldOfView ; 
Hdr.InPlanePhaseEncodingDirection  = Img.Hdr.InPlanePhaseEncodingDirection ; 

[rHat, cHat, sHat] = Img.getdirectioncosines( ) ;  

nSlices       = size( Img.img, 3 ) ;
nEchoes       = numel( Img.getechotime() ) ;
nAcquisitions = numel( Img.getacquisitiontime ) ;

nImg = nSlices*nEchoes*nAcquisitions ;

for iSlice = 1 : nSlices 
    for iEcho  = 1 : nEchoes 
        for iAcq = 1 : nAcquisitions  

            iImg = iSlice*iEcho*iAcq ;

            disp( [num2str(100*iImg/nImg) '% ...'] ) ;
    
            %-------
            % filename 
            sliceSuffix = '000000' ;
            sliceSuffix = [ sliceSuffix( length(iSlice) : end ) num2str(iSlice) '-' num2str(iEcho) '-' num2str(iAcq) ] ;
            sliceSuffix = ['-' sliceSuffix '.dcm'] ;
            filename    = [saveDirectory '/' Img.Hdr.PatientName.FamilyName sliceSuffix] ;

            %-------
            % image specific hdr info 
            Hdr.ImageNumber          = iSlice*iEcho*iAcq ;
            Hdr.InstanceNumber       = iSlice*iEcho*iAcq ;
            Hdr.AcquisitionTime      = Img.Hdrs{iSlice,iEcho,iAcq}.AcquisitionTime ;
            
            Hdr.ImagePositionPatient = [(X(1,1,iSlice)) (Y(1,1,iSlice)) (Z(1,1,iSlice))] ;

            Hdr.SliceLocation        = dot( Hdr.ImagePositionPatient, sHat ) ;
           
            dicomwrite( Img.img(:,:,iSlice,iEcho, iAcq) , filename, 'ObjectType', 'MR Image Storage', Hdr ) ;
            
            if( iSlice==1 )
                info                  = dicominfo( filename ) ;
                Hdr.StudyInstanceUID  = info.StudyInstanceUID ;
                Hdr.SeriesInstanceUID = info.SeriesInstanceUID ;
            end

        end
    end
end

%-------
if ~strcmp( imgFormat, 'dcm' )

    if isSavingSingleNiis
        list = MaRdI.findimages( saveDirectory ) ;
        for iImg = 1 : size( list, 1 )
            system(['dcm2niix -s y -f %f_%r ' list{iImg}]) ;
        end
    else
        system(['dcm2niix ' saveDirectory]) ;
    end

    if strcmp( imgFormat, 'nii' )
        delete( [saveDirectory '/*.dcm'] ) ;
    end
end

end
% =========================================================================

end
% =========================================================================
% =========================================================================
methods(Access=private)
% =========================================================================
function [] = rescaleimg( Img, isUndoing )
%RESCALEIMG
%
% []=RESCALEIMG( Img ) 
% []=RESCALEIMG( Img, isUndoing ) 
%
% Only alters phase images (+ Hdr) by rescaling to radians. Mag. is unaffected.

if ~Img.isphase() 
    return;
end

DEFAULT_ISUNDOING        = false ;
DEFAULT_RESCALEINTERCEPT = -(2^12) ;
DEFAULT_RESCALESLOPE     = 2 ;

if (nargin < 2) || isempty( isUndoing ) 
   isUndoing = DEFAULT_ISUNDOING ; 
end 

if ~isUndoing
    
    if  ~myisfield( Img.Hdr, 'RescaleIntercept' ) || isempty( Img.Hdr.RescaleIntercept )
        Img.Hdr.RescaleIntercept = DEFAULT_RESCALEINTERCEPT ;
    end

    if ~myisfield( Img.Hdr, 'RescaleSlope' ) || isempty( Img.Hdr.RescaleSlope )
        Img.Hdr.RescaleSlope = DEFAULT_RESCALESLOPE ;
    end

    % rescale to rad:
    Img.img = pi*(Img.Hdr.RescaleSlope .* double( Img.img ) ...
                + Img.Hdr.RescaleIntercept)/double(2^Img.Hdr.BitsStored) ;

    % update Hdr:
    Img.Hdr.ImageType  = 'ORIGINAL\SECONDARY\P\' ; 
    Img.Hdr.PixelRepresentation = uint8(1) ; % i.e. signed 
    Img.Hdr.PixelComponentPhysicalUnits = '0000H' ; % i.e. none

else
    
    Img.Hdr.RescaleIntercept = min( Img.img(:) ) ;
    Img.img                  = Img.img - Img.Hdr.RescaleIntercept ;
    Img.Hdr.RescaleSlope     = max( Img.img(:) )/double( 2^Img.Hdr.BitsStored ) ;
    Img.img                  = uint16( Img.img / Img.Hdr.RescaleSlope ) ;

end

end
% =========================================================================
function Img = reshapemosaic( Img )
%RESHAPEMOSAIC
% 
% Reshape Siemens mosaic into volume and remove padded zeros
%
% Adapted from dicm2nii by
% xiangrui.li@gmail.com 
% http://www.mathworks.com/matlabcentral/fileexchange/42997

assert( ~isempty(strfind( Img.Hdr.ImageType, 'MOSAIC' ) ), 'Corrupt image header?' ) ;       

nImgPerLine = ceil( sqrt( Img.getnumberofslices() ) ); % always nImgPerLine x nImgPerLine tiles

nRows    = size(Img.img, 1) / nImgPerLine; 
nColumns = size(Img.img, 2) / nImgPerLine; 
nEchoes  = size(Img.img, 4) ;  
nVolumes = size(Img.img, 5) ;

img = zeros([nRows nColumns Img.getnumberofslices() nEchoes nVolumes], class(Img.img));

for iImg = 1 : Img.getnumberofslices()

    % 2nd slice is tile(1,2)
    r = floor((iImg-1)/nImgPerLine) * nRows + (1:nRows);     
    c = mod(iImg-1, nImgPerLine) * nColumns + (1:nColumns);

    img(:, :, iImg, :) = Img.img(r, c, :);
end

Img.img = img ;


% -----
% Update header

% -----
% Correct Hdr.ImagePositionPatient 
%   see: http://nipy.org/nibabel/dicom/dicom_mosaic.html

%-------
% Rotation matrix: R
[r, c, s] = Img.getdirectioncosines ;
R = [r c s] ; 
                  
%-------
% Scaling matrix: S  
voxelSize = Img.getvoxelspacing() ;
S = diag(voxelSize); 

RS = R*S ;

Img.Hdr.ImagePositionPatient = Img.Hdr.ImagePositionPatient + ...
    RS(:,1)*(double(Img.Hdr.Rows) - nRows)/2 + ...
    RS(:,2)*(double(Img.Hdr.Columns) - nColumns)/2 ;


Img.Hdr.Rows    = uint16(nRows) ;
Img.Hdr.Columns = uint16(nColumns) ;

tmp = strfind( Img.Hdr.ImageType, '\MOSAIC' ) ;

if ~isempty(tmp) && (tmp > 1)
    Img.Hdr.ImageType = Img.Hdr.ImageType(1:tmp-1) ;
else
    Img.Hdr.ImageType = '' ;    
end

end
% =========================================================================
function [] = setslicenormalvector( Img )
%SETSLICENORMALVECTOR
%
% For determining voxel positions in 3d slice-stack
    
r = Img.Hdr.ImageOrientationPatient(4:6) ; 
c = Img.Hdr.ImageOrientationPatient(1:3) ; 
    
Img.Hdr.Img.SliceNormalVector = cross( c, r ) ;  

if size( Img.img, 3 ) > 1
% Determine: ascending or descending slices?

    % Estimate positions of last slice based on the 1st image:
    % 1. using
    [X1,Y1,Z1] = Img.getvoxelpositions() ;     

    % 2. using the reverse
    Img.Hdr.Img.SliceNormalVector = cross( r, c ) ;  
    [X2,Y2,Z2] = Img.getvoxelpositions() ; % estimate positions based on the 1st loaded image in the directory. 

    % Actual position corresponding to the slice direction can be increasing or
    % decreasing with slice/image number. So, which estimate is closer: 1 or 2? 
    if norm( Img.Hdrs{end,1,1}.ImagePositionPatient' - [ X1(1,1,end) Y1(1,1,end) Z1(1,1,end) ] ) < ...
            norm( Img.Hdrs{end,1,1}.ImagePositionPatient' - [ X2(1,1,end) Y2(1,1,end) Z2(1,1,end) ] ) 
        % if true, then 1. corresponds to the correct orientation
        Img.Hdr.Img.SliceNormalVector = cross( c, r ) ;
    end
end

end
% =========================================================================

end
% =========================================================================
% =========================================================================
methods(Hidden=true)
% =========================================================================
function [] = nii( Img )
%NII - Write MaRdI image to NiFtI file
%
% Wraps to NII( ) (which wraps to the NiFtI toolbox)   
%
%.....
%   Syntax
%
%   nii( Img ) 
%.....
%
% WARNING
%
%   nii() function is convenient for quickly writing a file to throw
%   into an external viewing application (e.g. ImageJ). 
%   The nifti Hdr info (i.e. orientation) is probably all wrong. 
%
%   To save NifTI's properly (takes longer) use Img.write() 

workingDir       = [ pwd '/' ] ;
Params.filename  = [ workingDir Img.Hdr.PatientName.FamilyName '_' num2str( Img.Hdr.SeriesNumber ) '_' Img.Hdr.SeriesDescription  ] ;
Params.voxelSize = Img.getvoxelspacing() ;

nii( squeeze(Img.img), Params )  ;

end
% =========================================================================

end
% =========================================================================
% =========================================================================
methods(Static)
% =========================================================================
function [mask, weights] = segmentspinalcanal_s( Params )
%SEGMENTSPINALCANAL_S
% 
% segment T2* multiecho data using the Spinal Cord Toolbox
%
% [ mask, weights ] = SEGMENTSPINALCANAL_S( Params )
%
% Params
%
%   .dataLoadDir 
%       DICOM folder
%   
%   .dataSaveDir
%       [default = './gre_seg/'] 
%
%   .isUsingPropsegCsf
%       [default = false]
%
%   .centerlineMethod
%       Method used to obtain the centerline: 
%           'midfov': mask is centered in the middle of the axial FOV
%           'spinalcord': mask follows the spinal cord centerline
%       [default = 'midfov']
%
% NOTE
%   The protocol is basically that of Topfer R, et al. Magn Reson Med, 2018. 
%   It hasn't been tested extensively for different acquisition protocols or systems 
%
% TODO
% SEGMENTSPINALCANAL_S
%   is the static form of MaRdI.segmentspinalcanal( Img, Params )
%
%   I (RT) was hoping Matlab would allow the 2 identically named methods (as in C)
%   given that the static form takes only 1 arg, and the other form requires 2...
%
%   --> either rectify this if pos. or change the method names or another alt.

%   .isForcingOverwrite
%       if .nii or .gz files exist already in dataSaveDir they will be deleted
%       [default = false]
%

mask = false ;

DEFAULT_DATALOADDIR        = [] ; % path to dicom folder
DEFAULT_DATASAVEDIR        = './gre_seg/'
DEFAULT_ISFORCINGOVERWRITE = false ;
DEFAULT_ISUSINGPROPSEGCSF  = true ; %use the propseg -CSF option
DEFAULT_CENTERLINEMETHOD   = 'midfov' ;
DEFAULT_CYLINDERSIZE       = 40 ;
DEFAULT_GAUSSIANSIZE       = 20 ;

if nargin < 1 || isempty(Params) || ~myisfield( Params, 'dataLoadDir' ) || isempty(Params.dataLoadDir)
    error('Function requires struct Params. with Params.dataLoadDir defined. See documentation.')
end

if  ~myisfield( Params, 'dataSaveDir' ) || isempty(Params.dataSaveDir)
    Params.dataSaveDir = DEFAULT_DATASAVEDIR ;
elseif ( Params.dataSaveDir(end) ~= '/' )
    Params.dataSaveDir(end+1) = '/';
end

if  ~myisfield( Params, 'isForcingOverwrite' ) || isempty(Params.isForcingOverwrite)
    Params.isForcingOverwrite = DEFAULT_ISFORCINGOVERWRITE ;
end

if  ~myisfield( Params, 'cylinderSize' ) || isempty(Params.cylinderSize)
    Params.cylinderSize = DEFAULT_CYLINDERSIZE ;
end

if  ~myisfield( Params, 'gaussianSize' ) || isempty(Params.gaussianSize)
    Params.gaussianSize = DEFAULT_GAUSSIANSIZE ;
end

if  ~myisfield( Params, 'centerlineMethod' ) || isempty(Params.centerlineMethod)
    Params.centerlineMethod = DEFAULT_CENTERLINEMETHOD ;
end

Params.tmpSaveDir = [ Params.dataSaveDir 'tmp_sct_' datestr(now, 30) '/'] ;
mkdir( Params.tmpSaveDir )

% if ~Params.isForcingOverwrite & exist( Params.dataSaveDir )
%     error('Params.dataSaveDir should not exist, or use input option Params.isForcingOverwrite == true')
% end

if  ~myisfield( Params, 'isUsingPropsegCsf' ) || isempty(Params.isUsingPropsegCsf)
    Params.isUsingPropsegCsf = DEFAULT_ISUSINGPROPSEGCSF ;
end

if ~exist( Params.dataSaveDir )
    mkdir( Params.dataSaveDir ) ;
end

dicm2nii( Params.dataLoadDir, Params.tmpSaveDir )

% rename
system( ['mv ' Params.tmpSaveDir '*.nii.gz ' Params.tmpSaveDir 't2s_allEchoes.nii.gz'] ) ;

% average across echoes
system( ['sct_maths -i ' Params.tmpSaveDir 't2s_allEchoes.nii.gz -mean t -o ' Params.tmpSaveDir 't2s.nii.gz'] ) ;

% switch between methods for obtaining a pixel location per slice
if Params.centerlineMethod == 'midfov'
    % create a vertical line centered in the axial FOV
    system( ['sct_create_mask -i ' Params.tmpSaveDir 't2s.nii.gz -p center -size 1 -f box -o ' Params.tmpSaveDir 't2s_centerline.nii.gz' ] ) ;
elseif Params.centerlineMethod == 'spinalcord'
    % get cord centerline
    % TODO: make the param -c an input Params.
    system( ['sct_get_centerline -i ' Params.tmpSaveDir 't2s.nii.gz -c t2 -o ' Params.tmpSaveDir 't2s_centerline'] ) ;
end

% create a binary cylindrical mask around the centerline
system( ['sct_create_mask -i ' Params.tmpSaveDir 't2s.nii.gz -p centerline,' ...
    Params.tmpSaveDir 't2s_centerline.nii.gz -size ' ...
    num2str(Params.cylinderSize) 'mm -f cylinder -o ' ...
    Params.tmpSaveDir 't2s_seg.nii.gz' ] ) ;

% create a soft gaussian mask around the centerline
system( ['sct_create_mask -i ' Params.tmpSaveDir 't2s.nii.gz -p centerline,' ...
    Params.tmpSaveDir 't2s_centerline.nii.gz -size ' ...
    num2str(Params.gaussianSize) 'mm -f gaussian -o ' ...
    Params.tmpSaveDir 't2s_weights.nii.gz' ] ) ;

% unzip the image volumes we wish to keep
system( ['gunzip ' Params.tmpSaveDir 't2s_seg.nii.gz -df'] ) ;
system( ['gunzip ' Params.tmpSaveDir 't2s_weights.nii.gz -df'] ) ;

% delete the other images
system( ['rm ' Params.tmpSaveDir '*.nii.gz'] ) ;

% move segmentation + weights
system( ['mv ' Params.tmpSaveDir 't2s_seg.nii ' Params.dataSaveDir 'gre_seg.nii'] ) ;  
system( ['mv ' Params.tmpSaveDir 't2s_weights.nii ' Params.dataSaveDir 'gre_weights.nii'] ) ;  
system( ['rm -r ' Params.tmpSaveDir] ) ;

Mask = load_untouch_nii( [ Params.dataSaveDir 'gre_seg.nii' ] );
mask = Mask.img ;
mask = logical(permute( mask, [2 1 3] )) ;
mask = flipdim( mask, 1 ) ;

Weights = load_untouch_nii( [ Params.dataSaveDir 'gre_weights.nii' ] );
weights = Weights.img ;
weights = double(permute( weights, [2 1 3] )) ;
weights = flipdim( weights, 1 ) ;
% normalize
weights = weights - min(weights(:)) ;
weights = weights/max(weights(:)) ;

end
% =========================================================================
function [ Hdr ] = dicominfosiemens( filename )
%DICOMINFOSIEMENS
%
% Hdr = DICOMINFOSIEMENS
%
% Wraps to Matlab's DICOMINFO to return the standard DICOM Hdr, but also
% add fields returned by PARSE_SIEMENS_SHADOW.
%
% (Also replaces original Hdr. fields with more precise values should they
% be available in Hdr.MrProt., e.g. In the case of Siemens socket-delivered
% images, some fields are missing, & others are simply of reduced precision, 
% e.g. Hdr.ImagingFrequency is rounded to kHz)
%
% See also DICOMINFO, PARSE_SIEMENS_SHADOW

Hdr = dicominfo( filename ) ;

% parses additional info (e.g. "MrProt" (re: protocol))
[ Hdr.Img, Hdr.Ser, Hdr.MrProt ] = parse_siemens_shadow( Hdr ) ;
fprintf('\n') ;

% parse_mrprot produces many, many warning messages (at least for DICOMs from the Prisma).
% This should suppress all but the 1st instance:
[ lastMsg, lastWarnId ] = lastwarn( ) ;
warning( 'off', lastWarnId ) ;

% for some reason (?) Img.Hdr.MrProt.sTXSPEC.asNucleusInfo.lFrequency appears
% to have multiple elements, only the first is non-empty and it contains the
% Larmor freq., but, explicitly/manually trying to access it produces an error.
% +for some reason (?) the error is avoided by first copying the value to the
% temporary variable (f0), hence:
f0 = Hdr.MrProt.sTXSPEC.asNucleusInfo.lFrequency  ;
Hdr.ImagingFrequency = f0 * 1E-6 ; % [units : MHz]

% absolute table position (.Hdr field doesn't exist for the socket-delivered .dcm images)
if ~myisfield( Hdr.Img, 'Private_0019_1013' ) && ~isempty( Hdr.Img.ImaAbsTablePosition )
    Img.Hdr.Private_0019_1013 = Hdr.Img.ImaAbsTablePosition ; % [units : mm?]
end

end
% =========================================================================

   
end
% =========================================================================
% =========================================================================


end
