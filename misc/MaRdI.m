classdef MaRdI < matlab.mixin.SetGet 
%MaRdI Ma(t)-R-dI(com)
%
%   Dicom iunto Matlab for Siemens MRI data
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
% RESLICEIMG()
%   griddata takes too long (possible solution: write interp function in cpp?)
%   rename to regridimg()?
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

properties(SetAccess=protected)
    Hdr ; % full Siemens DICOM header of 1st img (i.e. Img.img(:,:,1) )
    Hdrs ; % cell array of (truncated) DICOM headers courtesy of dicominfo()
end

properties(SetAccess=protected, Hidden = true)
    Ref ; % Reference properties - prior to manipulation
end

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
function [f0, g0, s0] = adjvalidateshim( Img )
%ADJVALIDATESHIM
% 
% [f0, g0, s0] = ADJVALIDATESHIM( Img )
% 
% ADJVALIDATESHIM is not a particularly revealing name for a function; however,
% it is based on the Siemens AdjValidate commandline tool, with a "-shim" input
% argument.
%
% f0 = scalar Larmor (Tx) freq. [ units : Hz ]
% g0 = 3-element vector of the linear gradient offsets (gX, gY, gZ) [units : DAC bits]
% s0 = 5-element vector of 2nd order shim currents [units : mA] 
%
% To convert to the units of the 3D Shim card on the Siemens (Prisma) console,
% see 
% ShimSpecs_IUGM_Prisma_fit.converttomultipole( )
% ShimSpecs_HGM_Prisma_fit.converttomultipole( )

% Larmor (Tx) freq.
f0 = Img.Hdr.MrProt.sTXSPEC.asNucleusInfo.lFrequency ; 

% linear gradients
g0 = [ Img.Hdr.MrProt.sGRADSPEC.asGPAData.lOffsetX ; ...
       Img.Hdr.MrProt.sGRADSPEC.asGPAData.lOffsetY ; ...
       Img.Hdr.MrProt.sGRADSPEC.asGPAData.lOffsetZ ] ;

% -------
% 2nd order shims (5-element vector)
s0 = Img.Hdr.MrProt.sGRADSPEC.alShimCurrent ;

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
        Img.Aux.Data.p       = interp1( Img.Aux.Data.t, Img.Aux.Data.p, tInterp, Params.interpolationMethod, 'extrap' )' ; 
        Img.Aux.Data.pRaw    = interp1( Img.Aux.Data.t, Img.Aux.Data.pRaw, tInterp, Params.interpolationMethod, 'extrap' )' ; 
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
% =========================================================================
function ImgCopy = copy(Img)
%COPY 
% 
% Make a copy of a MaRdI (i.e. handle) object.
% 
% ImgCopy = Copy( Img ) ;

ImgCopy     = MaRdI() ;

ImgCopy.img  = Img.img;
ImgCopy.Hdr  = Img.Hdr ;
ImgCopy.Hdrs = Img.Hdrs ;

if ~isempty( Img.Aux ) 
    if isa( Img.Aux, 'ProbeTracking' ) 
        ImgCopy.Aux = Img.Aux.copy() ;
    else
        error('Expected Img.Aux to be of type ProbeTracking') ;
    end
end

end
% =========================================================================
function [t0] = estimatekorigintime( Img ) 
%ESTIMATEKORIGINTIME
% 
% t0 = ESTIMATEKORIGINTIME( Img )
%
% Returns an estimate of when the k-space origin of an image was sampled
% relative to the AcquisitionTime (field in Siemens DICOM header) as a double
% in units of milliseconds. 
% 
% See also: MaRdI.getacquisitiontime()
% 
% NOTE: This is a crude estimate and only the case of Cartesian k-space
% sampling, beginning at the k_min periphery, has been considered in the
% current implementation!
 
nSlices = Img.getnumberofslices() ;
nEchoes = length( Img.getechotime() ) ;

if myisfield( Img.Hdr.MrProt, 'lRepetitions' ) 
    nMeasurements = ( Img.Hdr.MrProt.lRepetitions + 1) ;
else
    nMeasurements = 1;
end

assert( nMeasurements == size( Img.img, 5 ), 'Invalid .Hdr' ) ;

tAcq = Img.getacquisitiontime() ;

% Estimate time from excitation to k-space origin:
dt = 0; % [units: ms]

if Img.Hdr.MrProt.sKSpace.ucTrajectory == 1
    
    if ~strcmp( Img.Hdr.ScanningSequence, 'EPI' )
    % NOTE: The Siemens DICOM field SliceMeasurementDuration appears to be a misnomer:
    % Rather, it (usually?) refers to the duration of an entire volume. 
        dt = Img.Hdr.Img.SliceMeasurementDuration ; % [units: ms]
    else
    % *Exception*: Apparently EPI are handled differently as said DICOM field
    % is oddly large (probably includes dummy volumes?)
        dt = Img.Hdr.RepetitionTime/nSlices ; % [units: ms]
    end
    
    pf = Img.getpartialfourierfactors() ;
    
    if all( pf == 1 )
        dt = dt/2 ;
    else
        switch Img.Hdr.MRAcquisitionType 
            case '2D'
                dt = dt*( 3/2 - pf(2) ) ; 
            case '3D' % not sure what the function is here...
                error('Not implemented!') ;
        end
    end

else    
    display('Estimating time at which the k-space origin was sampled')
    warning(' Current strategy is simplistic and likely wrong for certain encoding trajectories.')
    dt = 0 ;
end

assert( dt >= 0, 'Unexpected value for time-to-origin since excitation (AcquisitionTime). See code + check DICOM entries are valid.' )

t0 = tAcq + dt ;

end
% =========================================================================
function [is] = exist( Img )
%EXIST
is = true ;
end
% =========================================================================
function [] = filter( Img, weights, Params )
%FILTER
%
% 3D low-pass (weighted or unweighted) or median filtering.
% Wraps to smooth3() or medfilt3() accordingly.
%
% [] = FILTER( Img )
% [] = FILTER( Img, weights )
% [] = FILTER( Img, weights, Params )
%
% Img 
%   the MaRdI-type image volume.
%
% weights
%   an array of data weights (>=0) to penalize (for smooth3) or exclude (for
%   medfilt3()) identifiable outliers.  Dimensions of weights must be the
%   same as the image volume (i.e. Img.getgridsize() )
%
% Params
%   an optional struct for which the following Params.fields are supported
%
%   .kernelSize
%       in number of voxels
%       default = [3 3 3] 
%
%   .method
%       'gaussian' OR 'box' OR 'median'
%       default = 'gaussian'
%
% 
% TODO
%   Add support for 2d (single slice) images 

DEFAULT_KERNELSIZE = [3 3 3] ;
DEFAULT_METHOD     = 'gaussian' ;

if nargin < 3 || isempty(Params)
    Params.dummy = [] ;
end

Params = assignifempty( Params, 'kernelSize', DEFAULT_KERNELSIZE ) ;
Params = assignifempty( Params, 'method', DEFAULT_METHOD ) ;

nVolumes = Img.getnumberofmeasurements() ;
nEchoes  = Img.getechotime() ;

if nargin < 2 || isempty( weights )
    weights = [ ones( Img.getgridsize() ) nEchoes nVolumes ] ;
else
    assert( all( size(weights) == size(Img.img) ), ...
        'Filter weights and image volume must possess the same dimensions' ) ;
end

for iVolume = 1 : nVolumes
    for iEcho = 1 : nEchoes

        img = Img.img(:,:,:, iEcho, iVolume ) ;
        w   = weights(:,:,:, iEcho, iVolume ) ;

        switch Params.method
            case 'median'
                img( ~w ) = NaN ; % medfilt3() will ignore these values
                img = medfilt3( img, Params.kernelSize ) ; 
                img( ~w ) = 0 ;
            otherwise 
                weightsSmoothed = smooth3( w, Params.method, Params.kernelSize ) ;
                Img.img = smooth3( w .* img, Params.method, Params.kernelSize ) ./ weightsSmoothed ; 
        end
        
        Img.img(:,:,:, iEcho, iVolume) = img ;
    end
end

end
% =========================================================================
function [t] = getacquisitiontime( Img ) 
% GETACQUISITIONTIME
% 
% t = GETACQUISITIONTIME( Img )
%
% Derives from the AcquisitionTime field in Siemens DICOM header to return
% an array of doubles describing the milliseconds elapsed since midnight
%
% dimensions of t: [ nSlices x nEchoes x nMeasurements ]
%
% t - t(1) yields the elapsed time since acquisition of the 1st k-space point
% in the series.
%
% For EPI MOSAIC (which has a single AcquisitionTime value for each volume),
% t( iSlice ) = AcquisitionTime (first slice) + iSlice*(volumeTR/nSlices)

nSlices = Img.getnumberofslices() ;
nEchoes = length( Img.getechotime() ) ;

if myisfield( Img.Hdr.MrProt, 'lRepetitions' ) 
    nMeasurements = ( Img.Hdr.MrProt.lRepetitions + 1) ;
else
    nMeasurements = 1;
end

%assert( nMeasurements == size( Img.img, 5 ), 'Invalid .Hdr' ) ; % EAO:
%commented this out so I could save a time series

t = zeros( nSlices, nEchoes, nMeasurements ) ;

if isempty(strfind( Img.Hdrs{1}.ImageType, 'MOSAIC' ) )
    
    for iMeasurement = 1 : nMeasurements
        for iEcho = 1 : nEchoes 
            for iSlice = 1 : nSlices
                % RT: commented 20191013, added convertacquisitiontime() -- which is probably more correct/accurate?
                %
                % t_iAcq = Img.Hdrs{iSlice,iEcho,iMeasurement}.AcquisitionTime ;
                % t_iAcq = str2double( t_iAcq(1:2) )*60*60 + str2double( t_iAcq(3:4) )*60 + str2double( t_iAcq(5:end) ) ;
                % t(iSlice, iEcho, iMeasurement) = t_iAcq ;
                
                t(iSlice, iEcho, iMeasurement) = convertacquisitiontime( Img.Hdrs{iSlice,iEcho,iMeasurement}.AcquisitionTime ) ;
            end
        end
    end

else % special case for EPI MOSAIC (single DICOM header for multiple slices)
   
    sliceOrder = Img.Hdr.MrProt.sSliceArray.alSliceAcqOrder ;
    sliceMeasurementDuration = Img.Hdr.RepetitionTime/nSlices ; % [units: ms]

    for iMeasurement = 1 : nMeasurements
        for iEcho = 1 : nEchoes 
            for iSlice = 1 : nSlices
                % See above comment 20191013
                % t_iAcq = Img.Hdrs{1,iEcho,iMeasurement}.AcquisitionTime  ;
                % t_iAcq = str2double( t_iAcq(1:2) )*60*60 + str2double( t_iAcq(3:4) )*60 + str2double( t_iAcq(5:end) ) ;
                t_iAcq = convertacquisitiontime( Img.Hdrs{1,iEcho,iMeasurement}.AcquisitionTime ) ;
                t_iAcq = t_iAcq + sliceOrder(iSlice)*sliceMeasurementDuration ;
                t(iSlice, iEcho, iMeasurement) = t_iAcq ;
            end
        end
    end

end

function t_s = convertacquisitiontime( AcquisitionTime )
%CONVERTACQUISITIONTIME
% 
% Converts the time-stamp string in the .AcquisitionTime DICOM header
%
% For details, see https://cfn.upenn.edu/aguirre/wiki/public:pulse-oximetry_during_fmri_scanning

% hours + minutes + seconds + fractional seconds :
t_s = str2double( AcquisitionTime(1:2) )*60*60 + str2double( AcquisitionTime(3:4) )*60 + ...
      str2double( AcquisitionTime(5:6) ) + (str2double( AcquisitionTime(7:10) )/10)/1000 ;
t_s = t_s * 1000 ; % [units: ms]

end

end
% =========================================================================
function [f0] = getimagingfrequency( Img ) 
%GETIMAGINGFREQUENCY    Returns Larmor frequency in Hz
% 
% f0 = GETIMAGINGFREQUENCY( Img ) ;  

f0 = Img.Hdr.MrProt.sTXSPEC.asNucleusInfo.lFrequency ;

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
function [nVolumes] = getnumberofmeasurements( Img )
%GETNUMBEROFMEASUREMENTS    Returns the number of measurements
%
% n = GETNUMBEROFMEASUREMENTS( Img )
%
% NOTE
%
% GETNUMBEROFMEASUREMENTS( Img ) is equivalent to n = size( Img.img, 5 ),
% however GETNUMBEROFMEASUREMENTS also checks the DICOM header in Img.Hdr and
% issues a warning if n differs from the expected value
% (Img.Hdr.MrProt.lRepetitions +1).

nVolumes = size( Img.img, 5 ) ;

if myisfield( Img.Hdr.MrProt, 'lRepetitions' ) 
    % number of measurements according to hdr:
    nVolumesHdr = Img.Hdr.MrProt.lRepetitions + 1 ;
    if nVolumes ~= nVolumesHdr
        warning([ num2str(nVolumes) ' image volumes have been loaded, however the DICOM Hdr indicates ' num2str(nVolumesHdr) ' were acquired.'] ) ;
    end 
end

end
% =========================================================================
function [nSlices] = getnumberofslices( Img ) 
%GETNUMBEROFSLICES  Returns number of acquired slices
% 
% nSlices = GETNUMBEROFSLICES( Img ) 
%
% NOTE: nSlices is not necessarily equal to size( Img.img, 3).
% e.g. For a 3d (slab) encoding, GETNUMBEROFSLICES returns 1. 

nSlices = Img.Hdr.MrProt.sSliceArray.lSize ;

end
% =========================================================================
function [pff] = getpartialfourierfactors( Img ) 
%GETPARTIALFOURIERFACTORS  Returns fraction of k-space coverage in each dim
% 
% pff = GETPARTIALFOURIERFACTORS( Img ) 
% 
% Returns 3-element vector of partial Fourier factors in read, phase (in-plane),
% and partition (slice) encoding directions.
%
% NOTE:
% Siemens uses an enumeration scheme to store Partial Fourier info in the DICOM header:
%
% pff --> Siemens DICOM Hdr entry 
% 4/8 --> 0x1  = 1
% 5/8 --> 0x2  = 2
% 6/8 --> 0x4  = 4
% 7/8 --> 0x8  = 8
% 8/8 --> 0x10 = 16 (i.e. no partial Fourier)
% 
% See: https://github.com/malaterre/GDCM/blob/master/Source/DataDictionary/CSAHeader.xml

dcmFields = { 'Img.Hdr.MrProt.sKSpace.ucReadoutPartialFourier' ;
              'Img.Hdr.MrProt.sKSpace.ucPhasePartialFourier' ;
              'Img.Hdr.MrProt.sKSpace.ucSlicePartialFourier' ; } ;

pfAsInt = [ Img.Hdr.MrProt.sKSpace.ucReadoutPartialFourier
            Img.Hdr.MrProt.sKSpace.ucPhasePartialFourier
            Img.Hdr.MrProt.sKSpace.ucSlicePartialFourier ]' ;

pff     = [ 1 1 1 ] ;

for iDim = 1 : 3
    switch pfAsInt(iDim)
        case 1
            pff(iDim) = 4/8 ;
        case 2
            pff(iDim) = 5/8 ;
        case 4
            pff(iDim) = 6/8 ;
        case 8
            pff(iDim) = 7/8 ;
        case 16
            pff(iDim) = 1 ;
        otherwise
            error( [ 'Unexpected value of ' num2str( pfAsInt(iDim) ) ' in ' dicomFields{iDim} ] ) ;
    end
end

end
% =========================================================================
function [echoTime] = getechotime( Img, iEcho ) 
%GETECHOTIME
% 
% TE = GETECHOTIME( Img )
% TE = GETECHOTIME( Img, iEcho )
%
% Returns vector of echo times in units of ms.
% If 2nd argument (echo index iEcho) is provided, GETECHOTIME returns the TE of
% the corresponding echo.

if nargin == 1
    if strcmp( Img.Hdr.SequenceName, '*fm2d2' ) && Img.isphase() 
        echoTime = ( Img.Hdr.MrProt.alTE(2) - Img.Hdr.MrProt.alTE(1) )/1000  ;
    else
        nEchoes  = size( Img.img, 4 ) ;
        echoTime = Img.Hdr.MrProt.alTE(1:nEchoes)/1000 ;
    end
else
    echoTime = Img.Hdr.MrProt.alTE(iEcho)/1000 ;
end 

end
% =========================================================================
function fieldOfView = getfieldofview( Img )
%GETFIELDOFVIEW
% 
% fov = GETFIELDOFVIEW( Img ) ;
%
% Returns field of view in units of mm : [Row Column Slice] dimensions
fieldOfView = [ Img.Hdr.PixelSpacing(1) * double( Img.Hdr.Rows ), ...
                Img.Hdr.PixelSpacing(2) * double( Img.Hdr.Columns ), ...
                Img.Hdr.SpacingBetweenSlices * size( Img.img, 3 ) ] ;
end
% =========================================================================
function gridSize = getgridsize( Img )
%GETGRIDSIZE    Image dimensions as 3-element vector (rows, columns, slices)
% 
% gridSize = GETGRIDSIZE( Img ) 

gridSize = size( Img.img ) ;

if numel( gridSize ) < 2
    error('Img.img should be at least 2d') ;
elseif numel( gridSize ) < 3
    gridSize = [gridSize 1] ;
else
    gridSize = gridSize(1:3) ;
end

end
% =========================================================================
function GYRO = getgyromagneticratio( Img )
%GETGYROMAGNETICRATIO
%
% Gyro = getgyromagneticratio( Img )
%
% Examines .Hdr of MaRdI-type Img for .ImagedNucleus and returns gyromagnetic
% ratio in units of rad/s/T.

switch Img.Hdr.ImagedNucleus 
    case '1H' 
        GYRO = 267.513E6 ; 
    otherwise
        error('Not implemented.') ;
end

end
% =========================================================================
function nVoxels = getnumberofvoxels( Img )
%GETNUMBEROFVOXELS
nVoxels = prod( Img.getgridsize( ) ) ;
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
function [X,Y,Z] = getvoxelpositions( Img )
% GETVOXELPOSITIONS
% 
% [X,Y,Z] = GETVOXELPOSITIONS( Img ) 
%
%   Returns three 3D arrays of doubles, each element containing the
%   location [units: mm] of the corresponding voxel with respect to 
%   DICOM's 'Reference Coordinate System'.

% from DICOM standard: https://www.dabsoft.ch/dicom/3/C.7.6.2.1.1/
%
% If Anatomical Orientation Type (0010,2210) is absent or has a value of
% BIPED, the x-axis is increasing to the left hand side of the patient. The
% y-axis is increasing to the posterior side of the patient. The z-axis is
% increasing toward the head of the patient.
%
% Arrays containing row, column, and slice indices of each voxel
assert( ~myisfield( Img.Hdr, 'AnatomicalOrientationType' ) || ...
        strcmp( Img.Hdr.AnatomicalOrientationType, 'BIPED' ), ...
        'Error: AnatomicalOrientationType not supported.' ) ;
        
[iRows,iColumns,iSlices] = ndgrid( [0:1:Img.Hdr.Rows-1], ...
                  [0:1:Img.Hdr.Columns-1], ...
                  [0:1:(size(Img.img,3)-1)] ) ; 

iRows    = double(iRows);
iColumns = double(iColumns);
iSlices  = double(iSlices);

%-------
% Rotation matrix: R
[r, c, s] = Img.getdirectioncosines ;
R = [r c s] ; 
                  
%-------
% Scaling matrix: S  
voxelSize = Img.getvoxelspacing() ;
S = diag(voxelSize); 

RS = R*S ;

%-------
% Scale and rotate to align row direction with x-axis, 
% column direction with y-axis, slice with z-axis
X1 = RS(1,1)*iRows + RS(1,2)*iColumns + RS(1,3)*iSlices;
Y1 = RS(2,1)*iRows + RS(2,2)*iColumns + RS(2,3)*iSlices;
Z1 = RS(3,1)*iRows + RS(3,2)*iColumns + RS(3,3)*iSlices;

%-------
% TRANSLATE w.r.t. origin (i.e. location of 1st element: .img(1,1,1))
X = Img.Hdr.ImagePositionPatient(1) + X1 ; 
Y = Img.Hdr.ImagePositionPatient(2) + Y1 ; 
Z = Img.Hdr.ImagePositionPatient(3) + Z1 ; 

end
% =========================================================================
function h = getvoxelspacing( Img )
%GETVOXELSPACING
% 
% h = GETVOXELSPACING( Img )
%       
% Returns 3-component grid-spacing vector [units: mm]:
%
%   h(1) : row spacing (between centers of adjacent rows, i.e. vertical spacing). 
%
%   h(2) : column spacing (between the centers of adjacent columns,
%   i.e. horizontal spacing).    
%
%   h(3) : slice spacing (between the centers of adjacent slices, i.e.
%   from the DICOM hdr, this is Hdr.SpacingBetweenSlices - for a 2D acquisition
%   this not necessarily the same as Hdr.SliceThickness).    

h = [ Img.Hdr.PixelSpacing(1) Img.Hdr.PixelSpacing(2) Img.Hdr.SpacingBetweenSlices ] ;

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
function [F] = resliceimg( Img, X_Ep, Y_Ep, Z_Ep, varargin ) 
%RESLICEIMG
% 
% Interpolate a MaRdI image (Img.img) and update Img.Hdr accordingly.
% 
% In general, RESLICEIMG() uses MATLAB's scatteredInterpolant class. 
% The exception is when the image input (Img.img) is 2d and the target
% output (prescribed by inputs X,Y,Z) is a volume. This scenario is
% incompatible with scatteredInterpolant, and nearest-neighbor substitution is
% used instead.
%
% -----   
% Basic Usage
%
% [] = RESLICEIMG( Img, X, Y, Z )
% [] = RESLICEIMG( Img, X, Y, Z, mask )
% 
% Inputs:
%
% X, Y, Z:  
%       2d or 3d arrays (size=output image grid) describing the X, Y, Z patient
%       coordinates (i.e. of the DICOM reference coordinate system) of the
%       target (output) voxels. In general, if one is interpolating from one
%       image grid (Img) to another (MaRdI-type object Img2), these arrays are
%       obtained by the call: [X,Y,Z] = Img2.getvoxelpositions()
% 
% mask: [Optional, default = true(size output image grid)]
%       A logical array (size=output image grid) specifying the subset of the
%       output voxels that are of interest. (i.e. voxels in the output image
%       with a corresponding mask entry == FALSE will simply be assigned zero).
%       Note: Specifying the region of interest for extrapolation with this
%       variable can greatly accelerate the interpolation!
%
% -----   
%
% Advanced Usage TODO
%
%   [F] = RESLICEIMG( Img, X, Y, Z, mask, F ) 
% 
%   case:
%       interpolationMethod [default='linear']
%       is a string supported by the scatteredInterpolant constructor.
%   F is the object of type 'scatteredInterpolant' used for interpolation.


%% -----
%NOTE: Terminology: 
% 'Ip' = interpolant/initial point
% 'Ep' = extrapolant/end point

F = [] ;

DEFAULT_INTERPOLATIONMETHOD  = 'linear' ;
DEFAULT_ISFORMINGINTERPOLANT = true ;

[X_Ip, Y_Ip, Z_Ip] = Img.getvoxelpositions( ) ;

%% -----
% Parse and check inputs
if MaRdI.compareimggrids( X_Ip, Y_Ip, Z_Ip, X_Ep, Y_Ep, Z_Ep ) 
    warning('Voxel positions are already identical. Not interpolating.');
    return ;
end
        
isFormingInterpolant = true ; 

if nargin < 5

    interpolationMethod  = DEFAULT_INTERPOLATIONMETHOD ;
    isFormingInterpolant = DEFAULT_ISFORMINGINTERPOLANT ;

elseif nargin >= 5 
    if islogical( varargin{1} ) ;
        maskEp = varargin{1} ;
    end

    if nargin == 6
        if ischar( varargin{2} )
            interpolationMethod = varargin{2} ;
        elseif isa( varargin{2}, 'scatteredInterpolant' ) ;
            F = varargin{2} ;
        else
            error( 'Unknown input. See HELP MaRdI.resliceimg' ) ;
        end
    end
end

isUsingScatteredInterpolant = [] ;

if (ndims(Img.img) > 1) && (ndims(Img.img) <= 5)

    gridSizeIp = Img.getgridsize() ;
    
    if gridSizeIp(3) > 1
        isUsingScatteredInterpolant = true ;
    else
        isUsingScatteredInterpolant = false ;
    end
else
    error('Dimensions of input Img.img must be >= 2, and <= 5') ;
end

if ndims( X_Ep ) == 2 % interpolating down to 2d single-slice
    gridSizeEp = [ size(X_Ep) 1 ] ;
elseif ndims( X_Ep ) == 3
    gridSizeEp = size(X_Ep) ;
else
    error('Expected 2d or 3d target interpolation grid')
end

%% -----
% Define variables
gridSizeIp = size( X_Ip ) ;

if myisfieldfilled( Img.Hdr, 'MaskingImage' ) 
    maskIp = logical( sum( sum( Img.Hdr.MaskingImage, 5 ), 4 ) ) ;
else
    maskIp = true( gridSizeIp ) ;
end

gridSizeEp = size( X_Ep ) ;
if ndims(gridSizeEp) == 2
    gridSizeEp = [gridSizeEp 1] ;
end

if ~exist('maskEp')
    if ~isUsingScatteredInterpolant 
        warning('No logical mask provided: For faster results, restrict the target/output voxels to those of interest by providing this mask!') ;
    end
    maskEp = true( gridSizeEp ) ;
end

iEp = find( maskEp(:) ) ; % indices of target voxels

iNearest = zeros( size(iEp) ) ;
nEp      = length( iEp ) ;

nImgVolumesDim4 = size(Img.img, 4 ) ; % nEchoes
nImgVolumesDim5 = size(Img.img, 5 ) ; % nMeasurements
nImgVolumes     = nImgVolumesDim4 * nImgVolumesDim5 ;

imgOut = zeros( [gridSizeEp nImgVolumesDim4 nImgVolumesDim5] ) ;

%% -------
if isFormingInterpolant 
    tic
    disp( 'Forming interpolant...' )
    disp( '(Computation time depends on input image size. This may take a few minutes.)' ) ;

    if isUsingScatteredInterpolant

        % The following avoids the error from scatteredInterpolant when one
        % attempts to form a 3d interpolant from a 2d input: 
        isValidDim0 = [ numel(unique(X_Ip(:))) numel(unique(Y_Ip(:))) numel(unique(Z_Ip(:))) ] > 1 ;
        r0          = [X_Ip(:) Y_Ip(:) Z_Ip(:)] ;
        
        isValidDim1 = [ numel(unique(X_Ep(:))) numel(unique(Y_Ep(:))) numel(unique(Z_Ep(:))) ] > 1 ;
        r1          = [X_Ep(:) Y_Ep(:) Z_Ep(:)] ;

        if nnz( isValidDim0 ) == 2
            
            assert( all( isValidDim1 == isValidDim0 ), ... 
                'Query points should sit within the same plane as the interpolant points' ) ;
                
            % coordinate of interpolant plane along normal dim:
            qn0 = unique( r0(:, ~isValidDim0) ) ;
            % coordinate of query plane along same dim:
            qn1 = unique( r1(:, ~isValidDim1) ) ;

            % This could instead be a warning? (e.g. if the 2d planes are indeed very close, should interp still be performed?)
            assert( qn0 == qn1, 'Query points should sit within the same plane as the interpolant points' ) ;

            % exclude the coordinate along the normal dim from the interpolation
            r1 = r1(:, isValidDim1) ; 

        end
        
        F                     = scatteredInterpolant() ;
        F.Method              = interpolationMethod ;
        F.ExtrapolationMethod = 'none' ;
        F.Points              = r0(:, isValidDim0) ;
    
    else % Map nearest neighbours
        
        % truncate IP voxels
        X_Ip = X_Ip(maskIp) ;
        Y_Ip = Y_Ip(maskIp) ;
        Z_Ip = Z_Ip(maskIp) ;

        for iR = 1 : nEp
            [~,iNearest(iR)] = min( sqrt( ( X_Ip - X_Ep( iEp(iR) ) ).^2 + ...
                                   ( Y_Ip - Y_Ep( iEp(iR) ) ).^2 + ...
                                   ( Z_Ip - Z_Ep( iEp(iR) ) ).^2 ) ) ;
        end
    end
    toc
end

%% -----
disp('Reslicing...')
tic
for iImgDim4 = 1 : nImgVolumesDim4
    for iImgDim5 = 1 : nImgVolumesDim5
        disp( ['Reslicing image volume...' num2str(iImgDim4*iImgDim5) ' of ' num2str(nImgVolumes) ]) ;
       
        imgIp = Img.img(:,:,:, iImgDim4, iImgDim5 ) ;
        imgEp = zeros( gridSizeEp ) ;
      
        if isUsingScatteredInterpolant  
            F.Values = imgIp(:) ;
            imgEp    = reshape( F( r1 ), gridSizeEp ) ;
        
        else % Nearest-neighbor substitution
            imgIp = imgIp(maskIp) ;
            for iR = 1 : nEp 
                imgEp( iEp(iR) ) = imgIp( iNearest(iR) ) ;
            end
        end

        imgOut(:,:,:, iImgDim4, iImgDim5 ) = imgEp ;
    end
end
toc

imgOut( isnan( imgOut ) ) = 0 ; 

Img.img = imgOut ; 

Img.Hdr.MaskingImage = Img.img ~= 0 ;
%% -----------------------------------------------------------------------

% ------------------------------------------------------------------------
% Update header
Img.Hdr.ImageType = 'DERIVED\SECONDARY\REFORMATTED' ;


Img.Hdr.ImagePositionPatient( 1 ) = X_Ep(1) ; 
Img.Hdr.ImagePositionPatient( 2 ) = Y_Ep(1) ;
Img.Hdr.ImagePositionPatient( 3 ) = Z_Ep(1) ;

%-------
% Rows 
Img.Hdr.Rows = size(Img.img, 1) ;

dx = X_Ep(2,1,1) - X_Ep(1,1,1) ;
dy = Y_Ep(2,1,1) - Y_Ep(1,1,1) ;
dz = Z_Ep(2,1,1) - Z_Ep(1,1,1) ;  

% vertical (row) spacing
Img.Hdr.PixelSpacing(1) = ( dx^2 + dy^2 + dz^2 )^0.5 ; 

% column direction cosine (expressing angle btw column direction and X,Y,Z axes)
Img.Hdr.ImageOrientationPatient(4) = dx/Img.Hdr.PixelSpacing(1) ;
Img.Hdr.ImageOrientationPatient(5) = dy/Img.Hdr.PixelSpacing(1) ;
Img.Hdr.ImageOrientationPatient(6) = dz/Img.Hdr.PixelSpacing(1) ;

%-------
% Columns 
Img.Hdr.Columns = size(Img.img, 2) ;       

dx = X_Ep(1,2,1) - X_Ep(1,1,1) ;
dy = Y_Ep(1,2,1) - Y_Ep(1,1,1) ;
dz = Z_Ep(1,2,1) - Z_Ep(1,1,1) ;  

% horizontal (column) spacing
Img.Hdr.PixelSpacing(2) = ( dx^2 + dy^2 + dz^2 )^0.5 ;

% row direction cosine (expressing angle btw column direction and X,Y,Z axes)
Img.Hdr.ImageOrientationPatient(1) = dx/Img.Hdr.PixelSpacing(2) ;
Img.Hdr.ImageOrientationPatient(2) = dy/Img.Hdr.PixelSpacing(2) ;
Img.Hdr.ImageOrientationPatient(3) = dz/Img.Hdr.PixelSpacing(2) ;

%-------
% Slices

if size( Img.img, 3 ) > 1
    Img.Hdr.SpacingBetweenSlices = ( (X_Ep(1,1,2) - X_Ep(1,1,1))^2 + ...
                                     (Y_Ep(1,1,2) - Y_Ep(1,1,1))^2 + ...
                                     (Z_Ep(1,1,2) - Z_Ep(1,1,1))^2 ) ^(0.5) ;
else
    Img.Hdr.SpacingBetweenSlices = 0 ;
end

Img.setslicenormalvector() ; % redefines sHat

[rHat, cHat, sHat] = Img.getdirectioncosines( ) ;  
Img.Hdr.SliceLocation = dot( Img.Hdr.ImagePositionPatient, sHat ) ;

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
function [] = setmaskingimage( Img, mask )
%SETMASKINGIMAGE
% 
% [] = SETMASKINGIMAGE( Img, mask ) 
%
% Copies valid mask (a logical array of 1's and 0's) to Img.Hdr.MaskingImage
%
% The purpose of this function is to specify the signal spatial support (e.g. 
% of mag, phase, field data) within the image grid. e.g. it is called during 
% phase unwrapping/field mapping, but it might also be called prior to regridding
% if the interpolation should exclude certain voxels.
% 
% To be valid, mask must either be the same size as Img.img OR the same size as
% Img.getgridsize() (i.e. the size of a single image volume of a
% multi-echo/multi-measurement stack), in which case, the single mask is simply
% copied such that the assigned Img.Hdr.MaskingImage always possesses the same
% dimensions as Img.img

assert( nargin == 2, 'Function requires at least 2 input arguments. See HELP MaRdI.setmaskingimage' ) ;

if ~islogical( mask )
    assert( numel(mask) == ( nnz( mask(:) == 0 ) + nnz( mask(:) == 1 ) ), ...
        'Input mask must consist exclusively of zeros and ones. See HELP MaRdI.setmaskingimage' ) ;
    mask = logical( mask ) ;
end

maskSize = size( mask ) ;

%% -----
% assign masking image 
if ndims( Img.img ) == ndims( mask ) && all( size(Img.img) == maskSize )
    Img.Hdr.MaskingImage = mask ;

elseif ndims( Img.img ) > ndims( mask ) 
    
    if numel( maskSize ) == 2
        maskSize = [ maskSize 1 ] ;
    end

    if ( numel( Img.getgridsize() ) == numel( maskSize ) ) && ... 
       ( all( Img.getgridsize() == maskSize ) ) % single mask provided
        Img.Hdr.MaskingImage = repmat( mask, [1 1 1 size(Img.img, 4) size(Img.img, 5)] ) ;
    else
        error('Input mask and Img.img should possess the same dimensions' ) ;
    end
else 
    error('Input mask and Img.img should possess the same dimensions' ) ;
end

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
function [] = unwrapphase( Phase, varargin )
%UNWRAPPHASE
%
% Interface to SUNWRAP, to FSL Prelude, or to Abdul-Rahman's path-based phase unwrapper
%
% .......
%   
% Usage
%
%   [] = UNWRAPPHASE( Phase )
%   [] = UNWRAPPHASE( Phase, Mag )
%   [] = UNWRAPPHASE( Phase, Mag, Options )
%   
%   Phase and Mag are objects of type MaRdI.
% 
%   Options is a struct that can contain the following fields:
%
%   .unwrapper 
%       == 'AbdulRahman_2007' [default if Phase.img is 3D] : 
%           calls UNWRAP3D. 
%           Not permitted if Phase.img is 2D (will default to SUNWRAP)
%           If Mag is supplied, Phase.Hdr.MaskingImage must
%           See HELP UNWRAP3D for description of permitted Options 
%
%       == 'Sunwrap' [default if Phase.img is 2D] : 
%           calls SUNWRAP. See HELP SUNWRAP for more details.
%
%       == 'FslPrelude' : calls PRELUDE. 
%           See HELP PRELUDE for description of permitted Options 
% 
%   .threshold  [default = 0.01] 
%       Relative threshold of Mag.img used to define the unwrapping region in
%       the SUNWRAP case and for the other 2 cases when Phase.Hdr.MaskingImage
%       is undefined (Mag.getreliabilitymask() is called).
%
%   NOTE: UNWRAP3D and PRELUDE support an Options.mask input to which
%   Phase.Hdr.MaskingImage will always be assigned. To assign the mask
%   manually, run Phase.setmaskingimage before calling Phase.unwrapphase. 


%% ------
% check inputs, assign defaults:
assert( strcmp( Phase.Hdr.PixelComponentPhysicalUnits, '0000H' ), 'SCALEPHASE2RAD() before unwrapping.' ) ;

if nargin > 1
    for iArg = 1 : length( varargin )
        switch class( varargin{iArg} )
            case 'MaRdI'
                Mag = varargin{iArg} ;
                assert( Mag.ismagnitude(), 'See HELP MaRdI.unwrapphase' ) ;
            case 'struct'
                Options = varargin{iArg} ;
            otherwise
                error('Unrecognized input. See HELP MaRdI.unwrapphase') ;
        end
    end
end

DEFAULTS.threshold     = 0.01 ;
DEFAULTS.isLinearizing = false ;

if ( size( Phase.img, 3 ) > 1 )
    is3d = true ;
    DEFAULTS.unwrapper = 'AbdulRahman_2007' ;    
else
    is3d = false ;
    DEFAULTS.unwrapper = 'Sunwrap' ;
end

Options.dummy = [];
Options       = assignifempty( Options, DEFAULTS ) ;

if ~any( strcmp( Options.unwrapper, {'Sunwrap','sunwrap','FslPrelude','AbdulRahman_2007','QGU'} ) )
    error('Unrecognized "Options.unwrapper" input. See HELP MaRdI.unwrapphase') ;
elseif strcmp( Options.unwrapper, 'AbdulRahman_2007' ) && ( size( Phase.img, 3 ) == 1 )
    warning('Options.unwrapper = AbdulRahman_2007 is incompatible with 2d images. Using Sunwrap method.') ;
    Options.unwrapper = 'Sunwrap' ;
end

%% ------
nVolumes = Phase.getnumberofmeasurements() ;
nEchoes  = size( Phase.img(), 4 ) ;

if ~exist( 'Mag' )
    assert( strcmp( Options.unwrapper, 'AbdulRahman_2007' ), ...
        ['Missed required input: corresponding Magnitude MaRdI-object. ' ...
        ' (Applies to Sunwrap and FslPrelude cases.) See HELP MaRdI.unwrapphase'] ) ;

    assert( myisfieldfilled( Phase.Hdr, 'MaskingImage' ), ...
        ['Either a logical masking array must be defined in Phase.Hdr.MaskingImage ' ...
         'OR the input must include the corresponding Magnitude image. See HELP MaRdI.unwrapphase'] ) ;

else
    assert( Phase.iscoincident( Mag ), ['Inputs Mag.img and Phase.img must correspond '...
        '(respective voxel positions and number of measurements should be identical'] )

    if ~myisfieldfilled( Phase.Hdr, 'MaskingImage' )
        mask = Mag.getreliabilitymask( Options.threshold ) ; 
        % indexing for potential bug, e.g. in Siemens dual echo field mapping, where nEchoes magnitude input > nEchoes phase input...
        % TODO: find a more elegant way of handling this...
        Phase.setmaskingimage( mask(:,:,:,1:nEchoes,1:nVolumes) ) ; 
    end
end

if ( nEchoes > 1 ) && Options.isLinearizing
    % phase evolution between 1st 2 echoes, to be used to correct 2*pi offsets
    % of individual echoes after unwrapping
    PhaseDiff      = Phase.copy() ;

    img            = Mag.img(:,:,:,1,:) .* exp( i*Phase.img(:,:,:,1,:) ) ;
    img(:,:,:,2,:) = Mag.img(:,:,:,2,:) .* exp( i*Phase.img(:,:,:,2,:) ) ;

    PhaseDiff.img  = angle( img(:,:,:,2,:) .* conj(img(:,:,:,1,:) ) ) ;

    PhaseDiff.Hdr.EchoTime = Phase.getechotime(2) - Phase.getechotime(1) ; % [units : ms]

    mask = Phase.Hdr.MaskingImage(:,:,:,1,:) & Phase.Hdr.MaskingImage(:,:,:,2,:) ;
    PhaseDiff.setmaskingimage( logical( mask ) & ~isnan( PhaseDiff.img ) ) ;

    PhaseDiff.unwrapphase( Mag, Options ) ;

    % if PhaseDiff is not centered around 0, shift it:
    mask      = PhaseDiff.Hdr.MaskingImage(:,:,:,1,:) ;
    phaseDiff = PhaseDiff.img(:,:,:,1,:) ; 
    n         = round( median( phaseDiff(mask) )/pi ) ;
    PhaseDiff.img(:,:,:,1,:) = PhaseDiff.img(:,:,:,1,:) - n*pi ;
end    


%% ------

for iVolume = 1 : nVolumes

    display(['Unwrapping volume ' num2str(iVolume) ' of ' num2str(nVolumes) ] ) ;    

    for iEcho = 1 : nEchoes
        
        if nEchoes > 1
            display(['Unwrapping echo ' num2str(iEcho) ' of ' num2str(nEchoes) ] ) ;    
        end

        switch Options.unwrapper

            case 'AbdulRahman_2007'
                
                Phase.img(:,:,:,iEcho, iVolume) = unwrap3d( Phase.img(:,:,:,iEcho, iVolume), logical(Phase.Hdr.MaskingImage(:,:,:,iEcho, iVolume)), Options ) ;

            case 'FslPrelude'

                if myisfield( Phase.Hdr, 'MaskingImage') && ~isempty(Phase.Hdr.MaskingImage)
                    Options.mask = single( Phase.Hdr.MaskingImage(:,:,:,iEcho,iVolume) ) ;
                end
                
                Options.voxelSize = Phase.getvoxelspacing() ;   

                Phase.img(:,:,:,iEcho, iVolume) = prelude( Phase.img(:,:,:,iEcho, iVolume), Mag.img(:,:,:,iEcho, iVolume), Options ) ;

            case {'Sunwrap', 'sunwrap'}
                
                iMag      = Mag.img(:,:,:,iEcho,iVolume) ;
                iMag      = iMag./max(iMag(:)) ;
                Phase.img(:,:,:,iEcho, iVolume) = sunwrap( iMag .* exp( 1i* Phase.img(:,:,:,iEcho,iVolume) ), Options.threshold ) ;
            
            case { 'QGU' }
                
                Options.mask = squeeze(single( Phase.Hdr.MaskingImage(:,:,:,iEcho,iVolume) ) ) ;
                if ((iVolume == 1) && (iEcho == 1))
                    [Phase.img(:,:,:,iEcho, iVolume),xpoint, ypoint] = QualityGuidedUnwrap2D_EAO(Phase.img(:,:,:,iEcho, iVolume), Options.mask);
                else
                    [Phase.img(:,:,:,iEcho, iVolume),xpoint, ypoint] = QualityGuidedUnwrap2D_EAO(Phase.img(:,:,:,iEcho, iVolume), Options.mask, xpoint, ypoint);
                end
                
        end

        if Options.isLinearizing && ( iEcho == 1 ) 
            % if phase(TE1) is not centered around 0, shift it:
            mask  = Phase.Hdr.MaskingImage(:,:,:,1,iVolume) ;
            phase = Phase.img(:,:,:,iEcho,iVolume) ; 
            n     = fix( median( phase(mask) )/pi ) ;
            if n ~= 0
                display( ['Recentering phase of 1st echo to sit between [-pi,pi] (global subtraction of ' num2str(n) 'pi)'] ) 
                Phase.img(:,:,:,iEcho,iVolume) = Phase.img(:,:,:,iEcho,iVolume) - n*pi*double(mask) ;
            end
        end

    end

    if ( nEchoes > 1 ) && Options.isLinearizing
    %% -----
    % Correct inter-echo wraps  
        for iEcho = 1 : nEchoes 
            % (assuming phase(TE=0) offset = 0):
            phaseEstimate = (Phase.getechotime(iEcho)/PhaseDiff.Hdr.EchoTime)*PhaseDiff.img ;

            err = phaseEstimate - Phase.img(:,:,:,iEcho,iVolume) ;
            err = median( err( Phase.Hdr.MaskingImage(:,:,:,iEcho,iVolume ) ) ) ; 

            % When the median deviation between estimate and unwrapped
            % measurement exceeds pi, correct the measurement by adding the
            % appropriate pi-multiple:
            n  = ( abs(err) > pi ) .* round( err/pi )  ;

            Phase.img(:,:,:,iEcho,iVolume) = Phase.img(:,:,:,iEcho,iVolume) + n*pi*Phase.Hdr.MaskingImage(:,:,:,iEcho,iVolume ) ;
        end
    end

end

Phase.img = double( Phase.img ) ;

%% -----
% if time series, correct for potential wraps between time points:
if nVolumes > 1
    display('Correcting for potential temporal phase wraps...')
    
    for iEcho = 1 : nEchoes
        % correct temporal wraps by comparing to phaseEstimate:
        phaseEstimate = median( Phase.img(:,:,:,iEcho,:),  5, 'omitnan' ) ;

        % normalize the estimate so its spatial median is within [-pi,pi]
        tmpMask = logical( prod( Phase.Hdr.MaskingImage(:,:,:,iEcho,:), 5 ) ) ;
        n       = fix( median( phaseEstimate(tmpMask) )/pi ) ;
        phaseEstimate = phaseEstimate - n*pi ;
        
        % Wherever the absolute deviation from the estimate exceeds pi,
        % correct the measurement by adding the appropriate pi-multiple:
        for iVolume = 1 : nVolumes
            dPhase = phaseEstimate - Phase.img( :,:,:,iEcho,iVolume )  ;
            n      = ( abs(dPhase) > pi ) .* fix( dPhase/pi ) ;
            Phase.img( :,:,:,iEcho,iVolume ) = Phase.img(:,:,:,iEcho,iVolume ) + n*pi ;
        end
    end
end

%% -----
% update header
Img.Hdr.ImageType         = 'DERIVED\SECONDARY\P\' ; 
Img.Hdr.SeriesDescription = ['phase_unwrapped_' Options.unwrapper ] ; 

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
methods(Access=protected)
% =========================================================================
function Interpolant = getinterpolant( Img, method, extrapolationMethod )
%GETINTERPOLANT     Return object of type scatteredInterpolant
% 
% GETINTERPOLANT returns an instance of Matlab's scatteredInterpolant class,
% useful for interpolating between different image grids (voxel positions).
% 
% Interpolant = GETINTERPOLANT( Img ) 
% Interpolant = GETINTERPOLANT( Img, method ) 
% Interpolant = GETINTERPOLANT( Img, method, extrapolationMethod ) 
%
% Default Interpolant property assignments:
%
%   .Method = 'linear' [i.e. the *interpolation* method]
%
%   .ExtrapolationMethod = 'none' 
%
%   .Values = [vectorized voxel values of 1st echo/measurement, i.e.: Img.img(:,:,:,1,1)]
%
% Note that all 3 properties can be reassigned at any point upon return.
%
% For info on the 2 optional arguments, see help scatteredInterpolant

DEFAULT_METHOD              = 'linear' ;
DEFAULT_EXTRAPOLATIONMETHOD = 'none' ;

if ( nargin < 2 ) || isempty(method)
    method = DEFAULT_METHOD ;
end

if ( nargin < 3 ) || isempty( extrapolationMethod )
    extrapolationMethod = DEFAULT_EXTRAPOLATIONMETHOD ;
end

disp( 'Forming interpolant...(Computation time is proportional to image size. This may take a few minutes.)' )

Interpolant = scatteredInterpolant() ;

Interpolant.Method              = method ;
Interpolant.ExtrapolationMethod = extrapolationMethod ;

[X,Y,Z] = Img.getvoxelpositions() ;

% The following avoids the error from scatteredInterpolant when one
% attempts to form a 3d interpolant from a 2d input: 
isValidDim0 = [ numel(unique(X(:))) numel(unique(Y(:))) numel(unique(Z(:))) ] > 1 ;
r0          = [X(:) Y(:) Z(:)] ;

Interpolant.Points = r0(:, isValidDim0) ;

v = Img.img(:,:,:,1,1) ;
Interpolant.Values = v(:) ;

end
% =========================================================================
function isSame = iscoincident( Img1, Img2 )
%ISCOINCIDENT   Check coincidence of 2 images 
% 
% isSame = ISCOINCIDENT( Img1, Img2 )
%
% Returns TRUE if Img1 and Img2 possess coincident voxel positions
% and number of measurements/volumes
%
% TODO: Check additional properties + add outputs for each?

assert( ( nargin == 2 ) && isa( Img2, 'MaRdI' ), 'Missed required input: 2 MaRdI-objects' )

isSame = false ;

if MaRdI.compareimggrids( Img1, Img2 )  
     
    isSame = true ;

    if ~strcmp( Img1.Hdr.SeriesDescription, Img2.Hdr.SeriesDescription )
        warning('A computation is being performed based on images acquired in seperate series.')
    end
end

end    
% =========================================================================
function [] = scalephasetofrequency( Img, undoFlag )
%SCALEPHASETOFREQUENCY
%
% Converts unwrapped phase [units:rad] to field [units: Hz]
% 
%   Field = scalephasetofrequency( UnwrappedPhase )
%   Phase = scalephasetofrequency( Field, -1 )
%   
%   The 'undo' mode with -1 as the 2nd argument scales from Hz back to rad

assert( ~Img.ismagnitude(), 'Input cannot be a magnitude image' ) ;

te = Img.getechotime()*(1E-3) ;

if (nargin < 2) || (undoFlag ~= -1)
    assert( Img.isphase() && strcmp( Img.Hdr.PixelComponentPhysicalUnits, '0000H' ), ...
        'Expected Phase Img input with voxel values in radians' ) ;
    
    scalingFactors = 1./(2*pi*te) ;
    Img.Hdr.PixelComponentPhysicalUnits = '0005H' ; % i.e. Hz

elseif (undoFlag == -1)
    assert( strcmp( Img.Hdr.PixelComponentPhysicalUnits, '0005H' ), ...
        'Expected Field Img input with voxel values in Hz' )
    
    scalingFactors = 2*pi*te ;
    Img.Hdr.PixelComponentPhysicalUnits = '0000H' ; % i.e. none
end

for iEcho = 1 : numel(te) 
    Img.img(:,:,:,iEcho,:) = scalingFactors(iEcho) * Img.img(:,:,:,iEcho,:) ;
end

end
% =========================================================================

end
% =========================================================================
% =========================================================================
methods(Access=private)
% =========================================================================
function [r,c,s] = getdirectioncosines( Img ) 
% GETDIRECTIONALCOSINES
% 
%   "The direction cosines of a vector are the cosines of the angles between
%   the vector and the three coordinate axes. Equivalently, they are the
%   contributions of each component of the basis to a unit vector in that
%   direction." 
%   https://en.wikipedia.org/wiki/Direction_cosine
%
% [r,c,s] = GETDIRECTIONALCOSINES( Img ) 
%   r: row *index* direction cosine
%   c: column index " "
%   s: slice index " "
%
%   NB: the *index* term. r & c may be defined obstrusely:
%   i.e. r is not the row direction cosine (c is!), it is the direction cosine
%   of the vector that points along the direction of increasing row indices
%   (i.e. it's in the column direction!)
      
% " To make a full rotation matrix, we can generate the last column from
% the cross product of the first two. However, Siemens defines, in its
% private CSA header, a SliceNormalVector which gives the third column, but
% possibly with a z flip, so that R is orthogonal, but not a rotation
% matrix (it has a determinant of < 0)."
%  -From http://nipy.org/nibabel/dicom/dicom_mosaic.html 
%
% NOTE: Hdr.Img.SliceNormalVector gets defined in the MaRdI constructor

c = Img.Hdr.ImageOrientationPatient(1:3) ; 
r = Img.Hdr.ImageOrientationPatient(4:6) ; 
s = Img.Hdr.Img.SliceNormalVector ;

end 
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
% NOTE
%   Temporarily placing cropimg(), nii(), and zeropad() here since these methods
%   1) might be deprecated
%   2) may or may not be useful to a typical user
% =========================================================================
function Img = cropimg( Img, gridSizeImgCropped, centralPoint )
%CROPIMG
%
% Img = CROPIMG( Img, croppedDims )
% Img = CROPIMG( Img, croppedDims, centralPoint )
% 
% centralPoint are the indices of the original img voxel to become the centralPoint of 
% the cropped array.
%
% default  (nargin == 2)
%   centralPoint = round( size(Img.img)/2 );
% 
% note: 
% if centralPoint +/- croppedDims/2 exceeds the original bounds, the array is cropped at the bound (as opposed to zero filling)
% -------
% *** TODO
%
%   make compatible for odd-sized arrays
    
gridSizeImgOriginal = size(Img.img) ;

if (nargin == 2) || isempty(centralPoint)
    centralPoint = round( gridSizeImgOriginal/2 ) ;
end

low  = centralPoint - round(gridSizeImgCropped/2) + [1 1 1] ;
high = low + gridSizeImgCropped - [1 1 1] ;  

for dim = 1: 3
   if low(dim) < 1
      low(dim) = 1 ;
   end
   if high(dim) > gridSizeImgOriginal(dim)
      high(dim) = gridSizeImgOriginal(dim) ;
  end
end
    
% Update header
[X, Y, Z] = Img.getvoxelpositions( ); 
x0        = X(low(1),low(2),low(3)) ;
y0        = Y(low(1),low(2),low(3)) ;
z0        = Z(low(1),low(2),low(3)) ;

Img.Hdr.ImagePositionPatient = [x0 y0 z0] ;

[rHat, cHat, sHat] = Img.getdirectioncosines( ) ;  

Img.Hdr.SliceLocation = dot( Img.Hdr.ImagePositionPatient, sHat ) ;

% crop img
Img.img = Img.img( low(1):high(1), low(2):high(2), low(3):high(3) ) ; 

if myisfield( Img.Hdr, 'MaskingImage' )  && ~isempty( Img.Hdr.MaskingImage ) 
    
   Img.Hdr.MaskingImage = Img.Hdr.MaskingImage( low(1):high(1), low(2):high(2), low(3):high(3) ) ; 

end

% Update header
Img.Hdr.Rows                 = size(Img.img, 1) ;
Img.Hdr.Columns              = size(Img.img, 2) ;       

end
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
function Img = zeropad( Img, padSize, padDirection )
%ZEROPAD
% Img = ZEROPAD( Img, padSize, padDirection )
%
% padSize = [ nZerosRows nZerosColumns nZerosSlices ] 
%
% padDirection == 'post' || 'pre' || 'both'
%
% -------
%   See also PADARRAY()
    
gridSizeOriginalImg = size(Img.img) ;

Img.img = padarray( Img.img, padSize, 0, padDirection ) ; 

if myisfield( Img.Hdr, 'MaskingImage' )  && ~isempty( Img.Hdr.MaskingImage ) 
    
   Img.Hdr.MaskingImage = padarray( Img.Hdr.MaskingImage, padSize, 0, padDirection ) ; 

end

% Update header
Img.Hdr.Rows                 = size(Img.img, 1) ;
Img.Hdr.Columns              = size(Img.img, 2) ;       

if ~strcmp( padDirection, 'post' )
% update image position 
% (i.e. location in DICOM RCS of 1st voxel in data array (.img))
    
    voxelSize = Img.getvoxelspacing() ;

    dr = voxelSize(1) ; % spacing btwn rows 
    dc = voxelSize(2) ; % spacing btwn columns
    ds = voxelSize(3) ;
    
    nr = padSize(1) ;
    nc = padSize(2) ;
    ns = padSize(3) ;
    
    [r, c, s] = Img.getdirectioncosines( ) ;

    % -1 because the zeros are padded before ('pre') 1st element (origin)        
    dx = -1*( r(1)*dr*nr + c(1)*dc*nc + s(1)*ds*ns ) ; 
    dy = -1*( r(2)*dr*nr + c(2)*dc*nc + s(2)*ds*ns ) ;
    dz = -1*( r(3)*dr*nr + c(3)*dc*nc + s(3)*ds*ns ) ;

    x1 = Img.Hdr.ImagePositionPatient(1) + dx ;
    y1 = Img.Hdr.ImagePositionPatient(2) + dy ;
    z1 = Img.Hdr.ImagePositionPatient(3) + dz ;

    Img.Hdr.ImagePositionPatient = [x1 y1 z1] ;

    Img.Hdr.SliceLocation = dot( Img.Hdr.ImagePositionPatient, s ) ;
end

end
% =========================================================================

end
% =========================================================================
% =========================================================================
methods(Static)
% =========================================================================
function isSame = compareimggrids( varargin )
%COMPAREIMGGRIDS 
%
% isSame = COMPAREIMGGRIDS( Img1, Img2 )
% isSame = COMPAREIMGGRIDS( X1, Y1, Z1, X2, Y2, Z2 )
%
% Returns TRUE if voxel positions of Img1 and Img2 are identical

%% -------
% check & parse inputs
assert( nargin >= 2, 'Function requires at least 2 input arguments. See HELP MaRdI.compareimggrids' ) 

if isa( varargin{1}, 'MaRdI' ) && isa( varargin{2}, 'MaRdI' )
    [X1, Y1, Z1] = varargin{1}.getvoxelpositions ;
    [X2, Y2, Z2] = varargin{2}.getvoxelpositions ;
elseif ( nargin ==6 ) && all( cellfun( @isnumeric, varargin ) )
    X1 = varargin{1} ;
    Y1 = varargin{2} ;
    Z1 = varargin{3} ;
    X2 = varargin{4} ;
    Y2 = varargin{5} ;
    Z2 = varargin{6} ;
else
    error('See HELP MaRdI.compareimggrids ')
end

%% -------
% compare grid positions
if ( numel(size(X1)) ~= numel(size(X2)) ) || any( size(X1) ~= size(X2) ) || any( X1(:) ~= X2(:) ) || any( Y1(:) ~= Y2(:) ) || any( Z1(:) ~= Z2(:) )
    isSame = false ;
elseif ( all(X1(:) == X2(:) ) && all( Y1(:) == Y2(:) ) && all( Z1(:) == Z2(:) ) ) 
    isSame = true ;
else
    error('Unexpected result: Check conditions apparently insufficient. Modify code')
end

end
% =========================================================================
function list = findimages( imgDir )
%FINDIMAGES
% 
% list = FINDIMAGES( imageDirectory ) 
%
% Calls dir() to return list of .dcm OR .IMA files in imageDirectory and its echo_* subdirectories

imgDir       = [ imgDir '/' ] ;

ListSubdirs  = dir( [ imgDir 'echo*'] );
nEchoSubdirs = length( ListSubdirs ) ;

if nEchoSubdirs > 0 
    
    List = finddcmorima( [imgDir ListSubdirs(1).name] ) ;
    nImgPerEcho = length(List) ;

    list = cell( nImgPerEcho*nEchoSubdirs, 1 ) ;

    for iImg = 1 : nImgPerEcho 
        list{iImg} = [ imgDir ListSubdirs(1).name '/' List(iImg).name ] ;
    end

    for iEcho = 2 : nEchoSubdirs
        ListiSubdir = finddcmorima( [imgDir ListSubdirs(iEcho).name] ) ;
        assert( length(ListiSubdir) == nImgPerEcho, 'Each echo subdirectory should contain the same number of images' ) ; 

        for iImg = 1 : nImgPerEcho
            list{ (iEcho-1)*nImgPerEcho + iImg} = [ imgDir ListSubdirs(iEcho).name '/' ListiSubdir(iImg).name ] ;
        end
    end

else
    List = finddcmorima( imgDir ) ;
    nImg = length(List) ;
    list = cell( nImg, 1 ) ;

    for iImg = 1 : nImg
        list{iImg} = [ imgDir List(iImg).name ] ;
    end
end


function List = finddcmorima( imgDir ) 
%Find .dcm OR .IMA files in imgDir
List   = dir( [ imgDir '/*.dcm'] );

if length(List) == 0 % try .IMA
    List = dir( [imgDir '/*.IMA'] ) ;
end

nImages = length( List ) ; 
assert( nImages ~= 0, 'No .dcm or .IMA files found in given directory' ) ;

end

end
% =========================================================================
function fullDir = getfulldir( dataLoadDir, iDir )
%GETFULLDIR
% 
% fullDir = GETFULLDIR( parentDir, index ) 
%
%   Searches parentDir[ectory] for subdirectory beginning with index- (e.g.
%   .../dataLoadDir/[index]-* ) to return its full path.
%
%   (Useful for rapidly initializing MaRdI with a dicom folder.)

if nargin < 2
    error('Function requires 2 input arguments: parent directory and index.')
end

if isnumeric(iDir)
    iDir = num2str( iDir ) ;
end

if length( iDir ) == 1
    iDir = ['0' iDir] ;
end

if ~strcmp( dataLoadDir(end), '/' )
    dataLoadDir(end+1) = '/' ;
end

Tmp = dir( [ dataLoadDir iDir '-*'] ) ;
fldrName = Tmp.name ;

[fullDir,~,~] = fileparts( [dataLoadDir fldrName '/'] ) ;

fullDir(end+1) = '/' ;

end
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
DEFAULT_CYLINDERSIZE       = 80 ;
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
if strcmp(Params.centerlineMethod, 'midfov')
    % create a vertical line centered in the axial FOV
    system( ['sct_create_mask -i ' Params.tmpSaveDir 't2s.nii.gz -p center -size 1 -f box -o ' Params.tmpSaveDir 't2s_centerline.nii.gz' ] ) ;
elseif strcmp(Params.centerlineMethod, 'spinalcord')
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
system( ['gunzip ' Params.tmpSaveDir 't2s_allEchoes.nii.gz -df'] ) ;

% delete the other images
system( ['rm ' Params.tmpSaveDir '*.nii.gz'] ) ;

% move segmentation, weights + t2* images
system( ['mv ' Params.tmpSaveDir 't2s_seg.nii ' Params.dataSaveDir 'gre_seg.nii'] ) ;  
system( ['mv ' Params.tmpSaveDir 't2s_weights.nii ' Params.dataSaveDir 'gre_weights.nii'] ) ; 
system( ['mv ' Params.tmpSaveDir 't2s_allEchoes.nii ' Params.dataSaveDir 'mgre.nii'] ) ; 
system( ['rm -r ' Params.tmpSaveDir] ) ;

% load FSLeyes to see segmentation overlaid on gre images
% this will help the user assess what field gradients are going to be
% averaged and compensated for 
system(['fsleyes ' Params.dataSaveDir 'mgre.nii ' Params.dataSaveDir 'gre_seg.nii -cm red -a 70.0 &' ]);

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
function [] = sortimages( unsortedDicomDir, sortedDicomDir, isCopying )
%SORTIMAGES 
%
% Arrange unsorted DICOMs (e.g. as delivered via Siemens socket) into organized 
% subdirectories.
% 
% [] = SORTDATA( unsortedDicomDir ) 
% [] = SORTDATA( unsortedDicomDir, sortedDicomDir ) 
% [] = SORTDATA( unsortedDicomDir, sortedDicomDir, isCopying ) 
% 
% If sortedDicomDir is unspecified, a subdirectory ('sorted') is created
% within unsortedDicomDir to contain the sorted images.
%
% isCopying (Boolean) is TRUE by default. Set to 0 to move the files instead
% of copying. Note that the files will still be renamed.

DEFAULT_SORTEDDICOMDIR = [ unsortedDicomDir '/sorted/' ] ;
DEFAULT_ISCOPYING      = true ;

%Dicom files inside the directory -----------------------------------------
listOfImages = dir( [ unsortedDicomDir '/*.dcm'] );

if length(listOfImages) == 0
    % try .IMA
    listOfImages = dir( [unsortedDicomDir '/*.IMA'] ) ;
end

nImages = length( listOfImages ) ;

assert( nImages~=0, 'No .dcm or .IMA files found in given directory' ) ;

if nargin < 2 || isempty( sortedDicomDir )
    sortedDicomDir = DEFAULT_SORTEDDICOMDIR ;
end

if ~exist( sortedDicomDir, 'dir') 
    mkdir( sortedDicomDir ); 
end

if nargin < 3 || isempty( isCopying )
    isCopying = DEFAULT_ISCOPYING ;
end

%Read Dicom files header and extract series names and numbers -------------
for iImage = 1 : nImages

    display( ['Sorting... ' num2str(100*iImage/nImages) ' %)' ])
    
    Hdr       = dicominfo( [unsortedDicomDir '/' listOfImages(iImage).name] );
    Hdr.Img   = parse_siemens_shadow( Hdr ) ;
    seriesDir = [ num2str(Hdr.SeriesNumber) '_' Hdr.SeriesDescription '/' ];

    if ( Hdr.SeriesNumber < 10)
        seriesDir = [ '0' seriesDir ] ;
    end
    
    % Create directories for each series
    if ~exist( [ sortedDicomDir seriesDir ], 'dir' )
        mkdir( [ sortedDicomDir seriesDir ] );
    end
    
    echoDir = [ sortedDicomDir seriesDir 'echo_' num2str( Hdr.EchoTime, 3 ) ] ;
         
    if ~exist( echoDir, 'dir' )
        mkdir( echoDir ) ;
    end
    
    iSlice   = Hdr.Img.ProtocolSliceNumber + 1 ;
    sliceStr = num2str( iSlice ) ;

    acqStr   = num2str( Hdr.AcquisitionNumber ) ;

    for ord = 3 : -1 : 1
        if iSlice < 10^ord
            sliceStr = [ '0' sliceStr ] ;
             acqStr  = [ '0' acqStr ] ;
        end
    end

    [~, ~, ext]    = fileparts( Hdr.Filename ) ;
    sortedFilename = fullfile( echoDir, strcat( Hdr.PatientName.FamilyName, '-', ...
                Hdr.Img.ImaCoilString, '-', sliceStr, '-', acqStr, ext ) ) ;
    
    if isCopying
        copyfile( Hdr.Filename, sortedFilename{1} ) ;
    else
        movefile( Hdr.Filename, sortedFilename{1} ) ;
    end

end
    
end
% =========================================================================
function [ studyDirs ] = tablestudy( sortedDicomDir )
%TABLESTUDY 
%
% Returns a cell array ( studyDirs ) pertaining to the study directory input
% ( sortedDicomDir ) where each element in the second column is a MaRdI-loadable 
% images series. (The 1st column is merely the row index.)
%
% e.g. Protocol to load MaRdI-object :
%
%   % omit semi-colon to display the full cell array (i.e. with the row indices)
%   [ studyDirs ] = MaRdI.tablestudy( sortedDicomDir ) 
%
%   % determine the row index of the series you want to load (e.g. 10):
%   Img = MaRdI( studyDirs{ 10, 2 } ) ;

assert( nargin == 1, 'Function requires sortedDicomDirectory as input argument.' ) ;

if ~strcmp( sortedDicomDir(end), '/' ) 
    sortedDicomDir(end+1) = '/' ;
end

studyDirs = cell( 0 ) ;

Tmp      = dir( [ sortedDicomDir ] );
Tmp      = Tmp( 3:end ) ; % ignore self ('.') and parent ('..') dirs
nEntries = length( Tmp ) ;

for iEntry = 1 : nEntries 

   if Tmp( iEntry ).isdir
   
       tmpSeriesSubdir = [ Tmp( iEntry ).name '/'] ;
    
        TmpEchoSubdirs = dir( [ sortedDicomDir tmpSeriesSubdir 'echo*' ] ) ;
        nEchoSubdirs   = length( TmpEchoSubdirs ) ;

        if nEchoSubdirs ~= 0

            for iEchoSubdir = 1 : nEchoSubdirs

                studyDirs{end+1, 2} = strcat( sortedDicomDir, tmpSeriesSubdir, TmpEchoSubdirs(iEchoSubdir).name )  ;
                iSeries = size( studyDirs, 1 ) ;
                studyDirs{ iSeries, 1 } = iSeries ;
            end

        % check if tmpSeriesSubdir itself contains images
        elseif length( dir( [ sortedDicomDir tmpSeriesSubdir '/*.dcm'] ) ) ~= 0 || ...
                length( dir( [ sortedDicomDir tmpSeriesSubdir '/*.IMA'] ) ) ~= 0 

           studyDirs{end+1, 2} = strcat( sortedDicomDir,tmpSeriesSubdir )  ;
            iSeries = size( studyDirs, 1 ) ;
            studyDirs{ iSeries, 1 } = iSeries ;

        end

   end
   
end

end
% =========================================================================
function [Params] = savefigure( img, Params )
%SAVEFIGURE
%
%   Description
%   
%   Write .png image to file using the 'export_fig' tool
%
%    .......................
%
%   Syntax
%
%   Parameters = SAVEFIGURE( img, Parameters )
%   
%   Returns employed Parameters struct.
%
%    .......................
%   
%   The following Parameter.fields are supported
%
%   .filename
%       default = './tmp'
%
%   .colormap
%       default = 'gray'
%
%   .scaling
%       default = [min(img) max(img)]
%
%   .magnification
%       default = 1
%
%   .isColorbar
%       default = false
%
%   .isBackgroundTransparent
%       default = false

DEFAULT_FILENAME      = './tmp.bin' ;
DEFAULT_COLORMAP      = 'gray' ;
DEFAULT_MAGNIFICATION = 1 ;
DEFAULT_ISCOLORBAR    = false ;
DEFAULT_ISBACKGROUNDTRANSPARENT = false ;
% =========================================================================
% Check inputs
% =========================================================================
if nargin < 1 || isempty(img)
    disp('Error: function requires at least 1 argument (2D image-matrix)')
    help(mfilename);
    return;  
end

if nargin < 2 || isempty(Params)
    disp('Default parameters will be used')
    Params.dummy = [] ;
end

if  ~myisfield( Params, 'filename' ) || isempty(Params.filename)
    Params.filename = DEFAULT_FILENAME ;
end

if  ~myisfield( Params, 'colormap' ) || isempty(Params.colormap)
    Params.colormap = DEFAULT_COLORMAP ;
end

if  ~myisfield( Params, 'magnification' ) || isempty(Params.magnification)
    Params.magnification = DEFAULT_MAGNIFICATION ;
end

if  ~myisfield( Params, 'isColorbar' ) || isempty(Params.isColorbar)
    Params.isColorbar = DEFAULT_ISCOLORBAR ;
end

if  ~myisfield( Params, 'isBackgroundTransparent' ) || isempty(Params.isBackgroundTransparent)
    Params.isBackgroundTransparent = DEFAULT_ISBACKGROUNDTRANSPARENT ;
end

if  ~myisfield( Params, 'scaling' ) || isempty(Params.scaling)
    Params.scaling = [min(img(:)) max(img(:))] ;
    if ( Params.scaling(2) == Params.scaling(1) )
        Params.scaling(2) = Inf ;
    end
end

[~,~,extension] = fileparts( Params.filename ) ;

if ~strcmp( extension, '.png' )
    Params.filename = [Params.filename '.png' ] ;
end

% =========================================================================
% Create figure
% =========================================================================

figure('units','normalized','outerposition',[0 0 1 1])

imagesc( img, Params.scaling ) ; 
colormap(Params.colormap); 
axis image ;

if Params.isColorbar
    colorbar;
    hcb=colorbar;
    set(hcb,'YTick',[])
end

set(gca,'XTick',[]) % Remove the ticks in the x axis
set(gca,'YTick',[]) % Remove the ticks in the y axis
set(gca,'Position',[0 0 1 1]) % Make the axes occupy the hole figure

if Params.isBackgroundTransparent
    export_fig(Params.filename, '-png', '-transparent', ['-m' num2str(Params.magnification)]) 
else
    export_fig(Params.filename, '-png', ['-m' num2str(Params.magnification)]) 
end

% close gcf

% crop out border voxels ?
img = imread( Params.filename, 'png' ) ;
imwrite( img( 5:end-3, 5:end-3, : ), Params.filename, 'png' ) ;

end
% =========================================================================

   
end
% =========================================================================
% =========================================================================


end
