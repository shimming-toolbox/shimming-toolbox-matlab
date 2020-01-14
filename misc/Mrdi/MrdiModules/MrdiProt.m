classdef (Sealed=true) MrdiProt 
%MrdiProt  MR DICOM Image Protocol: Image acquisition protocol properties

% =========================================================================
% Author::ryan.topfer@polymtl.ca
% =========================================================================

% =========================================================================
% =========================================================================
properties (Transient=true)
    Hdrs ;
end

properties (Dependent, SetAccess=private)

    % Image acquisition times [ms elapsed since midnight]. Array dims [ nSlices x nEchoes x nMeasurements ]
    % 
    % Derives from the AcquisitionTime field in Siemens DICOM header.
    %
    % acquisitionTime - acquisitionTime(1) yields the elapsed time since
    % acquisition of the 1st k-space point in the series.
    %
    % For EPI MOSAIC (which has a single AcquisitionTime value in the DICOM
    % header for each volume), 
    % acquisitionTime( iSlice ) = AcquisitionTime (first slice) + iSlice*(volumeTR/nSlices)
    acquisitionTime {mustBeNonnegative} = 0 ; 

    % Gyromagnetic ratio of imaged nucleus [units: rad/s/T] (Note: Only implemented for 1H protons!)
    gyro(1,1) {mustBeReal} = 0 ; 
    
    % Larmor Tx. frequency [units: Hz]
    imagingFrequency(1,1) {mustBeNonnegative} = 0 ; 
    
    % [x,y,z] coordinates of magnet isocenter in the patient coordinate system
    isocenter(1,3) {mustBeReal} = [0 0 0 ] ; 
    
    % Number of image measurements
    nMeasurements(1,1) {mustBeInteger} = 1 ; 
    
    % Number of acquired slices (=1 for 3d slab encoding)
    nSlices(1,1) {mustBeInteger} = 1 ;  

    % Fraction of k-space coverage in each dim [read phase partition]
    partialFourierFactors {mustBePositive} = [ 1 1 1 ] ; 
    
    % Slice [units: mm] : DICOM Hdr tag ( )
    sliceThickness(1,1) {mustBePositive} = 1 ; 
    
    % spacingBetweenSlices [units: mm] : DICOM Hdr tag ( )
    spacingBetweenSlices(1,1) {mustBePositive} = 1 ; 
    
    % Vector of echo times [units: ms]
    te {mustBeNonnegative} = 0 ; 
    
    % Repetition time [units: ms]
    tr(1,1) {mustBeNonnegative} = 0 ;

end
% =========================================================================
% =========================================================================

% =========================================================================
% =========================================================================    
methods
% =========================================================================    
function Prot = MrdiProt( Hdrs )
    
    if nargin == 1
        Prot.Hdrs = Hdrs ;
    end

end
% =========================================================================
function [] = associateaux( Img, Aux, Params )
%ASSOCIATEAUX - link image to corresponding auxiliary recording object
% 
% Method for MrdiProt or Mrdi??
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
function [t] = get.acquisitionTime( Prot ) 
% GETACQUISITIONTIME  Return array of image acquisition times [ms elapsed since midnight]
%
% Derives from the AcquisitionTime field in Siemens DICOM header.
%
% t = GETACQUISITIONTIME( Prot )
%
% dimensions of t: [ nSlices x nEchoes x nMeasurements ]
%
% t - t(1) yields the elapsed time since acquisition of the 1st k-space point
% in the series.
%
% For EPI MOSAIC (which has a single AcquisitionTime value for each volume),
% t( iSlice ) = AcquisitionTime (first slice) + iSlice*(volumeTR/nSlices)

nSlices       = Prot.nSlices ;
nEchoes       = length( Prot.te ) ;
nMeasurements = Prot.nMeasurements ;

t = zeros( nSlices, nEchoes, nMeasurements ) ;

if isempty(strfind( Prot.Hdrs{1}.ImageType, 'MOSAIC' ) )
    
    for iMeasurement = 1 : nMeasurements
        for iEcho = 1 : nEchoes 
            for iSlice = 1 : nSlices
                t(iSlice, iEcho, iMeasurement) = convertacquisitiontime( Prot.Hdrs{iSlice,iEcho,iMeasurement}.AcquisitionTime ) ;
            end
        end
    end

else % special case for EPI MOSAIC (single DICOM header for multiple slices)
   
    sliceOrder = Prot.Hdr.MrProt.sSliceArray.alSliceAcqOrder ;
    sliceMeasurementDuration = Prot.tr/nSlices ; % [units: ms]

    for iMeasurement = 1 : nMeasurements
        for iEcho = 1 : nEchoes 
            for iSlice = 1 : nSlices
                t_iAcq = convertacquisitiontime( Prot.Hdrs{1,iEcho,iMeasurement}.AcquisitionTime ) ;
                t_iAcq = t_iAcq + sliceOrder(iSlice)*sliceMeasurementDuration ;
                t(iSlice, iEcho, iMeasurement) = t_iAcq ;
            end
        end
    end

end

function t_s = convertacquisitiontime( AcquisitionTime )
%CONVERTACQUISITIONTIME  Converts the time-stamp string in the .AcquisitionTime DICOM header
% For details, see https://cfn.upenn.edu/aguirre/wiki/public:pulse-oximetry_during_fmri_scanning

    % hours + minutes + seconds + fractional seconds :
    t_s = str2double( AcquisitionTime(1:2) )*60*60 + str2double( AcquisitionTime(3:4) )*60 + ...
          str2double( AcquisitionTime(5:6) ) + (str2double( AcquisitionTime(7:10) )/10)/1000 ;
    t_s = t_s * 1000 ; % [units: ms]

end

end
% =========================================================================
function GYRO = get.gyro( Prot )
%GETGYRO  Return gyromagnetic ratio of imaged nucleus [units: rad/s/T]
%
% NOTE: Only supports 1H protons!

    switch Prot.Hdr.ImagedNucleus 
        case '1H' 
            GYRO = 267.513E6 ; 
        otherwise
            error('Not implemented.') ;
    end

end
% =========================================================================
function xyzIso = get.isocenter( Prot )
%ISOCENTER  Return magnet isocenter [x,y,z] coordinates w.r.t. patient coordinate system
% 
% xyzIso = ISOCENTER( Prot ) 

    xyzIso = Prot.Hdr.Prot.ImaRelTablePosition()' ;

    assert( xyzIso(1) == 0, 'Table shifted in L/R direction?' ) ;
    assert( xyzIso(2) == 0, 'Table shifted in A/P direction?' ) ;

end
% =========================================================================
function [f0] = get.imagingFrequency( Prot ) 
%GETIMAGINGFREQUENCY  Return Larmor frequency [units: Hz]

    f0 = Prot.Hdr.MrProt.sTXSPEC.asNucleusInfo.lFrequency ;

end
% =========================================================================
function [nVolumes] = get.nMeasurements( Prot )
%GETNMEASUREMENTS    Return the number of image measurements
%
% n = GETNMEASUREMENTS( Prot )
%
% NOTE: GETNMEASUREMENTS( Prot ) is equivalent to n = size( Prot.img, 5 ),
% however GETNUMBEROFMEASUREMENTS also checks the DICOM header in Prot.Hdr and
% issues a warning if n differs from the expected value
% (Prot.Hdr.MrProt.lRepetitions +1).

    nVolumes = size( Prot.img, 5 ) ;

    if myisfield( Prot.Hdr.MrProt, 'lRepetitions' ) 
        % number of measurements according to Hdr:
        nVolumesHdr = Prot.Hdr.MrProt.lRepetitions + 1 ;
        if nVolumes ~= nVolumesHdr
            warning([ num2str(nVolumes) ...
                ' image volumes have been loaded, however the DICOM Hdr indicates ' num2str(nVolumesHdr) ' were acquired.'] ) ;
        end 
    end

end
% =========================================================================
function [nSlices] = get.nSlices( Prot ) 
%GETNSLICES  Return number of acquired slices
%
% NOTE: nSlices is not necessarily equal to size( Prot.img, 3).
% e.g. For a 3d (slab) encoding, GETNUMBEROFSLICES returns 1. 

    nSlices = Prot.Hdr.MrProt.sSliceArray.lSize ;

end
% =========================================================================
function [pff] = get.partialFourierFactors( Prot ) 
%GETPARTIALFOURIERFACTORS  Return fraction of k-space coverage in each dim
% 
% pff = GETPARTIALFOURIERFACTORS( Prot ) 
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

    dcmFields = { 'Prot.Hdr.MrProt.sKSpace.ucReadoutPartialFourier' ;
                  'Prot.Hdr.MrProt.sKSpace.ucPhasePartialFourier' ;
                  'Prot.Hdr.MrProt.sKSpace.ucSlicePartialFourier' ; } ;

    pfAsInt = [ Prot.Hdr.MrProt.sKSpace.ucReadoutPartialFourier
                Prot.Hdr.MrProt.sKSpace.ucPhasePartialFourier
                Prot.Hdr.MrProt.sKSpace.ucSlicePartialFourier ]' ;

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
                error( [ 'Unexpected value of ' num2str( pfAsInt(iDim) ) ' in ' dcmFields{iDim} ] ) ;
        end
    end

end
% =========================================================================
function [te] = get.te( Prot ) 
%GETTE    Return vector of echo times [units: ms]
% 
% TE = GETECHOTIME( Prot )

    if strcmp( Prot.Hdr.SequenceName, '*fm2d2' ) && Prot.isphase() 
        te = ( Prot.Hdr.MrProt.alTE(2) - Prot.Hdr.MrProt.alTE(1) )/1000  ;
    else
        nEchoes  = size( Prot.img, 4 ) ;
        te = Prot.Hdr.MrProt.alTE(1:nEchoes)/1000 ;
    end

end
% =========================================================================
function tr = get.tr( Prot ) 
%GETTR repetition time [units: ms]

    tr = Prot.Hdr.RepetitionTime ;

end
% =========================================================================

end
% =========================================================================
% =========================================================================    
methods(Sealed=true)
% =========================================================================
function [f0, g0, s0] = adjvalidateshim( Prot )
%ADJVALIDATESHIM    Return protocol shim settings from Siemens private header
% 
% [f0, g0, s0] = ADJVALIDATESHIM( Prot )
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
f0 = Prot.imagingFrequency ; 

% linear gradients
g0 = [ Prot.Hdr.MrProt.sGRADSPEC.asGPAData.lOffsetX ; ...
       Prot.Hdr.MrProt.sGRADSPEC.asGPAData.lOffsetY ; ...
       Prot.Hdr.MrProt.sGRADSPEC.asGPAData.lOffsetZ ] ;

% -------
% 2nd order shims (5-element vector)
s0 = Prot.Hdr.MrProt.sGRADSPEC.alShimCurrent ;

end
% =========================================================================
function [t0] = estimatekorigintime( Prot ) 
%ESTIMATEKORIGINTIME  Return estimate of time when k-space origin was sampled
% 
% t0 = ESTIMATEKORIGINTIME( Prot )
%
% Returns an estimate of when the k-space origin of an image was sampled
% relative to the AcquisitionTime (field in Siemens DICOM header) as a double
% in units of milliseconds. 
% 
% See also: MrdiProt.getacquisitiontime()
% 
% NOTE: This is a crude estimate and only the case of Cartesian k-space
% sampling, beginning at the k_min periphery, has been considered in the
% current implementation!
 
nSlices       = Prot.nSlices ;
nEchoes       = length( Prot.te ) ;
nMeasurements = Prot.nMeasurements ;

tAcq = Prot.getacquisitiontime() ;

% Estimate time from excitation to k-space origin:
dt = 0; % [units: ms]

if Prot.Hdr.MrProt.sKSpace.ucTrajectory == 1
    
    if ~strcmp( Prot.Hdr.ScanningSequence, 'EPI' )
    % NOTE: The Siemens DICOM field SliceMeasurementDuration appears to be a misnomer:
    % Rather, it (usually?) refers to the duration of an entire volume. 
        dt = Prot.Hdr.Prot.SliceMeasurementDuration ; % [units: ms]
    else
    % *Exception*: Apparently EPI are handled differently as said DICOM field
    % is oddly large (probably includes dummy volumes?)
        dt = Prot.tr/nSlices ; % [units: ms]
    end
    
    pf = Prot.partialFourierFactors ;
    
    if all( pf == 1 )
        dt = dt/2 ;
    else
        switch Prot.Hdr.MRAcquisitionType 
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
function [Hdrs] = set.Hdrs( Hdrs )
%SETHDRS Update protocol with new cell array of DICOM Hdrs. 
    % Prot.Hdrs = Hdrs ; 
end
% =========================================================================

end
% =========================================================================
% =========================================================================
methods(Sealed=true)
% =========================================================================
% function isSame = iscoincident( Img1, Img2 )
% %ISCOINCIDENT   Check coincidence of 2 images 
% % 
% % isSame = ISCOINCIDENT( Img1, Img2 )
% %
% % Returns TRUE if Img1 and Img2 possess coincident voxel positions
% % and number of measurements/volumes
% %
% % TODO: Check additional properties + add outputs for each?
%
% assert( ( nargin == 2 ) && isa( Img2, 'MrdiUtil' ), 'Missed required input: 2 MrdiUtil-objects' )
%
% isSame = false ;
%
% if Mrdi.compareimggrids( Img1, Img2 )  
%      
%     isSame = true ;
%
%     if ~strcmp( Img1.Hdr.SeriesDescription, Img2.Hdr.SeriesDescription )
%         warning('A computation is being performed based on images acquired in seperate series.')
%     end
% end
%
% end    
% =========================================================================

end
% =========================================================================
% =========================================================================


end
