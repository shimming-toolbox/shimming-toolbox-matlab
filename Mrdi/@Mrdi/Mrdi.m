classdef Mrdi < MrdiUtil & dynamicprops
%Mrdi  MR DICOM Image: Generic data type for image representation and data handling.
%
% Incl. member methods for image processing 
% 
% e.g.
%   filter()
%   ...etc.

% =========================================================================
% Author::ryan.topfer@polymtl.ca
% =========================================================================

properties( Constant )
    
    % NOTE: Unsure re: best attributes (or type! as it could just be a
    % string...) for 'precision' property. Perhaps it could be configurable and
    % not Constant...  Matlab default of 'double' generally simplifies things,
    % however, since the images off the scanner are typically 12bit, avoiding
    % doubles seems reasonable to avoid memory issues 
    
    % Numeric type of the image data as a function handle (used to typecast img during construction)
    precision function_handle = @single ;
    
end

properties( SetAccess=protected )
    
    % Numeric array of image data (i.e. pixel or voxel values) 
    img {mustBeNumeric} ;

end


properties( Dependent )
    
    % Vector of image dimensions (i.e. =size( Obj.img ) ), ..... (See also differ from Obj.Grid.size)
    size(1,:) { mustBePositive, mustBeInteger } = [1 1 1] ;

    % Total number of 3d image volumes in stack (scalar integer) 
    %
    % nVolumes is equivalent to prod( Img.size(4:end) )
    nVolumes(1,1) { mustBePositive, mustBeInteger } = 1 ;

end

properties( AbortSet, SetObservable )

    % Logical array specifying signal support over the grid [default=true( Img.size )]
    %
    % .mask is a Boolean array (1's and 0's) with size = size(Img) which
    % demarcates the signal spatial support over the image grid (e.g. it may be
    % used in phase unwrapping, field mapping, or prior to regridding an image,
    % in order to exclude voxels suspected of being unreliable).
    % 
    % To set the property,
    %   Img.mask = mask ;
    %  
    % size of the input mask must be equal to Img.size OR Img.Grid.size 
    % (i.e. the size of a single image volume of, for instance, a
    % multi-echo/multi-measurement stack). In the latter case, the single mask
    % is simply copied across the extra dimensions such that the assigned
    % Img.mask always possesses the same dimensions as Img.img.
    mask {mustBeNumericOrLogical} ;
    

end

properties( SetAccess=immutable, Hidden ) 
    
    % NOTE: 
    % Best attributes for the .Hdrs property are open to debate:  
    % Set access should probably be kept protected so nothing dubious happens
    % to the metadata. However, re: GetAccess, though it is currently 'Hidden'
    % (readable from anywhere --> simplifies debugging), ultimately, it should
    % probably either be made visible to the user (if useful), or
    % concealed entirely (simpler interface) with GetAccess=private.

    % Struct array of metadata (e.g. DICOM headers)
    Hdrs struct = struct( [] ) ;
    
end


properties( SetAccess=protected ) 

    % % MrdiIo object (methods for writing/reading to/from file).
    % % See also MrdiIo
    % Io MrdiIo ;

    % MrdiGrid object describing image grid properties.
    % See also MrdiGrid
    Grid MrdiGrid ;

end
    
properties( AbortSet )

    % MrdiInterp object 
    Interpolator MrdiInterp ; 

end

% =========================================================================
% =========================================================================    
methods
% =========================================================================    
function Img = Mrdi( varargin )

    if nargin == 0
        return ;
    end

    if ischar( varargin{1} ) || isstring( varargin{1} )
        warning('Direct calls to the Mrdi constructor are discouraged. Instead, use MrdiIo.make( )')
        [imgs, Hdrs] = MrdiIo.loadandsortimages( varargin{:} ) ;
        Img          = Mrdi( imgs{1}, Hdrs{1} ) ;

    elseif isnumeric( varargin{1} ) && isstruct( varargin{2} )
        img          = varargin{1} ;
        Hdrs         = varargin{2} ;
        nSlicesTotal = prod( size( img, [3:ndims(img)] ) ) ;
        nHdrsTotal   = numel( Hdrs ) ;

        if isequal( nHdrsTotal , nSlicesTotal ) 
        % if isequal( size( varargin{2} ) , size( varargin{1}, [3:ndims(varargin{1})] ) ) 

            Img.img  = varargin{1} ;
            Img.Hdrs = varargin{2} ;

            % typecast img using default precision
            % NOTE: if useful, precision could otherwise be specified in a parameters struct, or perhaps determined from the input Hdrs?
            Img.img  = Mrdi.precision( Img.img ) ;

            Img.Grid = MrdiGrid( Img.Hdrs ) ;

            Img.mask = true( Img.size ) ; % set mask only after Grid initialization

        else
            error( ['Input Hdrs must possess an entry for each 2d image slice.\n ' ...
            '(i.e. size(Hdrs) == size( img, [3:ndims(img)] ) '] , '%s' ) ;
        end 

    end


end
% =========================================================================
function nVolumes = get.nVolumes( Img )
%NVOLUMES Return total number of 3d image volumes in stack (= prod( Img.size(4:end) ) )
    
    nVolumes = prod( Img.size(4:end) ) ; 

end
% =========================================================================
function imgSize = get.size( Img )
%SIZE Return total image size over each dimension (i.e. =size( Img.img ) )
    
    imgSize = size( Img.img ) ;
    
    if numel( imgSize ) == 2
        imgSize = [ imgSize 1 ] ;
    end

end
% =========================================================================
function [] = set.mask( Img, mask )
%SETMASK Assign valid mask (a logical array of 1's and 0's) to Img.mask

    if ~islogical( mask )
        if numel(mask) == ( nnz( mask(:) == 0 ) + nnz( mask(:) == 1 ) )
            mask = logical( mask ) ;
        else
            error( 'Input mask must consist only of zeros and ones.' ) ;
        end
    end

    maskSize = size( mask ) ;
    
    if numel( maskSize ) == 2
        % Img.size and Img.Grid.size will always return at least 3 elements
        % so, append 1 to compare sizes 
        maskSize = [maskSize 1] ; 
    end

    if isequal( maskSize, Img.size )
        Img.mask = mask ;
    elseif isequal( maskSize, Img.Grid.size )
        Img.mask = repmat( mask, [1 1 1 Img.size(4:end)] ) ;
    else 
        error('Input mask must possess the same size as the image object' ) ;
    end

end
% =========================================================================
function [] = set.Interpolator( Img, NewInterpolant )
%SET.INTERPOLATOR

DEFAULTS.class               = class( scatteredInterpolant ) ; 
DEFAULTS.constructor         = str2func( DEFAULTS.class ) ;

DEFAULTS.Method              = 'linear' ;
DEFAULTS.ExtrapolationMethod = 'none' ;

% disp( 'Forming interpolant...(Computation time is proportional to image size. This may take a few minutes.)' )
%
% Interpolant = scatteredInterpolant() ;
%
% Interpolant.Method              = method ;
% Interpolant.ExtrapolationMethod = extrapolationMethod ;
%
% [X,Y,Z] = Img.Grid.gridpositions() ;
%
%
%
% v = Img.img(:,:,:,1,1) ;
% Interpolant.Values = v(:) ;

    % if nargin == 1 
    %     Params = DEFAULTS ;
    % elseif ( nargin == 2 ) && isstruct( NewInterpolant )
    %     Params = assignifempty( NewInterpolant, DEFAULTS ) ;
    % else
    %     error( 'Invalid inputs' ) ;
    % end
    %
    % if isempty( Img.Interpolant ) 
    %     Img.Interpolant = Params.constructor() ;
    % end 
    %
    % switch Params.class
    %     case DEFAULTS.class
    %         Img.Interpolant.Method              = Params.Method ;
    %         Img.Interpolant.ExtrapolationMethod = Params.ExtrapolationMethod ;
    %     
    %     otherwise
    %        warning( ['The only class currently supported for Img.Interpolant is ' DEFAULTS.class ] );
    % end
    %
    % % The following avoids the error from scatteredInterpolant when one
    % % attempts to form a 3d interpolant from a 2d input: 
    % [X,Y,Z]     = Img.Grid.gridpositions() ;
    % isValidDim0 = [ numel(unique(X(:))) numel(unique(Y(:))) numel(unique(Z(:))) ] > 1 ;
    % r0          = [X(:) Y(:) Z(:)] ;
    % r0          = r0(:, isValidDim0) ;
    %
    % if isempty( Img.Interpolant )
    %     img     = Img.img(:,:,:,1) ; % initialize with 1st image volume (arbitrary)
    %     Img.Interpolant = scatteredInterpolant( r0, DEFAULTS.Method, DEFAULTS.ExtrapolationMethod ) ;
    % end

    %
    %     Img 
    % %     Interpolant.Points = [X(:) Y(:) Z(:)]
    %     Params = DEFAULTS ;
    %
    %     || 
    %    error('Invalid input') ;
    % end
    %
    % % if isa( Interpolant, 'scatteredInterpolant' )
    % %     Img.Interpolant = Interpolant ;
    %
    % if isstruct( Params )
    %     
    %     if isfield( Interpolant, 'Method' )
    %         Img.Interpolant.Method = Interpolant.Method ;
    %     end
    %     
    %     if isfield( Interpolant, 'ExtrapolationMethod' )
    %         Img.Interpolant.ExtrapolationMethod = Interpolant.ExtrapolationMethod ;
    %     end
    %
    % else

end
% =========================================================================    
end

% =========================================================================
% =========================================================================    
methods
    %.....
    [] = filter( Img, weights, Params )
    %.....
    [F] = regrid( Img, X_Ep, Y_Ep, Z_Ep, varargin )
end
% =========================================================================
% =========================================================================
methods(Access=protected)
    %.....
    Interpolant = getinterpolant( Img, method, extrapolationMethod )
end
% =========================================================================
% =========================================================================
methods(Hidden=true)
% NOTE
%   Temporarily placing cropimg() and zeropad() here since these methods
%   1) might be deprecated
%   2) may or may not be useful
    %.....
    Img = cropimg( Img, gridSizeImgCropped, centralPoint )
    %.....
    Img = zeropad( Img, padSize, padDirection )
end
% =========================================================================
% =========================================================================
methods
% =========================================================================
function timeAverage = timeaverage( Img )
%TIMEAVERAGE  Return mean( Img.img, 5 ) 
%
% Img = TIMEAVERAGE( Img ) 

    timeAverage = mean( Img.img, 5 ) ;

end
% =========================================================================
function timeStd = timestd( Img )
%TIMESTD Return std( Img.img, 0, 5 ) 
%
% standardDeviation = TIMESTD( Img ) 

    timeStd = std( Img.img, 0, 5 ) ;

end
% =========================================================================

% =========================================================================
% =========================================================================
end


% =========================================================================
% =========================================================================
methods
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
end
% =========================================================================
% =========================================================================

end
