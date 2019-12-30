classdef (Abstract) MrdiProc < handle
%MrdiProc   MR DICOM Image Processing 
%
% Member methods for image processing 
% 
% e.g.
%   filter()
%   resliceimg()
%   setmaskingimage()
%   ...etc.
% 
% For documentation, type doc MrdiProc
%
% =========================================================================
% Author::ryan.topfer@polymtl.ca
% =========================================================================

% =========================================================================
% =========================================================================    
methods
% =========================================================================    
function Proc = MrdiProc(  )

end
% =========================================================================
% =========================================================================
end

% =========================================================================
% =========================================================================    
methods
% =========================================================================    
function [] = filter( Img, weights, Params )
%FILTER     3D low-pass filtering
% 
% Description
%
%   Filtering can be weighted, unweighted, or median. Wraps to smooth3() or
%   medfilt3() accordingly.
%
% [] = FILTER( Img )
% [] = FILTER( Img, weights )
% [] = FILTER( Img, weights, Params )
%
% Img 
%   the Mrdi-type image volume.
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
% TODO
%   Add support for 2d (single slice) images 

DEFAULTS.kernelSize = [3 3 3] ;
DEFAULTS.method     = 'gaussian' ;

if nargin < 3 || isempty(Params)
    Params.dummy = [] ;
end

Params = assignifempty( Params, DEFAULTS ) ;

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
function [F] = resliceimg( Img, X_Ep, Y_Ep, Z_Ep, varargin ) 
%RESLICEIMG     Interpolate a Mrdi image object and update Img.Hdr accordingly.
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
%       image grid (Img) to another (Mrdi-type object Img2), these arrays are
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
if MrdiInfo.compareimggrids( X_Ip, Y_Ip, Z_Ip, X_Ep, Y_Ep, Z_Ep ) 
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
            error( 'Unknown input. See HELP MrdiProc.resliceimg' ) ;
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
function [] = setmaskingimage( Img, mask )
%SETMASKINGIMAGE Copies valid mask (a logical array of 1's and 0's) to Img.Hdr.MaskingImage
% 
% [] = SETMASKINGIMAGE( Img, mask ) 
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

assert( nargin == 2, 'Function requires at least 2 input arguments. See HELP MrdiProc.setmaskingimage' ) ;

if ~islogical( mask )
    assert( numel(mask) == ( nnz( mask(:) == 0 ) + nnz( mask(:) == 1 ) ), ...
        'Input mask must consist exclusively of zeros and ones. See HELP MrdiProc.setmaskingimage' ) ;
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

end
% =========================================================================
% =========================================================================
methods(Access=protected)
% =========================================================================
function Interpolant = getinterpolant( Img, method, extrapolationMethod )
%GETINTERPOLANT     Return scatteredInterpolant object
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


end
