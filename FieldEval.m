classdef FieldEval < MaRdI
%FIELDEVAL - [B0] Field Evaluation
%
% .......
% 
% Usage
%
% Field = FieldEval( )
%
% =========================================================================
% Notes
%
% Part of series of classes pertaining to shimming:
%
%    FieldEval
%    ProbeTracking
%    ShimCal
%    ShimCom
%    ShimOpt
%    ShimSpecs
%    ShimTest 
%    ShimUse
%
% FieldEval is a MaRdI subclass [FieldEval < MaRdI]
%     
% =========================================================================
% Updated::20170709::ryan.topfer@polymtl.ca
% =========================================================================

properties
end

% =========================================================================
% =========================================================================    
methods
% =========================================================================
function Field = FieldEval( Img )
%FIELDEVAL - [B0] Field Evaluation 

Field.Hdr = [] ;
Field.img = [] ; 

if nargin ~= 0
    Field.Hdr = Img.Hdr;
    Field.img = Img.img;
end

end
% =========================================================================
function ImgCopy = copy(Img)
%COPY 
% 
% Make a copy of a FieldEval (i.e. handle) object.
% 
% ImgCopy = Copy( Img ) ;

ImgCopy     = FieldEval() ;
ImgCopy.img = Img.img;
ImgCopy.Hdr = Img.Hdr ;

end
% =========================================================================
function [LocalField, BackgroundField] = extractharmonicfield( ...
        Field, Params )
%EXTRACTHARMONICFIELD
%
% Extract (smooth) harmonic field via RESHARP (Sun, H. Magn Res Med, 2014)
%
% ------
%
%   Syntax
%
%   [LocalField, BkgrField] = EXTRACTHARMONICFIELD( Field, Params )
%
%   Returns 2 FieldEval-type Field objects: 
%       LocalField : LocalField.img is the non-harmonic (high-pass) signal element
%       BkgrField  : BkgrField.img is the harmonic (low-pass) signal element
%
%   Inputs
%       
%   Field
%       FieldEval-type object containing GRE field data in Field.img
%
%  Params
%   .filterRadius 
%       scalar filter radius [units: mm] 
%           (default = 4)
%
%   .regularization
%       Tikhonov regularization parameter
%           (default = 0)
%
%   .maxIterations
%      max iterations of conjugate gradient solver 
%           (default = 500)
%
%   .tolerance
%      min acceptable discepancy (Ax - b)/|norm(b)| for conjugate gradient solver
%           (default = 1E-6)

DEFAULT_FILTERRADIUS   = 3;
DEFAULT_REGULARIZATION = 0 ;
DEFAULT_MAXITERATIONS = 500 ;
DEFAULT_TOLERANCE     = 1E-6 ;

if ~myisfield(Params, 'filterRadius') || isempty(Params.filterRadius)
    Params.filterRadius = DEFAULT_FILTERRADIUS ;
end

if ~myisfield(Params, 'regularization') || isempty(Params.regularization)
    Params.regularization = DEFAULT_REGULARIZATION ;
end

if  ~myisfield( Params, 'maxIterations' ) || isempty(Params.maxIterations)
    Params.maxIterations = DEFAULT_MAXITERATIONS ;
end

if  ~myisfield( Params, 'tolerance' ) || isempty(Params.tolerance)
    Params.tolerance = DEFAULT_TOLERANCE ;
end

voxelSize  = Field.getvoxelsize() ;
mask = Field.Hdr.MaskingImage ;

% make spherical/ellipsoidal convolution kernel (ker)
rx = round(Params.filterRadius/voxelSize(1)) ; 
ry = round(Params.filterRadius/voxelSize(2)) ;
rz = round(Params.filterRadius/voxelSize(3)) ;

[X,Y,Z] = ndgrid(-rx:rx,-ry:ry,-rz:rz);
h = (X.^2/rx^2 + Y.^2/ry^2 + Z.^2/rz^2 <= 1);
ker = h/sum(h(:));


% erode the mask by convolving with the kernel
reducedMask = shaver( Field.Hdr.MaskingImage, round(Params.filterRadius ./ voxelSize) ) ;  

% cvsize = Field.getgridsize + [2*rx+1, 2*ry+1, 2*rz+1] -1; % linear conv size
% mask_tmp = real(ifftn(fftn(mask,cvsize).*fftn(ker,cvsize)));
% mask_tmp = mask_tmp(rx+1:end-rx, ry+1:end-ry, rz+1:end-rz); % same size


% circularshift, linear conv to Fourier multiplication
csh = [rx,ry,rz]; % circularshift

% prepare convolution kernel: delta-ker
dker = -ker;
dker(rx+1,ry+1,rz+1) = 1-ker(rx+1,ry+1,rz+1);
DKER = fftn(dker,Field.getgridsize); % dker in Fourier domain

b = ifftn( conj(DKER).*fftn( circshift( ...
    reducedMask .* circshift( ifftn( DKER.*fftn(Field.img) ), -csh), csh)));
b = b(:);

localField = cgs(@Afun, b, Params.tolerance, Params.maxIterations );
localField = real( reshape( localField, Field.getgridsize)).*reducedMask;

LocalField                       = FieldEval();
LocalField.Hdr                   = Field.Hdr ;
LocalField.Hdr.MaskingImage      = reducedMask ;
LocalField.img                   = reducedMask .* localField;

BackgroundField                  = FieldEval() ;
BackgroundField.Hdr              = Field.Hdr ;
BackgroundField.img              = reducedMask .* (Field.img  - localField);
BackgroundField.Hdr.MaskingImage = reducedMask ;

function y = Afun(x)

    x = reshape( x, Field.getgridsize );
    y = ifftn( conj(DKER) .* fftn( circshift( reducedMask.* ...
        circshift(ifftn(DKER.*fftn(x)),-csh),csh))) + Params.regularization*x;
    y = y(:);
end

end
% =========================================================================
% =========================================================================
function Field = histogramfield( Field )
%HISTOGRAMFIELD



end
% =========================================================================
function Results = assessfielddistribution( Field, voi )
%ASSESSSHIM
%
% Results = ASSESSFIELDDISTRIBUTION( Field )
% Results = ASSESSFIELDDISTRIBUTION( Field, VOI )
% 
% VOI 
%    binary array the same size as Field.img indicating the region of
%    interest over which field calculations are made. 
%    default: Field.Hdr.MaskingImage
%
% Results contains fields
%
%   .volume
%       volume of region of interest (VOI) [units: cm^3]
%
%   .std
%       standard deviation of Field.img over the volume of interest
%   
%   .norm
%       L2 norm of the field over the VOI 

if nargin < 2 || isempty(voi)
    voi = Field.Hdr.MaskingImage ;
end

voi = logical( voi ) ;

Results.volume = nnz( voi ) .* prod( 0.1*Field.getvoxelsize() )  ; % [units: cm^3]
Results.mean   = mean( Field.img( voi ) ) ;
Results.median = median( Field.img( voi ) ) ;
Results.std    = std( Field.img( voi ) ) ;
Results.norm   = norm( Field.img( voi ), 2 ) ;

Results.meanAbs   = mean( abs( Field.img( voi ) ) ) ;
Results.medianAbs = median( abs( Field.img( voi ) ) ) ;

end
% =========================================================================
function [] = plotfieldhistogram( Field, Params )
%PLOTFIELDHISTOGRAM
% 
% [] = PLOTFIELDHISTOGRAM( Field )
% [] = PLOTFIELDHISTOGRAM( Field, voi )
% [] = PLOTFIELDHISTOGRAM( Field, voi, Params )
%
% Plots histogram of the distribution in Field.img() using the Matlab 
% HISTOGRAM function.
% 
% Inputs
% 
% voi
%   volume of interest binary mask (same size as Field.img) indicating
%   region over which the distribution is to be binned.
%   [default : Field.Hdr.MaskingImage]
%
% Params
%
%   .xLimits
%        horizontal limits in Hz 
%        [default : max(abs(field(:))) .* [ -1 1 ] ]
%
%   .yLimits
%        vertical limits as count #
%
%   .fontSize
%       [default : 16]
%   .yAxisLocation
%       [default : 'left'] 
%
%   .textContent
%        string appears in textbox (e.g. 'Standard Deviation = 5 Hz')
%        [default : '']
%
%   .textPosition
%        [X Y]
%        position of textbox relative to graph's origin 


DEFAULT_ISTEXTBOX     = false ;
DEFAULT_YAXISLOCATION = 'left' ;
DEFAULT_FONTSIZE      = 16 ;

% -------
% Check inputs
% -------
if nargin < 1 || isempty( Field ) 
    disp('Error: function require 1 argument (field vector, i.e. delta B0)')
    help(mfilename);
    return;  
end

if nargin < 2 || isempty(voi)
    voi = Field.Hdr.MaskingImage ;
else
    voi = logical( voi ) ;
end

field = Field.img( voi ) ;

if nargin < 3 || isempty(Params)
    disp('Default parameters will be used')
    Params.dummy = [] ;
end

if  ~myisfield( Params, 'xLimits' ) || isempty( Params.xLimits )
    Params.xLimits = max(abs(field)) .* [ -1 1 ] ;
end

assert( length(Params.xLimits) == 2 ) ;


if  ~myisfield( Params, 'xTick' ) 
    midPoint = (max(Params.xLimits) + min(Params.xLimits))/2 ;
    Params.xTick = [round(Params.xLimits(1)/2) midPoint round(Params.xLimits(2)/2)] ;
end


if  ~myisfield( Params, 'yAxisLocation' ) || isempty(Params.yAxisLocation)
    Params.yAxisLocation = DEFAULT_YAXISLOCATION ;
end

if  ~myisfield( Params, 'fontSize' ) || isempty(Params.fontSize)
    Params.fontSize = DEFAULT_FONTSIZE ;
end

if  ~myisfield( Params, 'textPosition' ) || isempty(Params.textPosition) || ...
    ~myisfield( Params, 'textContent' ) || isempty(Params.textContent)
    Params.isTextBox = DEFAULT_ISTEXTBOX ;
else
    Params.isTextBox = true;
end

scnsize = get(0,'ScreenSize');
figure('Color',[1 1 1],...
    'Position',scnsize-[-30 -40 60 120]);

if ~myisfield( Params, 'binWidth') || isempty( Params.binWidth )
    Histo = histogram( field, 'BinMethod', 'sqrt' ) ;
    Params.nBins = ceil( sqrt(numel(field(:)))) ;
else
    Histo = histogram( field, Params.nBins, 'BinWidth', Params.binWidth );
end

if ~myisfield( Params, 'yLimits' )
    Params.yLimits = [ min( Histo.Values ) round(1.05*max( Histo.Values ) ) ] ;
end

if  ~myisfield( Params, 'yTick' ) 
    midPoint = round((max(Params.yLimits) + min(Params.yLimits))/2) ;
    Params.yTick = [Params.yLimits(1) midPoint Params.yLimits(2)] ;
end


disp('Plot Params :')
disp(' ')
disp(' .xLimits')
disp([' ' num2str(Params.xLimits)])
disp(' .yLimits')
disp([' ' num2str(Params.yLimits)])
disp(' .nBins')
disp([' ' num2str(Params.nBins)])

Histo

ylabel('Counts');
xlabel('Field (Hz)')

axis([Params.xLimits Params.yLimits]) ;

ax = gca ;

ax.XTick = Params.xTick ; 
ax.YTick = Params.yTick ;

ax.FontSize = Params.fontSize ;
ax.AmbientLightColor =[1 1 1];

ax.YAxisLocation = Params.yAxisLocation ; 

if Params.isTextBox
    text(Params.textPosition(1),Params.textPosition(2),...
        Params.textContent,...
        'FontSize', Params.fontSize ) ;
end

end
% =========================================================================
function VoxelShiftMap = computevoxelshiftmap( Field, varargin )
%COMPUTEVOXELSHIFTMAP
% 
% Returns a copy of Field as VoxelShiftMap wherein the Field.img data has been
% transformed to a voxel shift map (with units of pixels as opposed to Hz) 
% according to the formula (generalized from Jezzard & Balaban, 1995) :
%
% VoxelShiftMap.img = Field.img / pixelBandwidthPe ;
%
% VoxelShiftMap = COMPUTEVOXELSHIFTMAP( Field, pixelBandwidthPe ) 
% VoxelShiftMap = COMPUTEVOXELSHIFTMAP( Field, R, nPhaseEncode, echoSpacing ) 
% VoxelShiftMap = COMPUTEVOXELSHIFTMAP( Field, R, nPhaseEncode, echoSpacing, nInterleaves ) 
% 
%
% According to https://lcni.uoregon.edu/kb-articles/kb-0003
%
%   "You can either use the actual echo spacing and the actual number of phase
%   encodes for the number of echoes, or the effective echo spacing and the
%   number of reconstructed phase lines (easier to get from the DICOM)."
% 
%
% The pixelBandwidth in the phase encode direction accounts for these terms:
%
% pixelBandwidthPe = ( R * nInterleaves ) / ( nPhaseEncode * echoSpacing ) 
% 
% pixelBandwidthPe
%   Bandwidth per pixel (Hz/pixel) in phase encode direction. For EPI this may
%   be stored in the DICOM field (0019, 1028). This is equivalent to the
%   inverse of the "effective echo spacing", divided by the size of the image 
%   matrix in the phase encode dimension.
% 
% R
%   The in-plane acceleration factor.
% 
% nPhaseEncode
%   The actual size (e.g. # rows or columns) of the *reconstructed* image
%   in the phase encode dimension (i.e. this will be greater or equal to
%   the *acquired* number of phase encode lines, which depends on partial
%   Fourier + parallel imaging)
%
% echoSpacing
%   The actual inter-echo spacing (time between the k_x = 0 crossing of 
%   one phase encode line, to the zero-crossing on the next line). Note
%   that the actual echoSpacing = the "effective echo spacing" divided by
%   the in-plane acceleration factor R.
%
% nInterleaves
%   number of interleaves for multi-shot
%   [default: 1, single-shot]
%
%
% Another useful website with explanations for terms like phase-encode bandwidth
%
% http://support.brainvoyager.com/functional-analysis-preparation/27-pre-processing/459-epi-distortion-correction-echo-spacing.html
%
% NOTE 
%   
%   I've only assumed that the definition of pixelBandwidthPe above is consistent
%   with that of Siemens (ie. DICOM field (0019,1028)). In particular,
%   I haven't verified that nInterleaves figures into their definition or not.
%   Should be OK for single-shot.

DEFAULT_NINTERLEAVES = 1 ;

assert( (nargin >= 2) && ~isempty( Field ) ) ; 

if nargin < 4
    
    pixelBandwidthPe = varargin{1} ;

    if nargin > 2
        display('Warning: Unexpected number of input arguments. Assuming 2nd argument is the phase encode pixel BW and ignoring the arguments that follow it.');
    end

elseif nargin == 4
   
    nInterleaves = DEFAULT_NINTERLEAVES ;

   % pixelBandwidthPe = ( R * nInterleaves ) / ( nPhaseEncode * echoSpacing ) 
   pixelBandwidthPe = (varargin{1}*nInterleaves)/( varargin{2}*varargin{3} ) ;

elseif nargin == 5

   % pixelBandwidthPe = ( R * nInterleaves ) / ( nPhaseEncode * echoSpacing ) 
   pixelBandwidthPe = (varargin{1}*varargin{4})/( varargin{2}*varargin{3} ) ;

end

VoxelShiftMap     = Field.copy() ;
VoxelShiftMap.img = Field.img ./ ( pixelBandwidthPe ) ;

VoxelShiftMap.Hdr.PixelComponentPhysicalUnits = '0000H' ; % i.e. none

end
% =========================================================================

end
% =========================================================================
% =========================================================================
methods(Static)
% =========================================================================
function [Field, Extras] = mapfield( ImgArray, Params, ObjectiveImg )
% MAPFIELD
%
% Field = MAPFIELD( ImgArray ) 
% Field = MAPFIELD( ImgArray, Params ) 
% Field = MAPFIELD( ImgArray, Params, ObjectiveImg ) 
% 
% ImgArray --- a cell array where the 1st dimension (i.e. row) corresponds to
% echo number and the index along the second dimension (i=1,2) corresponds to
% Magnitude & Wrapped Phase respectively (each a MaRdI-type object) 
%
% if nEchoes == 1 
%   (i.e. ImgArray = { Mag, Phase } )
%
%   It is assumed that Phase is in fact a phase-difference image.
%
% if nEchoes == 2
%   (i.e. ImgArray = { MagEcho1, PhaseEcho1 ; MagEcho2, PhaseEcho2 } )
%
% Params --- *must* contain the following fields
%
% if nEchoes == 1
%   .echoTimeDifference [units: ms]
%       
%
% Params --- may contain the following fields
%
%   .mask
%       binary array indicating phase region to be unwrapped 
%       [default: formed by thresholding magnitude images > Params.threshold)
%
%   .threshold  
%   (as a fraction (<1) of max measured magnitude intensity)
%   Determines the phase region to be unwrapped (i.e. areas of low signal are
%   ignored) 
%       [default: 0.01] 
%
%   .isFilteringField 
%       Applies harmonic ('RESHARP', i.e. low-pass) filtering & returns filtered field map.
%       [default: false]
%   
%   .filteringMask
%       binary array indicating local ROI to retain + RESHARP filter 
%       [default: uses the region retained after unwrapping, which may generally be = Params.mask)
%
%   .unwrapper
%       'AbdulRahman_2007' [default], calls unwrap3d( ), which wraps to the Abdul-Rahman binary
%       'FSL_Prelude', calls prelude( ), which wraps to FSL-prelude 

% TODO 
%
% ...Untested features
%
% ObjectiveImg --- a MaRdI-type object *example* image (e.g. EPI volume!) 
% that, most importantly, has its voxels positioned precisely where the field
% information/shim is desired. That is, if a spherical harmonic fitting of
% the field is performed (see below) then the fitted-field will be interpolated
% at the voxel positions from ObjectiveImg. 
%
% What is the immediate purpose of this? To enabled gapped/partial gre_field_map
% acquisitions in order to reduce the acquisition time for real-time shim training
% i.e. these generally involve the subject holding their breath, so, for
% feasibility + patient comfort, the acquisitions should be as brief as
% possible.
%
%   .isFittingSphericalHarmonics 
%       Fits field map to spherical harmonic basis set & returns the fitted field map.
%       See doc FieldEval.extractharmonicfield() for relevant Params.
%       [default: false]
%
%   .ordersToGenerate
%       Orders of SH incorporated into fitting 
%       [default: [0:1:8]]

DEFAULT_ISCORRECTINGPHASEOFFSET        = true ; % deprecated
DEFAULT_ISUNWRAPPINGECHOESINDIVIDUALLY = false ; % deprecated
DEFAULT_ISFILTERINGFIELD               = false ;
DEFAULT_THRESHOLD                      = 0.01 ;

DEFAULT_ISFITTINGSPHERICALHARMONICS    = false ;
DEFAULT_ORDERSTOGENERATE               = [0:8] ;

DEFAULT_UNWRAPPER                      = 'AbdulRahman_2007' ;

assert( (nargin >= 1) && ~isempty( ImgArray ) ) ;

if ~myisfield( Params, 'isCorrectingPhaseOffset' ) || isempty( Params.isCorrectingPhaseOffset ) 
    Params.isCorrectingPhaseOffset = DEFAULT_ISCORRECTINGPHASEOFFSET ;
end

if ~myisfield( Params, 'isUnwrappingEchoesIndividually' ) || isempty( Params.isUnwrappingEchoesIndividually ) 
    Params.isUnwrappingEchoesIndividually = DEFAULT_ISUNWRAPPINGECHOESINDIVIDUALLY ;
end

if ~myisfield( Params, 'isFilteringField' ) || isempty( Params.isFilteringField ) 
    Params.isFilteringField = DEFAULT_ISFILTERINGFIELD ;
end

if ~myisfield( Params, 'threshold' ) || isempty( Params.threshold ) 
    Params.threshold = DEFAULT_THRESHOLD ;
end

if ~myisfield( Params, 'isFittingSphericalHarmonics' ) || isempty( Params.isFittingSphericalHarmonics ) 
    Params.isFittingSphericalHarmonics = DEFAULT_ISFITTINGSPHERICALHARMONICS ;
end

if ~myisfield( Params, 'ordersToGenerate' ) || isempty( Params.ordersToGenerate ) 
    Params.ordersToGenerate = DEFAULT_ORDERSTOGENERATE ;
end

if ~myisfield( Params, 'unwrapper' ) || isempty( Params.unwrapper ) 
    Params.unwrapper = DEFAULT_UNWRAPPER ;
end

Extras = [] ;

% .......
nEchoes = size( ImgArray, 1 ) ;
% -------
% define spatial support for unwrapping
if ~myisfield( Params, 'mask' ) || isempty( Params.mask )

    Params.mask = ones( ImgArray{1,1}.getgridsize ) ;

    for iEcho = 1 : nEchoes 

        Params.mask = Params.mask .* ( ImgArray{ iEcho, 1 }.img > Params.threshold ) ;

    end

end

if nEchoes == 1
    
    PhaseDiff = ImgArray{ 1, 2 } ;
    PhaseDiff.Hdr.EchoTime = Params.echoTimeDifference ;

else
    
    % -------
    % phase difference image via complex division
    PhaseDiff      = ImgArray{ 1, 2}.copy() ;

    img            = ImgArray{ 1, 1 }.img .* exp(-i*ImgArray{ 1, 2 }.img) ;
    img(:,:,:,2)   = ImgArray{ 2, 1 }.img .* exp(-i*ImgArray{ 2, 2 }.img) ;

    PhaseDiff.img = angle( img(:,:,:,2) ./ img(:,:,:,1) ) ;

    PhaseDiff.Hdr.EchoTime = ( ImgArray{ 1, 2 }.Hdr.EchoTime - ImgArray{ 2, 2 }.Hdr.EchoTime );
end

% -------
% 3d path-based unwrapping
PhaseDiff.Hdr.MaskingImage = Params.mask ;

if strcmp( Params.unwrapper, 'AbdulRahman_2007' )

    PhaseDiff = PhaseDiff.unwrapphase(  ) ;

elseif strcmp( Params.unwrapper, 'FSL_Prelude' )
    
    Options.voxelSize = ImgArray{1,1}.getvoxelsize() ;
    PhaseDiff = PhaseDiff.unwrapphase( ImgArray{1,1}, Options ) ;

end

Field     = FieldEval( PhaseDiff.scalephasetofrequency( ) ) ;

if Params.isFilteringField

    if myisfield( Params, 'filteringMask' ) && ~isempty( Params.filteringMask )
        
        % intersection of unwrapped region with desired local region
        Field.Hdr.MaskingImage = logical(Field.Hdr.MaskingImage) & logical(Params.filteringMask) ;    

    end
    
    [Extras.LocalField,Field] = Field.extractharmonicfield( Params ) ;

end

if Params.isFittingSphericalHarmonics % UNTESTED
    % -------
    % fit spherical harmonic basis set to input Field 

    % generates basis set, with field positions same as those of input Field 
    Shims = ShimOptSHarmonics( Params, Field ) ;

    % calculate fitting coefficients ('currents')
    Shims = Shims.optimizeshimcurrents( Params ) ;
    
    Extras.FieldResidual = Field.img + Shims.Model.field ;

    Field.img = -Shims.Model.field ;
    
    if (nargin == 3) & ~isempty( ObjectiveImg )
        % Interpolate the field @ VoxelPositions
        %
        % Main purpose: to enable gapped slices in the field map acquisitions
        % for real-time shim training --- by reducing nSlices, acq. time is
        % reduced, & therefore, the necessary duration of the breath hold.

        [X0, Y0, Z0] = Field.getvoxelpositions( ) ; % original
        [X, Y, Z]    = ObjectiveImg.getvoxelpositions( ) ; % final

        % recalculate the set of harmonics at the given voxel positions      
        basisFields = ShimOptSHarmonics.generatebasisfields( Params.ordersToGenerate, X, Y, Z ) ;
        
        % scale each harmonic by the fitted 'currents' (coefficients)
        for iHarmonic = 1 : size( basisFields, 4 ) 
            basisFields(:,:,:, iHarmonic) = Shims.Model.currents(iHarmonic) * basisFields(:,:,:, iHarmonic) ;
        end
        
        Field.img = sum( -basisFields, 4 ) ;
        

        disp( ['Interpolating phase/field mask...' ]) ;
            
        Field.Hdr.MaskingImage = griddata( X0, Y0, Z0, Field.Hdr.MaskingImage, X, Y, Z, 'nearest' ) ;

        % if new positions are outside the range of the original, 
        % interp3/griddata replaces array entries with NaN
        Field.Hdr.MaskingImage( isnan( Field.Hdr.MaskingImage ) ) = 0 ; 

        % -------
        % Update Hdr 
        %
        % Note: the Hdr could probably simply be copied from ObjectiveImg but recalculating the entries 
        % is more general ('extensible') should the future user not have a
        % fully-formed 'ObjectiveImg' set of dicoms, but merely the target
        % voxel positions [X,Y,Z]
        % 
        % That said, the way SliceLocation is updated below may not always be correct.
        % (borrowed from MaRdI.resliceimg() )
        
        Field.Hdr.ImagePositionPatient( 1 ) = X(1) ; 
        Field.Hdr.ImagePositionPatient( 2 ) = Y(1) ;
        Field.Hdr.ImagePositionPatient( 3 ) = Z(1) ;

        %-------
        % Rows 
        Field.Hdr.Rows = size(Field.img, 1) ;

        dx = X(2,1,1) - X(1,1,1) ;
        dy = Y(2,1,1) - Y(1,1,1) ;
        dz = Z(2,1,1) - Z(1,1,1) ;  

        % vertical (row) spacing
        Field.Hdr.PixelSpacing(1) = ( dx^2 + dy^2 + dz^2 )^0.5 ; 

        % column direction cosine (expressing angle btw column direction and X,Y,Z axes)
        Field.Hdr.ImageOrientationPatient(4) = dx/Field.Hdr.PixelSpacing(1) ;
        Field.Hdr.ImageOrientationPatient(5) = dy/Field.Hdr.PixelSpacing(1) ;
        Field.Hdr.ImageOrientationPatient(6) = dz/Field.Hdr.PixelSpacing(1) ;

        %-------
        % Columns 
        Field.Hdr.Columns = size(Field.img, 2) ;       

        dx = X(1,2,1) - X(1,1,1) ;
        dy = Y(1,2,1) - Y(1,1,1) ;
        dz = Z(1,2,1) - Z(1,1,1) ;  

        % horizontal (column) spacing
        Field.Hdr.PixelSpacing(2) = ( dx^2 + dy^2 + dz^2 )^0.5 ;

        % row direction cosine (expressing angle btw column direction and X,Y,Z axes)
        Field.Hdr.ImageOrientationPatient(1) = dx/Field.Hdr.PixelSpacing(2) ;
        Field.Hdr.ImageOrientationPatient(2) = dy/Field.Hdr.PixelSpacing(2) ;
        Field.Hdr.ImageOrientationPatient(3) = dz/Field.Hdr.PixelSpacing(2) ;

        %-------
        % Slices
        Field.Hdr.NumberOfSlices       = size(Field.img, 3) ;
        Field.Hdr.SpacingBetweenSlices = ( (X(1,1,2) - X(1,1,1))^2 + ...
                                           (Y(1,1,2) - Y(1,1,1))^2 + ...
                                           (Z(1,1,2) - Z(1,1,1))^2 ) ^(0.5) ;

        [~, ~, sHat] = Field.getdirectioncosines( ) ;  
        Field.Hdr.SliceLocation = dot( Field.Hdr.ImagePositionPatient, sHat ) ;
    end

end

end
% =========================================================================

end
% =========================================================================
% =========================================================================

end

