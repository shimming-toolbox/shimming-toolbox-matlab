classdef FieldEval < MaRdI
%FIELDEVAL - [B0] Field Evaluation
%
% .......
% 
% Usage
%
% Field = FieldEval( Mag, Phase ) 
% Field = FieldEval( Mag, Phase, Params ) 
%
% Mag and Phase should either be paths to the respective DICOM directories, 
% OR, instantiated MaRdI-type image objects (e.g. Mag = MaRdI( path_to_Mag_DICOMs ) )
%
% Params may contain the following fields
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
%   .unwrapper
%       'Sunwrap' [default for 2d image (single slice)], calls sunwrap( ) (Maier, et al. MRM 2015)
%       'AbdulRahman_2007' [default for 3d image volume], calls unwrap3d( ), which wraps to the Abdul-Rahman binary
%       'FslPrelude', calls prelude( ), which wraps to FSL-prelude 
%
% .......
%
% NOTE
%
% FieldEval is a MaRdI subclass [FieldEval < MaRdI]
%     
% =========================================================================
% Author::ryan.topfer@polymtl.ca
% =========================================================================

properties
    Model ;
end

% =========================================================================
% =========================================================================    
methods
% =========================================================================
function Field = FieldEval( varargin )

Field.img   = [] ; 
Field.Hdr   = [] ;
Field.Hdrs  = [] ;
Field.Aux   = [] ;
Field.Model = [] ;

if nargin ~= 0 
    
    if nargin == 1
        Params.dummy = [] ;
    end

    if isa( varargin{1}, 'MaRdI' )
        % convert MaRdI-type Img object to FieldEval
        Img        = varargin{1} ;
        Field.img  = Img.img ;
        Field.Hdr  = Img.Hdr ;
        Field.Hdrs = Img.Hdrs ;
        Field.Aux  = Img.Aux ;
    
    elseif isa( varargin{1}, 'cell' ) 
        % if Img input is a cell array of MaRdI-images, this is presumed to be a set of GRE images destined for B0-mapping 
        Img = varargin{1} ; 

        for iImg = 1 : numel( Img ) 
            assert( isa( Img{ iImg }, 'MaRdI' ) ) ;
        end
        
        if nargin == 2
           Params = varargin{2} ;
        end
        
        Field = FieldEval.mapfield( Img, Params ) ;
    
    elseif ischar( varargin{1} ) 
       % input should be the path-string to the magnitude and phase directories of the GRE images 

       assert( nargin >= 2, 'Expected pathToMagnitude and pathToPhase (DICOM directories) as comma-separated arguments' )
       imgPath1 = varargin{1} ;
       imgPath2 = varargin{2} ; 
        
       if ( nargin == 3 ) 
           assert( isstruct( varargin{3} ), 'Expected Parameters struct as third input argument' ) ;
           Params = varargin{3} ;
       end
       
       Field = FieldEval.mapfield( imgPath1, imgPath2, Params ) ;

    end
end

end
% =========================================================================
function ImgCopy = copy(Img)
%COPY 
% 
% Make a copy of a FieldEval (i.e. handle) object.
% 
% ImgCopy = Copy( Img ) ;

%% these 2 lines are the only parts that differs from MaRdI.copy():
ImgCopy       = FieldEval() ;
ImgCopy.Model = Img.Model ;
%% 

ImgCopy.img   = Img.img;
ImgCopy.Hdr   = Img.Hdr ;
ImgCopy.Hdrs  = Img.Hdrs ;

if ~isempty( Img.Aux ) && myisfield( Img.Aux, 'Tracker' ) 
    ImgCopy.Aux.Tracker = Img.Aux.Tracker.copy() ;
end

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

voxelSize  = Field.getvoxelspacing() ;
mask = Field.Hdr.MaskingImage ;

% make spherical/ellipsoidal convolution kernel (ker)
rx = round(Params.filterRadius/voxelSize(1)) ; 
ry = round(Params.filterRadius/voxelSize(2)) ;
rz = round(Params.filterRadius/voxelSize(3)) ;

[X,Y,Z] = ndgrid(-rx:rx,-ry:ry,-rz:rz);
h = (X.^2/rx^2 + Y.^2/ry^2 + Z.^2/rz^2 <= 1);
ker = h/sum(h(:));


% cvsize = Field.getgridsize + [2*rx+1, 2*ry+1, 2*rz+1] -1; % linear conv size
% mask_tmp = real(ifftn(fftn(mask,cvsize).*fftn(ker,cvsize)));
% mask_tmp = mask_tmp(rx+1:end-rx, ry+1:end-ry, rz+1:end-rz); % same size


% circularshift, linear conv to Fourier multiplication
csh = [rx,ry,rz]; % circularshift

tfs = Field.img ;
% tfs  = padarray( Field.img, csh ) ;
% mask = padarray( mask, csh ) ;
%
imsize = size( tfs );

% erode the mask by convolving with the kernel
reducedMask = shaver( mask, round(Params.filterRadius ./ voxelSize) ) ;  

% prepare convolution kernel: delta-ker
dker = -ker;
dker(rx+1,ry+1,rz+1) = dker(rx+1,ry+1,rz+1) + 1 ;
DKER = fftn( dker, imsize ); % dker in Fourier domain

b = ifftn( conj(DKER).*fftn( circshift( ...
    reducedMask .* circshift( ifftn( DKER.*fftn( tfs ) ), -csh), csh) ) );
b = b(:) ;

localField = cgs( @Afun, b, Params.tolerance, Params.maxIterations ) ;
localField = real( reshape( localField, imsize ) ) .* reducedMask;

LocalField                       = Field.copy();
LocalField.Hdr.MaskingImage      = reducedMask ;
LocalField.img                   = reducedMask .* localField;

BackgroundField                  = Field.copy() ;
BackgroundField.Hdr.MaskingImage = reducedMask ;
BackgroundField.img              = reducedMask .* (Field.img  - localField);

function y = Afun(x)

    x = reshape( x, imsize );
    y = ifftn( conj(DKER) .* fftn( circshift( reducedMask.* ...
        circshift(ifftn(DKER.*fftn(x)),-csh),csh))) + Params.regularization*x;
    y = y(:);
end

end
% =========================================================================
function Stats = assessfielddistribution( Field, voi, filename )
%ASSESSSHIM
%
% Stats = ASSESSFIELDDISTRIBUTION( Field )
% Stats = ASSESSFIELDDISTRIBUTION( Field, VOI )
% Stats = ASSESSFIELDDISTRIBUTION( Field, VOI, filename )
% 
% VOI 
%    binary array the same size as Field.img indicating the region of
%    interest over which field calculations are made. 
%    default: Field.Hdr.MaskingImage
%
% filename
%   output to text file using writetable()
%
% Stats contains fields
%
%   .volume
%       volume of region of interest (VOI) [units: cm^3]
%
%   .mean
%       mean value of the field over the VOI
%
%   .median
%       median value of the field over the VOI
%
%   .std
%       standard deviation of Field.img over the VOI
%   
%   .rmsePerCm3
%       L2 norm of the field (i.e. residual) over the VOI normalized by the volume
%
%   .meanAbs
%       mean absolute value of the field over the VOI.
%
%   .medianAbs
%       median absolute value of the field over the VOI
%
%   .min
%   .max
assert( ndims( Field.img ) <=3, 'Multiple volumes not supported. TODO' )

if nargin < 2 || isempty(voi)
    voi = Field.Hdr.MaskingImage ;
end

voi = logical( voi ) ;

Stats.volume     = nnz( voi ) .* prod( 0.1*Field.getvoxelspacing() )  ; % [units: cm^3]
Stats.mean       = mean( Field.img( voi ) ) ;
Stats.median     = median( Field.img( voi ) ) ;
Stats.std        = std( Field.img( voi ) ) ;
Stats.rmsePerCm3 = norm( Field.img( voi ), 2 )/Stats.volume ;
Stats.meanAbs    = mean( abs( Field.img( voi ) ) ) ;
Stats.stdAbs     = std( abs( Field.img( voi ) ) ) ;
Stats.min        = min( ( Field.img( voi ) ) ) ;
Stats.max        = max( ( Field.img( voi ) ) ) ;

if nargin == 3 && ischar( filename ) 
    measure = {'Volume (cm^3)'; 'Mean (Hz)' ; 'Median (Hz)' ; 'St. dev. (Hz)' ; 'RMSE/Volume (Hz/cm^3)' ; 'Mean[abs.] (Hz)'; 'Std[abs.] (Hz)'; 'Min (Hz)'; 'Max (Hz)'} ;
    value   = num2str([ Stats.volume ; Stats.mean ; Stats.median ; Stats.std ; Stats.rmsePerCm3 ; Stats.meanAbs ; Stats.stdAbs ; Stats.min ; Stats.max ;], 4 ) ;
    writetable( table( measure, value ), filename ) ;
end

end
% =========================================================================
function [mask] = getvaliditymask( Field, maxAbsField )
%GETVALIDITYMASK 
%
% Returns binary mask - TRUE where field values are well defined and within 
% the expected range
%
% mask = GETVALIDITYMASK( Field )
% mask = GETVALIDITYMASK( Field, maxAbsField ) 
%
% .......................
%   
%
% maxAbsField 
%   maximum absolute voxel value assumed to represent an accurate field
%   measurement. Voxels with abs-values greater than this might stem from
%   errors in the unwrapping.  [default: 500 Hz]
%
% (Set to Inf to ignore the criterion)

DEFAULT_MAXABSFIELD        = 500 ;

if nargin < 2
    maxAbsField = DEFAULT_MAXABSFIELD ;
end

% NOTE : This assumes a valid measurement will never be exactly zero. 
% Perhaps better/additional conditions could be used ?
mask = ~( (Field.img==0) & (-Field.img==0) ) ; 

mask = mask & logical( Field.Hdr.MaskingImage ) ;

mask = mask & ( abs(Field.img) <= maxAbsField ) ;

end
% =========================================================================
function [] = plotfieldhistogram( Field, voi, Params )
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
function [Riro] = getriro( Field, p )
%GETRIRO
%
% Riro = GETRIRO( Field )
% Riro = GETRIRO( Field, p )
%
% Returns estimate of the respiration-induced resonance offset corresponding
% to tracker measurement p assuming the linear field model 
% (i.e. Riro[ p(t) ] ).
%
% If nargin == 1, Riro is a copy of Field.Model.Riro 
% (e.g. inspired-expired field difference);

assert( myisfield( Field.Model, 'Riro' ) && ~isempty( Field.Model.Riro ) ) ;

Riro = Field.Model.Riro.copy() ;

if nargin == 2
    % return instantaneous RIRO corresponding to tracker value p
    assert( isscalar(p) ) ;

    dfdp = Field.Model.Riro.img()/Field.Model.Riro.Aux.Tracker.Data.p ;
    dp   = Field.Aux.Tracker.debias( p ) ; 

    Riro.img = dfdp*dp ;
    Riro.Aux.Tracker.Data.p = p ;
end

end
% =========================================================================
function [] = chartriro( Field, p )
%CHARTRIRO
%
% [] = CHARTRIRO( Field )
% [] = CHARTRIRO( Field, p )
%
assert( myisfield( Field.Model, 'Riro' ) && ~isempty( Field.Model.Riro ) ) ;

dbstop in FieldEval at 599

pR = resample( p, 30, 100 ) ; % resample from 100 to 30 Hz
pR = medfilt1( pR, 5 ) ;

nFrames = length( pR ) ;
riro = zeros( [ Field.Model.Riro.getgridsize() nFrames ] ) ;

for iFrame = 1 : nFrames
    iFrame/nFrames
    iRiro = Field.getriro( pR( iFrame ) ) ;
    riro(:,:,:, iFrame) = iRiro.img ;    
end

end
% =========================================================================
function [] = scalefieldstrength( Field, B01 )
%SCALEFIELDSTRENGTH
%
% [] = SCALEFIELDSTRENGTH( Field, B01 )
%
% Scales Field (and, if present, FieldEval objects within Field.Model) to
% new (scanner) "main" field strength B01.

B00 = Field.Hdr.MrProt.sProtConsistencyInfo.flNominalB0 ; % original field [units: Tesla]

Field.img = Field.img * (B01/B00) ; % update field [units: Hz] 
Field.Hdr.MrProt.sProtConsistencyInfo = B01 ; % update header

if myisfield( Field.Model, 'Riro' )
    Field.Model.Riro.scalefieldstrength( B01 ) ;
end
if myisfield( Field.Model, 'Zero' )
    Field.Model.Zero.scalefieldstrength( B01 ) ;
end

end
% =========================================================================

end
% =========================================================================
% =========================================================================
methods(Static)
% =========================================================================
function [Field] = mapfield( varargin )
%MAPFIELD
%
% Field = MAPFIELD( Mag, Phase ) 
% Field = MAPFIELD( Mag, Phase, Params ) 
%
% Mag and Phase should either be paths to the respective DICOM directories, 
% OR, instantiated MaRdI-type image objects (e.g. Mag = MaRdI( path_to_Mag_DICOMs ) )
%
% Params may contain the following fields
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
%   .unwrapper
%       'Sunwrap' [default for 2d image (single slice)], calls sunwrap( ) (Maier, et al. MRM 2015)
%       'AbdulRahman_2007' [default for 3d image volume], calls unwrap3d( ), which wraps to the Abdul-Rahman binary
%       'FslPrelude', calls prelude( ), which wraps to FSL-prelude 

DEFAULT_THRESHOLD    = 0.01 ;
DEFAULT_UNWRAPPER_2D = 'Sunwrap' ; 
DEFAULT_UNWRAPPER_3D = 'AbdulRahman_2007' ;

assert( ( nargin >= 2 ) && ...
        isa( varargin{1}, 'MaRdI' ) || ischar( varargin{1} ) && ...
        isa( varargin{2}, 'MaRdI' ) || ischar( varargin{2} ), ...
        'Invalid input.' ) ;

for iImg = 1 : 2
    if isa( varargin{iImg}, 'MaRdI' )
        Img = varargin{iImg} ;
    elseif ischar( varargin{iImg} )
        display( ['Loading img ' num2str(iImg) ' of 2'] );
        Img = MaRdI( varargin{iImg} ) ;
    end

    if ~isempty( strfind( Img.Hdr.ImageType, '\M\' ) )
        Mag = Img.copy() ;
    elseif ~isempty( strfind( Img.Hdr.ImageType, '\P\' ) )
        Phase = Img.copy() ;
    else
        error('Unexpected image type: Neither magnitude nor phase?') ;
    end
end

if isstruct( varargin{end} )
    Params = varargin{end} ;
else
    Params.dummy = [];
end

if ~myisfield( Params, 'threshold' ) || isempty( Params.threshold ) 
    Params.threshold = DEFAULT_THRESHOLD ;
end

if ~myisfield( Params, 'unwrapper' ) || isempty( Params.unwrapper ) 
    nSlices = size( Phase.img, 3 ) ;
    if nSlices == 1
        Params.unwrapper = DEFAULT_UNWRAPPER_2D ;
    else
        Params.unwrapper = DEFAULT_UNWRAPPER_3D ;
    end
end
    
nEchoes       = length( Phase.echotime ) ;
nMeasurements = size( Phase.img, 5 ) ;

% -------
% define spatial support for unwrapping
if ~myisfield( Params, 'mask' ) || isempty( Params.mask )

    Params.mask = true( [Phase.getgridsize() nEchoes nMeasurements] ) ;

    for iMeasurement = 1 : nMeasurements
        for iEcho = 1 : nEchoes 
            mag_i = Mag.img(:,:,:,iEcho,iMeasurement) ;
        
            Params.mask(:,:,:,iEcho,iMeasurement) = Params.mask(:,:,:,iEcho,iMeasurement) .* ...
                ( mag_i ./ max( mag_i(:) ) ) >= Params.threshold ; 
        end
    end

    Params.mask = prod( Params.mask, 4 ) ;
end

Params.mask = logical( Params.mask ) ;

% % ------- TODO Add support for separate excitations for each TE
% if ImgArray{ 1, 2}.Hdr.MrProt.lContrasts == 1 
%     
%     if nEchoes == 1 % single echo phase image
%         % in this case only, PhaseDiff refers to the difference between t=TE and t=0
%         PhaseDiff = ImgArray{ 1, 2 }.copy() ; 
%     else % echoes acquired using separate excitations
%         %TODO:implement check/assertion that the echoes were indeed acquired using separate excitations 
%         
%         % correct for possible frequency difference between excitations using f0_Tx of 1st TE as reference
%         f0_1 = ImgArray{1,2}.Hdr.ImagingFrequency*1E6 ; % [units: Hz]
%
%         for iEcho = 2 : nEchoes
%             
%             delta_f0 = ImgArray{iEcho, 2}.Hdr.ImagingFrequency*1E6 - f0_1 ;
%             
%             if ( delta_f0 ~= 0 )
%                 warning('Modifying the MaRdI inputs...')
%                 % TODO: work-around that doesn't modify the inputs in this way
%                 ImgArray{ iEcho, 2 }.img = ImgArray{ iEcho, 2 }.img + 2*pi*delta_f0*ImgArray{ iEcho, 2 }.Hdr.EchoTime/1000 ;
%                 ImgArray{ iEcho, 2 }.Hdr.ImagingFrequency = ImgArray{1,2}.Hdr.ImagingFrequency ; 
%             end 
%
%         end
%
%         PhaseDiff      = ImgArray{ 1, 2}.copy() ;
%
%         img            = ImgArray{ 1, 1 }.img .* exp(i*ImgArray{ 1, 2 }.img) ;
%         img(:,:,:,2)   = ImgArray{ 2, 1 }.img .* exp(i*ImgArray{ 2, 2 }.img) ;
%
%         PhaseDiff.img = angle( img(:,:,:,2) ./ img(:,:,:,1) ) ;
%
%         PhaseDiff.Hdr.EchoTime = ImgArray{ 2, 2 }.Hdr.EchoTime - ImgArray{ 1, 2 }.Hdr.EchoTime ; % [units : ms]
%         % Cornell Complex fitting...        
%         M  = zeros( [ImgArray{1,1}.getgridsize() nEchoes] ) ;
%         TE = zeros( nEchoes, 1 ) ;
%
%         for iEcho = 1 : nEchoes
%             M( :,:,:, iEcho ) = ImgArray{ iEcho, 1}.img .* exp(i*-ImgArray{ iEcho, 2}.img) ;
%             TE(iEcho) = ImgArray{iEcho,1}.Hdr.EchoTime/1000 ;
%         end
%
%         [p1, dp1, relres, p0, iter]=Fit_ppm_complex_TE(M,TE);
%
%         P1 = ImgArray{1,2}.copy() ;
%         P1.img = p1 ;
%         P1.Hdr.EchoTime = TE(2) - TE(1) ;
%         Field3 = P1.copy() ;
%         Field3.scalephasetofrequency() ;
%         nii((p1/(1000*(TE(2)-TE(1)))/(2*pi)));
%     end
%     
% else
%
%
% end

if nEchoes == 1 
    % Input phase should be a phase *difference* image; alternatively, PhaseDiff will refer
    % to the difference between t=TE and t=0
    if( numel( Phase.Hdr.MrProt.alTE ) < 2 )
        warning( ['Expected Siemens phase difference image with at least 2 TEs specified in the DICOM header.' ...
           'Output field estimate may retain an offset term'] )
    end

    PhaseDiff = Phase.copy() ;
    PhaseDiff.Hdr.EchoTime = Phase.echotime ;

else
    % -------
    % phase difference image via complex division of first 2 echoes
    PhaseDiff      = Phase.copy() ;

    img            = Mag.img(:,:,:,1,:) .* exp( i*Phase.img(:,:,:,1,:) ) ;
    img(:,:,:,2,:) = Mag.img(:,:,:,2,:) .* exp( i*Phase.img(:,:,:,2,:) ) ;

    PhaseDiff.img  = angle( img(:,:,:,2,:) .* conj(img(:,:,:,1,:) ) ) ;

    PhaseDiff.Hdr.EchoTime = Phase.echotime(2) - Phase.echotime(1) ; % [units : ms]

end

PhaseDiff.Hdr.MaskingImage = logical( Params.mask ) & ~isnan( PhaseDiff.img ) ;

% -------
% 3d path-based unwrapping
PhaseDiff = PhaseDiff.unwrapphase( Mag, Params ) ;

% if time series, correct for potential 2pi wraps between time points:
nImg = size( PhaseDiff.img, 5 ) ;

if nImg > 1
    display('Correcting for temporal wraps in field time-series')

    % correct temporal wraps by comparing to phaseEstimate:
    phaseEstimate = median(PhaseDiff.img,  5, 'omitnan' ) ;
    % normalize the estimate so its spatial median is within [-pi,pi]
    tmpMask = logical( prod( Params.mask, 5 ) ) ;
    n       = round(median(phaseEstimate(tmpMask))/pi) ;
    phaseEstimate = phaseEstimate - n*pi ;

    % Wherever the absolute deviation from the estimate exceeds pi,
    % correct the measurement by adding the appropriate pi-multiple:
    for iImg = 1 : nImg
        dPhase = phaseEstimate - PhaseDiff.img( :,:,:,1,iImg )  ;
        n      = ( abs(dPhase) > pi ) .* round( dPhase/pi ) ;
        PhaseDiff.img( :,:,:,1,iImg ) = PhaseDiff.img(:,:,:,1,iImg ) + n*pi ;
    end
    
end

PhaseDiffInRad = PhaseDiff.copy() ; % PhaseDiffInRad.img : [units: rad]

PhaseDiff.scalephasetofrequency( ) ; % PhaseDiff.img : [units: Hz]
Field     = FieldEval( PhaseDiff ) ;

if nEchoes > 2
    error('TODO: field-fitting over multiple echoes');
    % approx filter (smoothing kernel) diameter in mm:
    filterDiameter = 11 ;
    
    % filter size in number of voxels 
    Params.kernelSize = round( filterDiameter ./ ImgArray{ 1, 1 }.getvoxelspacing() ) ;
    % if size is an even number of voxels along any dimension, increment up:
    evenDims = ~logical( mod( Params.kernelSize, 2 ) ) ;
    Params.kernelSize( evenDims ) = Params.kernelSize( evenDims ) + ones( [ 1 nnz( evenDims ) ] ) ; 
    Params.method     = 'gaussian' ; 

    % Unwrap later echoes using the phase difference/TE estimate from the 1st 2
    % echoes, assuming a linear model of phase evolution with TE
    % As described by Juchem C., Proc. Intl. Soc. Mag. Reson. Med. 21 (2013):
    % See recorded presentation, around 13:50:
    % https://cds.ismrm.org/protected/13MPresentations/E071/
    % https://cds.ismrm.org/protected/13MPresentations/abstracts/7298.pdf

    TE      = zeros( nEchoes, 1 ) ;
    PhaseTE = cell( nEchoes, 1 ) ;

    for iEcho = 1 : nEchoes
        TE(iEcho) = ImgArray{ iEcho, 2 }.Hdr.EchoTime ; % [units: ms]
    end
    
    % Correction for phase offsets following Sun H, et al. ISMRM 2016, Abstract # 2987:
    %
    %   Estimate phase @ TE = 0 from 1st 2 echoes and subtract this offset out from
    %   raw phase of subsequent echoes:

    % -----
    % Filter phase diff estimate before estimating offset + later phases:
    PhaseDiffInRad.filter( ImgArray{2,1}.img, Params ) ;
    
    PhaseOffset     = PhaseDiffInRad.copy() ;
    
    PhaseOffset.img = angle( ( ImgArray{ 1, 1 }.img .* exp(i*ImgArray{ 1, 2 }.img) ) .* ...
        conj( exp( i*TE(1)*(PhaseDiffInRad.img/PhaseDiffInRad.Hdr.EchoTime) ) ) ) ;
    
    PhaseOffset = PhaseOffset.unwrapphase( ImgArray{2,1}, Params ) ;

    % -----
    % lowpass filter unwrapped phase offset before correction of individual echoes:
    PhaseOffsetC = PhaseOffset.copy() ;
    PhaseOffset.filter( ImgArray{2,1}.img, Params ) ;

    % -----
    % Perform correction and unwrap phases :
    for iEcho = 1 : nEchoes
        
        img            = ImgArray{ iEcho, 1 }.img .* exp(i*ImgArray{ iEcho, 2 }.img) ;
        
        PhaseTE{iEcho} = ImgArray{ iEcho, 2}.copy() ; % copies .Hdr
        
        PhaseTE{iEcho}.img = angle( img .* conj( exp( i*PhaseOffset.img ) ) ) ;

        PhaseTE{iEcho}.Hdr.MaskingImage = PhaseDiff.Hdr.MaskingImage ;
        PhaseTE{iEcho} = PhaseTE{iEcho}.unwrapphase( ImgArray{iEcho, 1}, Params ) ;
        tmp(:,:,:,iEcho) = PhaseTE{iEcho}.img ;

    end
    
    % -----
    % Correct for temporal wraps between echoes (assuming phase offset = 0):
    PhaseCorrectedTE = cell( nEchoes, 1 ) ;
    
    for iEcho = 1 : nEchoes 

        phaseEstimate = (TE(iEcho)/PhaseDiffInRad.Hdr.EchoTime)*PhaseDiffInRad.img ;

        err            = phaseEstimate - PhaseTE{iEcho}.img ;
        % Wherever the deviation between estimate and unwrapped measurement exceeds pi,
        % correct the measurement by adding the appropriate pi-multiple:
        n              = ( abs(err) > pi ) .* round( err/pi ) ;

        PhaseCorrectedTE{iEcho}.img = PhaseTE{iEcho}.img + n*pi ;
        
        tmp2(:,:,:,iEcho) = PhaseCorrectedTE{iEcho}.img ;

    end
    
    % -----
    % magnitude sum of squares
    mag = 0;

    for iEcho = 1 : nEchoes
     mag = mag + ( ImgArray{ iEcho, 1 }.img ).^2 ;
    end

    mag = sqrt( mag ) ;  

    % if is fitting with zero phase offset:
    % -----
    % Magnitude-weighted sum of squares fit of corrected phase to echo time:
    num = 0;
    den = 0;

    for iEcho = 1 : nEchoes
     num = num + TE(iEcho)* ( ImgArray{ iEcho, 1 }.img .^2 ) .* PhaseCorrectedTE{iEcho}.img ;
     den = den + ( TE(iEcho)* ImgArray{ iEcho, 1 }.img ) .^2 ;
    end

    PhaseDiffAvg     = PhaseTE{1}.copy() ;
    PhaseDiffAvg.img = num./den ;
    PhaseDiffAvg.Hdr.EchoTime = 1 ; % [units: ms]
    PhaseDiffAvg.scalephasetofrequency( ) ; % scales to Hz

    Field2     = FieldEval( PhaseDiffAvg ) ;
    
    Nii.filename = './Field1' ;
    nii(mask.*Field.img, Nii);
    
    Nii.filename = './Field2' ;
    nii(mask.*Field2.img, Nii);
    
    Nii.filename = './Field3' ;
    nii(mask.*Field3.img, Nii);

    
    % dPhasedt = num./den ;    

    % tmp11 = mag .* dPhasedt ;
    % tmp12 = smooth3( mag, 'gaussian', Params.filterSize ) ;
    % tmp13 = smooth3( tmp11, 'gaussian', Params.filterSize ) ;
    % tmp14 = tmp13./tmp12 ;
    % nii(dPhasedt - tmp14);
    
    % Fit PhaseCorrectedTE to TE
    %
    mask = Params.mask ;

    for iEcho = 1 : nEchoes 
        mask = mask & ( PhaseTE{iEcho}.img ~=0 ) & ( ImgArray{iEcho,1}.img ~=0 );
    end

    nVoxelsVoi = nnz( mask ) ;

    nVoxelsImg = numel( PhaseTE{1}.img(:) ) ;
    indicesVoi = find( mask(:) ) ;

    % % To speed things up, construct masking/truncation operator M to limit the fit to 
    % % voxels within the reliable region of the image (i.e. region of sufficient SNR)
    Mvoi = sparse( [1:nVoxelsVoi], indicesVoi, ones([nVoxelsVoi 1]), nVoxelsVoi, nVoxelsImg ) ;
   
    % This aggrandized Mvoi is to be applied to the concatenated phase-offset + phase-change (dPh/dTE) vector 
    M = [ Mvoi sparse( size(Mvoi,1), size(Mvoi,2) ) ;
          sparse( size(Mvoi,1), size(Mvoi,2) ) Mvoi ] ;

    % % At this point, Mvoi would work for a single echo, but instead we have multiple echoes...
    %
    % % Full truncation operator M: 
    % %   1 extra row for every extra echo,
    % %   nVoxelsImg extra columns for every extra echo 
    % M  = sparse( nVoxelsVoi*nEchoes, nVoxelsImg * nEchoes ) ; 

    % for iEcho = 1 : nEchoes
    %     M( iEcho:nEchoes:(end-nEchoes+iEcho), iEcho:nEchoes:(end-nEchoes+iEcho) ) = Mvoi ;
    % end

    % phase looks like, ph:
    %
    % phase( 1st voxel: 1st echo )
    % phase( 1st voxel: 2nd echo )
    % ...
    % phase( 1st voxel: Nth echo )
    % phase( 2nd voxel: 1st echo )
    % ...
    % ...
    % phase( last voxel:last echo )
    %
    % data weights W follow the same alternating pattern but on the diagonal of a sparse matrix 
    
    ph = zeros( nVoxelsVoi*nEchoes, 1 ) ;
    W  = zeros( nVoxelsVoi*nEchoes, 1 ) ;

    for iEcho = 1 : nEchoes
        ph(iEcho:nEchoes:(end-nEchoes+iEcho)) = PhaseCorrectedTE{iEcho}.img( mask ) ;
        W(iEcho:nEchoes:(end-nEchoes+iEcho))  = ImgArray{iEcho, 1}.img( mask ) ;
        % ph(iEcho:nEchoes:(end-nEchoes+iEcho)) = PhaseCorrectedTE{iEcho}.img ;
        % W(iEcho:nEchoes:(end-nEchoes+iEcho))  = ImgArray{iEcho, 1}.img  ;
    end

    W = spdiags( W, 0, nEchoes*nVoxelsVoi, nEchoes*nVoxelsVoi ) ;
    W = W/(max(diag(W))) ; 

    % Linear operator A should be the following:
    %
    % A = sparse( nEchoes*nVoxels, 2*nVoxels ) ;
    %
    %   for iVoxel = 1 : nVoxels 
    %         iCol = iVoxel*2 - 1 ;
    %       iRow = nEchoes*(iVoxel-1) + 1 ;
    %       A( iRow:(iRow+nEchoes-1), iCol:(iCol+1) ) = A0 ;
    %   end
    %
    % ... but that takes too long.
    % Here is a shortcut:

    A0 = sparse([ ones( nEchoes, 1 ) TE ]) ;
    A  = A0 ;

    while size( A, 2 ) < 2*nVoxelsVoi

        A  = [ A  sparse( size(A,1), size(A,2) ) ; 
               sparse( size(A,1), size(A,2) ) A ] ;
    end

    % crop to correct size
    A = A(1:nVoxelsVoi*nEchoes, 1:2*nVoxelsVoi) ;
    
    % x = [M0; M1]*x, rearranges the values of x:
    % rather than alternating elements of x being phase offset (odd indices)
    % and phase change (even indices), the top half of x will now be the offset
    % and the bottom will be the phase change
    % (making it easier to apply specific regularization penalties to each phase component)
    M0 = sparse([ 1 0 ]) ;
    M1 = sparse([ 0 1 ]) ;

    while size( M0, 2 ) < 2*nVoxelsImg

        M0  = [ M0  sparse( size(M0,1), size(M0,2) ) ; 
               sparse( size(M0,1), size(M0,2) ) M0 ] ;
        M1  = [ M1  sparse( size(M1,1), size(M1,2) ) ; 
               sparse( size(M1,1), size(M1,2) ) M1 ] ;
    end

    % crop to correct size
    M0 = M0(1:nVoxelsImg, 1:2*nVoxelsImg) ;
    M1 = M1(1:nVoxelsImg, 1:2*nVoxelsImg) ;
    
    % Params for conjugate-gradient optimization
    CgParams.tolerance     = 1E-4 ;
    CgParams.maxIterations = 500 ;    
    CgParams.isDisplayingProgress = true ;    
    
    tic
    x = cgls( A'*A, ... % least squares operator
              A'*ph, ... % effective solution vector
              zeros( [2*nVoxelsVoi 1] ) , ... % % initial 'guess' solution vector
           CgParams ) ;
          
    toc
    
    phaseOffsetUw = zeros( ImgArray{1,1}.getgridsize() ) ;
    phaseOffsetUw( mask ) = x(1:2:end-1) ;  % [units: rad]
    % phaseOffsetUw( mask ) = M0*x ;  % [units: rad]

    dPhdtUw = zeros( ImgArray{1,1}.getgridsize() ) ;
    dPhdtUw( mask ) = x(2:2:end) ; % [units: rad/ms]
    % dPhdtUw( mask ) = M1*x ; % [units: rad/ms]

    % TODO
    %   Using truncation mask and/or weighting matrix takes way too long :(
    %   optimize somehow?
    %
    % e.g use unweighted solution as starting guess for weighted:
    %

    [Dx,Dy,Dz] = createdifferenceoperators( ImgArray{1,1}.getgridsize(), ImgArray{1,1}.getvoxelspacing, 2 ) ;

    L = Dx + Dy + Dz ;
    
    T = sparse( [ [1:2:2*nVoxelsVoi-1] [2:2:2*nVoxelsVoi] ], ...
                [ [1:nVoxelsVoi] [(nVoxelsVoi+1):(2*nVoxelsVoi) ] ], ...
                ones( 2*nVoxelsVoi, 1 ) ) ;

    % Augmented system:
    b = [ W*ph; zeros(2*nVoxelsVoi, 1) ] ;
    % b = [ ph; zeros(2*nVoxelsVoi, 1) ] ;
    % b = [ W*ph; ] ;
    % b = [ ph; ] ;
    
    mu0 = 1E-5 ;
    mu1 = 1E-5 ; 
    % penalty/regularizationg operator
    G = [ mu0*Mvoi*speye( nVoxelsImg, nVoxelsImg ), sparse( nVoxelsVoi, nVoxelsImg ) ; 
          sparse( nVoxelsVoi, nVoxelsImg ), mu1*Mvoi*L ] ;

    % A0 = W*A*T*M*[ M0; M1 ];
    % A0 = W*A*M ;
    
    A0 = W*A*T*M*[ M0; M1 ];
    A1 = [ A0 ; G*[ M0; M1 ] ] ;
    % A1 = A0 ;
    tic
    x = cgls( A1'*A1, ... % least squares operator
              A1'*b, ... % effective solution vector
              zeros( [2*nVoxelsImg 1] ), ...  %x, ...%... % initial 'guess' solution vector
              CgParams  ) ;     
    toc
  
    xx = [M0;M1]*x ;
    phaseOffset = reshape( xx(1:(length(xx)/2)) , ImgArray{1,1}.getgridsize() ) ;;  % [units: rad]
    dPhdt       = reshape( xx((length(xx)/2)+1:end), ImgArray{1,1}.getgridsize() ) ; % [units: rad/ms]


    % tic
    % x = cgls( A'*W'*W*A, ... % least squares operator
    %           A'*W'*W*ph, ... % effective solution vector
    %           x, ...  %x, ...%... % initial 'guess' solution vector
    %           CgParams  ) ;     
    % toc
  


    %
    % tic
    % x = cgls( A'*W'*W*A, ... % least squares operator
    %           A'*W'*W*ph, ... % effective solution vector
    %           zeros( [2*nVoxelsVoi 1] ), ...  %x, ...%... % initial 'guess' solution vector
    %           CgParams  ) ;     
    % toc
    %
    % tic
    % x = cgls( A'*W'*M'*M*W*A, ... % least squares operator
    %           A'*W'*M'*M*W*ph, ... % effective solution vector
    %           x, ...
    %           CgParams  ); %... % initial 'guess' solution vector
    % toc
    
    % phaseOffset = zeros( ImgArray{1,1}.getgridsize() ) ;
    % phaseOffset = reshape( x(1:2:end-1) , ImgArray{1,1}.getgridsize() ) ;;  % [units: rad]
    % dPhdt = zeros( ImgArray{1,1}.getgridsize() ) ;
    % dPhdt = reshape( x(2:2:end), ImgArray{1,1}.getgridsize() ) ; % [units: rad/ms]
    xx = [M0;M1]*x ;
    phaseOffset = reshape( xx(1:(length(xx)/2)) , ImgArray{1,1}.getgridsize() ) ;;  % [units: rad]
    dPhdt       = reshape( xx((length(xx)/2)+1:end), ImgArray{1,1}.getgridsize() ) ; % [units: rad/ms]

    PhaseDiffFit     = PhaseTE{1}.copy() ;
    % PhaseDiffFit.img = dPhasedt - tmp14 ;
    % PhaseDiffFit.img = dPhdt ;
    PhaseDiffFit.Hdr.EchoTime = 1 ; % [units: ms]
    PhaseDiffFit.scalephasetofrequency( ) ; % scales PhaseDiffFit.img to Hz

    Field     = FieldEval( PhaseDiffFit ) ;

end

Field.Hdr.SeriesDescription = [ 'B0Field_measured_' Field.Hdr.SeriesDescription  ] ;

end
% =========================================================================
function [Field] = modelfield( Fields, Params )
% MODELFIELD
%
% Map respiration-induced resonance offset (RIRO) assuming a linear model
% of field variation w/breath-amplitude (i.e. Topfer et al. MRM 2018)
%
% [Field] = MODELFIELD( Fields ) 
% [Field] = MODELFIELD( Fields, Params ) 
%
% Fields is a linear cell array, e.g.
%
%   Fields{1} = FieldInspired;
%   Fields{2} = FieldExpired;
%       where each Fields entry is a FieldEval-type object


% Optional input
%   pDc : DC auxiliary pressure reading corresponding to the mean respiratory
%       state 
%       [default: (FieldInspired.Aux.Data.p + FieldExpired.Aux.Data.p ) /2 ]
%
% Returns FieldEval-type objects
%   
%   Riro : Respiration induced resonance-offset (i.e. Field shift from respiration (FieldInspired - FieldExpired) )
%   Field : The constant (DC) field
%
% If input both fields have Field.Aux.Data.p defined as single scalar
% values (e.g. associated pressure measurements: pIn, pEx) then 
%
%    Field.img        = ( pIn*FieldExpired.img - pEx*FieldInspired.img )/(pIn - pEx) ;
%
% else
%    pIn = 1 ; 
%    pEx = -1 ;
%    such that Field.img as defined ^ becomes the avg. of the 2 input Fields

DEFAULT_MAXABSFIELD        = 600 ;
DEFAULT_MAXFIELDDIFFERENCE = 150 ;

if nargin < 1
    error('Function requires at least 1 input (e.g. linear cell array of FieldEval-type objects)') ;
elseif nargin == 1
    Params.dummy = [] ;
end

if ~myisfield( Params, 'maxAbsField' ) || isempty( Params.maxAbsField ) 
    Params.maxAbsField = DEFAULT_MAXABSFIELD ;
end

if ~myisfield( Params, 'maxFieldDifference' ) || isempty( Params.maxFieldDifference ) 
    Params.maxFieldDifference = DEFAULT_MAXFIELDDIFFERENCE ;
end

if iscell( Fields )
    % if difference in Larmor frequency >= 1 Hz, issue an error:
    assert( abs( 1000*double( Fields{1}.Hdr.ImagingFrequency - Fields{2}.Hdr.ImagingFrequency ) ) < 1.0, ... 
        'Expected the 2 given field maps to have been acquired at the same Larmor frequency. Unimplemented feature.' )

    if ~isempty( Fields{1}.Aux.Tracker.Data.p ) 
        assert( isscalar( Fields{1}.Aux.Tracker.Data.p), ...
            'Expected single scalar value for Fields{1}.Aux.Data.p' ) ;

        pIn = Fields{1}.Aux.Tracker.Data.p ;
    else 
        pIn = 1 ;
    end

    if ~isempty( Fields{2}.Aux.Tracker.Data.p ) 
        assert( isscalar( Fields{2}.Aux.Tracker.Data.p), ...
            'Expected single scalar value for Fields{2}.Aux.Data.p' ) ;

        pEx = Fields{2}.Aux.Tracker.Data.p ;
    else 
        pEx = -1 ;
    end

    if ~myisfield( Params, 'pBreathing' ) || isempty( Params.pBreathing ) 
        Params.pBreathing = [ pIn pEx ] ;
    end

    pDc            = median( Params.pBreathing ) ;
    pShiftTraining = pIn - pEx ;

    Field     = Fields{1}.copy() ; 
    Riro      = Fields{1}.copy() ; 
    FieldZero = Fields{1}.copy() ; 

    FieldZero.img                = ( pIn*Fields{2}.img - pEx*Fields{1}.img )/(pShiftTraining) ;
    FieldZero.Aux.Tracker.Data.p = 0 ;
    Field.Model.Zero             = FieldZero ;
    % no way of knowing what values might be reasonable for this 'Zero' field, so
    % there is no call to .getvaliditymask()

    Riro.img                = Fields{1}.img - Fields{2}.img ;
    Riro.Hdr.MaskingImage   = Riro.getvaliditymask( Params.maxFieldDifference ) ;
    % The training fields themselves (e.g. inspired/expired breath-holds) might not
    % be representative of the typical field shift if it's defined for any
    % given phase of the respiratory cycle as the deviation from the expected
    % (mean/DC) value.
    %
    % Scale the shift by the RMS pressure deviation observed during regular breathing:
    pShiftRms = rms( Params.pBreathing - pDc ) ;

    Riro.img  = ( pShiftRms/pShiftTraining ) * Riro.img ;
    Riro.Aux.Tracker.Data.p = pShiftRms ;

    Field.Model.Riro = Riro ;

    Field.img                = pDc*(Field.Model.Riro.img/pShiftRms) + Field.Model.Zero.img ;
    Field.Aux.Tracker.Data.p = pDc ;

    Field.Hdr.MaskingImage   = Field.getvaliditymask( Params.maxAbsField ) ;

elseif isa( Fields, 'FieldEval' )
    % fit Fields.img time series (along 4th dimension) to tracker time series
    nAcq    = size( Fields.img, 4 ) ;
    assert( nAcq == length( Fields.Aux.Tracker.Data.p ), 'Expected a single tracker measurement for each image' )

    mask = ( sum( Fields.Hdr.MaskingImage, 4 ) == nAcq ) ;
    
    nVoxels = prod( Fields.getgridsize() ) ;
    nVoxelsVoi = nnz(mask) ;
    % I       = speye( nVoxels, nVoxels ) ;
    I       = speye( nVoxelsVoi, nVoxelsVoi ) ;
    
    % linear operator: A
    A       = [ I Fields.Aux.Tracker.Data.p(1)*I ] ;
    
    % solution vector: Bt
    Bt  = Fields.img( :, :, :, 1) ;
    Bt  = Bt( mask ) ;
    % Bt      = squeeze( Fields.img( :, :, :, 1) ) ;
    % Bt      = Bt(:) ;

    for iT = 2 : nAcq 
        disp( [ num2str(100*iT/nAcq, 3) '%'])
        
        A    = [ A ; I Fields.Aux.Tracker.Data.p(iT)*I ] ;
        
        tmpB = Fields.img( :, :, :, iT) ;
        tmpB = tmpB( mask ) ;
        % tmpB = squeeze( Fields.img( :, :, :, iT) ) ;
        Bt   = [ Bt ; tmpB(:) ] ;
    end
    Bt(isnan(Bt)) = 0;
    
    % isRegularizing = false ; % possible TODO: add option?
    %     % lamda(1): regularization parameter for the static field component
    %     % lamda(2): " " for the RIRO component
    %     lamda = [0 0];
    %
    % if isRegularizing
    % A0 = A;
    %     % construct masking/truncation operator M to limit the fit to 
    %     % voxels within the reliable region of the image (i.e. region of sufficient SNR)
    %     nVoxelsVoi = nnz( Fields.Hdr.MaskingImage(:,:,:,1) ) ;
    %     indicesVoi = find( Fields.Hdr.MaskingImage(:,:,:,1)~= 0 ) ;
    %     Mvoi = sparse( [1:nVoxelsVoi], indicesVoi, ones([nVoxelsVoi 1]), nVoxelsVoi, nVoxels ) ;
    %     
    %     [Dx, Dy, Dz] = createdifferenceoperators( Fields.getgridsize(), Fields.getvoxelspacing(), 2) ;
    %
    %     L = Dx + Dy + Dz ;
    %
    %     R  = [lamda(1)*Mvoi*L sparse( size(Mvoi*L,1),size(Mvoi*L,2) ) ; 
    %           sparse( size(Mvoi*L,1),size(Mvoi*L,2) ) lamda(2)*Mvoi*L ] ; 
    %
    %     A  = [A; R] ;
    %     Bt = [Bt ; zeros(2*nVoxelsVoi, 1 ) ] ;
    % end

    % Params for conjugate-gradient optimization
    CgParams.tolerance     = 1E-4 ;
    CgParams.maxIterations = 500 ;    
    CgParams.isDisplayingProgress = true ;    

    x = cgls( A'*A, ... % least squares operator
              A'*Bt, ... % effective solution vector
              zeros( [2*nVoxelsVoi 1] ) , ... % % initial 'guess' solution vector
              CgParams ) ;
    
    Field                    = Fields.copy() ; 
    Field.img                = Field.img(:,:,:,1) ;
    Field.img( mask )        = x(1:nVoxelsVoi) ;
    Field.img( ~mask )       = 0 ;
    Field.Hdr.MaskingImage   = mask ; 
    Field.Hdr.MaskingImage   = Field.getvaliditymask( Params.maxAbsField ) ;
    Field.Aux.Tracker.Data.p = 0 ; 

    Riro                    = Fields.copy() ; 
    Riro.img                = Riro.img(:,:,:,1) ;
    % scale RIRO by RMSE of physio signal
    pShiftRms               = rms( Fields.Aux.Tracker.Data.p - mean(Fields.Aux.Tracker.Data.p) ) ;
    Riro.img( mask )        = pShiftRms .* x(nVoxelsVoi+1:end) ;
    Riro.img( ~mask )       = 0 ;
    Riro.Hdr.MaskingImage   = mask ; 
    Riro.Hdr.MaskingImage   = Riro.getvaliditymask( Params.maxFieldDifference ) ;
    Riro.Aux.Tracker.Data.p = pShiftRms ; 
    Field.Model.Riro        = Riro ;

end

end
% =========================================================================

end
% =========================================================================
% =========================================================================

end
