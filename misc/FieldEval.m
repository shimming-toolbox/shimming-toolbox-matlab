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
% See MaRdI documentation for additional features.
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

Params.dummy = [] ;

if nargin ~= 0 
    
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

if ~isempty( Img.Aux ) 
    if isa( Img.Aux, 'ProbeTracking' ) 
        ImgCopy.Aux = Img.Aux.copy() ;
    else
        error('Expected Img.Aux to be of type ProbeTracking') ;
    end
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
% maxAbsField 
%   maximum absolute voxel value assumed to represent an accurate field
%   measurement. Voxels with abs-values greater than this might stem from
%   errors in the unwrapping.  [default: 500 Hz]
%   (Set to Inf to ignore the criterion)

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

    dfdp = Field.Model.Riro.img()/Field.Model.Riro.Aux.Data.p ;
    dp   = Field.Aux.debias( p ) ; 

    Riro.img = dfdp*dp ;
    Riro.Aux.Data.p = p ;
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

%% -----
% parse inputs
assert( nargin >= 2, 'Not enough input arguments.' ) ;

for iImg = 1 : 2
    if isa( varargin{iImg}, 'MaRdI' )
        Img = varargin{iImg} ;
    elseif ischar( varargin{iImg} )
        display( ['Loading img ' num2str(iImg) ' of 2'] );
        Img = MaRdI( varargin{iImg} ) ;
    else 
        error('Invalid input') ;
    end

    if Img.ismagnitude() 
        Mag = Img.copy() ;
    elseif Img.isphase() 
        Phase = Img.copy() ;
    else
        error('Unexpected image type: Neither magnitude nor phase?') ;
    end
end

if ~exist('Phase') 
    error('Input did not include Phase image (MaRdI object). See help FieldEval.mapfield' ) 
elseif ~exist('Mag') 
    error('Input did not include Magnitude image (MaRdI object). See help FieldEval.mapfield' ) 
end

if isstruct( varargin{end} )
    Params = varargin{end} ;
else
    Params.dummy = [];
end

%% -----
% assign defaults to unassigned parameters
DEFAULTS.threshold = 0.01 ;

Params = assignifempty( Params, DEFAULTS ) ;
    
% -------
% define spatial support for unwrapping
if ~myisfieldfilled( Params, 'mask' ) 
    Params.mask = Mag.getreliabilitymask( Params.threshold ) ;
    Params.mask = prod( Params.mask, 4 ) ;
end

Params.mask = logical( Params.mask ) ;

nEchoes       = length( Phase.getechotime() ) ;
nMeasurements = size( Phase.img, 5 ) ;

if nEchoes == 1 
    % Input phase should be a phase *difference* image; 
    % alternatively, PhaseDiff will refer to the difference between t=TE and t=0
    if( numel( Phase.Hdr.MrProt.alTE ) < 2 )
        warning( ['Expected Siemens phase difference image with at least 2 TEs specified in the DICOM header.' ...
           'Output field estimate may retain an offset term'] )
    end

    PhaseDiff = Phase.copy() ;
    PhaseDiff.Hdr.EchoTime = Phase.getechotime() ;

elseif nEchoes == 2
    % -------
    % phase difference image via complex division of first 2 echoes
    PhaseDiff      = Phase.copy() ;

    img            = Mag.img(:,:,:,1,:) .* exp( i*Phase.img(:,:,:,1,:) ) ;
    img(:,:,:,2,:) = Mag.img(:,:,:,2,:) .* exp( i*Phase.img(:,:,:,2,:) ) ;

    PhaseDiff.img  = angle( img(:,:,:,2,:) .* conj(img(:,:,:,1,:) ) ) ;

    PhaseDiff.Hdr.EchoTime = Phase.getechotime(2) - Phase.getechotime(1) ; % [units : ms]
else
    error('nEchoes > 2 not supported. TODO') ;
end

PhaseDiff.setmaskingimage( logical( Params.mask ) & ~isnan( PhaseDiff.img ) ) ;

PhaseDiff.unwrapphase( Mag, Params ) ;

%NOTE: Don't call PhaseDiff.scalephasetofrequency since it will use the wrong EchoTime in the multi-echo case...
PhaseDiff.img  = PhaseDiff.img/(2*pi*PhaseDiff.Hdr.EchoTime/1000) ;  
PhaseDiff.Hdr.PixelComponentPhysicalUnits = '0005H' ; % i.e. Hz

Field = FieldEval( PhaseDiff ) ;

Field.Hdr.SeriesDescription = [ 'B0Field_measured_' Field.Hdr.SeriesDescription  ] ;

end
% =========================================================================
function [Field] = modelfield( Fields, Params )
%MODELFIELD     Fit B0 maps to auxiliary respiratory recording
%
% .......
% 
% Description 
%
% Maps static B0 and respiration-induced resonance offset (RIRO) assuming a
% linear model of field variation w/breath-amplitude (Ref: Topfer et al. MRM 2018)
%
% [FieldFit] = MODELFIELD( Fields ) 
% [FieldFit] = MODELFIELD( Fields, Params ) 
% 
% .......
% 
% Usage
%
% Returns FieldEval-type object FieldFit, where FieldFit.img is the respiration-independent static field estimate,
% and FieldFit.Model.Riro.img is the respiration-dependent component.
% 
% Inputs:
%
% Fields:
%
%   Case 1: Fields pertains to a field map time-series:
%   
%       Fields should be a single object of type FieldEval, with Fields.Aux containing
%       the corresponding respiratory (e.g. bellows) recording in the form of a
%       ProbeTracking object.
%
%   Case 2: Fields pertains to separate inspired and expired field maps: 
%
%       Fields should be a cell array, with Fields{1} and Fields{2}
%       respectively containing the FieldEval objects corresponding to the
%       'inspired' and 'expired' acquisitions.
%
% Params: 
%   Optional parameters struct can possess any of the following entries:
%
%   .maxAbsField [default = 600]
%       maximum absolute field value a voxel can possess in units of Hz to be deemed reliable
%
%   .maxAbsFieldDifference [default = 150]
%       maximum absolute field value a voxel can possess in units of Hz to be deemed reliable
%
%   pDc : DC auxiliary pressure reading corresponding to the mean respiratory state 
%   
%       
%       [default in Case 2: (FieldInspired.Aux.Data.p + FieldExpired.Aux.Data.p ) /2 ]



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

DEFAULTS.maxAbsField        = 600 ;
DEFAULTS.maxFieldDifference = 150 ;

if nargin < 1
    error('Function requires at least 1 input (e.g. linear cell array of FieldEval-type objects)') ;
elseif nargin == 1
    Params.dummy = [] ;
end

Params = assignifempty( Params, DEFAULTS ) ;

if isa( Fields, 'FieldEval' )
    % fit Fields.img time series (along 5th dimension) to tracker time series
    
    nAcq = size( Fields.img, 5 ) ;
    
    if ~myisfield( Fields.Aux, 'Data' ) || ~myisfieldfilled( Fields.Aux.Data, 'p' ) 
        error('Nothing to fit: Fields.Aux was empty but should contain a valid Auxiliary recording (ProbeTracking object).') ;
    elseif nAcq ~= length( Fields.Aux.Data.p )
        error( 'Expected a single Aux value for each image. See HELP MaRdI.associateaux()' )
    else
        disp( ['Modeling B0 field (fitting field time-series to Aux recording)'] )
    end

    mask = ( sum( Fields.Hdr.MaskingImage, 5 ) == nAcq ) ;
    
    nVoxels    = prod( Fields.getgridsize() ) ;
    nVoxelsVoi = nnz(mask) ;
    I          = speye( nVoxelsVoi, nVoxelsVoi ) ;
   
    p = Fields.Aux.Data.p - mean( Fields.Aux.Data.p ) ;
    pShiftRms = rms( p ) ;
    
    % linear operator: A
    A       = [ I p(1)*I ] ;
    
    % solution vector: Bt
    Bt  = Fields.img( :, :, :, 1, 1) ;
    Bt  = Bt( mask ) ;

    disp( ['Preparing linear fitting operator...'] )
    
    for iT = 2 : nAcq 
        disp( [ num2str(100*iT/nAcq, 3) '%'])
        
        A    = [ A ; I p(iT)*I ] ;
        
        tmpB = Fields.img( :, :, :, 1, iT) ;
        tmpB = tmpB( mask ) ;
        Bt   = [ Bt ; tmpB(:) ] ;
    end
    
    Bt(isnan(Bt)) = 0;

    % possible TODO: add fit-regulatization option?
    %
    % e.g.
    % isRegularizing = false ; 
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
    
    disp( ['Performing fit...'] )

    % Params for conjugate-gradient optimization
    CgParams.tolerance     = 1E-4 ;
    CgParams.maxIterations = 500 ;    
    CgParams.isDisplayingProgress = true ;    

    x = cgls( A'*A, ... % least squares operator
              A'*Bt, ... % effective solution vector
              zeros( [2*nVoxelsVoi 1] ) , ... % initial 'guess' solution vector
              CgParams ) ;
    
    Field                    = Fields.copy() ; 
    Field.img                = Field.img(:,:,:,1,1) ;
    Riro                     = Fields.copy() ; 
    Riro.img                 = Riro.img(:,:,:,1,1) ;
    
    Field.img( mask )        = x(1:nVoxelsVoi) ;
    Field.img( ~mask )       = 0 ;
    Field.Hdr.MaskingImage   = mask ; 
    Field.Hdr.MaskingImage   = Field.getvaliditymask( Params.maxAbsField ) ;
    Field.Aux.Data.p         = 0 ; 
    
    % scale RIRO by RMSE of physio signal
    Riro.img( mask )        = pShiftRms .* x(nVoxelsVoi+1:end) ;
    Riro.img( ~mask )       = 0 ;
    Riro.Hdr.MaskingImage   = mask ; % not redundant: mask used in following call to .getvaliditymask()
    Riro.Hdr.MaskingImage   = Riro.getvaliditymask( Params.maxFieldDifference ) ;
    Riro.Aux.Data.p         = pShiftRms ; 
    Riro.Aux.Specs.limits   = [ min(Fields.Aux.Data.p) max(Fields.Aux.Data.p) ] ;

    Field.Model.Riro        = Riro ;

    return ;

elseif iscell( Fields ) % Breath-hold case:
    % if difference in Larmor frequency >= 1 Hz, issue an error:
    assert( abs( 1000*double( Fields{1}.Hdr.ImagingFrequency - Fields{2}.Hdr.ImagingFrequency ) ) < 1.0, ... 
        'Expected the 2 given field maps to have been acquired at the same Larmor frequency. Unimplemented feature.' )

    if ~isempty( Fields{1}.Aux.Data.p ) 
        assert( isscalar( Fields{1}.Aux.Data.p), 'Expected single scalar value for Fields{1}.Aux.Data.p' ) ;

        pIn = Fields{1}.Aux.Data.p ;
    else 
        pIn = 1 ;
    end

    if ~isempty( Fields{2}.Aux.Data.p ) 
        assert( isscalar( Fields{2}.Aux.Data.p), 'Expected single scalar value for Fields{2}.Aux.Data.p' ) ;

        pEx = Fields{2}.Aux.Data.p ;
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

    FieldZero.img        = ( pIn*Fields{2}.img - pEx*Fields{1}.img )/(pShiftTraining) ;
    FieldZero.Aux.Data.p = 0 ;
    Field.Model.Zero     = FieldZero ;
    % no way of knowing what values might be reasonable for this 'Zero' field, so
    % there is no call to .getvaliditymask()

    Riro.img              = Fields{1}.img - Fields{2}.img ;
    Riro.Hdr.MaskingImage = Riro.getvaliditymask( Params.maxFieldDifference ) ;
    
    % The training fields themselves (e.g. inspired/expired breath-holds) might not
    % be representative of the typical field shift if it's defined for any
    % given phase of the respiratory cycle as the deviation from the expected
    % (mean/DC) value.
    %
    % Hence, scale the shift by the RMS pressure deviation observed during regular breathing:
    pShiftRms        = rms( Params.pBreathing - pDc ) ;
    Riro.img         = ( pShiftRms/pShiftTraining ) * Riro.img ;
    Riro.Aux.Data.p  = pShiftRms ;

    Field.Model.Riro = Riro ;

    Field.img        = pDc*(Field.Model.Riro.img/pShiftRms) + Field.Model.Zero.img ;
    Field.Aux.Data.p = pDc ;

    Field.Hdr.MaskingImage = Field.getvaliditymask( Params.maxAbsField ) ;

end

end
% =========================================================================

end
% =========================================================================
% =========================================================================

end
