classdef MaRdI < matlab.mixin.SetGet 
%MaRdI Ma(t)-R-dI(com)
%
% .......
%   
% Description
%
%   Dicom into Matlab for Siemens MRI data
%
% .......
%   
% Usage
%
% Img = MaRdI( dataLoadDirectory )
% Img = MaRdI( dataLoadDirectory, Params )
%
%   Img contains fields
%
%       .img    (array of images - 3D if there are multiple DICOMs in directory, 
%           4D for time series data (WIP))
%
%       .Hdr    (header of the first DICOM file read by dir( dataLoadDirectory ) ) 
%
% =========================================================================
% Updated::20160921::ryan.topfer@polymtl.ca
% =========================================================================

% =========================================================================
% =========================================================================
properties
    img ;
    Hdr ;
end

% =========================================================================
% =========================================================================    
methods
% =========================================================================    
function Img = MaRdI( dataLoadDirectory )

Img.img = [] ;
Img.Hdr = [] ;
    
if nargin == 1 && ~isempty(dataLoadDirectory)

    listOfDicoms = dir( [ dataLoadDirectory '/*.dcm'] );

    if length(listOfDicoms) == 0
        % try .IMA
        isImaDir = MaRdI.changefilesuffix( dataLoadDirectory ) ;
        
        if isImaDir == true
            display('Converted .IMA to .dcm');
        else
            error('No .dcm or .IMA files found in given directory') ;
        end

        listOfDicoms = dir( [ dataLoadDirectory '/*.dcm'] );
    end

    Img.Hdr      = dicominfo( [ dataLoadDirectory '/' listOfDicoms(1).name] ) ;
    
    Img.img      = double( dicomread( [ dataLoadDirectory '/' listOfDicoms(1).name] ) )  ;

    Img.Hdr.NumberOfSlices = uint16( length(listOfDicoms) ) ;
    
    for sliceIndex = 2 : Img.Hdr.NumberOfSlices
        Img.img(:,:,sliceIndex) = dicomread([ dataLoadDirectory '/' listOfDicoms(sliceIndex).name]) ;
    end

    Params = [] ;

    if ~isempty( strfind( Img.Hdr.ImageType, 'MOSAIC' ) ) 
        
        Params.isNormalizingMagnitude = false ;        
        Img = reshapemosaic( Img ) ;

    end

    Img = Img.scaleimgtophysical( Params ) ;   

    if ~myisfield( Img.Hdr, 'SpacingBetweenSlices' ) 
        Img.Hdr.SpacingBetweenSlices = Img.Hdr.SliceThickness ;
    end
end

end
% =========================================================================
% *** TODO
% 
% ..... 
% MaRdI( ) [constructor] 
%
%   Option to call Img = MaRdI( img, Hdr ) 
%   such that an object can be instantiated using an array of doubles and 
%   an appropriate dicom Hdr.  
% ..... 
% CROPIMG
%   make compatible for odd-sized arrays
% ..... 
% RESLICEIMG()
%   griddata takes too f-ing long.
%   write interp function in cpp
%   Clean up 'resliceimg()'
%   
% ..... 
% Changing .ImageType:
%   Invalid replacement occurs in a few places.
%   eg. in scaleimgtophysical()
%   Original entry is replaced by
%    Img.Hdr.ImageType  = 'ORIGINAL\SECONDARY\P\' ; 
%   -> Actually, things like DIS2D may have been applied and that appears after
%   the \P\ (or \M\). 
%   Solution: replace 'PRIMARY' with 'SECONDARY' and leave the rest.
%
% ..... 
% Saving FIELD as DICOM
%      -> (do not save image type as phase)
%      -> must save proper ImagePositionPatient info
% 
% ..... 
% Write static comparegrid() function or something to ensure operations involving
% multiple images are valid (e.g. voxel positions are the same) 
%
% ..... 
% GETDIRECTIONCOSINES()
%   Weird flip-conditional for slice direction? Is this OK??
%   Siemens private header apprently contains field SliceNormalVector.
%
% ..... 
% WRITE()
%   Neither 'as dcm' nor 'as nii' functionalities work properly.
% 
% ..... 
% EXTRACTHARMONICFIELD()
%  Clean up 
%     
% =========================================================================

function ImgCopy = copy(Img)
%COPY 
% 
% Make a copy of a MaRdI (i.e. handle) object.
% 
% ImgCopy = Copy( Img ) ;

ImgCopy     = MaRdI() ;
ImgCopy.img = Img.img;
ImgCopy.Hdr = Img.Hdr ;

% from https://www.mathworks.com/matlabcentral/newsreader/view_thread/257925
% % Instantiate new object of the same class.
% ImgCopy = feval(class(Img));
% Copy all non-hidden properties.
% ImgProperties = properties(Img)
% for iProperty = 1 : length(ImgProperties)
%     ImgCopy.(ImgProperties{i}) = Img.(ImgProperties{i});
% end
end
% =========================================================================
function Img3 = minus(Img1, Img2)
%MINUS 
% 
% Img3 = MINUS( Img1, Img2 ) 
%
% Returns MaRdI-type object Img3
% where 
% .Hdr is a copy of Img1.Hdr 
% and 
% Img3.img = Img1.img - Img2.img 

Img3 = Img1.copy( ) ;

Img3.img = Img1.img - Img2.img ;

end
% =========================================================================
function Img3 = plus(Img1, Img2)
%PLUS 
% 
% Img3 = PLUS( Img1, Img2 ) 
%
% Returns MaRdI-type object Img3
% where 
% .Hdr is a copy of Img1.Hdr 
% and 
% Img3.img = Img1.img + Img2.img 

Img3 = Img1.copy( ) ;

Img3.img = Img1.img - Img2.img ;

end
% =========================================================================
function Img = scaleimgtophysical( Img, Params )
%SCALEIMGTOPHYSICAL
%
% Img = SCALEIMGTOPHYSICAL( Img ) 
% Img = SCALEIMGTOPHYSICAL( Img, Params ) 

DEFAULT_ISUNDOING              = false ;

if (nargin < 2) || isempty( Params ) 
    Params.dummy = [] ;
end 

if  ~myisfield( Params, 'isUndoing' ) || isempty( Params.isUndoing )
    Params.isUndoing = DEFAULT_ISUNDOING ;
end


if ~Params.isUndoing
    
    if ~isempty( strfind( Img.Hdr.ImageType, '\M\' ) ) % is raw SIEMENS mag
        
        DEFAULT_ISNORMALIZINGMAGNITUDE = true ;
        
        if  ~myisfield( Params, 'isNormalizingMagnitude' ) || isempty( Params.isNormalizingMagnitude )
            Params.isNormalizingMagnitude = DEFAULT_ISNORMALIZINGMAGNITUDE ;
        end
            
        if Params.isNormalizingMagnitude
            Img = Img.normalizeimg() ;
        end

        Img.Hdr.ImageType  = 'ORIGINAL\SECONDARY\M\' ; 
        

    elseif ~isempty( strfind( Img.Hdr.ImageType, '\P\' ) ) % is raw SIEMENS phase
    
        DEFAULT_RESCALEINTERCEPT = -(2^12) ;
        DEFAULT_RESCALESLOPE     = 2 ;
        
        if  ~myisfield( Img.Hdr, 'RescaleIntercept' ) || isempty( Img.Hdr.RescaleIntercept )
            Img.Hdr.RescaleIntercept = DEFAULT_RESCALEINTERCEPT ;
        end

        if ~myisfield( Img.Hdr, 'RescaleSlope' ) || isempty( Img.Hdr.RescaleSlope )
            Img.Hdr.RescaleSlope = DEFAULT_RESCALESLOPE ;
        end

        Img.img = (Img.Hdr.RescaleSlope .* double( Img.img ) ...
                    + Img.Hdr.RescaleIntercept)/double(2^Img.Hdr.BitsStored) ;
    

        Img.img            = pi*Img.img ; % scale to rad
        Img.Hdr.ImageType  = 'ORIGINAL\SECONDARY\P\' ; 
        Img.Hdr.PixelRepresentation = uint8(1) ; % i.e. signed 
    
    else 
        error('Unknown Hdr.ImageType or invalid .Hdr') ;
    end

    Img.Hdr.PixelComponentPhysicalUnits = '0000H' ; % i.e. none

elseif Params.isUndoing
    
    Img.Hdr.RescaleIntercept = min( Img.img(:) ) ;
    Img.img                  = Img.img - Img.Hdr.RescaleIntercept ;
    Img.Hdr.RescaleSlope     = max( Img.img(:) )/double(2^Img.Hdr.BitsStored ) ;
    Img.img                  = uint16( Img.img / Img.Hdr.RescaleSlope ) ;

end

end
% =========================================================================
function Img = normalizeimg( Img )
%NORMALIZEIMG
%
% Img = NORMALIZEIMG( Img )

Img.img = double(Img.img)/max(Img.img(:)) ; % normalize 

end
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
Img.Hdr.NumberOfSlices       = size(Img.img, 3) ;


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
Img.Hdr.NumberOfSlices       = size(Img.img, 3) ;

if ~strcmp( padDirection, 'post' )
% update image position 
% (i.e. location in DICOM RCS of 1st voxel in data array (.img))
    
    voxelSize = Img.getvoxelsize() ;

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
function Img = scalephasetofrequency( Img, undoFlag )
%SCALEPHASETOFREQUENCY
%
% Converts unwrapped phase [units:rad] to field [units: Hz]
% 
%   Field = scalephasetofrequency( UnwrappedPhase )
%
% .....
%
% UnwrappedPhase.Hdr.EchoTime              [units: ms]

scalingFactor = 1/( 2*pi*Img.Hdr.EchoTime*(1E-3)  ) ;

if (nargin < 2) || (undoFlag ~= -1)
    assert( strcmp( Img.Hdr.PixelComponentPhysicalUnits, '0000H' ) )

    Img.img       = scalingFactor * Img.img ;
    Img.Hdr.PixelComponentPhysicalUnits = '0005H' ; % i.e. Hz

elseif (undoFlag == -1)
    assert( strcmp( Img.Hdr.PixelComponentPhysicalUnits, '0005H' ) )

    Img.img       = (scalingFactor^-1 )* Img.img ;
    Img.Hdr.PixelComponentPhysicalUnits = '0000H' ; % i.e. none
end

end
% =========================================================================
function GYRO = getgyromagneticratio( Img )
%GETGYROMAGNETICRATIO
%
% Examines .Hdr of MaRdI-type Img for .ImagedNucleus and returns gyromagnetic
% ratio in units of rad/s/T 
%
% .....
%
% Gyro = getgyromagneticratio( Img )

switch Img.Hdr.ImagedNucleus 
    case '1H' 
        GYRO = 267.513E6 ; 
    otherwise
        error('Unknown nucleus. Do something useful: add case and corresponding gyromagnetic ratio.') ;
end

end
% =========================================================================
function Phase = unwrapphase( Phase, varargin )
%UNWRAPPHASE
% Interface to FSL prelude or to Abdul-Rahman's path-based phase unwrapper
%
%    .......................
%    To call FSL Prelude   
% 
%       Phase = UNWRAPPHASE( Phase, Mag, PreludeOptions )
%       
%       Phase and Mag are objects of type MaRdI
%       See HELP PRELUDE for description of PreludeOptions 
%   
%    .......................
%    To call the Abdul-Rahman unwrapper  
% 
%	    Phase = UNWRAPPHASE( Phase, Options )
%       Phase is an object of type MaRdI and has defined Phase.Hdr.MaskingImage
%       See HELP UNWRAP3D for description of Options 
%
%    .......................
% 
% Abdul-Rahman H, et al. Fast and robust three-dimensional best path phase
% unwrapping algorithm. Applied Optics, Vol. 46, No. 26, pp. 6623-6635, 2007.

    assert( strcmp( Phase.Hdr.PixelComponentPhysicalUnits, '0000H' ), 'SCALEPHASE2RAD() beforehand.' ) ;
    
    isCallingPrelude = false ;
    Options.dummy    = [];
     
    if nargin == 1
        varargin{2} = Options ;
    end

   if nargin < 3 
    
        assert( myisfield( Phase.Hdr, 'MaskingImage') && ~isempty(Phase.Hdr.MaskingImage), ...
            'Logical masking array must be defined in Phase.Hdr.MaskingImage ' ) ;
        
        seriesDescriptionUnwrapper = 'AbdulRahman_2007' ;
        
        if nargin == 2
            Options = varargin{1} ;
        end
        Phase.img = unwrap3d( Phase.img, logical(Phase.Hdr.MaskingImage), Options ) ;

    elseif nargin == 3 
        isCallingPrelude = true ;
        seriesDescriptionUnwrapper = 'FSL_Prelude' ;
        
        Mag     = varargin{1} ;
        Options = varargin{2} ;

        if myisfield( Phase.Hdr, 'MaskingImage') && ~isempty(Phase.Hdr.MaskingImage)
            Options.mask = single( Phase.Hdr.MaskingImage ) ;
        end
        
        Options.voxelSize = Phase.getvoxelsize() ;   

        Phase.img = prelude( Phase.img, Mag.img, Options ) ;

    end

    Phase.img = double( Phase.img ) ;

    % update header
    Img.Hdr.ImageType         = 'DERIVED\SECONDARY' ; 
    Img.Hdr.SeriesDescription = ['phase_unwrapped_' seriesDescriptionUnwrapper ] ; 
end
% =========================================================================
function Field = mapfrequencydifference( Phase1, Phase2 )
%MAPFREQUENCYDIFFERENCE
% 
% Field = MAPFREQUENCYDIFFERENCE( UnwrappedEcho1, UnwrappedEcho2 ) 

Field = MaRdI() ;

assert( strcmp(Phase1.Hdr.PixelComponentPhysicalUnits, '0000H') && ...
        strcmp(Phase2.Hdr.PixelComponentPhysicalUnits, '0000H'), ...
'Inputs should be 2 MaRdI objects where .img is the unwrapped phase.' ) ;

%-------
% mask
if myisfield( Phase1.Hdr, 'MaskingImage' ) || ~isempty( Phase1.Hdr.MaskingImage ) ; 
    mask = Phase1.Hdr.MaskingImage ;
end

if myisfield( Phase2.Hdr, 'MaskingImage' ) || ~isempty( Phase2.Hdr.MaskingImage ) ; 
    mask = mask .* Phase2.Hdr.MaskingImage ;
end

%-------
% field map
Field.img = (Phase2.img - Phase1.img)/(Phase2.Hdr.EchoTime - Phase1.Hdr.EchoTime) ;
Field.img = mask .* Field.img ;
Field.img = (1000/(2*pi)) * Field.img ;


Field.Hdr = Phase1.Hdr ;

%-------    
% update header
Field.Hdr.MaskingImage = mask ;
Field.Hdr.EchoTime = abs( Phase2.Hdr.EchoTime - Phase1.Hdr.EchoTime ) ;

Field.Hdr.PixelComponentPhysicalUnits = '0005H' ; % Hz
Field.Hdr.ImageType         = 'DERIVED\SECONDARY\FREQUENCY' ; 
Field.Hdr.SeriesDescription = ['frequency_' Field.Hdr.SeriesDescription] ; 

end
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
%   r: row index direction cosine
%   c: column index " "
%   s: slice index " "
      
    r = Img.Hdr.ImageOrientationPatient(4:6) ; 
    c = Img.Hdr.ImageOrientationPatient(1:3) ; 
    s = cross( c, r ) ;  
    
    % " To make a full rotation matrix, we can generate the last column from
    % the cross product of the first two. However, Siemens defines, in its
    % private CSA header, a SliceNormalVector which gives the third column, but
    % possibly with a z flip, so that R is orthogonal, but not a rotation
    % matrix (it has a determinant of < 0)."
    %   -From http://nipy.org/nibabel/dicom/dicom_mosaic.html 
    %
    % In general, (I think) the coordinate is increasing along slice dimension
    % with slice index:
    [~,iS] = max( abs(s) ) ;
    if s(iS) < 0
        s = -s ;
    end

end 
% =========================================================================
function fieldOfView = getfieldofview( Img )
%GETFIELDOFVIEW
% 
% fov = getfieldofview( Img ) ;
%
% Returns field of view in units of mm : [Row Column Slice] dimensions
fieldOfView = [ Img.Hdr.PixelSpacing(1) * double( Img.Hdr.Rows ), ...
                Img.Hdr.PixelSpacing(2) * double( Img.Hdr.Columns ), ...
                Img.Hdr.SpacingBetweenSlices * double( Img.Hdr.NumberOfSlices ) ] ;
end
% =========================================================================
function gridSize = getgridsize( Img )
%GETGRIDSIZE
gridSize = double( [ Img.Hdr.Rows, Img.Hdr.Columns, Img.Hdr.NumberOfSlices ] ) ;
end
% =========================================================================
function nVoxels = getnumberofvoxels( Img )
%GETNUMBEROFVOXELS
nVoxels = prod( Img.getgridsize( ) ) ;
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

    % form DICOM standard: https://www.dabsoft.ch/dicom/3/C.7.6.2.1.1/
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
                      [0:1:(Img.Hdr.NumberOfSlices-1)] ) ; 
    
    iRows    = double(iRows);
    iColumns = double(iColumns);
    iSlices  = double(iSlices);
    
    %-------
    % Rotation matrix: R
    [r, c, s] = Img.getdirectioncosines ;
    R = [r c s] ; 
                      
    %-------
    % Scaling matrix: S  
    voxelSize = Img.getvoxelsize() ;
    S = diag(voxelSize); 

    RS = R*S ;

    %-------
    % Scale and rotate to align row direction with x-axis, 
    % column direction with y-axis, slice with z-axis
    X1 = RS(1,1)*iRows + RS(1,2)*iColumns + RS(1,3)*iSlices;
    Y1 = RS(2,1)*iRows + RS(2,2)*iColumns + RS(2,3)*iSlices;
    Z1 = RS(3,1)*iRows + RS(3,2)*iColumns + RS(3,3)*iSlices;

    %-------
    % TRANSLATE w.r.t. origin 
    % (i.e. location of 1st element: .img(1,1,1))
    X = Img.Hdr.ImagePositionPatient(1) + X1 ; 
    Y = Img.Hdr.ImagePositionPatient(2) + Y1 ; 
    Z = Img.Hdr.ImagePositionPatient(3) + Z1 ; 

end
% =========================================================================
function voxelSize = getvoxelsize( Img )
%GETVOXELSIZE
% 
% voxelSize = GETVOXELSIZE( Img )
%       
% Returns 3 component vector of voxel dimensions :
%
%   voxelSize(1) : row spacing in mm (spacing between the centers of adjacent
%   rows, i.e. vertical spacing). 
%
%   voxelSize(2) : column spacing in mm (spacing between the centers of
%   adjacent columns, i.e. horizontal spacing).    
%
%   voxelSize(3) : slice spacing in mm (spacing between the centers of
%   adjacent slices, i.e. from the DICOM hdr, this is Hdr.SpacingBetweenSlices 
%   - for a 2D acquisition this NOT necessarily the same as Hdr.SliceThickness).    

voxelSize = [ Img.Hdr.PixelSpacing(1) Img.Hdr.PixelSpacing(2) ...
        Img.Hdr.SpacingBetweenSlices ] ;
end
% =========================================================================
function Img = resliceimg( Img, X_1, Y_1, Z_1, interpolationMethod ) 
%RESLICEIMG
%
%   Img = RESLICEIMG( Img, X, Y, Z )
%   Img = RESLICEIMG( Img, X, Y, Z, interpolationMethod ) 
%
%   X, Y, Z MUST refer to X, Y, Z patient coordinates (i.e. of the DICOM
%   reference coordinate system)
%   
%   Optional interpolationMethod is a string supported by griddata().
%
%   See: help griddata  

%------ 
% Reslice to new resolution
if nargin < 5
    interpolationMethod = 'linear' ;
end

nImgStacks = size(Img.img, 4) ;

[X_0, Y_0, Z_0] = Img.getvoxelpositions( ) ;

tmp = zeros( [ size(X_1) nImgStacks ] ) ;

for iImgStack = 1 : nImgStacks     
    disp( ['Reslicing image stack...' num2str(iImgStack) ' of ' num2str(nImgStacks) ]) ;
    
    tmp(:,:,:,iImgStack) = ...
        griddata( X_0, Y_0, Z_0, Img.img(:,:,:,iImgStack), X_1, Y_1, Z_1, interpolationMethod ) ;
end

Img.img = tmp ; 
% if new positions are outside the range of the original, 
% interp3/griddata replaces array entries with NaN
Img.img( isnan( Img.img ) ) = 0 ; 

% ------------------------------------------------------------------------

% ------------------------------------------------------------------------
% Update header
Img.Hdr.ImageType = 'DERIVED\SECONDARY\REFORMATTED' ;


Img.Hdr.ImagePositionPatient( 1 ) = X_1(1) ; 
Img.Hdr.ImagePositionPatient( 2 ) = Y_1(1) ;
Img.Hdr.ImagePositionPatient( 3 ) = Z_1(1) ;

%-------
% Rows 
Img.Hdr.Rows = size(Img.img, 1) ;

dx = X_1(2,1,1) - X_1(1,1,1) ;
dy = Y_1(2,1,1) - Y_1(1,1,1) ;
dz = Z_1(2,1,1) - Z_1(1,1,1) ;  

% vertical (row) spacing
Img.Hdr.PixelSpacing(1) = ( dx^2 + dy^2 + dz^2 )^0.5 ; 

% column direction cosine (expressing angle btw column direction and X,Y,Z axes)
Img.Hdr.ImageOrientationPatient(4) = dx/Img.Hdr.PixelSpacing(1) ;
Img.Hdr.ImageOrientationPatient(5) = dy/Img.Hdr.PixelSpacing(1) ;
Img.Hdr.ImageOrientationPatient(6) = dz/Img.Hdr.PixelSpacing(1) ;

%-------
% Columns 
Img.Hdr.Columns = size(Img.img, 2) ;       

dx = X_1(1,2,1) - X_1(1,1,1) ;
dy = Y_1(1,2,1) - Y_1(1,1,1) ;
dz = Z_1(1,2,1) - Z_1(1,1,1) ;  

% horizontal (column) spacing
Img.Hdr.PixelSpacing(2) = ( dx^2 + dy^2 + dz^2 )^0.5 ;

% row direction cosine (expressing angle btw column direction and X,Y,Z axes)
Img.Hdr.ImageOrientationPatient(1) = dx/Img.Hdr.PixelSpacing(2) ;
Img.Hdr.ImageOrientationPatient(2) = dy/Img.Hdr.PixelSpacing(2) ;
Img.Hdr.ImageOrientationPatient(3) = dz/Img.Hdr.PixelSpacing(2) ;

%-------
% Slices
Img.Hdr.NumberOfSlices       = size(Img.img, 3) ;

if Img.Hdr.NumberOfSlices > 1
    Img.Hdr.SpacingBetweenSlices = ( (X_1(1,1,2) - X_1(1,1,1))^2 + ...
                                     (Y_1(1,1,2) - Y_1(1,1,1))^2 + ...
                                     (Z_1(1,1,2) - Z_1(1,1,1))^2 ) ^(0.5) ;
else
    Img.Hdr.SpacingBetweenSlices = 0 ;
end

[rHat, cHat, sHat] = Img.getdirectioncosines( ) ;  

Img.Hdr.SliceLocation = dot( Img.Hdr.ImagePositionPatient, sHat ) ;
% Img.Hdr.SliceLocation = dot( [X_1(1) Y_1(1) Z_1(1)],  ... 
%     cross( Img.Hdr.ImageOrientationPatient(1:3), Img.Hdr.ImageOrientationPatient(4:6) ) ) ;

end
% =========================================================================
function [] = write( Img, saveDirectory, imgFormat )
%WRITE Ma(t)-R-dI(com)
%
%.....
%   Syntax
%
%   WRITE( Img )
%   WRITE( Img, saveDirectory )
%   WRITE( Img, saveDirectory, imgFormat ) 
%
%   default saveDirectory = './tmp'
%   
%   imgFormat may be 'dcm' (default) or 'nii'
%
%.....
%
% Adapted from dicom_write_volume.m (D.Kroon, University of Twente, 2009)

DEFAULT_SAVEDIRECTORY = './tmp' ;
DEFAULT_IMGFORMAT     = 'dcm' ;

if nargin < 2 || isempty(saveDirectory)
    saveDirectory = DEFAULT_SAVEDIRECTORY ;
end

if nargin < 3 || isempty(imgFormat)
    imgFormat = DEFAULT_IMGFORMAT ;
end

fprintf(['\n Writing images to: ' saveDirectory ' ... \n'])

[~,~,~] = mkdir( saveDirectory ) ;

% if strcmp( imgFormat, 'dcm' ) 
%     Img = scaleimgtophysical( Img, -1 ) ;
% end

[X,Y,Z] = Img.getvoxelpositions() ;

%-------
% write Hdr

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
Hdr.ImageType               = Img.Hdr.ImageType ;
Hdr.StudyDescription        = ['StudyMAT' num2str(SN)];
Hdr.SeriesDescription       = ['StudyMAT' num2str(SN)];
Hdr.Manufacturer            = Img.Hdr.Manufacturer ;
Hdr.ScanningSequence        = Img.Hdr.ScanningSequence ;
Hdr.SequenceVariant         = Img.Hdr.SequenceVariant ;
Hdr.ScanOptions             = Img.Hdr.ScanOptions ;
Hdr.MRAcquisitionType       = Img.Hdr.MRAcquisitionType ;
% Hdr.StudyInstanceUID        = Img.Hdr.StudyInstanceUID ;
Hdr.SliceThickness          = Img.Hdr.SliceThickness ;
Hdr.SpacingBetweenSlices    = Img.Hdr.SpacingBetweenSlices ;
Hdr.PatientPosition         = Img.Hdr.PatientPosition ;
Hdr.PixelSpacing            = Img.Hdr.PixelSpacing ;

Hdr.ImageOrientationPatient = Img.Hdr.ImageOrientationPatient ;
Hdr.SliceLocation           = Img.Hdr.SliceLocation ; 
Hdr.NumberOfSlices          = Img.Hdr.NumberOfSlices ; 

[rHat, cHat, sHat] = Img.getdirectioncosines( ) ;  

for sliceIndex = 1 : Img.Hdr.NumberOfSlices 
    
    %-------
    % filename 
    sliceSuffix = '000000' ;
    sliceSuffix = [sliceSuffix( length(sliceIndex) : end ) num2str(sliceIndex) ] ;
    sliceSuffix = ['-' sliceSuffix '.dcm'] ;
    filename    = [saveDirectory '/' Img.Hdr.PatientName.FamilyName sliceSuffix] ;

    %-------
    % slice specific hdr info 
    Hdr.ImageNumber          = sliceIndex ;
    Hdr.InstanceNumber       = sliceIndex ;
    Hdr.ImagePositionPatient = [(X(1,1,sliceIndex)) (Y(1,1,sliceIndex)) (Z(1,1,sliceIndex))] ;

    Hdr.SliceLocation        = dot( Hdr.ImagePositionPatient, sHat ) ;
   
    dicomwrite( Img.img(:,:,sliceIndex) , filename, ...; % Hdr ) ;
        'ObjectType', 'MR Image Storage',Hdr ) ;
    
    if( sliceIndex==1 )
        info                  = dicominfo( filename ) ;
        Hdr.StudyInstanceUID  = info.StudyInstanceUID ;
        Hdr.SeriesInstanceUID = info.SeriesInstanceUID ;
    end

end

%-------
if strcmp( imgFormat, 'nii' ) 
    
    tmpSaveDirectory = [saveDirectory '_nii'] ;
    [~,~,~] = mkdir( tmpSaveDirectory ) ;
    
    dicm2nii( saveDirectory, tmpSaveDirectory ) ;
    rmdir( saveDirectory, 's' ) ;

    copyfile( tmpSaveDirectory, saveDirectory ) ;
    rmdir( tmpSaveDirectory, 's' ) ;
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

Img.Hdr.NumberOfSlices = Img.Hdr.Private_0019_100a ;

nImgPerLine = ceil( sqrt( double(Img.Hdr.NumberOfSlices) ) ); % always nImgPerLine x nImgPerLine tiles

nRows    = size(Img.img, 1) / nImgPerLine; 
nColumns = size(Img.img, 2) / nImgPerLine; 
nVolumes = size(Img.img, 3 ) ;

img = zeros([nRows nColumns Img.Hdr.NumberOfSlices nVolumes], class(Img.img));


for iImg = 1 : double(Img.Hdr.NumberOfSlices)

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
voxelSize = Img.getvoxelsize() ;
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
function timeStd = timestd( Img )
%TIMESTD
% 
% standardDeviation = TIMESTD( Img ) 
% 
% Assumes 4th dimension of Img.img corresponds to time:
%   
%   standardDeviation = std( Img.img, 0, 4 ) ;

timeStd = std( Img.img, 0, 4 ) ;

end
% =========================================================================
function timeAverage = timeaverage( Img )
%TIMEAVERAGE
% 
% Img = TIMEAVERAGE( Img) 
% 
% Assumes 4th dimension of Img.img corresponds to time:
%   
%   timeAverage = mean( Img.img, 4 ) ;

timeAverage = mean( Img.img, 4 ) ;

end
% =========================================================================

end
% =========================================================================
% =========================================================================
methods(Static)
% =========================================================================
function [isImaDir] = changefilesuffix( imgDirectory )
%CHANGEFILESUFFIX
% 
% Renames all .IMA files in imgDirectory to .dcm
% Returns true if directory contains .IMA files, false otherwise.
% 
% isImaDir = CHANGEFILESUFFIX( imgDirectory )

listOfImages = dir( [imgDirectory '/*.IMA'] ) ;

listOfImages.name ;

nImg = length(listOfImages) ;

if nImg == 0
    isImaDir = false ;
    display('Given directory does not contain .IMA files') ;

else 
    isImaDir = true ;

    for img = 1 : length(listOfImages)
        [PATHSTR,NAME,EXT] = fileparts( listOfImages(img).name ) ;
        [s,msg,~] = movefile( [ imgDirectory '/' NAME EXT ], [imgDirectory '/' NAME '.dcm'] ) ;
    end
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
function [] = writeimg( img, Params )
%WRITEIMG
%
%   Description
%   
%   Write .png image to file using the 'export_fig' tool
%
%    .......................
%
%   Syntax
%
%   WRITEIMG( img, Parameters )
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
end

[~,~,extension] = fileparts( Params.filename ) ;

if ~strcmp( extension, '.png' )
    Params.filename = [Params.filename '.png' ] ;
end

% =========================================================================
% Create figure
% =========================================================================


figure

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
% function [] = separateechoesintosubdirectories( dataLoadDir ) 
% %separateechoesintosubdirectories
% %
% %   [] = SEPARATEECHOESINTOSUBDIRECTORIES( dataLoadDir ) 
% %
% % Siemens img files on Martinos 7T do not separate echoes into different 
% % directories. As it stands, this means that instantiation of a MaRdI-type
% % object of these files will load them all in a single 3d array.
% % This function reorganizes them into subdirectories 
%
% nImg = size( Img.img, 3 ) ;
%
% echoTimes = [] ;
% img       = Img.img(:,:,1) ;
%
% % for iImg = 2 : nImg
% %     
% %     Img.Hdr.
%
%
% end

   
end
% =========================================================================
% =========================================================================




end
