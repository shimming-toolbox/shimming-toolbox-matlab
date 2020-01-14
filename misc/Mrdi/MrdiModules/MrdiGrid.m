classdef (Sealed = true) MrdiGrid 
%MrdiGrid  MR DICOM Image Grid: Image grid properties.

% =========================================================================
% Author::ryan.topfer@polymtl.ca
% =========================================================================


% =========================================================================
% =========================================================================    

% NB: Properties that can be calculated on-the-fly from other properties are 
% 'Dependent' (effectively, return values from inline functions).
%
% Properties unlikely to be of direct use to the user are kept Hidden.
%
% Properties format:
%   % Decription (appears in doc display)
%   name(size) type {validationFunctions} = [default] ;

properties( Constant = true ) % constant until there is reason to change

    % Units of fov, grid positions, grid spacing
    units char = ['mm'] ; 

end

properties ( Dependent, SetAccess=private )

    % Field of view dimensions: [Row, Column, Slice]
    fov(1,3) double {mustBeNonnegative} ; 
    
    % Number of grid points (i.e. number of image pixels or voxels)
    nPoints(1,1) uint64 {mustBeInteger} ; 

    % [Column, Row, Slice] orientations as 3x3 mtx of direction cosine column vectors 
    % 
    % From <a href="matlab: web('https://en.wikipedia.org/wiki/Direction_cosine')">Wikipedia</a>:
    %   "The direction cosines of a vector are the cosines of the angles between
    %   the vector and the three coordinate axes. Equivalently, they are the
    %   contributions of each component of the basis to a unit vector in that
    %   direction."     
    rotationMatrix(3,3) double {mustBeReal} = [eye(3,3)] ;
    
    % 4D array of grid positions: 
    % xCoordinates=xyz(:,:,:,1) 3-column mtx: rows=[x y z]-patient coordinates of the vectorized grid positions. 
    % See also MrdiGrid.gridpositions
    xyz(:,:,:,3) double {mustBeReal} ; 
    
end

properties 

    % Grid coordinate system as string. Default = 'PCS' (patient coordinate system)
    coordinateSystem char = ['PCS'] ;

end

properties ( SetAccess=private )

    % Image array dimensions: [#Rows, #Columns, #Slices]
    size(1,3) uint64 {mustBeInteger, mustBePositive} = [1 1 1] ; 
    
    % Grid-spacing vector: [Between-Rows, "-Columns, "-Slices]
    spacing(1,3) double {mustBeNonnegative} = [1 1 1] ;
    
end

properties ( Dependent, SetAccess=private, Hidden=true )

    % Slice [X;Y,Z] positions. DICOM Hdr Tag: (0020,1041)
    sliceLocation(3,:) double {mustBeReal} ; 
    
    % 3rd column vector of rotationMatrix: Unit vector indicating slice
    % direction (ascending vs. descending)
    sliceNormalVector(3,1) ; 

end

properties ( SetAccess=private, Hidden=true )

    % Column and Row direction cosines. DICOM Hdr Tag: (0020,0037)
    imageOrientationPatient(6,1) double {mustBeReal} ; %

    % [X Y Z] coordinates of 1st voxel,img(1,1,1). DICOM Hdr Tag: (0020,0032)
    imagePositionPatient(3,:) {mustBeReal} ; 

end

% =========================================================================
% =========================================================================

% =========================================================================
% =========================================================================    
methods
% =========================================================================    
function Grid = MrdiGrid( Hdrs )
    
    if nargin == 1
        Grid = Grid.initializefromdcm( Hdrs ) ;
    end

end
% =========================================================================
function fov = get.fov( Grid )
%GETFOV  Return field of view dimensions [units: mm]: [Row, Column, Slice]  

    fov = Grid.spacing .* double( Grid.size ) ;

end
% =========================================================================
function nPoints = get.nPoints( Grid )
%GETNVOXELS  Return # image voxels of a single image volume as scalar double 
    
    nPoints = prod( Grid.size ) ;

end
% =========================================================================
function [x,y,z] = gridpositions( Grid )
%GRIDPOSITIONS  Return [x,y,z]: Three 2- or 3-D arrays of voxel positions 
% 
% [x,y,z] = GRIDPOSITIONS( Grid ) 
%
% Returns three 3D arrays of doubles, each element containing the
% location [units: mm] of the corresponding voxel with respect to 
% DICOM's 'Reference Coordinate System'.

    % Arrays of voxel row, column, and slice indices
    [iRows,iColumns,iSlices] = ndgrid( [0:1:Grid.size(1)-1], ...
                      [0:1:Grid.size(2)-1], ...
                      [0:1:Grid.size(3)-1] ) ; 

    % Rotation and Scaling matrix: RS  
    RS = Grid.rotationMatrix * diag(Grid.spacing) ;

    % Scale and rotate to align row direction with x-axis, column direction
    % with y-axis, slice with z-axis; then translate w.r.t origin (via ImagePositionPatient)
    x = ( RS(1,1)*iRows + RS(1,2)*iColumns + RS(1,3)*iSlices ) + Grid.imagePositionPatient(1) ;
    y = ( RS(2,1)*iRows + RS(2,2)*iColumns + RS(2,3)*iSlices ) + Grid.imagePositionPatient(2) ;
    z = ( RS(3,1)*iRows + RS(3,2)*iColumns + RS(3,3)*iSlices ) + Grid.imagePositionPatient(3) ;

end
% =========================================================================
function [R] = get.rotationMatrix( Grid ) 
%GETROTATIONMATRIX
% 
% R = GETROTATIONMATRIX( Img ) 
%   
%   i.e.  R = [r c s], where
%   r: row *index* direction cosine
%   c: column index " "
%   s: slice cosine vector: ==  cross( c, r ) OR -cross( c, r ) (i.e. cross( r, c ))
%   depending on whether slices are in acending or descending order,
%
%   Note the term *index*: r is not the row direction cosine (c is), rather, it
%   is the direction cosine of the vector that points along the direction of
%   increasing row indices (i.e. it is in the column direction!)
%   From Wikipedia:  
%   "The direction cosines of a vector are the cosines of the angles between
%   the vector and the three coordinate axes. Equivalently, they are the
%   contributions of each component of the basis to a unit vector in that
%   direction." 
% https://en.wikipedia.org/wiki/Direction_cosine
 
    c = Grid.imageOrientationPatient(1:3) ; 
    r = Grid.imageOrientationPatient(4:6) ; 
    s = Grid.sliceNormalVector ;

    R = [r c s] ;

end 
% =========================================================================
function [sliceLocation] = get.sliceLocation( Grid )
%SLICELOCATION Slice [X;Y;Z] positions as defined by the SliceLocation field 
% of the DICOM header

    for iSlice = 1 : Grid.size(3)
        sliceLocation(:,iSlice) = dot( Grid.imagePositionPatient(:,iSlice), Grid.sliceNormalVector ) ;
    end 

end
% =========================================================================
function [s] = get.sliceNormalVector( Grid )
%GETSLICENORMALVECTOR Unit vector indicating slice direction (ascending vs. descending)
%
% sliceNormalVector = 3rd column vector of rotationMatrix
%
% Defining:
%   c = Grid.imageOrientationPatient(1:3) 
%   r = Grid.imageOrientationPatient(4:6)
%   then, s ==  cross( c, r ) OR -cross( c, r ) (i.e. cross( r, c ))
% 
% For multiple slices, the calculation in get.sliceNormalVector() determines
% the correct sign, and defaults arbitrarily to cross( r, c ) in the single
% slice case.
%
% For more info see discussion at http://nipy.org/nibabel/dicom/dicom_mosaic.html 
    
if Grid.size(3) == 0 % arbitrary orientation
    s = cross( Grid.imageOrientationPatient(4:6), Grid.imageOrientationPatient(1:3) ) ;
else
    ds = Grid.imagePositionPatient(:,end) - Grid.imagePositionPatient(:,1) ;
    s  = ds/( Grid.spacing(3) * ( double(Grid.size(3)) - 1 ) ) ;
end

end
% =========================================================================
function [xyz] = get.xyz( Grid )
%GETXYZ Return grid coordinates as 3-column matrix [ x(:) y(:) z(:) ] 

    [x,y,z] = Grid.gridpositions() ;
    
    xyz          = zeros( [size(x) 3] ) ;
    xyz(:,:,:,1) = x ;
    xyz(:,:,:,2) = y ;
    xyz(:,:,:,3) = z ;

end
% =========================================================================
function [xyz] = set.xyz( Grid, xyz4d )
%SETXYZ

if nargin == 1 
    return ;
end

[X,Y,Z] = deal( xyz4d(:,:,:,1), xyz4d(:,:,:,2), xyz4d(:,:,:,3) ) ;

%% -----
sizeNew = size(X) ;

%% ------
% update spacing

% Displacements [x; y; z;] between rows, columns, and slices, respectively: dR, dC, dS 
dR = [ X(2,1,1) - X(1,1,1) ; Y(2,1,1) - Y(1,1,1) ; Z(2,1,1) - Z(1,1,1) ] ;
dC = [ X(1,2,1) - X(1,1,1) ; Y(1,2,1) - Y(1,1,1) ; Z(1,2,1) - Z(1,1,1) ] ;
    
if sizeNew(3) > 1
    dS = [ X(1,1,2) - X(1,1,1) ; Y(1,1,1) - Y(1,1,2) ; Z(1,1,2) - Z(1,1,1) ] ;
else
    dS = 0 ;
end

spacingNew = sqrt( [ sum(dR.^2) sum(dC.^2) sum(dS.^2) ]  ) ; 

%% -----
% Direction cosines:

% column (expressing angle between column direction and X,Y,Z axes)
imageOrientationPatientNew(4:6) = dR/spacingNew(1) ;
    
% row (expressing angle btw row direction and X,Y,Z axes)
imageOrientationPatientNew(1:3) = dC/spacingNew(2) ;

%% -----
% update imagePositionPatient:
imagePositionPatientNew = [ X(1,1,:) ; Y(1,1,:) ; Z(1,1,:) ] ;


end
% =========================================================================

end
% =========================================================================
% =========================================================================
methods ( Access=private )
% =========================================================================
function [Grid] = initializefromdcm( Grid, Hdrs ) 

% "If Anatomical Orientation Type (0010,2210) is absent or has a value of
% BIPED, the x-axis is increasing to the left hand side of the patient. The
% y-axis is increasing to the posterior side of the patient. The z-axis is
% increasing toward the head of the patient."
%
% from DICOM standard: https://www.dabsoft.ch/dicom/3/C.7.6.2.1.1/
% assert( ~myisfield( Grid.Img.Hdr, 'AnatomicalOrientationType' ) || ...
%             strcmp( Grid.Img.Hdr.AnatomicalOrientationType, 'BIPED' ), ...
%             'Error: AnatomicalOrientationType not supported.' ) ;

Grid.size = uint64( [Hdrs{1}.Rows Hdrs{2}.Columns size(Hdrs,1) ] ) ;

Grid.imageOrientationPatient = Hdrs{1}.ImageOrientationPatient ;

if myisfield( Hdrs{1}, 'SpacingBetweenSlices' )
    Grid.spacing = [ Hdrs{1}.PixelSpacing(1) Hdrs{1}.PixelSpacing(2) Hdrs{1}.SpacingBetweenSlices ] ;
else
    Grid.spacing = [ Hdrs{1}.PixelSpacing(1) Hdrs{1}.PixelSpacing(2) Hdrs{1}.SliceThickness ] ;
end

for iSlice = 1 : Grid.size(3)
    Grid.imagePositionPatient(:,iSlice) = Hdrs{iSlice,1}.ImagePositionPatient ;
end


% if size( Grid.Img.img, 3 ) > 1
% % Determine: ascending or descending slices?
%
%     % Estimate positions of last slice using the position of the 1st image 
%     % and the 2 possibilities for the SliceNormalVector:
%
%     % 1. using cross( c, r ) ;
%     % [X1,Y1,Z1] = [Grid.voxelPositions.X, Grid.voxelPositions.Y, Grid.voxelPositions.Z] ;     
%     [X1,Y1,Z1] = deal( Grid.voxelPositions.X,  Grid.voxelPositions.Y, Grid.voxelPositions.Z ) ;
%     
%     % 2. using the reverse
%     Grid.SliceNormalVector = cross( r, c ) ;
%     [X2,Y2,Z2] = deal( Grid.voxelPositions.X,  Grid.voxelPositions.Y, Grid.voxelPositions.Z ) ;
%     
%     % Actual position corresponding to the slice direction can be increasing or
%     % decreasing with slice/image number. So, which estimate is closer: 1 or 2? 
%     if norm( Grid.Img.Hdrs{end,1,1}.ImagePositionPatient' - [ X1(1,1,end) Y1(1,1,end) Z1(1,1,end) ] ) < ...
%             norm( Grid.Img.Hdrs{end,1,1}.ImagePositionPatient' - [ X2(1,1,end) Y2(1,1,end) Z2(1,1,end) ] ) 
%         % if true, then 1. corresponds to the correct orientation
%         Grid.SliceNormalVector = cross( c, r ) ;
%     end
% end
%     
% function [X,Y,Z] = calculategridpositions( size, spacing, rotationMatrix, imagePositionPatient ) 
%
%     % Arrays of voxel row, column, and slice indices
%     [iRows,iColumns,iSlices] = ndgrid( [0:1:Grid.size(1)-1], ...
%                       [0:1:Grid.size(2)-1], ...
%                       [0:1:Grid.size(3)-1] ) ; 
%
%     % Rotation and Scaling matrix: RS  
%     RS = Grid.rotationMatrix * diag(Grid.spacing) ;
%
%     % Scale and rotate to align row direction with x-axis, column direction
%     % with y-axis, slice with z-axis; then translate w.r.t origin (via ImagePositionPatient)
%     X = ( RS(1,1)*iRows + RS(1,2)*iColumns + RS(1,3)*iSlices ) + Grid.Img.Hdr.ImagePositionPatient(1) ;
%     Y = ( RS(2,1)*iRows + RS(2,2)*iColumns + RS(2,3)*iSlices ) + Grid.Img.Hdr.ImagePositionPatient(2) ;
%     Z = ( RS(3,1)*iRows + RS(3,2)*iColumns + RS(3,3)*iSlices ) + Grid.Img.Hdr.ImagePositionPatient(3) ;
%
% end

end
% =========================================================================
function [] = updategrid( Grid, X, Y, Z )
%UPDATEGRID
 
if ( nargin ~= 4 ) || ~isequal( size(X), size(Y), size(Z) ) 
    error( ['Function requires 4 arguments: MrdiGrid object followed by ' ...
        '3 identically-sized position arrays X,Y,Z with regular spacing'] ) ;
end

%% -----
sizeNew = size(X) ;

%% ------
% update spacing

% Displacements [x; y; z;] between rows, columns, and slices, respectively: dR, dC, dS 
dR = [ X(2,1,1) - X(1,1,1) ; Y(2,1,1) - Y(1,1,1) ; Z(2,1,1) - Z(1,1,1) ] ;
dC = [ X(1,2,1) - X(1,1,1) ; Y(1,2,1) - Y(1,1,1) ; Z(1,2,1) - Z(1,1,1) ] ;
    
if sizeNew(3) > 1
    dS = [ X(1,1,2) - X(1,1,1) ; Y(1,1,1) - Y(1,1,2) ; Z(1,1,2) - Z(1,1,1) ] ;
else
    dS = 0 ;
end

spacingNew = sqrt( [ sum(dR.^2) sum(dC.^2) sum(dS.^2) ]  ) ; 

%% -----
% Direction cosines:

% column (expressing angle between column direction and X,Y,Z axes)
imageOrientationPatientNew(4:6) = dR/spacingNew(1) ;
    
% row (expressing angle btw row direction and X,Y,Z axes)
imageOrientationPatientNew(1:3) = dC/spacingNew(2) ;

%% -----
% update imagePositionPatient:
imagePositionPatientNew = [ X(1,1,:) ; Y(1,1,:) ; Z(1,1,:) ] ;

end
% =========================================================================

end
% =========================================================================
% =========================================================================    
methods ( Hidden = true )
% Hidden for simplicity (overloaded operators should be self-explanatory)
% =========================================================================
function isEqual = eq( Grid1, Grid2 )
%EQ  Return TRUE if positions of 2 MrdiGrid objects coincide 

    if ( ( nargin ~= 2 ) || ~isa( Grid2, 'MrdiGrid' ) )
        error('Function requires 2 MrdiGrid objects as arguments.') ; 
    elseif ~strcmp( Grid1.coordinateSystem, Grid2.coordinateSystem ) 
        error( 'Input grid coordinate systems must match.' ) ;
    elseif ~strcmp( Grid1.units, Grid2.units ) 
        error( 'Input grid units must match.' ) ;
    else 
        isEqual = MrdiGrid.comparegridpositions( Grid1.xyz, Grid2.xyz ) ;    
    end

end
% =========================================================================
function isNotEqual = neq( Grid1, Grid2 ) 
%NEQ  Return FALSE if positions of 2 MrdiGrid objects coincide 

    isNotEqual = ~eq( Grid1, Grid2 ) ;

end
% =========================================================================    

end
% =========================================================================    
% =========================================================================    
methods(Static)
% =========================================================================
function isSame = comparegrids( varargin )
%COMPAREGRIDS  Return TRUE if voxel positions coincide 
%
% isSame = COMPAREGRIDS( Grid1, Grid2 )
% isSame = COMPAREGRIDS( xyz1, xyz2 )
% isSame = COMPAREGRIDS( X1, Y1, Z1, X2, Y2, Z2 )

assert( nargin >= 2, 'Function requires at least 2 input arguments. See HELP MrdiGrid.comparegridpositions' ) 

if nargin == 2
    if isa( varargin{1}, 'MrdiGrid' ) && isa( varargin{2}, 'MrdiGrid' )
        isSame = ( varargin{1} == varargin{2} ) ;   
        return ;
    elseif isnumeric( varargin{1} ) && isnumeric( varargin{2} )
        isSame = isequal( varargin{1}, varargin{2} ) ; 
        return ;
    else
        error('See HELP MrdiGrid.comparegridpositions ')
    end
elseif ( nargin ==6 ) && all( cellfun( @isnumeric, varargin ) )
    isSame  = isequal( varargin{1}, varargin{4} ) ... % compare X-coordinates
           && isequal( varargin{2}, varargin{5} ) ... % compare Y-coordinates 
           && isequal( varargin{3}, varargin{6} ) ;   % compare Z-coordinates 
    return ;
else
    error('See HELP MrdiGrid.comparegridpositions ')
end

end
% =========================================================================

end
% =========================================================================
% =========================================================================


end
