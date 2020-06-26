classdef Grid 
%%Grid An image grid 
%

%% ========================================================================
% Properties
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

properties 

    % Grid coordinate system as string. Default = 'PCS' (patient coordinate system)
    coordinateSystem char = ['PCS'] ;

end

properties( Dependent, SetAccess=private )

    % Field of view dimensions: [Row, Column, Slice]
    fov(1,3) double {mustBeNonnegative} ; 
    
    % Number of grid points (i.e. number of image pixels or voxels)
    nPoints(1,1) {mustBeInteger} ; 

    % [Column, Row, Slice] orientations as 3x3 mtx of direction cosine column vectors 
    % 
    %> The direction cosines of a vector are the cosines of the angles between
    %> the vector and the three coordinate axes. Equivalently, they are the
    %> contributions of each component of the basis to a unit vector in that
    %> direction.
    % 
    % https://en.wikipedia.org/wiki/Direction_cosine/
    rotation(3,3) double {mustBeReal} = [eye(3,3)] ;
    
    % 4D array of grid positions: 
    % xCoordinates=xyz(:,:,:,1) 3-column mtx: rows=[x y z]-patient coordinates of the vectorized grid positions. 
    %
    % See also img.Grid.gridpositions
    xyz(:,:,:,3) double {mustBeReal} ; 
    
end

properties( SetAccess=private )

    % Image array dimensions: [#Rows, #Columns, #Slices]
    size(1,3) {mustBeInteger, mustBePositive} = [1 1 1] ; 
    
    % Grid-spacing vector: [Between-Rows, "-Columns, "-Slices]
    spacing(1,3) double {mustBeNonnegative} = [1 1 1] ;
    
end

properties( Dependent, SetAccess=private, Hidden=true )

    % Slice [X;Y,Z] positions. 
    %
    % DICOM Hdr Tag: (0020,1041)
    sliceLocation(3,:) double {mustBeReal} ; 
    
    % 3rd column vector of `rotation` matrix 
    %
    % Unit vector indicating slice direction (ascending vs. descending)
    sliceNormalVector(3,1) ; 

end

properties( SetAccess=private, Hidden=true )

    % Column and Row direction cosines. 
    %
    % DICOM Hdr Tag: (0020,0037)
    imageOrientationPatient(6,1) double {mustBeReal} ; 

    % [X Y Z] coordinates of 1st voxel,img(1,1,1). 
    %
    % DICOM Hdr Tag: (0020,0032)
    imagePositionPatient(3,:) {mustBeReal} ; 

end

%% ========================================================================
% Methods 
% =========================================================================    

% =========================================================================
% =========================================================================    
methods
% =========================================================================    
function self = Grid( Hdrs )
    
    if nargin == 0
        return ;  
    elseif ( nargin == 1 ) && isstruct( Hdrs )
        self = self.initializefromdicom( Hdrs ) ;
    else
        error('Invalid input') ;
    end

end
% =========================================================================
function fov = get.fov( self )
%GETFOV  Return field of view dimensions [units: mm]: [Row, Column, Slice]  

    fov = self.spacing .* double( self.size ) ;

end
% =========================================================================
function nPoints = get.nPoints( self )
%GET.NVOXELS  Return # image voxels of a single image volume as scalar double 
    
    nPoints = prod( self.size ) ;

end
% =========================================================================
function [R] = get.rotation( self ) 
%GET.ROTATION
% 
% R = GET.ROTATION( self ) 
%   
% i.e.  `R = [r c s];` where
% r: row *index* direction cosine
% c: column index " "
% s: slice cosine vector: ==  cross( c, r ) OR -cross( c, r ) (i.e. cross( r, c ))
% depending on whether slices are in acending or descending order,
%
% Note the term *index*: r is not the row direction cosine (c is), rather, it
% is the direction cosine of the vector that points along the direction of
% increasing row indices (i.e. it is in the column direction!)
%
% __ETC__
% See also [Wikipedia](https://en.wikipedia.org/wiki/Direction_cosine):
%
%> The direction cosines of a vector are the cosines of the angles between
%> the vector and the three coordinate axes. Equivalently, they are the
%> contributions of each component of the basis to a unit vector in that
%> direction.
 
    c = self.imageOrientationPatient(1:3) ; 
    r = self.imageOrientationPatient(4:6) ; 
    s = self.sliceNormalVector ;

    R = [r c s] ;

end 
% =========================================================================
function [sliceLocation] = get.sliceLocation( self )
%GET.SLICELOCATION Slice [X;Y;Z] positions as defined by the SliceLocation field of the DICOM hdr 

    for iSlice = 1 : self.size(3)
        sliceLocation(:,iSlice) = ...
            dot( self.imagePositionPatient(:,iSlice), self.sliceNormalVector ) ;
    end 

end
% =========================================================================
function [s] = get.sliceNormalVector( self )
%GETSLICENORMALVECTOR Unit vector indicating slice direction (ascending vs. descending)
%
% sliceNormalVector = 3rd column vector of rotation
%
% Defining:
%   c = self.imageOrientationPatient(1:3) 
%   r = self.imageOrientationPatient(4:6)
%   then, s ==  cross( c, r ) OR -cross( c, r ) (i.e. cross( r, c ))
% 
% For multiple slices, the calculation in get.sliceNormalVector() determines
% the correct sign, and defaults arbitrarily to cross( r, c ) in the single
% slice case.
%
% For more info see discussion at 
% http://nipy.org/nibabel/dicom/dicom_mosaic.html 
    
    if self.size(3) == 0 % arbitrary orientation
        s = cross( self.imageOrientationPatient(4:6), self.imageOrientationPatient(1:3) ) ;
    else
        ds = self.imagePositionPatient(:,end) - self.imagePositionPatient(:,1) ;
        s  = ds/( self.spacing(3) * ( double(self.size(3)) - 1 ) ) ;
    end

end
% =========================================================================
function [xyz] = get.xyz( self )
%GET.XYZ Return grid coordinates as 3-column matrix [ x(:) y(:) z(:) ] 

    [x,y,z] = self.gridpositions() ;
    
    xyz          = zeros( [size(x) 3] ) ;
    xyz(:,:,:,1) = x ;
    xyz(:,:,:,2) = y ;
    xyz(:,:,:,3) = z ;

end
% =========================================================================
function [xyz] = set.xyz( self, xyz4d )
%SET.XYZ

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
methods( Hidden = true )
% Hidden for simplicity (overloaded operators should be self-explanatory)
% =========================================================================
function isEqual = eq( self, Grid2 )
%EQ  Return `true` if positions of two Grid objects coincide 

    if ( ( nargin ~= 2 ) || ~isa( Grid2, 'img.Grid' ) )
        error('Function requires 2 img.Grid objects as arguments.') ; 
    elseif ~strcmp( self.coordinateSystem, Grid2.coordinateSystem ) 
        error( 'Input grid coordinate systems must match.' ) ;
    elseif ~strcmp( self.units, Grid2.units ) 
        error( 'Input grid units must match.' ) ;
    else 
        isEqual =img.Grid.comparegridpositions( self.xyz, Grid2.xyz ) ;    
    end

end
% =========================================================================
function isNotEqual = neq( Grid1, Grid2 ) 
%NEQ  Return `false` if positions of two Grid objects coincide 

    isNotEqual = ~eq( Grid1, Grid2 ) ;

end
% =========================================================================    

end
% =========================================================================
% Larger methods in separate files
% =========================================================================
methods
    %.....
    [x,y,z] = gridpositions( self )
end

methods( Static )
    %.....
    isSame = comparegrids( varargin )
end

methods( Access=private )
    %.....
    [self] = initializefromdicom( self, Hdrs )
    %.....
    []     = update( self, X, Y, Z )
end
% =========================================================================
%
% =========================================================================


end
