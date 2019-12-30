classdef (Abstract) MrdiInfo < handle
%MrdiInfo   MR DICOM Image Info
% 
% Member methods pertain to querying image attributes
% 
% e.g.
%   getfieldofview()
%   getgridsize()
%   getnumberofvoxels()
%   ...etc.
% 
% To list all methods, type methods('MrdiInfo')
% For documentation, type doc MrdiInfo
%
% =========================================================================
% Author::ryan.topfer@polymtl.ca
% =========================================================================

% =========================================================================
% =========================================================================    
methods
% =========================================================================    
function Info = MrdiInfo(  )

end
% =========================================================================
% =========================================================================
end

% =========================================================================
% =========================================================================    
methods(Sealed=true)
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
function nVoxels = getnumberofvoxels( Img )
%GETNUMBEROFVOXELS  Return # image voxels of a single image volume 
%
% nVoxels = GETNUMBEROFVOXELS( Img )
%
% Returns a double precision scalar.

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
%GETVOXELSPACING Return 3-component grid-spacing vector [units: mm]
%
% h = GETVOXELSPACING( Img )
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

end
% =========================================================================
% =========================================================================
methods(Access=protected)
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
% NOTE: Hdr.Img.SliceNormalVector gets defined in the Mrdi constructor

c = Img.Hdr.ImageOrientationPatient(1:3) ; 
r = Img.Hdr.ImageOrientationPatient(4:6) ; 
s = Img.Hdr.Img.SliceNormalVector ;

end 
% =========================================================================

end
% =========================================================================
% =========================================================================


end
