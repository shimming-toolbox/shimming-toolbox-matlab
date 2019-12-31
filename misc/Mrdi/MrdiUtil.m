classdef (Abstract) MrdiUtil < handle
%MrdiUtil   MR DICOM Image Utility  
%
% Member methods for handling Mrdi objects 
% 
% e.g. public methods:
%
%   copy()
%   getproperties()
%
% e.g. protected methods:
%
%   copyproperties()
%   exist()
%
%   ...etc.
% 
% For documentation, type doc MrdiUtil
%
% =========================================================================
% Author::ryan.topfer@polymtl.ca
% =========================================================================

% =========================================================================
% =========================================================================    
methods
% =========================================================================    
function Proc = MrdiUtil(  )

end
% =========================================================================
% =========================================================================
end

% =========================================================================
% =========================================================================    
methods
% =========================================================================
function ImgCopy = copy( Img )
%COPY   Return copy of a Mrdi-object.
% 
% ImgCopy = COPY( Img ) ;

constructorHandle = str2func( class( Img ) ) ;
ImgCopy = constructorHandle() ;

ImgCopy.copyproperties( Img ) ;

end
% =========================================================================
function Props = getproperties( Img )
%GETPROPERTIES  Return struct containing copies of a Mrdi-object's properties
% 
% Props = GETPROPERTIES( Img ) ;

C   = metaclass( Img ) ;
Tmp = C.Properties ;

for iProp = 1: length( Tmp )
    if ~Tmp{iProp}.Dependent
       Props.( Tmp{iProp}.Name ) = Img.( Tmp{iProp}.Name ) ;
    end
end

end
% =========================================================================

end
% =========================================================================
% =========================================================================
methods(Access=protected)
% =========================================================================
function [] = copyproperties( ImgCopy, ImgOriginal )
%COPYPROPERTIES  Copies properties from one MaRdI-object to another
%
% COPYPROPERTIES( Copy, Original ) 
%
% Copies the properties of Mrdi-object Original to Mrdi-object Copy 
% (overwriting its initial values)

assert( isa( ImgOriginal, 'MrdiUtil' ) ) ;
Props = ImgOriginal.getproperties() ;
names = fieldnames( Props ) ;

for iProp = 1: length( names )
    ImgCopy.( names{iProp} ) = ImgOriginal.( names{iProp} ) ;
end

end
% =========================================================================
function isSame = iscoincident( Img1, Img2 )
%ISCOINCIDENT   Check coincidence of 2 images 
% 
% isSame = ISCOINCIDENT( Img1, Img2 )
%
% Returns TRUE if Img1 and Img2 possess coincident voxel positions
% and number of measurements/volumes
%
% TODO: Check additional properties + add outputs for each?

assert( ( nargin == 2 ) && isa( Img2, 'MrdiUtil' ), 'Missed required input: 2 MrdiUtil-objects' )

isSame = false ;

if Mrdi.compareimggrids( Img1, Img2 )  
     
    isSame = true ;

    if ~strcmp( Img1.Hdr.SeriesDescription, Img2.Hdr.SeriesDescription )
        warning('A computation is being performed based on images acquired in seperate series.')
    end
end

end    
% =========================================================================

end
% =========================================================================
% =========================================================================
methods(Static)
% =========================================================================
function isSame = compareimggrids( varargin )
%COMPAREIMGGRIDS  Return TRUE if voxel positions of 2 Mrdi image objects coincide 
%
% isSame = COMPAREIMGGRIDS( Img1, Img2 )
% isSame = COMPAREIMGGRIDS( X1, Y1, Z1, X2, Y2, Z2 )

%% -------
% check & parse inputs
assert( nargin >= 2, 'Function requires at least 2 input arguments. See HELP MaRdI.compareimggrids' ) 

if isa( varargin{1}, 'MrdiUtil' ) && isa( varargin{2}, 'MrdiUtil' )
    [X1, Y1, Z1] = varargin{1}.getvoxelpositions ;
    [X2, Y2, Z2] = varargin{2}.getvoxelpositions ;
elseif ( nargin ==6 ) && all( cellfun( @isnumeric, varargin ) )
    X1 = varargin{1} ;
    Y1 = varargin{2} ;
    Z1 = varargin{3} ;
    X2 = varargin{4} ;
    Y2 = varargin{5} ;
    Z2 = varargin{6} ;
else
    error('See HELP MrdiUtil.compareimggrids ')
end

%% -------
% compare grid positions
if ( numel(size(X1)) ~= numel(size(X2)) ) || any( size(X1) ~= size(X2) ) || any( X1(:) ~= X2(:) ) || any( Y1(:) ~= Y2(:) ) || any( Z1(:) ~= Z2(:) )
    isSame = false ;
elseif ( all(X1(:) == X2(:) ) && all( Y1(:) == Y2(:) ) && all( Z1(:) == Z2(:) ) ) 
    isSame = true ;
else
    error('Unexpected result: Check conditions apparently insufficient. Modify code')
end

end
% =========================================================================

end
% =========================================================================
% =========================================================================
methods(Access=protected, Hidden=true)
% =========================================================================
function [is] = exist( Img )
%EXIST  Return true (object exists)

is = true ;

end

end
% =========================================================================
% =========================================================================


end
