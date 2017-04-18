function[ dataArray ] = croparray( dataArray, gridSizeOut ) 
%CROPARRAY returns cropped/central portion of 3D array
%
% Syntax
%   B = croparray( A, desiredSize )
%   
%   returns central portion of 3D array A: data outside 'desiredSize' is
%   trimmed away.
%
% see also PADARRAY( )
%
% =========================================================================
% Updated::20170210::ryan.topfer@polymtl.ca
% =========================================================================

% TODO
% Arbitrary dimensions

gridSizeIn = size( dataArray ) ;
assert( all( gridSizeIn >= gridSizeOut ), 'desired grid size must be <= original' ) ; 

midPoint  = round( gridSizeIn / 2 ) ; 
low       = midPoint - round( gridSizeOut / 2 ) + 1 ; 
high      = low + gridSizeOut - 1;

dataArray = dataArray( low(1):high(1), low(2):high(2), low(3):high(3) ) ;


end

