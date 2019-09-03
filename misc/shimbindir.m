function shimDir = shimbindir( )
%SHIMBINDIR  Returns binary directory of shim reference maps as a string
%
% shimDir = SHIMBINDIR( ) ;
%
% by default, shimDir = '~/Matlab/shimReferenceMaps/'
%
% If you want to store the binaries in a different directory, replace 
% the corresponding entry in this function with the appropriate directory.

% change this if the reference map binaries are to be stored elsewhere:
shimDir = '~/Matlab/shimming/shimReferenceMaps/' ;

assert( exist(shimDir) == 7, 'Shim binary directory not found.' ) ;

end
