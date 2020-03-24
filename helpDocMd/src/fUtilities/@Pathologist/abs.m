function [ pathOut, pathType ] = abs( pathIn )
%abs Validate and convert file system paths
%    
%    [pathOut, pathType] = abs( pathIn )  
% 
% Wraps to the Matlab function `fileattrib` to check each element of the input
% array `pathIn` for valid file system paths (which can be relative or
% abs/absolute). It returns two arrays of the same type (string, char, or
% cellstr) and size (arbitrary) as the input.
%
% The possibilities for a given return element `(i)` are:
% 
% `pathIn(i)`  | `pathOut(i)`| `pathType(i)`                       
% -------------|-------------|-------------------------------------
%  is a file   | abs path   | "file"                              
%  is a folder | abs path   | "directory"                         
%  else        |   ""        | the error message from [fileattrib] 
%
% **Note** for invalid input paths, the error message returned from
% `fileattrib( pathIn(i) )` will presumably be 'No such file or directory.'
% However the MATLAB documentation is unclear whether other possibilities exist
% for the case where the function is called with a single argument, as is done
% here. For more info, refer to the MATLAB documentation: 
% [fleattrib](https://www.mathworks.com/help/matlab/ref/fileattrib.html/)
%
% See also 
% FILEATTRIB 
    arguments
        pathIn { mustBeStringOrCharOrCellstr } ;
    end

    P = Pathologist( pathIn ) ;
    [ pathOut, pathType ] = callfileattrib( P.data ) ;
    [ pathOut, pathType ] = P.returnasinput( pathOut, pathType ) ;

end % abs()

function [ pathOut, pathType ] = callfileattrib( pathInStr )
%CALLFILEATTRIB
%
% NOTE: Returns 1&2 from arrayfun-call @fileattrib are cell arrays respectively containing:
% 1. scalar logicals == true for valid paths.
% 2. structs of file/folder attributes for valid paths; char vectors of error messages otherwise.
[far1, far2] = arrayfun( @fileattrib, pathInStr, 'UniformOutput', false ) ;
isPath       = reshape( [ far1{:} ], size( pathInStr ) ) ;

% initialize returns 
pathOut      = repmat( "", size( pathInStr ) ) ;
pathType     = repmat( "", size( pathInStr ) ) ;

if any( ~isPath ) 
    pathType( ~isPath ) = string( far2( ~isPath ) ) ;
end 

if any( isPath )

    PathAttributes          = [ far2{isPath} ] ;
    pathOut( isPath )       = string( { PathAttributes.Name } ) ;

    validPathType           = strings( size( PathAttributes ) ) ;
    isDir                   = [ PathAttributes.directory ] ;
    validPathType( isDir )  = "directory" ;
    validPathType( ~isDir ) = "file" ;

    pathType( isPath )      = validPathType ;

end 
 
end % callfileattrib
