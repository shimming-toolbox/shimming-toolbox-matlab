function [ pathOut ] = abspath( pathIn )
% ABSPATH Return absolute path from relative path input 
%
% ### Syntax ###
%
% pathOut = ABSPATH( pathIn ) 
%
% #### Inputs ####
% 
% - pathIn : relative file system path(s) as a string array, 'cell-string', or
% character array of any dimension. (If an element of pathIn is already an
% absolute path, it will be returned as the output.)
% 
% #### Outputs ####
%
% - pathOut: absolute path(s), sized and typed according to the input. 
%
% ### References ###
%
% See also dir, fileattrib
    arguments
        pathIn {mustBeStringOrCharOrCellstr, mustBeFileOrFolder} ;
    end

pathInStr = strip( string( pathIn ) ) ;

[~, Values] = arrayfun( @fileattrib, pathInStr ) ;

%% Match return value class/type and size to those of the input:
pathOut = reshape( { Values.Name }, size( pathInStr ) ) ; % cellstr

switch class( pathIn )
    case 'string'
        pathOut = string( pathOut ) ;
    case 'char'
        pathOut = char( pathOut ) ;
    otherwise % double-check...
        assert( iscellstr( pathIn ) && iscellstr( pathOut ), 'Unexpected input type' ) ;
end

end
