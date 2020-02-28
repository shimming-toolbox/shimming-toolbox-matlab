function [] = mustBeFileOrFolder( A )
%% MUSTBEFILEORFOLDER Validate value is a path (or set of paths) to existing files or folders
%
% MUSTBEFILEORFOLDER is a validation function which wraps to isfile() and
% isfolder() and throws an error if the input argument is not comprised solely
% of existing file system paths. (The input can be a string array, character
% array, or cell array of character vectors (aka "cell-string") of any size.)
%
% ---
% ### Usage ###
%
% [] = MUSTBEFILEORFOLDER( A ) 
%
% ### References ###
%
% See also 
%
% <https://www.mathworks.com/help/matlab/ref/isfolder.html isfolder>
%
% <https://www.mathworks.com/help/matlab/ref/isfile.html isfile>
%
% <https://www.mathworks.com/help/matlab/matlab_prog/argument-validation-functions.html validation functions> 
%%
    arguments
        A {mustBeStringOrCharOrCellstr} ;
    end

A      = strip( string( A ) ) ;
A      = A(:) ;
isPath = [ isfile( A ) | isfolder( A ) ] ;

if any( ~isPath )
    error( strcat( 'Not an existing file or directory path:\n', join( A(~isPath),"\n") ), '%s' ) ;
end

end



