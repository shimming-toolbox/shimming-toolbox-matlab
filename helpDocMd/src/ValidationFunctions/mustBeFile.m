function [] = mustBeFile( A )
%% MUSTBEFILE Validate value is a filename (or set of filenames) on the path 
%
% MUSTBEFILE is a validation function
% which wraps to isfile() and issues an error if the input argument is not comprised solely of 
% existing files on the path. (The input can be a string array, character array, or
% cell array of character vectors (aka "cell-string") of any size.)
% 
% ### Usage ###
%
% [] = MUSTBEFILE( A ) 
%
% ### References ###
%
% See also 
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
isPath = isfile( A ) ;

if any( ~isPath )
    error( strcat( 'Not an existing file path:\n', join( A(~isPath),"\n") ), '%s' ) ;
end

end


