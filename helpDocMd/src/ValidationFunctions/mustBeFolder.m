function [] = mustBeFolder( A )
%% MUSTBEFOLDER Validate value is a folder (or set of folders) on the path 
%
% MUSTBEFOLDER is a validation function
% which wraps to isfolder() and issues an error if the input argument is not comprised solely of 
% existing folders on the path. (The input can be a string array, character array, or
% cell array of character vectors (aka "cell-string") of any size.)
% 
% ### Usage ###
%
% [] = MUSTBEFOLDER( A ) 
%
% ### References ###
%
% See also 
%
% <https://www.mathworks.com/help/matlab/ref/isfolder.html isfolder>
%
% <https://www.mathworks.com/help/matlab/matlab_prog/argument-validation-functions.html validation functions> 
%%
    arguments
        A {mustBeStringOrCharOrCellstr} ;
    end

A      = strip( string( A ) ) ;
A      = A(:) ;
isPath = isfolder( A ) ;

if any( ~isPath )
    error( strcat( 'Not an existing directory path:\n', join( A(~isPath),"\n") ), '%s' ) ;
end

end

