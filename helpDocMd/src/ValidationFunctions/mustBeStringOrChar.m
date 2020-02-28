function [] = mustBeStringOrChar( A )
%% MUSTBESTRINGORCHAR Validate value is a string or character array
%
% MUSTBESTRINGORCHAR is a validation function which issues an error if the
% input argument is not a string or character array.
% 
% Usage
%
% [] = MUSTBESTRINGORCHAR( A ) 
% 
%
% ### References ###
%
% See also
%
% <https://www.mathworks.com/help/matlab/matlab_prog/argument-validation-functions.html validation functions>
%%

if ~strcmp( class(A), {'string', 'char'} )
    error('Value must be a string or char array.') ;
end

end
