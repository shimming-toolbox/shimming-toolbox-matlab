function [] = mustBeStringScalarOrCharVector( A )
%% MUSTBESTRINGSCALARORCHARVECTOR Validate value is a string scalar or character vector
%
% MUSTBESTRINGSCALARORCHARVECTOR is a validation function which issues an error if the
% input argument is not a string scalar or character vector.
% 
% Usage
%
% [] = MUSTBESTRINGSCALARORCHARVECTOR( A ) 
%
% ### References ###
%
% See also 
%
% <https://www.mathworks.com/help/matlab/matlab_prog/argument-validation-functions.html validation functions>

switch class( A )
    case 'string'
        if isscalar( A )
            return ;
        end
    
    case 'char'
        if prod( size( A ) ) == numel( A )
            return ;
        end
end
    
error('Value must be a string scalar or character vector.') ;

end

