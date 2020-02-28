function [] = mustBeString( A )
%% MUSTBESTRING Validate value is a string 
%
% MUSTBESTRING is a validation function which issues an error if the
% input argument is not a string 
% 
% Usage
%
% [] = MUSTBESTRING( A ) 
%
% ### References ###
%
% See also 
%
% <https://www.mathworks.com/help/matlab/matlab_prog/argument-validation-functions.html validation functions>
%%

if ~isstring( A )
    error('Value must be a string.') ;
end

end
