function [] = mustBeBoolean( A )
%% MUSTBEBOOLEAN Validate value is a 1 or 0, true or false
%
% MUSTBEBOOLEAN is a
% <https://www.mathworks.com/help/matlab/matlab_prog/argument-validation-functions.html validation function> 
% which issues an error if the input argument is not comprised solely of 1's or
% 0's, true or false. (The input can be an array of any size.)
% 
% The wrapper function is merely shorthand for: mustBeMember( A, [ 0 1 ] )
% 
% ### Usage ###
%
% [] = MUSTBEBOOLEAN( A ) 
%
% ### References ###
%
% See also 
%
% <https://www.mathworks.com/help/matlab/ref/mustbemember.html mustBeMember>
%
% <https://www.mathworks.com/help/matlab/matlab_prog/argument-validation-functions.html validation functions> 
%%

mustBeMember( A, [ 0 1 ] ) ;

end
