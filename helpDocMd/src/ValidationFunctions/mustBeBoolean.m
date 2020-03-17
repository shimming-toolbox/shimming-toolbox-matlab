function [] = mustBeBoolean( A )
%MUSTBEBOOLEAN Validate value equals 1 or 0, true or false
%
% Issues an error when the input is not comprised solely of 1's or 0's, true or false. 
%    
%    [] = mustBeBoolean( A ) 
%
% Input `A` can be of any size.
%
% mustBeBoolean is a *validation* function:
% <https://www.mathworks.com/help/matlab/matlab_prog/argument-validation-functions.html> 

if all( islogical(A) )
    return ;
elseif isnumeric(A) && all( [A==0] | [A==1] )
    return ;
else
    error( 'Value must equal 1 or 0.' ) ;
end

end
