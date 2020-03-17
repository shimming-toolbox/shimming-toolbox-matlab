function [] = mustBeA( A, permittedTypes )
%MUSTBEA Validate variable data-type/class 
% 
%   [] = mustBeA( A, permittedTypes )
%
% Issues an error if the class/datatype of `A` is not included among those
% listed in `permittedTypes`. 
%
% __NOTE__
% `mustBeA` merely asserts that `ismember( class(A), permittedTypes )`
% evaluates to `true`. However, since nested function calls are not allowed
% for function argument validation (See 1.) `mustBeA` can be used instead. 
%
% __SEE__
% 1. [Arg. validation](https://www.mathworks.com/help/matlab/matlab_prog/function-argument-validation-1.html)
% 2. [class](https://www.mathworks.com/help/matlab/ref/class.html)
% 3. [mustBeMember](https://www.mathworks.com/help/matlab/ref/mustbemember.html)
% 4. [isa](https://www.mathworks.com/help/matlab/ref/isa.html)
%
% SEE ALSO
% class, mustBeMember, isa

if ismember( class( A ), permittedTypes )
    return;
end

error( 'mustBeA:invalidType', ['Input type/class must be == [' ...
       char( "'" + strjoin( permittedTypes, "'|'") +"'") '].\n'], '%' ) ;  ;

end


