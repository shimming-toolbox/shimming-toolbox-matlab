function [] = mustBeA( A, allowedTypes )
%MUSTBEA Validate variable data-type/class 
% 
%   [] = mustBeA( A, allowedTypes )
%
% Throws an error if the datatype (aka "class") of `A` is not included among
% the specified `allowedTypes` (a string array of permitted classes).
%
% __DESCRIPTION__
% `mustBeA` is equivalent to `assert( 1==ismember( class(A), allowedTypes ) );`
% (which hardly deserves its own file) However, as nested function calls are not
% permitted for argument [validation], this function can be used instead.
%
% __SEE__
% 1. [validation](https://www.mathworks.com/help/matlab/matlab_prog/function-argument-validation-1.html)
% 2. [class](https://www.mathworks.com/help/matlab/ref/class.html)
% 3. [mustBeMember](https://www.mathworks.com/help/matlab/ref/mustbemember.html)
% 4. [isa](https://www.mathworks.com/help/matlab/ref/isa.html)
%
% SEE ALSO
% class, mustBeMember, isa

if ismember( class( A ), allowedTypes )
    return 
end

error( 'mustBeA:invalidType', ['Input type/class must be == [' ...
       char( "'" + strjoin( allowedTypes, "'|'") +"'") '].\n'], '%' ) ;  ;

end


