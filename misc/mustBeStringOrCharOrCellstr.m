function [] = mustBeStringOrCharOrCellstr( A )
%MUSTBESTRINGORCHARORCELLSTR Validate value is a string, char, or a cell containing strings/chars.
%
% MUSTBESTRINGORCHARORCELLSTR( A ) issues an error if A contains non-characters.
%
% See also: mustBeNumeric    

switch class( A )

    case {'string', 'char'}
        return ;

    case 'cell'
        for i = 1 : numel( A )
            mustBeStringOrChar( A{i} ) ;
        end 

    otherwise 
        error('Value must be a string, char array, or a cell containing strings or chars') ;

end

end
