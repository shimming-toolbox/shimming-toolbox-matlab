function [] = mustBeStringOrCharOrCellstr( A )
%MUSTBESTRINGORCHARORCELLSTR Validate value is string, char, or cell of chars
%
% [] = mustBeStringOrCharOrCellstr( A )
%
% MUSTBESTRINGORCHARORCELLSTR is a validation function which issues an error if
% the input argument is not a string or character array, i.e.:
%  if ~( isstring( A ) || ischar( A ) || iscellstr( A ) )

if ~( isstring( A ) || ischar( A ) || iscellstr( A ) )
    error('Value must be a string, char, or cell array containing chars exclusively') ;
end

end
