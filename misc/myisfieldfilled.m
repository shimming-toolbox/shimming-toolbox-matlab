function isFieldFilled = myisfieldfilled( StructObjIn, fieldName )
%MYISFIELDFILLED    Check if field of struct or object exists AND is non-empty
%
% isFieldFilled = MYISFIELDFILLED( StructObjIn, fieldName )
% 
% StructObjIn is a structure, object, or an array of structures to search
% fieldName is the name of the field for which the function searches
%
% -> Returns TRUE if fieldName exists *AND* is the field is filled (non-empty)
% -> Returns FALSE otherwise
% 
% MYISFIELDFILLED wraps to MYISFIELD( ) and is equivalent to
%
% isFieldFilled = myisfield( StructObjIn, fieldName ) && ~isempty( StructObjIn(1).(fieldName) ) 
%
% See also: MYISFIELD, ASSIGNIFEMPTY

if ( nargin ~= 2 ) || ~ischar( fieldName ) || ...
        ~( isempty( StructObjIn) || isstruct( StructObjIn ) || isobject( StructObjIn ) )
    help(mfilename); 
    return; 
else
    isFieldFilled = false ;
    if isempty( StructObjIn )
        return ;
    elseif myisfield( StructObjIn, fieldName ) ;
        if ~isempty( StructObjIn(1).(fieldName) )
            isFieldFilled = true ;
            return ;
        end
    end
end

end
