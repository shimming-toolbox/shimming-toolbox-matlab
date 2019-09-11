function [StructObjOut] = assignifempty( StructObjIn, fieldName, x )
%ASSIGNIFEMPTY  Assigns a value to an object or struct field if it is empty
% 
% [StructObjOut] = ASSIGNIFEMPTY( StructObjIn, fieldName, x )
%
% ASSIGNIFEMPTY calls MYISFIELDFILLED to check if the data struct or object
% StructObjIn has a non-empty field called fieldName and returns StructObjOut:
% a copy of the input StructObjIn, wherein fieldName is assigned the value of x
% if and only if the entry was non-existent or empty.

if ( nargin ~= 3 ) || ~( isstruct( StructObjIn ) || isobject( StructObjIn ) ) || ~ischar( fieldName )
    help(mfilename); 
    return;
end
    
StructObjOut = StructObjIn ;

if ~myisfieldfilled( StructObjIn, fieldName )
    StructObjOut(1).( fieldName ) = x ;
end

end
