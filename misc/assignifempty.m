function [StructObjOut] = assignifempty( StructObjIn, varargin )
%ASSIGNIFEMPTY  Assigns values to object or struct fields when empty
% 
% [StructObjOut] = ASSIGNIFEMPTY( StructObjIn, fieldName, x )
% [StructObjOut] = ASSIGNIFEMPTY( StructObjIn, Assignments )
%
% ASSIGNIFEMPTY( StructObjIn, fieldName, x )
%
%   Calls MYISFIELDFILLED to check if the data struct or object StructObjIn has
%   a non-empty field called fieldName and returns StructObjOut: a copy of the
%   input StructObjIn, wherein fieldName is assigned the value of x if and only
%   if the entry was non-existent or empty.
%
% [StructObjOut] = ASSIGNIFEMPTY( StructObjIn, Assignments )
%
%   Assignments is a struct or object possessing all the fields/properties that
%   are to be assigned (copied) to StructObjOut. This is equivalent to making
%   repeated calls to ASSIGNIFEMPTY() with multiple fieldNames and x-values.
%
%   e.g. 
%       StructObjIn.w1 = 'Hello' ;
%   
%       Assignments.w2 = 'World' ;
%       Assignments.w3 = '!' ;
%
%       StructObjOut = assignifempty( StructObjIn, Assignments ) ;
%       display( StructObjOut ) ;
%
%       StructObjOut =
%
%           w1: 'Hello'
%           w2: 'World'
%           w3: '!'

if ( nargin < 2 ) || ( nargin > 3 ) || ~( isstruct( StructObjIn ) || isobject( StructObjIn ) ) 
    help(mfilename); return;

elseif nargin == 2
    Assignments = varargin{1} ;
    if ~( isstruct( Assignments ) || isobject( Assignments ) ) 
        help(mfilename); return;
    else
        fieldNames   = fieldnames( Assignments ) ;
        StructObjOut = StructObjIn ;
        
        for iField = 1 : length( fieldNames )
            StructObjOut = assignifempty( StructObjOut, fieldNames{iField}, Assignments(1).( fieldNames{iField} ) ) ;
        end
    end

elseif nargin == 3
    fieldName = varargin{1} ;
    x         = varargin{2} ;

    if ~ischar( fieldName )
        help(mfilename); return;
    else    
        StructObjOut = StructObjIn ;

        if ~myisfieldfilled( StructObjIn, fieldName )
            StructObjOut(1).( fieldName ) = x ;
        end
    end
end

end
