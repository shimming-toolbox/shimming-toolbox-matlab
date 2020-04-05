function [docStr] = documentbasic( Dr, Att, headingLevel )
%DOCUMENTBASIC Return string vector of rudimentary documentation
%    
%    [docStr] = documentbasic( Self )
%    [docStr] = documentbasic( Self, Att )
%    [docStr] = documentbasic( Self, Att, headingLevel )
%
% When called without a second argument, `documentbasic` derives the following details
% from `Self.Info.Attributes` to return the documentation string vector `docStr`:
% - Name : of the script, function, or class
% - mType: script, function, or class file type
% - Description: header line from the help/documentation
% - DetailedDescription: body of the help/documentation
% 
% When called with a second argument, `documentbasic` works similarly to the
% above case, however, details are instead derived from attributes-struct `Att`
% (mType may be omitted in this case, but the other field names must be
% present.)
%
% Optional 3rd argument is a scalar integer (= 0,1,2,3,4,5, or 6) [default=1]
% indicating the number of '#' signs to precede the 'Name' value (docStr's
% first element). for markdown syntax, this controls the heading level.
    arguments
        Dr Documentor ;
        Att struct = Dr.Info.Attributes ;
        headingLevel(1,1) {mustBeInteger, mustBeNonnegative, mustBeLessThan(headingLevel, 7) } = 1 ;
    end

assert( all( isfield( Att, ["Name" "Description" "DetailedDescription"] ) ) ) ;

if headingLevel == 0
    docStr = [ ""; Att.Name ; "" ] ;
else
    docStr = [ string( repmat( '#', [1 headingLevel] ) ) + strcat(" ", Att.Name ) ; "" ] ;
end

if isfield( Att, 'mType' )
    assert( ~strcmp( Att.mType, "NA" ), 'Invalid mFile' ) ;
    docStr = [ docStr ; strcat( "**Filetype:** _MATLAB&reg; ", Att.mType, "_" ) ; "" ] ;
end

if ( Att.Description ~= "" )
    docStr = [ docStr ; strcat( "**Synopsis:** _",  Att.Description, "_" ) ; "" ] ;
end

if any( Att.DetailedDescription ~= "" ) && ~isequal( Att.DetailedDescription, Att.Description )
    docStr = [ docStr ; Att.DetailedDescription ; "" ] ;
end 

end
