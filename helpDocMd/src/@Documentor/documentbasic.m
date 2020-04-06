function [docStr] = documentbasic( Info, headingLevel )
%DOCUMENTBASIC Return string vector of rudimentary documentation
%    
%    [docStr] = documentbasic( Info )
%    [docStr] = documentbasic( Info, headingLevel )
    arguments
        Info struct ;
        headingLevel(1,1) {mustBeInteger, mustBeNonnegative, mustBeLessThan(headingLevel, 7) } = 1 ;
    end

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

assert( all( isfield( Info, ["Name" "Description" "DetailedDescription"] ) ) ) ;

if headingLevel == 0
    docStr = [ ""; Info.Name ; "" ] ;
else
    docStr = [ string( repmat( '#', [1 headingLevel] ) ) + strcat(" ", Info.Name ) ; "" ] ;
end

if isfield( Info, 'mType' )
    assert( ~strcmp( Info.mType, "NA" ), 'Invalid mFile' ) ;
    docStr = [ docStr ; strcat( "**Filetype:** _MATLAB&reg; ", Info.mType, "_" ) ; "" ] ;
end

if ( Info.Description ~= "" )
    docStr = [ docStr ; strcat( "**Synopsis:** _",  Info.Description, "_" ) ; "" ] ;
end

if any( Info.DetailedDescription ~= "" ) && ~isequal( Info.DetailedDescription, Info.Description )
    docStr = [ docStr ; Info.DetailedDescription ; "" ] ;
end 

end
