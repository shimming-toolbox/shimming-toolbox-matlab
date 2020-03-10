function [docStr] = documentbasic( Dr )
%DOCUMENTBASIC Return string vector of rudimentary documentation for all .m file types
%
% DOCUMENTBASIC Documents the following .m file details:
% - Name : of the script, function, or class
% - Type: script, function, or class
% - Description: header line from the help/documentation
% - DetailedDescription: body of the help/documentation

Info = Dr.Info.Attributes ;

assert( ~strcmp(Info.mType, "NA"), 'Invalid mFile' ) ;

docStr = [ strcat( "# ", Info.Name ) ; "" ];
docStr = [ docStr ;  strcat(  "**Filetype:** _MATLAB&reg; ", Info.mType, "_" ) ; "" ] ;
docStr = [ docStr ;  strcat(  "**Synopsis:** _",  Info.Description, "_" ) ; "" ] ;
docStr = [ docStr ; Info.DetailedDescription ; "" ] ;

end 

