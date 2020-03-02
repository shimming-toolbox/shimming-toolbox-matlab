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

docStr    = strings( 6,1 ) ;
docStr(1) = strcat( "# ", Info.Name, " #", " (a MATLAB(R) ", Info.mType, ")" ) ;
docStr(2) = "" ;
docStr(3) = strcat( "_", Info.Description, "_" ) ;
docStr(4) = "" ;
docStr(5) = "### Description ###"
docStr    = [docStr ; Info.DetailedDescription ] ;

end 

