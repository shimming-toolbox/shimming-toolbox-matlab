function [docStr] = documentfunction( Dr )
%DOCUMENTFUNCTION adds function-specific info to documentation  
% (NOTE: for now, this is just nArgin/nArgout but this should be elaborated
% in Informer.m -- e.g. by parsing the function arguments block when it exists)

assert( Dr.Info.mType=="function", 'mFile is not a function' ) ;

Info   = Dr.Info.Attributes ;

docStr = Dr.documentbasic( ) ;

% remove fields included in documentbasic 
Info   = rmfield( Info, ["mType" ; "Name" ; "Description" ; "DetailedDescription"] ) ;
fields = string( fieldnames( Info ) ) ;

docStr = [docStr ; "" ; "### Attributes" ; "" ] ;

for iField = 1 : numel(fields)
    field = char( fields( iField ) ) ;
    docStr(end+1) = "" ;
    docStr(end+1) = strcat( "- ", fields(iField), " : ", string( Info.( field ) ) ) ;
end

end
