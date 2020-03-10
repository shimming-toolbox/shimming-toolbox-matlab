function [docStr] = documentfunction( Dr )
%DOCUMENTFUNCTION adds function-specific info to documentation  
% (NOTE: for now, this is just nArgin/nArgout but this should be elaborated
% in Informer.m -- e.g. by parsing the function arguments block when it exists)

Info = Dr.Info.Attributes ;

assert( Info.mType=="function", 'mFile is not a function' ) ;

docStr = documentbasic( Dr ) ;

% remove fields included in documentbasic 
Info   = rmfield( Info, ["mType" ; "Name" ; "Description" ; "DetailedDescription"] ) ;
fields = string( fieldnames( Info ) ) ;

docStr = [docStr ; "" ; "### Attributes  "] ;

for iField = 1 : numel(fields)
    field = char( fields( iField ) ) ;
    docStr(end+1) = "" ;
    docStr(end+1) = strcat( "- ", fields(iField), " : ", string( Info.( field ) ) ) ;
end

end
