function [ docFiles ] = printall( Dr )  
%PRINTALL Prints documentation for the complete list of `mFiles` 
%
% ### Syntax
%    
%    [docFiles] = PRINTALL( Self ) 
%
% Documentation is written to `Dr.docFiles`

for iM = 1 : length( Dr.mFiles )

    Dr.iM        = iM ;
    docFiles(iM) = Dr.printdoc ;

end

end
