function [ docFile ] = printdoc( Dr )  
%PRINTDOC Write documentation to file
% 
% ### Syntax
%    
%    [docFile] = PRINTDOC( Self ) 
%
% Prints the current contents of `Self.mdDoc` to `docFile` (i.e. to `Dr.docFiles(Dr.iM)`)
% To overwrite an existing file, first set `Self.isOverwriting = true`

docFile = Dr.docFiles( Dr.iM ) ;

assert( Dr.isOverwriting || ~exist( docFile ), ...
    ['Doc file already exists. Assign a different file path for the output,' ...
     'or set ' inputname(1) '.isOverwriting =true to force overwrite'], '%s' );

[fid, errMsg] = fopen( docFile, 'w+' ) ;
assert( fid~=-1, ['Write failed: ' errMsg], '%s' ) ;

fprintf( strcat("Writing doc: ", docFile, "\n") ) ;
fprintf( fid, '%s\n', Dr.mdDoc ) ;

fclose(fid);

end
