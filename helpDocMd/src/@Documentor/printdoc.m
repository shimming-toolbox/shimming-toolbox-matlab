function [ filename ] = printdoc( Dr )  
%PRINTDOC Write documentation to file
% 
% ### Syntax
%    
%    [filename] = PRINTDOC( Self ) 
%
% Prints the contents of `Self.mdDoc` to `filename`. 
% To overwrite an existing file, first set `Self.isOverwriting = true`
Dr.dirInTop
Dr.dirOutTop
filename = [ Dr.dirOutTop + filesep + Dr.Info.Attributes.Name + Dr.extOut ] ;

assert( Dr.isOverwriting || ~exist( filename ), ...
    ['Doc file already exists. Assign a different file path for the output,' ...
     'or set ' inputname(1) '.isOverwriting =true to force overwrite'], '%s' );

[fid, errMsg] = fopen( filename, 'w+') ;
assert( fid~=-1, ['Write failed: ' errMsg], '%s' ) ;

fprintf( strcat("Writing doc: ", filename, "\n") ) ;
fprintf( fid, '%s\n', Dr.mdDoc);

fclose(fid);

end
