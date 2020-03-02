function [ pathOut ] = write( Dr )  
%WRITE Write documentation to file
% 
% WRITE creates (or, optionally, overwrites) a .md file and writes to it the
% contents of Self.mdDoc.
%
% ### Syntax ###
%
% [pathOut] = write( Self ) 
%
% To enable overwriting of existing files, set Self.isOverwriting = true 

pathOut = [Dr.dirOutTop + filesep + Dr.nameIn + Dr.extOut] ;

assert( Dr.isOverwriting || ~exist( pathOut ), ...
    ['Doc file already exists. Assign a different file path for the output,' ...
     'or set ' inputname(1) '.isOverwriting =true to force overwrite'], '%s' );

[fid, errMsg] = fopen( pathOut, 'w+') ;
assert( fid~=-1, ['Write failed:' errMsg] ) ;

fprintf( strcat("Writing doc: ", pathOut, "\n") ) ;
fprintf( fid, '%s\n', Dr.mdDoc);

fclose(fid);

end

