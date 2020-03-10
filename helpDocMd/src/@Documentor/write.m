function [ pathOut ] = write( Dr )  
%WRITE Write documentation to file
% 
% ### Syntax
%    
%    [pathOut] = write( Self ) 
%
% writes the contents of Self.mdDoc to file. To _overwrite_ an existing file,
% before calling WRITE, set 
%     Self.isOverwriting = true 

pathOut = [ Dr.dirOutTop + filesep + Dr.Info.Attributes.Name + Dr.extOut ] ;

assert( Dr.isOverwriting || ~exist( pathOut ), ...
    ['Doc file already exists. Assign a different file path for the output,' ...
     'or set ' inputname(1) '.isOverwriting =true to force overwrite'], '%s' );

[fid, errMsg] = fopen( pathOut, 'w+') ;
assert( fid~=-1, ['Write failed:' errMsg] ) ;

fprintf( strcat("Writing doc: ", pathOut, "\n") ) ;
fprintf( fid, '%s\n', Dr.mdDoc);

fclose(fid);

end

