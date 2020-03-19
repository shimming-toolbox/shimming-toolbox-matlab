function [ dFile ] = printdoc( Dr, iM )
%PRINTDOC Write documentation to file(s)
% 
% ### Syntax
%    
%    [dFiles] = PRINTDOC( Dr ) 
%    [dFile]  = PRINTDOC( Dr, iM ) 
%
% When called with a single argument (Documentor instance `Dr`), PRINTDOC 
% iteratively descends the list of .m files in `Dr.mFiles`, printing default
% documentation to the output file paths listed in `Dr.dFiles`.

%
% This is equivalent to calling
% If printing completes without error, 
%
% PRINTDOC prints the
% current contents of `Dr.mdDoc` to `docFile` (i.e. to `Dr.dFiles(Dr.iM)`).

% If string scalar "all" is given as the second argument, 
%
%
%
% If scalar index `iM` is given as the second argument, PRINTDOC updates `Dr`
% to print `Dr.dFiles(iM)` (short-hand for `Dr.iM = iM ; Dr.printdoc() ;`).
%
% If printing completes without errors, PRINTDOC returns the path(s) to the
% printed documentation as a string vector.
% 
% **Note 1** 
% When called with a second argument, documentation will be printed
% using the default formating. To bypass it (one file at a time): 
%
% 1. Set `Dr.iM` to the element of interest in `Dr.mFiles` to update the
% `Dr.mdDoc` string vector.
%
% 2. Edit or replace the default content of `Dr.mdDoc` as desired. 
%
% 3. Print by calling `Dr.printdoc()` without the second argument.
% 
% **Note 2** 
% To overwrite existing files, set `Dr.isOverwriting = true` prior to calling PRINTDOC.

% if nargin == 1
%     for iM = 1 : numel( Dr.mFiles )
%         docFile(iM) = Dr.printdoc( iM ) ;
%     end
% else
%     try
%         docFile(iM) = Dr.printdoc( iM ) ;
%     catch Me
%         display( help( mfilename('fullfile') ) ) ;
%         display('Print failed.')
%         Me.rethrow() ;
%     end
% end

Dr.iM = iM ;

assert( Dr.isOverwriting || ~exist( Dr.dFiles( Dr.iM ) ), ...
    ['Doc file already exists. Assign a different file path for the output,' ...
     'or set ' inputname(1) '.isOverwriting =true to force overwrite'], '%s' );

[fid, errMsg] = fopen( Dr.dFiles( Dr.iM ), 'w+' ) ;
assert( fid~=-1, ['Print abort-fail: ' errMsg], '%s' ) ;

fprintf( strcat("Writing doc: ", Dr.dFiles( Dr.iM ), "\n") ) ;
fprintf( fid, '%s\n', Dr.mdDoc ) ;

fclose(fid);

docFile = Dr.dFiles( Dr.iM ) ;

end

% errMsg = [ 'Optional 2nd input must be either:\n ' ...
%            '1. A valid scalar index to the vector of .m file paths ' ...
%            '( i.e. ' inputname(1) '.mFiles ); Or,\n ' ...
%            '2. A string = "all" (to iteratively document all .m files).\n ' ] ;
