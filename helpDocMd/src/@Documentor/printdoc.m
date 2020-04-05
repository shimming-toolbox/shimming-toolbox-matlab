function [ dFile ] = printdoc( Dr, iM )
%PRINTDOC Write documentation to file(s)
%    
%    [dFiles] = PRINTDOC( Dr ) 
%    [dFile]  = PRINTDOC( Dr, iM ) 
%
% When called with a single argument (Documentor instance `Dr`), PRINTDOC 
% iteratively descends the list of .m files in `Dr.mFiles`, printing default
% documentation to the output file paths listed in `Dr.dFiles`.


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

if nargin == 1
    for iM = 1 : numel( Dr.mFiles )
        fprintf( strcat("Preparing doc ", num2str(iM), "/", num2str(numel( Dr.mFiles )), "\n" ) );
        dFile(iM) = Dr.printdoc( iM ) ;
    end

    return ;
end

try
    Dr.iM   = iM ;
    dFile   = Dr.dFiles( Dr.iM ) ;
    dFolder = fileparts( dFile ) ;

    assert( Dr.isOverwriting || ~exist( dFile ), ...
        ['Doc file already exists. Assign a different output file path (`dFiles` property),' ...
         'or change the `isOverwriting` property to `true` to force overwrite'], '%s' );

    if ~isfolder( dFolder )
        [isMade, msg, msgId] = mkdir( dFolder ) ;
        assert( isMade, msgId, ['Directory creation failed: ' msg ] )
    end    

    [fid, errMsg] = fopen( dFile, 'w+' ) ;
    assert( fid~=-1, ['Print abort-fail for :' char(dFile) '\n' errMsg], '%s' ) ;

    fprintf( fid, '%s\n', Dr.mdDoc ) ;
    fprintf( strcat("Written to: ", dFile, "\n") ) ;

    fclose(fid);

catch Me
    Me.rethrow() ;
end

end
