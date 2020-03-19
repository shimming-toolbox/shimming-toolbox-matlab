function [mFiles] = findfiles( src, isSearchRecursive )
%FINDFILES Return list of documentable .m files
%
% ### Syntax
%    
%    mFiles = Documentor.findfiles( src )
%    mFiles = Documentor.findfiles( src, isSearchRecursive )
%
% Documentor.findfiles returns a list of *documentable* .m files `mFiles`
% based on the given set of file system paths `src`.
%
% -When `src` path elements point to .m files, Documentor.findfiles verifies that 
% the Documentor class will indeed be able to document them.
%
% -When `src` path elements point to directories, Documentor.findfiles first
% searches the directories for .m files (see **Note 3.** below), followed by the
% compatibliy check described above.
% 
% -When called with a single argument, subfolders of any directory elements in
% `src` are recursively searched for files (i.e. by default,
% `isSearchRecursive=true`). To restrict the search depth to directories
% explicitly included in `src`, the function can be called with `0` as the
% second argument.
%
% ### Inputs
%
% -`src`: a string-, char-, or cellstr-array of file system paths
%
% Optional:
% 
% -`isSearchRecursive [default=1]`: a scalar logical
%  
% ### Outputs 
%
% `mFiles`: String-vector list of *Documentor-compatible* .m files 
%
% ### Implementation notes: 
%
% **1. Warnings and errors:** 
%
% If identified, invalid paths or incompatible .m
% files suggested by the given inputs will elicit warning messages and, if no
% documentable .m files are found, an error is issued.
% 
% **2. Re:'Documentability':** 
%
% *Documentability* of a given .m file (i.e. whether it is included among the
% returned `mFiles` list) is defined according to the file-type returned by
% `mType = Informer.mfiletype( mFile )`:
%
% -`if mType == "NA"`: 
% The file is considered invalid and thereby omitted from `mFiles`. 
% This will be the case for nominal .m files (non-MATLAB files with
% '.m' file extensions) as well as invalid MATLAB files (source files with
% buggy implementations that preclude assessment with `Informer.mfiletype`).
%
% -`if mType == "method"`: 
% Standalone source files pertaining to class-methods are, likewise, omitted from
% `mFiles` as methods are included as part of the overall class documentation,
% rather than documented separately.
%
% **3. .m file search:** 
%
% When a directory figures among the entries of the input `src` paths,
% `Documentor.findfiles` searches for any .m files contained therein by
% wrapping to the standalone function *findfiles.m*, effectively calling:
%     
%    [~,mFiles] = findfiles( src( isfolder(src) ), "*.m", isSearchRecursive ) ;
%
% (If `src` contains multiple directories, findfiles.m will be called iteratively,
% with `mFiles` accordingly appended.)
%
% ### References
%
% See also
%
% - findfiles.m (standalone function)
% - Informer.mfiletype
    arguments
        src {mustBeStringOrCharOrCellstr} ;
        isSearchRecursive(1,1) {mustBeBoolean} = true ;
    end


%% Find .m files:
mFiles = src( isfile( src ) ) ;

if any( isfolder( src ) ) % dir search

    dirsIn = src( isfolder( src ) ) ;    

    for iDir = 1 : length( dirsIn )
        
        [~, mFilesiDir] = findfiles( dirsIn(iDir), "*.m", isSearchRecursive ) ;
        
        if ~isempty(mFilesiDir) && ~isequal( mFilesiDir, "")
            mFiles = [ mFiles ; mFilesiDir ] ;
        end
    end

end

if isempty(mFiles) || isequal( mFiles, "")
    error('*** Failed to find any .m files using the given inputs ***') ;
end

%% Validate .m files:
mTypes       = Informer.mfiletype( mFiles ) ;
invalidFiles = mFiles( mTypes=="NA" ) ;

if numel( invalidFiles ) > 0 
    display( invalidFiles ) ; % in case the list is long, print it above the warning 
    warning( '*** The above .m files are invalid (possibly incomplete) and cannot be documented  ***' ) ;
    
    % Remove invalid files
    mFiles( mTypes=="NA" ) = [] ;
    mTypes( mTypes=="NA" ) = [] ;
end

%% Remove any standalone class-method files
mFiles( mTypes=="method" ) = [] ;

assert( numel( mFiles ) > 0, 'Failed to find any Documentor-compatible .m files.') ;

end

function [ srcValid ] = validateinputpaths( srcIn ) 
%Validate input paths 

srcIn    = strip( string( srcIn ) ) ; 
srcIn    = srcIn(:) ; 
isPath = [ isfile( srcIn ) | isfolder( srcIn ) ] ;

if all( isPath ) % Easiest case: srcIn consists entirely of valid paths already accessible to MATLAB:
    srcValid = abspath( srcIn ) ; 
    return ;

else 
   ************
    % 2nd easiest case: srcIn consists entirely of valid paths already accessible to MATLAB:
    % src paths (if complete) could be good but not on the current MATLAB path, quickly:

    error( strcat( 'Not an existing file or directory path:\n', join( A(~isPath),"\n") ), '%s' ) ;
end

if all( isfile( src ) | isfolder( src ) )
    srcValidated = abspath( src ) ;
else
    try
        srcValidated = abspath( src ) ;
    end

% In case relative path(s) were provided, convert to absolute:
% NOTE: abspath() will throw an error if called with an invalid path---
% likely undesirable (e.g. if *most* paths are files, but a few are missing).
% Instead, use try/catch: log indices of any invalid paths and issue a warning if present
iInvalid = [] ;

for iPath = 1 : numel( src )
    try
        src( iPath ) = abspath( src(iPath) ) ;  
    catch Me
        iInvalid( end+1 ) = iPath ; 
    end
end

if numel( iInvalid ) == numel( src ) 
    error( 'Failed to validate any of the input paths.' ) ;

elseif iInvalid > 0
    display( src(iInvalid) ) ; % in case the list is long, print it above the warning 
    warning( '*** Input paths listed above do not exist or could otherwise not be validated ***' ) ;

badPath   =  ;
warning( strcat(, src(iPath) ) ) ;
if nBadPaths > 0
    src( iBadPaths ) = [] ;
end

end
