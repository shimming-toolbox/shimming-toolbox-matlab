function [mFiles] = findfiles( src, isRecursive )
% FINDFILES Return list of .m files to document from directory search
%    
%    mFiles = Documentor.findfiles( src )
%    mFiles = Documentor.findfiles( src, isRecursive )
%
% Returns a list of *documentable* .m files `mFiles` based on the set of file
% system paths specified in `src`.
%
% `src` 
% :    a [string-, char-, or cellstr-] array of paths to .m source files,
%      and/or directories in which to search for them.
% 
% `isRecursive`
% :    a Boolean toggle to include subdirectories in the search. To restrict
%      the search directories explicitly included in `src`, the function
%      should be called with `0` as the second argument.
%      [default=`true`]
%
% `mFiles`
% :     String-vector list of *Documentor-compatible* .m files 
%
% When source files are included in the input, the function verifies that the
% Documentor class will indeed be able to document them (see Note 2. below).
%
% When `src` elements point to directories, the function first searches the
% directories for .m files, followed by the compatibility check (see Note 3.
% below).
%
% __NOTES__
%
% **1. Warnings and errors:** 
%
% If identified, invalid paths or incompatible .m files suggested by the given
% inputs will elicit warning messages and, if no documentable .m files are
% found whatsoever, an error is issued.
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
% For more info, refer to the documentation for Informer.mfiletype
%
% -`if mType == "method"`: 
% Standalone source files pertaining to class-methods are, likewise, omitted
% from `mFiles` as methods are included as part of the overall class
% documentation, rather than documented separately.
%
% **3. .m file search:** 
%
% When a directory figures among the entries of the input `src` paths,
% `Documentor.findfiles` searches for any .m files contained therein by
% wrapping to the standalone function *findfiles.m*, effectively calling:
%     
%    [~,mFiles] = findfiles( src( isfolder(src) ), "*.m", isRecursive ) ;
%
% (If `src` contains multiple directories, findfiles.m will be called iteratively,
% with `mFiles` accordingly appended.)
%
% NOTE/TODO: Handling the case where multiple directories are included and
% `isRecursive == true`. There is the potential here to include directories
% multiple times (can be easily filtered out ahead of time using
% Pathologist.subfolders -- just need to implement the filter). More
% importantly, a folder the user didn't really want included might be included
% anyway if it happens to be a sub of one that was specified... Should a
% warning be issued?
% 
% __ETC__
%
% - standalone function: findfiles.m 
% - Informer.mfiletype
%
% See also
% FINDFILES
    arguments
        src {mustBeStringOrCharOrCellstr} ;
        isRecursive(1,1) {mustBeBoolean} = true ;
    end

% TODO: make Pathologist initialize `data` w/full paths when possible...
src    = abspath( src ) ;
P      = Pathologist(src) ; 
mFiles = P.files ;
dirsIn = P.folders ;

%% Validate input paths exist 
if all( ~P.isvalid )
    error( [ 'helpDocMd:Documentor.findfiles:invalidPath', ...
             'Failed to validate any of the given input paths..'] ) ;

elseif any( ~P.isvalid )
    % in case the list is long, print it above the warning
    display( src( ~P.isvalid) ) ;  
    warning( '*** Input paths listed above do not exist or could otherwise not be validated ***' ) ;
end

if isRecursive && numel( dirsIn ) > 1
% TODO filter potential repeat-dirs & verify it might actually be a problem
% before issuing a warning...
% if any( Dirs.folders == Dirs.baseDir ) 
    warning( [ 'Mutliple folders included with recursive option:'  ...
        'Search may be slow and/or inadvertently include undesired directories...' ] ) ;
end

%% Directory search .m :
if ~isempty( dirsIn ) 
    for iDir = 1 : length( dirsIn )
        mFilesiDir = findfiles( dirsIn(iDir), "*.m", isRecursive ) ;
        if ~isempty(mFilesiDir) && ~isequal( mFilesiDir, "")
            mFiles = [ mFiles ; mFilesiDir ] ;
        end
    end
end

if isempty(mFiles) || isequal( mFiles, "" )
    error( [ 'helpDocMd:Documentor.findfiles:invalidPath', ...
             'Failed to find any .m files using the given inputs.'] ) ;
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
