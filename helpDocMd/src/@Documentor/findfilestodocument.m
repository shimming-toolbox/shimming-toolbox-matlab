function [mFiles] = findfilestodocument( Dr, pathIn )
%FINDFILESTODOCUMENT Return list of .m files to document from directory search
% 
% FINDFILESTODOCUMENT searches a directory for .m files and then removes any
% class method files from the list (methods are included as part of the
% overall class documentation).
%
% ### Syntax
%    
%    mFiles = Dr.FINDFILESTODOCUMENT( pathIn )
%
% ### Implementation details
%
% FINDFILESTODOCUMENT wraps to findfiles with the function call:
%     
%    [~,mFiles] = findfiles( folder, "*.m", Dr.isSearchRecursive ) ;
%
% .m file types are then determined using Informer.mfiletype( mFiles ) and, if
% present in the list, methods .m files are removed.

errMsg = ['Input argument to Documentor must be a file or directory path string ' ...
            ' or a 1-D array of file path strings'] ;

%% check input path(s)
if all( ischar( pathIn ) ) || iscellstr( pathIn )
    pathIn = string( pathIn ) ;
elseif ~isstring( pathIn ) || ndims( pathIn ) > 2
   error( errMsg  ) ;
end

pathIn    = strtrim( pathIn ) ;
badPaths  = "" ;
iBadPaths = 0 ;
nBadPaths = 0 ;

for iPath = 1 : numel( pathIn )
% NOTE: the call to abspath will throw an error if these are not valid paths.
% This may not be desirable (e.g. if *most* paths are files, but a few are
% missing). Therefore use try/catch to issue a warning about the missing files.
    try
        % convert to fullpath in case the relative path was provided
        pathIn(iPath) = abspath( pathIn(iPath) ) ;  
    catch Me
        iBadPaths = iPath ; % copy index
        nBadPaths = nBadPaths + 1 ;
        badPaths( nBadPaths ) = pathIn(iPath) ;
        warning( strcat("Input path does not exist: ", badPaths( nBadPaths ) ) ) ;
    end
end

if nBadPaths > 0
   pathIn(iBadPaths) = [] ;
end

assert( numel( pathIn ) > 0, "No valid paths found" ) ;

isInputFolder = any( isfolder( pathIn ) ) ;

if isInputFolder
    assert( numel( pathIn ) == 1, errMsg ) ;
    [~,mFiles] = findfiles( pathIn, "*.m", Dr.isSearchRecursive ) ;
else
    mFiles = pathIn ;
end

%% if present in mFiles list, remove standalone class methods files:
mTypes = Informer.mfiletype( mFiles ) ;

methodFiles = mFiles( mTypes=="method" ) ;

mFiles( mTypes=="method" ) = [] ;

% return a column vector
if size( mFiles, 2 ) > size( mFiles, 1 )
    mFiles = mFiles' ;
end

if numel( mFiles ) < 1
    warning('Failed to find any .m files to document') ;
end

end
