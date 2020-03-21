function [paths, List] = findfiles( sFolder, sFilePattern, isRecursive, isExcludingHidden )
% FINDFILES Search a directory for filenames matching a pattern
% 
%    [paths, List] = findfiles( sFolder, sFilePattern, isRecursive, isExcludingHidden )
% 
% Looks for files in `sFolder` with filenames matching `sFilePattern` by
% calling Matlab function [dir] and returns the file paths as elements of a
% string column vector `paths`. Structs output by `dir()` are arrayed and
% returned as `List.
% 
% INPUTS
%     
%   sFolder=["."]
%     The base directory of the search.    
%
%   sFilePattern=["*.*"] 
%     The searchPattern of interest. If provided as a string array, patterns
%     are searched successively. (The default corresponds to including  all
%     files with explicit with explicit file extensions.)
%
%   isRecursive=[true|1]
%     Toggle to include (1) or exclude (0) subdirectories in the search.
%
%   isExcludingHidden=[true|1]
%     Toggle to include (1) or exclude (0) hidden files (i.e. for
%     Unix: filenames beginning with "." )  
%
% ETC 
%
%   For more info, refer to the documentation for
%   [dir](https://www.mathworks.com/help/matlab/ref/dir.html)
%
% See also DIR
    arguments
        sFolder(1,:) { mustBeStringScalarOrCharVector, mustBeFolder } = "." ;
        sFilePattern  {mustBeStringOrCharOrCellstr} = "*.*" ;    
        isRecursive(1,1) {mustBeBoolean} = true ;
        isExcludingHidden(1,1) {mustBeBoolean} = true ;
    end

searchStr = strcat( sFolder, filesep )  ;

if isRecursive
    searchStr = strcat( searchStr, '**', filesep )  ;
end

sFilePattern = string( searchPattern ) ;

List = dir( strcat(searchStr,sFilePattern(1)) ) ;

for iPattern = 2 : numel( sFilePattern )
    List = [ List ; dir( strcat( searchStr, sFilePattern(iPattern) ) ) ] ;
end

%% remove folders
List = List( ~[ List(:).isdir ] ) ;

%% remove 'hidden' files
if isExcludingHidden 
    iHidden = startsWith( string( { List(:).name } ), "." ) ;
    List    = List( ~iHidden ) ;
end

paths = "" ;
for iFile = 1 : length( List ) 
    paths( iFile ) = string( fullfile( List(iFile).folder, List(iFile).name ) ) ;
end

paths = paths' ;

end
