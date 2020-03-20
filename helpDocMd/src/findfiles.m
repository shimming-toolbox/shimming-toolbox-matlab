function [List, paths] = findfiles( searchDir, searchPattern, isSearchRecursive, isExcludingHidden )
% FINDFILES Search a directory for filenames matching a pattern
%     
%    [List, paths] = findfiles( searchDir, searchPattern, isSearchRecursive, isExcludingHidden )
%
% Returns the list of files in `searchDir` and its subdirectories with file
% names matching `searchPattern`. 
%
% __Description__
% Wraps to Matlab function `dir` [1] and returns List: 1-D struct array of
% Matlab files. Full file paths from the List.folder and List.name fields are
% output as a string vector in the second return argument.
% 
% Each of the input arguments is optional, so `findfiles( )` will yield the
% same result as `findfiles( ".", "*.*", 1, 1  )`
% __OPTIONS__                                           
%     
% `searchDir` The base directory of the search. [default="."]   
%
% `searchPattern` The searchPattern of interest. If provided as a string array,
% patterns searched successively
% [default, all files with explicit file extensions: "*.*"]  
%
% `isSearchRecursive` Toggle to include subdirectories in search when true
% [default=1]  
%
% `isExcludingHidden` Toggle to exclude hidden files when true.
% (i.e. for Unix: filenames beggining with "." )
% [default=1]  
%
% Inputs:
%  
%| Name              |`Default`   | {type} (size)   |  description |
%| ----------        | ---------- | ----------------| ----------  -------- |
%|`searchDir `       |`"."`       |{string}(1,1)    | Base directory of search
%|`searchPattern`    |`"*.*"      | {string}(any)   | e.g "*.m" for Matlab source code|
%|`isSearchRecursive`|`true`      | {logical} (1,1) |                         |
%|`isExcludingHidden`|`true`      | {logical} (1,1) |                         |
%
%
% __SEE__
% [1]: https://www.mathworks.com/help/matlab/ref/dir.html
%
% See also 
% DIR
    arguments
        searchDir(1,:) { mustBeStringScalarOrCharVector, mustBeFolder } = "." ;
        searchPattern  {mustBeStringOrCharOrCellstr} = "*.*" ;    
        isSearchRecursive(1,1) {mustBeBoolean} = true ;
        isExcludingHidden(1,1) {mustBeBoolean} = true ;
    end

searchStr = strcat( searchDir, filesep )  ;

if isSearchRecursive
    searchStr = strcat( searchStr, '**', filesep )  ;
end

searchPattern = string( searchPattern ) ;

List = dir( strcat(searchStr,searchPattern(1)) ) ;

for iPattern = 2 : numel( searchPattern )
    List = [ List ; dir( strcat( searchStr, searchPattern(iPattern) ) ) ] ;
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
