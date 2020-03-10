function [List, paths] = findfiles( searchDir, searchPattern, isSearchRecursive, isExcludingHidden )
% FINDFILES Returns list of files from directory search
%
% ### Syntax
%    
%    [List, paths] = FINDFILES( searchDir, searchPattern, isSearchRecursive )
%
% Wraps to dir() and returns List: 1-D struct array of Matlab files. Full file
% paths from the List.folder and List.name fields are output as a string vector
% in the second return argument.
% 
% ### Inputs
%      
% - searchDir [default="./"]  
%   Search directory
%
% - searchPattern [default, all files with explicit file extensions: "*.*"]  
%   searchPattern of interest (can be an array, in which case the patterns are search successively)
%
% - isSearchRecursive [default=1]  
%   Boolean toggle to include subdirectories in search. 
%
% - isExcludingHidden [default=1]  
%   Boolean toggle to exclude hidden files (i.e. for Unix: filenames beggining with "." )
%
% ### References ###
%
% See also
%
% <https://www.mathworks.com/help/matlab/ref/dir.html dir>
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
