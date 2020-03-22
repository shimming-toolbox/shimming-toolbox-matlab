function [paths, List] = findfiles( sFolder, sPattern, isRecursive, isExcludingHidden, returnType )
% FINDFILES Search a directory for filenames matching a pattern
%     
%    [paths, List] = findfiles( sFolder, sPattern, isRecursive, isExcludingHidden, returnType )
% 
% Looks for files and/or subfolders in `sFolder` with names matching `sPattern`
% by calling Matlab function [dir] and returns the file paths as elements of a
% string column vector `paths`. The single-element structs output by `dir()`
% are arrayed and returned as `List.
% 
% __INPUTS__ 
%     
%   sFolder=["."]
%     The base directory of the search.    
%
%   sPattern=["*.*"] 
%     The searchPattern of interest. If provided as a string array, patterns
%     are searched successively. (The default corresponds to including  all
%     files with explicit with explicit file extensions.)
%
%   isRecursive=[true|1]
%     Toggle to include (1) or exclude (0) subdirectories in the search.
%
%   isExcludingHidden=[true|1]
%     Toggle to include (1) or exclude (0) hidden files (i.e. for
%     Unix: filenames beginning with ".")
%
%   returnType=["files"]
%     Selects what type of path elements are retained in the two outputs : 
%     Options are: "files", "folders", or "both". 
% 
% ETC 
%
%   For more info, refer to the documentation for
%   [dir](https://www.mathworks.com/help/matlab/ref/dir.html)
%
% See also DIR
    arguments
        sFolder(1,:) { mustBeStringScalarOrCharVector, mustBeFolder } = "." ;
        sPattern  {mustBeStringOrCharOrCellstr} = "*.*" ;    
        isRecursive(1,1) {mustBeBoolean} = true ;
        isExcludingHidden(1,1) {mustBeBoolean} = true ;
        returnType(1,:) { mustBeStringScalarOrCharVector, ...
            mustBeMember( returnType, ["files" "folders" "both"] ) } = "files" ;
    end

%% Cast as string
sPattern  = string( sPattern ) ;
searchStr = string( fullfile( sFolder, filesep ) ) ;

if isRecursive
    searchStr = fullfile( searchStr, "**", filesep )  ;
end

%% Search
List = dir( strcat(searchStr, sPattern(1) ) ) ;

for iPattern = 2 : numel( sPattern )
    List = [ List ; dir( strcat( searchStr, sPattern(iPattern) ) ) ] ;
end

[paths, iUnique] = unique( fullfile( string({List(:).folder})', string({List(:).name})' ) ) ;
List             = List(iUnique) ; 

if isExcludingHidden 
    if isunix
        iHidden = contains( paths, filesep + "." ) ;
    else 
        [~, Attributes] = arrayfun( @fileattrib, paths ) ;
        iHiddenFolders  = [ Attributes.hidden(:) & Attributes(:).directory ] ;
        iHiddenFiles    = [ Attributes.hidden(:) & ~Attributes(:).directory ] ;
        iInclude        = ~[ iHiddenFiles | iHiddenFolders ] ; 
        % ^redundant if 'hidden' is attributed recursively. not using a PC so idk.)
    end
    iInclude = ~iHidden ;
else
    iInclude = true( size( paths ) ) ; 
end

switch returnType
    case "files"
        iInclude = [ iInclude & ~[ List(:).isdir ]' ] ;
    case "folders"
        iInclude = [ iInclude & [ List(:).isdir ]' ] ; 
    % otherwise % (both)
end

paths = unique( paths( iInclude ) ) ;
List  = List( iInclude ) ;

%% Recast to input type
switch class( sFolder )
    case 'char'
        paths = char( paths ) ;
    case 'cell'
        paths = cellstr( paths ) ;
    otherwise
        paths = string( paths ) ; % should already be true
end

end
