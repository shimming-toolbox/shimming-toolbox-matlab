function [List] = findimagefiles( sFolder, fileExt, isRecursive )
%FINDIMAGEFILES Returns list of image files from `dir()` search
%      
%      List = findimagefiles( sFolder, fileExt, isRecursive )
% 
% Looks for image files in `sFolder` with extensions `fileExt` and returns
% the file matches as elements of a 1-D struct array `List`.
%
% __INPUTS__
%   
%   sFolder=["."]  
%     Base directory of the search as a string scalar or char vector. 
%
%   fileExt=[img.Maker.supportedInputs]  
%     File extension(s) of interest as a string vector (e.g. `[".dcm", ".IMA" ]`)
%
%   isRecursive=[false|0]  
%     Toggle to include (1) or exclude (0) subdirectories.
%
% __ETC__ 
%
% `findimagefiles` wraps to MATLAB function `dir()`
%
% See also 
% DIR — https://www.mathworks.com/help/matlab/ref/dir.html
% img.Maker.loadandsortimages
    arguments
        sFolder(1,:) {mustBeFileOrFolder} = "." ;
        fileExt string  = img.Maker.supportedInputs ;
        isRecursive(1,1) {mustBeBoolean} = false ;
    end

%% Cast as string
searchStr = string( fullfile( sFolder, filesep ) ) ;

if isRecursive
    searchStr = fullfile( searchStr, "**", filesep )  ;
end

searchStr = searchStr + "*" ;

%% Search
List = dir( strcat( searchStr, fileExt(1) ) ) ;

for iPattern = 2 : numel( fileExt )
    List = [ List ; dir( strcat( searchStr, fileExt(iPattern) ) ) ] ;
end

[paths, iUnique] = unique( fullfile( string({List(:).folder})', string({List(:).name})' ) ) ;
List             = List(iUnique) ; 

%% Exclude directories and hidden files
if isunix
    List = List( ~[ List(:).isdir ]' & ~contains( paths, filesep + "." ) ) ;
else 
    [~, Attributes] = arrayfun( @fileattrib, paths ) ;
    List = List( ~[ Attributes.directory(:) ] & ~[ Attributes.hidden(:) ] ) ;
end

display( [ num2str(length(List)) ' images found [totaling ' ...
           num2str( sum([ List(:).bytes ]) ) ' bytes].' ] ) ; 

end
