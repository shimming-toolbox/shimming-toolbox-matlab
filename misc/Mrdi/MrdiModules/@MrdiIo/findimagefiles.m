function [List] = findimagefiles( imgPath, fileExt, isRecursive )
%FINDIMAGEFILES Return list of image files from dir() search
% 
% Wraps to dir() and returns List: 1-D struct array of image files
% 
% List = FINDIMAGEFILES( imgPath, fileExt )
% 
% Inputs
%   
%   imgPath [default: './']
%       search directory
%
%   fileExt [default: MrdiIo.supportedInputs]
%       image file extension of interest
%
%   isRecursive [default: TRUE]
%       also searches subdirectories if true
%
% See also dir, MrdiIo.loadandsortimages()
    arguments
        imgPath(1,:) char = pwd ;
        fileExt {mustBeMember( fileExt, { '.dcm' ; '.IMA' } )} = MrdiIo.supportedInputs ;
        isRecursive {mustBeNumericOrLogical} = true ;
    end

if isRecursive
    searchStr = [imgPath filesep '**' filesep '*'] ;
    display( [ 'Searching recursively for images in ' imgPath ] ) ; 
else
    searchStr = [imgPath filesep '*'] ;
    display( [ 'Searching for images in ' imgPath ] ) ; 
end

if ischar( fileExt )
    fileExt = { fileExt } ;
end

List = dir( [searchStr fileExt{1} ] ) ;

if length( fileExt ) > 1 % look for all file extensions
    for iExt = 2 : length( fileExt )
        List_iExt = dir( [searchStr fileExt{iExt} ] ) ;
        List(end+1:end+length(List_iExt)) = List_iExt ;
    end
end

nImg   = length( List ) ;
nBytes = sum( [ List(:).bytes ] ) ;

display( [ num2str(nImg) ' images found [totaling ' num2str( nBytes ) ' bytes].' ] ) ; 

end

