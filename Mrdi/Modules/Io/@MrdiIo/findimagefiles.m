%%  Returns list of image files from dir() search
% 
% Wraps to dir() and returns List: 1-D struct array of image files
% 
% List = MrdiIo.findimagefiles( searchDir, searchExt, isSearchRecursive )
% 
% Inputs
%   
%   searchDir [default: './']
%       search directory
%
%   searchExt [default: [".dcm", ".IMA" ]]
%       image file extension of interest
%
%   isSearchRecursive [default: TRUE]
%       also searches subdirectories if true
%
% See also dir, MrdiIo.loadandsortimages()

function [List] = findimagefiles( searchDir, searchExt, isSearchRecursive )
%FINDIMAGEFILES Returns list of image files from dir() search
% 
% Wraps to dir() and returns List: 1-D struct array of image files
% 
% List = MrdiIo.findimagefiles( searchDir, searchExt, isSearchRecursive )
% 
% Inputs
%   
%   searchDir [default: './']
%       search directory
%
%   searchExt [default: [".dcm", ".IMA" ]]
%       image file extension of interest
%
%   isSearchRecursive [default: TRUE]
%       also searches subdirectories if true
%
% See also dir, MrdiIo.loadandsortimages()
    arguments
        searchDir(1,:) = './' ;
        searchExt {mustBeMember( searchExt, [".dcm", ".IMA" ] ) } = MrdiIo.FileExt.supportedInputs ;
        isSearchRecursive {mustBeNumericOrLogical} = true ;
    end

%% Check inputs
narginchk(0,3) ;
validateattributes( searchDir, {'string', 'char'}, {'scalartext'}, mfilename, 'searchDir', 1 ) ;

if isSearchRecursive
    searchStr = [searchDir filesep '**' filesep '*'] ;
    display( [ 'Searching recursively for images in ' searchDir ] ) ; 
else
    searchStr = [searchDir filesep '*'] ;
    display( [ 'Searching for images in ' searchDir ] ) ; 
end

if ischar( searchExt )
    searchExt = { searchExt } ;
end

List = dir( [searchStr searchExt{1} ] ) ;

if length( searchExt ) > 1 % look for all file extensions
    for iExt = 2 : length( searchExt )
        List_iExt = dir( [searchStr searchExt{iExt} ] ) ;
        List(end+1:end+length(List_iExt)) = List_iExt ;
    end
end

nImg   = length( List ) ;
nBytes = sum( [ List(:).bytes ] ) ;

display( [ num2str(nImg) ' images found [totaling ' num2str( nBytes ) ' bytes].' ] ) ; 

end

