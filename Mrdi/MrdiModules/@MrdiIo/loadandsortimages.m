function [Imgs, Hdrs] = loadandsortimages( imgPath, fileExt, nBytesMax )
%LOADANDSORTIMAGES Load and sort image files
%
% LOADANDSORTIMAGES searches a directory for supported image files and sorts the
% images+headers into 2 cell arrays which respective contain numeric arrays of
% images, and struct arrays of headers. 
%
% Should the directory contain images from multiple acquisition series, each
% series is a separate element of the returned cells.
%
% Sorting of the returned arrays occurs according to: slice, echo time,
% measurement number, channel number.
% For the image arrays, these dimensions are respectively 3, 4, 5, and 6;
% For the headers, they are respectively 1,2,3, and 4.
% 
% Usage
%
%   [Imgs, Hdrs] = LOADANDSORTIMAGES( imgPath, fileExt, nBytesMax ) 
%
% Inputs
%
%   imgPath [default: './' ]
%       directory to search for images (char)
%
%   fileExt [default: { '.dcm' ; '.IMA' }]
%       file extensions of images (char or char-containing cell array)
%       (fileExt must be a member of the given default values)
%
%   nBytesMax [default: 2 * 1E9]
%       byte limit on load. Function aborts if limit is exceeded.
%       (Note: matlab 'memory' function exists to query available memory but
%       currently only works on Windows)
    arguments
        imgPath(1,:) char = pwd ;
        fileExt {mustBeMember( fileExt, { '.dcm' ; '.IMA' } )} = MrdiIo.supportedInputs ;
        nBytesMax(1,1) {mustBePositive} = 2*( 1E9 ) ;
    end

Imgs = [] ;
Hdrs = [] ;

%% ----- 
% Find files
List = MrdiIo.findimagefiles( imgPath, fileExt ) ;

if isempty(List)
   return ;
end

assert( sum( [ List(:).bytes ] )  < nBytesMax, [ 'Byte limit exceeded [%s GB].\n ' ...
    'Select a new image path or byte limit and run again.' ], num2str( nBytesMax/1E9 ) ) ;

%% ----- 
% Extract specific file types (e.g. Dicom) as they will need to be treated differently
filenames = string( [ { ( List(:).name ) } ] );

%% ----- 
DicomList = List( endsWith( filenames, MrdiIo.dicomExt ) ) ;

if ~isempty( DicomList )
    [Imgs, Hdrs] = MrdiIo.loadandsortdicoms( DicomList ) ;
end

%% ----- 
% Etc. (if other supported file types are added)

end
