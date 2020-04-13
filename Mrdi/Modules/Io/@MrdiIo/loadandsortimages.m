function [Imgs, Hdrs] = loadandsortimages( List, nBytesMax )
%LOADANDSORTIMAGES Load and sort image files
%
%     [Imgs, Hdrs] = MrdiIo.loadandsortimages( List ) 
%     [Imgs, Hdrs] = MrdiIo.loadandsortimages( List, nBytesMax ) 
%
% MrdiIo.loadandsortimages accepts a list of image files to read (i.e. returned from
% MrdiIo.findimagefiles() ), loads the supported file types (), .. and sorts
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
% __OPTIONS__ 
%
%   nBytesMax [default: 2 * 1E9]
%       byte limit on load. Function aborts if limit is exceeded.
%       (Note: matlab 'memory' function exists to query available memory but
%       currently only works on Windows)
% 
% __ETC__ 
% See also 
% MrdiIo.findimagefiles, MrdiIo.make
    arguments
        List(:,1) struct ;
        nBytesMax(1,1) {mustBePositive} = MrdiIo.defaults('nBytesMax') ;
    end

%% Check inputs
narginchk(1,2) ;

Imgs = [] ;
Hdrs = [] ;

if isempty(List)
    warning('List of image files was empty.') ;
    return ;
elseif ~isfield( List, 'bytes' ) || ~isfield( List, 'name' )
    error('Invalid input: List should be a struct returned from MrdiIo.findimagefiles()')
end

%% Check total file size 
assert( sum( [ List(:).bytes ] )  < nBytesMax, [ 'Byte limit exceeded [%s GB].\n ' ...
    'Select a new image path or byte limit and run again.' ], num2str( nBytesMax/1E9 ) ) ;

%% Extract specific file types (e.g. Dicom) as they will need to be handled differently
filenames = string( [ { ( List(:).name ) } ] );

%% Dicoms: 
if any( ismember( fileExt, MrdiIo.FileExt.dicom ) ) ; 

    DicomList = List( endsWith( filenames, MrdiIo.FileExt.dicom ) ) ;

    if ~isempty( DicomList )
        [Imgs, Hdrs] = MrdiIo.loadandsortdicoms( DicomList ) ;
    end
end

%% Etc: (if other supported file types are added)

end
